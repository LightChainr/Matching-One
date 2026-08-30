#!/usr/bin/env python3
"""Eliminate a one-gain ambient-H1 explanation of the N650 mixed residual."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import mpmath as mp

from analyze_p200_n650_context_morphism import zero_test
from score_p200_n650_mixed_join import PRIMARY, SECONDARY, load_inputs


ROOT = Path(__file__).resolve().parents[1]
JOINT_ORDER = tuple(PRIMARY) + tuple(SECONDARY)
DERIVED_ORDER = (
    "beta_ES",
    "H1_orthogonal_OS",
    "H1_residual_ED",
    "H1_residual_OD",
    "S_plane_determinant",
)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def repository_path(path: Path) -> str:
    return path.resolve().relative_to(ROOT).as_posix()


def jackknife_covariance(values: list[list[mp.mpf]]) -> list[list[mp.mpf]]:
    count = len(values)
    center = [mp.fsum(row[index] for row in values) / count for index in range(len(values[0]))]
    factor = mp.mpf(count - 1) / count
    return [
        [
            factor
            * mp.fsum(
                (row[first] - center[first]) * (row[second] - center[second])
                for row in values
            )
            for second in range(len(center))
        ]
        for first in range(len(center))
    ]


def _means(rows: list[dict]) -> tuple[list[mp.mpf], list[mp.mpf]]:
    total_samples = sum(row["samples"] for row in rows)
    primary = [
        mp.fsum(row["primary"][index] for row in rows) / (2 * total_samples)
        for index in range(4)
    ]
    ambient = [
        mp.fsum(row["secondary"][index] for row in rows) / (2 * total_samples)
        for index in range(4)
    ]
    return primary, ambient


def _derived(primary: list[mp.mpf], ambient: list[mp.mpf]) -> list[mp.mpf]:
    if ambient[0] == 0:
        raise ValueError("ambient_ES must be nonzero to fix the common gain")
    beta = primary[0] / ambient[0]
    return [
        beta,
        primary[2] - beta * ambient[2],
        primary[1] - beta * ambient[1],
        primary[3] - beta * ambient[3],
        primary[2] * ambient[0] - primary[0] * ambient[2],
    ]


def _float_matrix(values: list[list[mp.mpf]]) -> list[list[float]]:
    return [[float(value) for value in row] for row in values]


def _subset_test(
    names: list[str], mean: list[mp.mpf], covariance: list[list[mp.mpf]], indices: list[int]
) -> dict:
    submean = [float(mean[index]) for index in indices]
    subcovariance = [
        [float(covariance[first][second]) for second in indices] for first in indices
    ]
    return zero_test(names, submean, subcovariance)


def analyze(batch_path: Path, metadata_path: Path, prediction_path: Path) -> dict:
    mp.mp.dps = 80
    metadata, prediction_hash, rows = load_inputs(batch_path, metadata_path, prediction_path)
    if metadata["samples"] != 20000 or metadata["batches"] != 100:
        raise ValueError("this result is frozen to the revealed N650 20k archive")

    primary, ambient = _means(rows)
    joint = primary + ambient
    derived = _derived(primary, ambient)

    joint_leave_one: list[list[mp.mpf]] = []
    derived_leave_one: list[list[mp.mpf]] = []
    for omitted in range(len(rows)):
        leave = rows[:omitted] + rows[omitted + 1 :]
        p_leave, h_leave = _means(leave)
        joint_leave_one.append(p_leave + h_leave)
        derived_leave_one.append(_derived(p_leave, h_leave))

    joint_covariance = jackknife_covariance(joint_leave_one)
    derived_covariance = jackknife_covariance(derived_leave_one)

    # A scalar endpoint-H1 mechanism predicts collinearity of the common
    # (matching-even, matching-odd) vectors.  The determinant is a gain-free
    # certificate of the component perpendicular to ambient H1.
    determinant = derived[4]
    determinant_se = mp.sqrt(derived_covariance[4][4])
    determinant_z = determinant / determinant_se
    determinant_p = mp.erfc(abs(determinant_z) / mp.sqrt(2))

    split_half = []
    for training_parity in (0, 1):
        training = [row for row in rows if row["batch"] % 2 == training_parity]
        heldout = [row for row in rows if row["batch"] % 2 != training_parity]
        p_train, h_train = _means(training)
        p_test, h_test = _means(heldout)
        beta_train = p_train[0] / h_train[0]
        split_half.append(
            {
                "training_batch_parity": training_parity,
                "heldout_batch_parity": 1 - training_parity,
                "beta_ES_training": float(beta_train),
                "heldout_residual": {
                    "H1_orthogonal_OS": float(p_test[2] - beta_train * h_test[2]),
                    "H1_residual_ED": float(p_test[1] - beta_train * h_test[1]),
                    "H1_residual_OD": float(p_test[3] - beta_train * h_test[3]),
                },
            }
        )

    prediction = json.loads(prediction_path.read_text(encoding="utf-8"))
    local = prediction["toy_exact_normalization"]["p_ref"]
    return {
        "schema": "matching-one.p200-n650-h1-orthogonal-residual.v1",
        "issues": [200, 249, 255],
        "status": "post_reveal_model_elimination_no_new_production",
        "input": {
            "batch_csv": repository_path(batch_path),
            "metadata": repository_path(metadata_path),
            "prediction": repository_path(prediction_path),
            "samples": metadata["samples"],
            "batches": metadata["batches"],
            "sha256": {
                "batch_csv": file_sha256(batch_path),
                "metadata": file_sha256(metadata_path),
                "prediction": prediction_hash,
            },
        },
        "observable_boundary": {
            "local_incidence_subtraction": "already configurationwise in R_c=J_full,c-J_local,c",
            "exact_local_baseline_per_source_fiber": {
                "black": float(local["mean_decimal"][0]),
                "white": float(local["mean_decimal"][1]),
            },
            "ambient_H1": "same-stream convention-labelled endpoint rank diagnostic",
            "path_or_state_memory": "not identified; no ordered intermediate state is stored",
        },
        "joint_same_stream_state": {
            "state_order": list(JOINT_ORDER),
            "mean": [float(value) for value in joint],
            "delete_one_covariance": _float_matrix(joint_covariance),
        },
        "one_gain_H1_projection": {
            "hypothesis": "(primary_ES,primary_OS)=beta*(ambient_ES,ambient_OS)",
            "gain_fixed_by": "beta=primary_ES/ambient_ES",
            "derived_state_order": list(DERIVED_ORDER),
            "derived_mean": [float(value) for value in derived],
            "derived_delete_one_covariance": _float_matrix(derived_covariance),
            "gain_free_collinearity_certificate": {
                "formula": "primary_OS*ambient_ES-primary_ES*ambient_OS",
                "estimate": float(determinant),
                "se": float(determinant_se),
                "z": float(determinant_z),
                "two_sided_p_value": float(determinant_p),
                "log10_p_value": float(mp.log10(determinant_p)),
            },
            "orthogonal_common_mode": _subset_test(
                ["H1_orthogonal_OS"], derived, derived_covariance, [1]
            ),
            "geometry_difference_after_same_gain": _subset_test(
                ["H1_residual_ED", "H1_residual_OD"],
                derived,
                derived_covariance,
                [2, 3],
            ),
            "all_unfitted_directions": _subset_test(
                ["H1_orthogonal_OS", "H1_residual_ED", "H1_residual_OD"],
                derived,
                derived_covariance,
                [1, 2, 3],
            ),
            "split_half_direction_check": split_half,
            "decision": (
                "reject one shared endpoint-H1 gain; a large matching-odd static mixed-factor "
                "direction remains, while ED/OD remain null"
            ),
            "saturated_alternative": (
                "separate unconstrained even/odd gains are not identifiable with only the two common S rows"
            ),
        },
        "interpretation": {
            "changed_mechanism_space": (
                "the N650 static interaction is not exhausted by one ambient-H1 endpoint direction"
            ),
            "does_not_show": (
                "chronological memory, noncommuting joins, a Jordan block, or a new continuum field"
            ),
            "minimal_next_acquisition": {
                "design": "semantics-matched Gaussian-by-annulus 2x2 rectangle",
                "randomness": "one common random block for all four cells",
                "observer_tuple": [
                    "typed colour layer",
                    "ambient H1 before/after each context",
                    "local incidence J_local",
                    "connected R=J_full-J_local",
                    "ED/OD geometry label",
                ],
                "covariance": "archive all cell means and one complete delete-one covariance",
                "reason": (
                    "the N650 archive has no synonymous annulus context or ordered intermediate state"
                ),
            },
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batches", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--prediction", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = analyze(args.batches, args.metadata, args.prediction)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

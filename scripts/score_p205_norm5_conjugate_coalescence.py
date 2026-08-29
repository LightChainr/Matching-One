#!/usr/bin/env python3
"""Score the frozen Issue #205 same-parent norm-5 coalescence experiment.

The scorer keeps the preregistered fixed-p observable and the frozen affine
H4/H12/H8 residuals.  Within each size, C/A/B are evaluated under the common
priority field and their full 3x3 delete-one covariance is propagated into the
one-dimensional affine residual.  The two sizes use independent seeds.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Mapping, Sequence

import mpmath as mp
import yaml

from analyze_matching_parity_derivatives_fast import H, combine, obs, read, remove
from analyze_threshold_rank_orientation import read_histograms, validate_moments


SIZES = (325, 425)
GEOMETRY_ORDER = ("C", "A", "B")
MODEL_ORDER = ("H4", "H12", "H8")
EXPECTED = {
    325: {
        "C": (15, 10), "A": (17, 6), "B": (18, 1),
        "seed": 2026105501, "smith_C": (5, 65),
    },
    425: {
        "C": (20, 5), "A": (16, 13), "B": (19, 8),
        "seed": 2026105502, "smith_C": (5, 85),
    },
}


@dataclass(frozen=True)
class PairRun:
    n: int
    partner: str
    histogram_path: Path
    moments_path: Path
    metadata_path: Path
    metadata: Mapping[str, object]
    data: Mapping[tuple[int, str, int], H]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def orientation_rows_sha256(path: Path, orientation: str) -> str:
    """Hash the exact CSV data-row bytes for one orientation."""
    marker = f",{orientation},".encode()
    selected = b"".join(
        line for line in path.read_bytes().splitlines(keepends=True)
        if marker in line
    )
    if not selected:
        raise ValueError(f"{path}: no {orientation} rows")
    return hashlib.sha256(selected).hexdigest()


def load_contract(prediction_path: Path, experiment_path: Path) -> dict[str, object]:
    prediction = yaml.safe_load(prediction_path.read_text(encoding="utf-8"))
    experiment = yaml.safe_load(experiment_path.read_text(encoding="utf-8"))
    if prediction.get("status") != "prospective_before_noncyclic_C_target_reveal":
        raise ValueError("prediction is not the frozen prospective artifact")
    if experiment.get("status") != "pilot_ready_frozen_before_C_reveal":
        raise ValueError("experiment is not the frozen pilot contract")
    if tuple(experiment["score_order"][:3]) != (
        "joint_H4_fixed_p_ref_affine_residual",
        "fixed_H12_affine_residual",
        "fixed_H8_affine_residual",
    ):
        raise ValueError("frozen score order changed")
    order = tuple(prediction["fixed_angular_adversaries"]["frozen_order"])
    if ("H4",) + order != MODEL_ORDER:
        raise ValueError("frozen harmonic order changed")

    weights: dict[str, dict[int, tuple[int, int, int]]] = {
        name: {} for name in MODEL_ORDER
    }
    for n in SIZES:
        key = f"N{n}"
        weights["H4"][n] = tuple(
            int(value) for value in prediction["primary_H4"][key]["integer_residual_C_A_B"]
        )
        for name in ("H12", "H8"):
            weights[name][n] = tuple(
                int(value)
                for value in prediction["fixed_angular_adversaries"][name][key][
                    "integer_residual_C_A_B"
                ]
            )
        for name in MODEL_ORDER:
            if sum(weights[name][n]) != 0 or weights[name][n][0] == 0:
                raise ValueError(f"{name} N={n}: residual does not cancel H0")
    return {
        "prediction": prediction,
        "experiment": experiment,
        "p_ref": mp.mpf(str(prediction["observable"]["p_ref"])),
        "weights": weights,
    }


def load_pair(
    n: int, partner: str, histogram: Path, moments: Path, metadata_path: Path,
    contract: Mapping[str, object],
) -> PairRun:
    if n not in SIZES or partner not in ("A", "B"):
        raise ValueError("pair must be one frozen N and partner A/B")
    data = read(histogram)
    validate_moments(moments, read_histograms(histogram))
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    expected = EXPECTED[n]
    experiment = contract["experiment"][f"N{n}"]
    pair_index = 0 if partner == "A" else 1
    pair_contract = experiment["pairs"][pair_index]
    design = metadata["designs"][0]
    required = {
        "git_commit", "seed", "replica_counter_first",
        "replica_counter_last_exclusive", "samples_per_pair", "batches",
        "threads_requested",
    }
    missing = required - set(metadata)
    if missing:
        raise ValueError(f"N={n} C-{partner}: metadata lacks {sorted(missing)}")
    if tuple(design["first"]) != expected["C"] or tuple(design["second"]) != expected[partner]:
        raise ValueError(f"N={n} C-{partner}: geometry representation changed")
    if tuple(design["first_smith_invariants"]) != expected["smith_C"]:
        raise ValueError(f"N={n} C-{partner}: C Smith class changed")
    if int(metadata["seed"]) != expected["seed"]:
        raise ValueError(f"N={n} C-{partner}: seed changed")
    if [int(metadata["replica_counter_first"]), int(metadata["replica_counter_last_exclusive"])] != [
        int(value) for value in pair_contract["replica_counter"]
    ]:
        raise ValueError(f"N={n} C-{partner}: replica interval changed")
    pilot = contract["experiment"]["pilot"]
    if (
        int(metadata["samples_per_pair"]) != int(pilot["samples_per_pair"])
        or int(metadata["batches"]) != int(pilot["batches"])
        or int(metadata["threads_requested"]) != int(pilot["threads_per_job"])
    ):
        raise ValueError(f"N={n} C-{partner}: production allocation changed")
    if len(str(metadata["git_commit"])) != 40:
        raise ValueError(f"N={n} C-{partner}: incomplete commit provenance")
    batches = int(metadata["batches"])
    if sorted({key[2] for key in data}) != list(range(batches)):
        raise ValueError(f"N={n} C-{partner}: incomplete batch grid")
    for orientation, representation in (("first", expected["C"]), ("second", expected[partner])):
        rows = [row for key, row in data.items() if key[1] == orientation]
        if len(rows) != batches or {(row.a, row.b) for row in rows} != {representation}:
            raise ValueError(f"N={n} C-{partner}: {orientation} descriptor changed")
        if sum(row.samples for row in rows) != int(metadata["samples_per_pair"]):
            raise ValueError(f"N={n} C-{partner}: sample count mismatch")
    return PairRun(n, partner, histogram, moments, metadata_path, metadata, data)


def aligned_rows(run: PairRun, orientation: str) -> list[H]:
    return [run.data[key] for key in sorted(run.data) if key[1] == orientation]


def jackknife_covariance(rows: Sequence[Sequence[mp.mpf]]) -> list[list[mp.mpf]]:
    if len(rows) < 2 or not rows or not rows[0]:
        raise ValueError("jackknife covariance needs at least two nonempty rows")
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise ValueError("jackknife rows have inconsistent widths")
    count = len(rows)
    means = [mp.fsum(row[j] for row in rows) / count for j in range(width)]
    scale = mp.mpf(count - 1) / count
    return [[
        scale * mp.fsum(
            (row[i] - means[i]) * (row[j] - means[j]) for row in rows
        )
        for j in range(width)] for i in range(width)]


def score_size(ca: PairRun, cb: PairRun, p_ref: mp.mpf) -> dict[str, object]:
    if ca.n != cb.n or ca.partner != "A" or cb.partner != "B":
        raise ValueError("score_size requires aligned C-A and C-B runs")
    signature_fields = (
        "git_commit", "seed", "replica_counter_first",
        "replica_counter_last_exclusive", "samples_per_pair", "batches",
    )
    if any(ca.metadata[field] != cb.metadata[field] for field in signature_fields):
        raise ValueError(f"N={ca.n}: C-A/C-B common-field signature differs")
    c_hist_a = orientation_rows_sha256(ca.histogram_path, "first")
    c_hist_b = orientation_rows_sha256(cb.histogram_path, "first")
    c_mom_a = orientation_rows_sha256(ca.moments_path, "first")
    c_mom_b = orientation_rows_sha256(cb.moments_path, "first")
    if c_hist_a != c_hist_b or c_mom_a != c_mom_b:
        raise ValueError(f"N={ca.n}: duplicated C stream is not byte-identical")

    geometry_rows = {
        "C": aligned_rows(ca, "first"),
        "A": aligned_rows(ca, "second"),
        "B": aligned_rows(cb, "second"),
    }
    combined = {name: combine(rows) for name, rows in geometry_rows.items()}
    point = {name: obs(combined[name], p_ref)["M"] for name in GEOMETRY_ORDER}
    deleted = []
    for batch in range(len(geometry_rows["C"])):
        deleted.append([
            obs(remove(combined[name], geometry_rows[name][batch]), p_ref)["M"]
            for name in GEOMETRY_ORDER
        ])
    covariance = jackknife_covariance(deleted)
    return {
        "N": ca.n,
        "point": point,
        "covariance": covariance,
        "common_C_histogram_rows_sha256": c_hist_a,
        "common_C_moments_rows_sha256": c_mom_a,
        "pair_runs": (ca, cb),
    }


def residual_score(
    point: Mapping[str, mp.mpf], covariance: Sequence[Sequence[mp.mpf]],
    integer_weights: Sequence[int],
) -> dict[str, object]:
    scale = Fraction(1, int(integer_weights[0]))
    weights = [mp.mpf(weight * scale.numerator) / scale.denominator for weight in integer_weights]
    residual = mp.fsum(weights[i] * point[name] for i, name in enumerate(GEOMETRY_ORDER))
    variance = mp.fsum(
        weights[i] * covariance[i][j] * weights[j]
        for i in range(3) for j in range(3)
    )
    if variance <= 0:
        raise ValueError("affine residual variance is not positive")
    z = residual / mp.sqrt(variance)
    return {
        "integer_weights_C_A_B": list(integer_weights),
        "normalized_weights_C_A_B": [mp.nstr(value, 30) for value in weights],
        "residual": mp.nstr(residual, 30),
        "variance": mp.nstr(variance, 25),
        "standard_error": mp.nstr(mp.sqrt(variance), 20),
        "signed_z": mp.nstr(z, 18),
    }


def render(
    pairs: Mapping[int, tuple[PairRun, PairRun]], contract: Mapping[str, object],
    prediction_path: Path, experiment_path: Path,
) -> dict[str, object]:
    if int(pairs[325][0].metadata["seed"]) == int(pairs[425][0].metadata["seed"]):
        raise ValueError("the two frozen sizes must use independent seeds")
    by_size = {
        n: score_size(pairs[n][0], pairs[n][1], contract["p_ref"])
        for n in SIZES
    }
    models = []
    for name in MODEL_ORDER:
        size_scores = {}
        chi_square = mp.mpf(0)
        for n in SIZES:
            row = residual_score(
                by_size[n]["point"], by_size[n]["covariance"],
                contract["weights"][name][n],
            )
            size_scores[str(n)] = row
            chi_square += mp.mpf(row["signed_z"]) ** 2
        models.append({
            "name": name,
            "by_size": size_scores,
            "joint_chi_square": mp.nstr(chi_square, 20),
            "degrees_of_freedom": 2,
            "chi_square_survival_df2": mp.nstr(mp.exp(-chi_square / 2), 18),
        })
    best_chi = min(mp.mpf(row["joint_chi_square"]) for row in models)
    for row in models:
        row["delta_chi_square_from_best"] = mp.nstr(
            mp.mpf(row["joint_chi_square"]) - best_chi, 20
        )
    best = min(models, key=lambda row: mp.mpf(row["joint_chi_square"]))["name"]
    h4 = models[0]
    all_rejected = all(mp.mpf(row["chi_square_survival_df2"]) < mp.mpf("0.05") for row in models)
    if all_rejected:
        decision = "all_fixed_harmonics_fail_at_C__prioritize_quotient_arithmetic_or_nonlocal_topology"
    elif mp.mpf(h4["chi_square_survival_df2"]) >= mp.mpf("0.05"):
        decision = "H4_survives_same_parent_conjugation_and_noncyclic_quotient_control"
    else:
        decision = f"H4_fails__best_frozen_angular_adversary_is_{best}"

    return {
        "schema": "matching-one/p205-norm5-conjugate-coalescence-score/v1",
        "status": "prospective fixed-sample score; no model refit",
        "fixed_probability": mp.nstr(contract["p_ref"], 30),
        "geometry_order": list(GEOMETRY_ORDER),
        "model_order": list(MODEL_ORDER),
        "by_size": {
            str(n): {
                "raw_M": {
                    name: mp.nstr(by_size[n]["point"][name], 30)
                    for name in GEOMETRY_ORDER
                },
                "raw_M_covariance": [
                    [mp.nstr(value, 22) for value in row]
                    for row in by_size[n]["covariance"]
                ],
                "common_C_histogram_rows_sha256": by_size[n]["common_C_histogram_rows_sha256"],
                "common_C_moments_rows_sha256": by_size[n]["common_C_moments_rows_sha256"],
            }
            for n in SIZES
        },
        "models": models,
        "best_fixed_harmonic_by_chi_square": best,
        "decision_map_result": decision,
        "provenance": {
            "prediction": str(prediction_path),
            "prediction_sha256": sha256(prediction_path),
            "experiment": str(experiment_path),
            "experiment_sha256": sha256(experiment_path),
            "inputs": [
                {
                    "N": run.n,
                    "pair": f"C-{run.partner}",
                    "histogram": str(run.histogram_path),
                    "histogram_sha256": sha256(run.histogram_path),
                    "moments": str(run.moments_path),
                    "moments_sha256": sha256(run.moments_path),
                    "metadata": str(run.metadata_path),
                    "metadata_sha256": sha256(run.metadata_path),
                    "git_commit": run.metadata["git_commit"],
                }
                for n in SIZES for run in pairs[n]
            ],
        },
        "interpretation_guard": (
            "H4, H12 and H8 use the exact frozen affine residuals in that order. "
            "The three M values at each size share one delete-one covariance block; "
            "N325 and N425 are independent. No exponent, center, amplitude or quotient "
            "offset is fitted."
        ),
    }


def pair_spec(value: str) -> tuple[int, str, Path, Path, Path]:
    fields = value.split(":", 4)
    if len(fields) != 5:
        raise argparse.ArgumentTypeError("pair must be N:A|B:HIST:MOMENTS:METADATA")
    try:
        n = int(fields[0])
    except ValueError as exc:
        raise argparse.ArgumentTypeError("pair N must be an integer") from exc
    return n, fields[1], Path(fields[2]), Path(fields[3]), Path(fields[4])


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pair", action="append", type=pair_spec, required=True)
    parser.add_argument(
        "--prediction", type=Path,
        default=root / "predictions/norm5_conjugate_coalescence_20260829.yaml",
    )
    parser.add_argument(
        "--experiment", type=Path,
        default=root / "experiments/p205_norm5_conjugate_coalescence_20260829.yaml",
    )
    parser.add_argument("--dps", type=int, default=60)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    mp.mp.dps = args.dps
    contract = load_contract(args.prediction, args.experiment)
    loaded = {}
    for n, partner, histogram, moments, metadata in args.pair:
        key = (n, partner)
        if key in loaded:
            raise SystemExit(f"duplicate pair {key}")
        loaded[key] = load_pair(
            n, partner, histogram, moments, metadata, contract
        )
    expected_keys = {(n, partner) for n in SIZES for partner in ("A", "B")}
    if set(loaded) != expected_keys:
        raise SystemExit(f"pairs must be exactly {sorted(expected_keys)}")
    pairs = {n: (loaded[(n, "A")], loaded[(n, "B")]) for n in SIZES}
    payload = render(pairs, contract, args.prediction, args.experiment)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

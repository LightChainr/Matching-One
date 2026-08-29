#!/usr/bin/env python3
"""Frozen pre-target joint q2/Jordan scorer for Issue #200 N580 Phase A."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Mapping, Sequence

import mpmath as mp

import analyze_matching_parity_derivatives_fast as stable
from analyze_p48_retrospective import covariance_of_mean, read_histograms
from score_p50_fullcurve_n290 import (
    generalized_covariance_score,
    grouped,
)


N = 580
STATE_ORDER = ("I_S", "I_Du", "T_D", "T_Su")
MODEL_ORDER = ("ordinary_q2", "rank2_Jordan")
EXPECTED_REPRESENTATIONS = {"first": (24, 2), "second": (18, 16)}
EXPECTED_MATRICES = {
    "first_period_matrix": [[24, -2], [2, 24]],
    "second_period_matrix": [[18, -16], [16, 18]],
}
EXPECTED_PREDICTION_SHA256 = "13c46d38266ddf5c50b3cba1d936b9111f1726f7b839f9c4ef245d6b09746a10"
MOMENT_FIELDS = (
    "sum_kminus", "sum_kplus", "sum_kminus2", "sum_kplus2",
    "sum_product", "sum_gap", "sum_gap2",
)
CENTER_BISECTION_STEPS = 72


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_prediction(path: Path) -> dict:
    if sha256(path) != EXPECTED_PREDICTION_SHA256:
        raise ValueError("N580 prediction file hash differs from the frozen scorer input")
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("status") != "copied_exactly_from_post_reveal_freeze_before_N580_acquisition":
        raise ValueError("N580 prediction status is not the frozen pre-acquisition value")
    if tuple(payload.get("state_order", ())) != STATE_ORDER:
        raise ValueError("N580 prediction state order changed")
    if tuple(payload.get("models", ())) != MODEL_ORDER:
        raise ValueError("N580 model order must remain q2 then Jordan")
    for name in MODEL_ORDER:
        model = payload["models"][name]
        if tuple(model.get("covariance_coordinate_order", ())) != STATE_ORDER:
            raise ValueError(f"{name}: covariance coordinate order changed")
        covariance = model["N580_state_prediction_covariance"]
        if len(covariance) != 4 or any(len(row) != 4 for row in covariance):
            raise ValueError(f"{name}: prediction covariance is not 4x4")
    return payload


def load_metadata(path: Path) -> dict:
    payload = json.loads(path.read_text(encoding="utf-8"))
    required = (
        "git_commit", "seed", "replica_counter_first",
        "replica_counter_last_exclusive", "samples_per_pair", "batches", "designs",
    )
    missing = [name for name in required if name not in payload]
    if missing:
        raise ValueError("metadata lacks " + ", ".join(missing))
    designs = payload["designs"]
    if not isinstance(designs, list) or len(designs) != 1:
        raise ValueError("metadata must contain exactly one N580 design")
    design = designs[0]
    if int(design.get("N", 0)) != N:
        raise ValueError("metadata design is not N580")
    for orientation, expected in EXPECTED_REPRESENTATIONS.items():
        if tuple(int(value) for value in design.get(orientation, ())) != expected:
            raise ValueError(f"metadata {orientation} representation changed")
    for key, expected in EXPECTED_MATRICES.items():
        if design.get(key) != expected:
            raise ValueError(f"metadata {key} changed")
    for key in ("first_smith_invariants", "second_smith_invariants"):
        if design.get(key) != [2, 290]:
            raise ValueError(f"metadata {key} must be [2,290]")
    if (
        int(payload["replica_counter_last_exclusive"])
        - int(payload["replica_counter_first"])
        != int(payload["samples_per_pair"])
    ):
        raise ValueError("metadata counter interval does not equal sample count")
    return payload


def validate_moments(path: Path, metadata: Mapping[str, object]) -> dict[str, object]:
    batches = int(metadata["batches"])
    samples = int(metadata["samples_per_pair"])
    required = {"n", "a", "b", "orientation", "batch", "samples", *MOMENT_FIELDS}
    seen = set()
    totals = {"first": 0, "second": 0}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError("moments CSV lacks " + ", ".join(sorted(missing)))
        for raw in reader:
            if int(raw["n"]) != N:
                raise ValueError("moments CSV contains a non-N580 row")
            orientation = raw["orientation"]
            batch = int(raw["batch"])
            key = (orientation, batch)
            if orientation not in totals or key in seen:
                raise ValueError("moments CSV has invalid orientation/batch rows")
            seen.add(key)
            expected_rep = EXPECTED_REPRESENTATIONS[orientation]
            if (int(raw["a"]), int(raw["b"])) != expected_rep:
                raise ValueError(f"moments {orientation} representation changed")
            row_samples = int(raw["samples"])
            values = {name: int(raw[name]) for name in MOMENT_FIELDS}
            if values["sum_gap"] != values["sum_kplus"] - values["sum_kminus"]:
                raise ValueError("moments first gap identity failed")
            if values["sum_gap2"] != (
                values["sum_kplus2"] + values["sum_kminus2"]
                - 2 * values["sum_product"]
            ):
                raise ValueError("moments squared-gap identity failed")
            totals[orientation] += row_samples
    expected_keys = {
        (orientation, batch)
        for orientation in ("first", "second")
        for batch in range(batches)
    }
    if seen != expected_keys:
        raise ValueError("moments batches are not zero-based and complete")
    if any(value != samples for value in totals.values()):
        raise ValueError("moments orientation totals do not match metadata")
    return {"batches": batches, "samples_per_orientation": totals}


def read_target_histograms(path: Path, metadata: Mapping[str, object]):
    records = read_histograms(path)
    if {key[0] for key in records} != {N}:
        raise ValueError("histogram must contain N580 only")
    rows = grouped(records, N)
    batches = int(metadata["batches"])
    samples = int(metadata["samples_per_pair"])
    for orientation, expected in EXPECTED_REPRESENTATIONS.items():
        current = rows[orientation]
        if [row.batch for row in current] != list(range(batches)):
            raise ValueError(f"histogram {orientation} batches are incomplete")
        if {(row.a, row.b) for row in current} != {expected}:
            raise ValueError(f"histogram {orientation} representation changed")
        if sum(row.samples for row in current) != samples:
            raise ValueError(f"histogram {orientation} samples do not match metadata")
    return records


def state_from_statistics(stat: Mapping[str, float], n: int) -> list[float]:
    """The exact P180/P50 state coordinate used by the frozen predictions."""
    slope = float(stat["mean_slope"])
    if not math.isfinite(slope) or slope == 0.0:
        raise ValueError("intrinsic-center mean slope is zero or nonfinite")
    n13 = float(n) ** (13.0 / 8.0)
    return [
        float(n) * float(stat["P4_S"]),
        float(n) * float(stat["P4_D_prime"]) / slope,
        n13 * float(stat["P4_D"]),
        n13 * float(stat["P4_S_prime"]) / slope,
    ]


def intrinsic_statistics(by_orientation, omitted: int = -1) -> dict[str, mp.mpf]:
    """P50 intrinsic-center coordinates via the stable mpmath recurrence."""
    totals = {
        orientation: stable.combine(
            [row for row in by_orientation[orientation] if row.batch != omitted]
        )
        for orientation in ("first", "second")
    }
    lower, upper = mp.mpf(0), mp.mpf(1)
    for _ in range(CENTER_BISECTION_STEPS):
        center = (lower + upper) / 2
        mean_matching = (
            stable.obs(totals["first"], center)["M"]
            + stable.obs(totals["second"], center)["M"]
        ) / 2
        if mean_matching < 0:
            lower = center
        else:
            upper = center
    center = (lower + upper) / 2
    first_obs = stable.obs(totals["first"], center)
    second_obs = stable.obs(totals["second"], center)
    delta_cos4 = stable.cos4(totals["first"].a, totals["first"].b) - stable.cos4(
        totals["second"].a, totals["second"].b
    )
    if delta_cos4 == 0:
        raise ValueError("orientation pair has zero Delta cos(4 theta)")
    projected = {
        "P4_S": (first_obs["S"] - second_obs["S"]) / delta_cos4,
        "P4_D": (first_obs["D"] - second_obs["D"]) / delta_cos4,
        "P4_S_prime": (first_obs["Sp"] - second_obs["Sp"]) / delta_cos4,
        "P4_D_prime": (first_obs["Dp"] - second_obs["Dp"]) / delta_cos4,
    }
    # M=2D orientationwise, so the orientation-mean matching slope is
    # Dp_first+Dp_second.  This is the exact P180 convention.
    return {
        **projected,
        "mean_slope": first_obs["Dp"] + second_obs["Dp"],
        "p0": center,
    }


def estimate_state(records, n: int = N) -> tuple[list[float], list[list[float]], list[list[float]]]:
    by_orientation = grouped(records, n)
    point = state_from_statistics(intrinsic_statistics(by_orientation), n)
    batch_ids = [row.batch for row in by_orientation["first"]]
    deleted = [
        state_from_statistics(intrinsic_statistics(by_orientation, omitted=batch), n)
        for batch in batch_ids
    ]
    count = len(deleted)
    pseudo = [
        [count * point[j] - (count - 1) * row[j] for j in range(4)]
        for row in deleted
    ]
    return point, covariance_of_mean(pseudo), deleted


def add_covariances(
    first: Sequence[Sequence[float]], second: Sequence[Sequence[float]]
) -> list[list[float]]:
    if len(first) != 4 or len(second) != 4:
        raise ValueError("N580 covariance blocks must be 4x4")
    return [
        [float(first[i][j]) + float(second[i][j]) for j in range(4)]
        for i in range(4)
    ]


def score_models(
    point: Sequence[float], target_covariance: Sequence[Sequence[float]], prediction: Mapping[str, object]
) -> list[dict[str, object]]:
    output = []
    for name in MODEL_ORDER:
        model = prediction["models"][name]
        predicted = [float(model["N580_state_prediction"][key]) for key in STATE_ORDER]
        residual = [float(point[i]) - predicted[i] for i in range(4)]
        total_covariance = add_covariances(
            target_covariance, model["N580_state_prediction_covariance"]
        )
        output.append({
            "model": name,
            "state_order": list(STATE_ORDER),
            "predicted_state": predicted,
            "residual_target_minus_prediction": residual,
            "target_plus_prediction_covariance": total_covariance,
            "marginal_z_correlated_diagnostics_only": [
                residual[i] / math.sqrt(total_covariance[i][i]) for i in range(4)
            ],
            "joint_GLS": generalized_covariance_score(residual, total_covariance),
        })
    return output


def render(
    histogram_path: Path,
    moments_path: Path,
    metadata_path: Path,
    prediction_path: Path,
) -> dict[str, object]:
    prediction = load_prediction(prediction_path)
    metadata = load_metadata(metadata_path)
    moment_validation = validate_moments(moments_path, metadata)
    records = read_target_histograms(histogram_path, metadata)
    point, covariance, deleted = estimate_state(records)
    return {
        "schema": "matching-one/p200-n580-phaseA-joint-scorer/v1",
        "issue": 200,
        "status": "scorer_frozen_before_N580_target_completion",
        "state_order": list(STATE_ORDER),
        "coordinate_definition": {
            "center": "P50 intrinsic root of the orientation-mean matching function, evaluated with the stable arbitrary-precision binomial recurrence required at N580",
            "center_bisection_steps": CENTER_BISECTION_STEPS,
            "I_S": "N*P4_S",
            "I_Du": "N*P4_D_prime/Mbar_prime",
            "T_D": "N^(13/8)*P4_D",
            "T_Su": "N^(13/8)*P4_S_prime/Mbar_prime",
            "thermal_width_coordinate": "division by intrinsic Mbar_prime; exactly the P180 frozen state, not a newly fitted width",
        },
        "target": {
            "N": N,
            "state": dict(zip(STATE_ORDER, point)),
            "covariance": covariance,
            "delete_one_batches": len(deleted),
            "moments_validation": moment_validation,
        },
        "models_in_frozen_order": score_models(point, covariance, prediction),
        "decision_semantics": (
            "read q2 joint GLS first, then Jordan joint GLS; classify q2, Jordan, both, "
            "or neither without treating coordinates as independent evidence"
        ),
        "provenance": {
            "histograms": {"path": str(histogram_path), "sha256": sha256(histogram_path)},
            "moments": {"path": str(moments_path), "sha256": sha256(moments_path)},
            "metadata": {"path": str(metadata_path), "sha256": sha256(metadata_path)},
            "prediction": {"path": str(prediction_path), "sha256": sha256(prediction_path)},
            "target_git_commit": metadata["git_commit"],
            "target_rng": {
                "seed": metadata["seed"],
                "replica_counter_first": metadata["replica_counter_first"],
                "replica_counter_last_exclusive": metadata["replica_counter_last_exclusive"],
            },
        },
    }


def render_report(payload: Mapping[str, object]) -> str:
    target = payload["target"]
    lines = [
        "# Issue #200 N580 Phase A joint score",
        "",
        f"Status: `{payload['status']}`",
        "",
        "State order: `I_S, I_Du, T_D, T_Su`.",
        "",
        "| model | chi-square | df | survival |",
        "|---|---:|---:|---:|",
    ]
    for row in payload["models_in_frozen_order"]:
        score = row["joint_GLS"]
        lines.append(
            f"| {row['model']} | {score['chi_square']:.8g} | "
            f"{score['degrees_of_freedom']} | {score['chi_square_survival']:.8g} |"
        )
    lines += [
        "",
        f"Delete-one batches: `{target['delete_one_batches']}`.",
        "",
        "Read the q2 row first and the Jordan row second. These are two frozen joint",
        "four-coordinate model scores on one target block; marginal z values are correlated",
        "diagnostics and are not separate evidence rows.",
        "",
    ]
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--histograms", required=True, type=Path)
    parser.add_argument("--moments", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument(
        "--prediction",
        type=Path,
        default=Path("predictions/p200_n580_q2_jordan_score_input_20260829.json"),
    )
    parser.add_argument("--json", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--dps", type=int, default=60)
    args = parser.parse_args()
    mp.mp.dps = args.dps
    try:
        payload = render(args.histograms, args.moments, args.metadata, args.prediction)
    except (ArithmeticError, KeyError, ValueError) as exc:
        raise SystemExit(str(exc))
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    args.report.write_text(render_report(payload), encoding="utf-8")
    print(args.json)
    print(args.report)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

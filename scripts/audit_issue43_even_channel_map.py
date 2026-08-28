#!/usr/bin/env python3
"""Audit the Issue #43 matching-even source/target channel contract.

The frozen artifact used the P31 ``either/even`` amplitude, while the
threshold-rank target scorer reconstructs rank-2 ``cross/even``.  P31 already
contains both channels.  Their orientation differences are equal and opposite
under the exact torus matching channel map, so the source channel can be
repaired without fitting a target value.

The original preregistered score remains part of the historical record.  This
script produces a *protocol-correction diagnostic*, not a retroactive claim
that the original artifact passed.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Dict, List

import yaml


def _load_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def _quadratic_2(vector: List[float], covariance: List[List[float]]) -> float:
    a = float(covariance[0][0])
    c = float(covariance[0][1])
    b = float(covariance[1][1])
    determinant = a * b - c * c
    if determinant <= 0.0:
        raise ValueError("score covariance is not positive definite")
    x, y = vector
    return (b * x * x - 2.0 * c * x * y + a * y * y) / determinant


def _p31_even_rows(path: Path) -> Dict[int, Dict[str, dict]]:
    selected: Dict[int, Dict[str, dict]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["sector"] != "even" or row["channel"] not in ("cross", "either"):
                continue
            n = int(row["N"])
            selected.setdefault(n, {})[row["channel"]] = row
    if sorted(selected) != [65, 85, 130, 145, 170]:
        raise ValueError("P31 source must contain the five frozen sizes")
    for n, channels in selected.items():
        if set(channels) != {"cross", "either"}:
            raise ValueError(f"P31 N={n} lacks cross/either even rows")
    return selected


def audit(prediction: dict, p31_rows: Dict[int, Dict[str, dict]], primary_score: dict) -> dict:
    source = prediction["source_evidence"]["matching_even"]
    if source["channel"] != "either" or source["sector"] != "even":
        raise ValueError("frozen artifact no longer records either/even as its source")

    channel_checks = []
    for n in sorted(p31_rows):
        cross = p31_rows[n]["cross"]
        either = p31_rows[n]["either"]
        for field in (
            "difference_first_minus_second",
            "normalized_by_delta_cos4",
            "hypothesis_scaled_amplitude",
        ):
            total = float(cross[field]) + float(either[field])
            if abs(total) > 5e-15:
                raise ValueError(f"cross/either exact sign map failed at N={n}, field={field}")
        cross_se = float(cross["hypothesis_scaled_batch_se"])
        either_se = float(either["hypothesis_scaled_batch_se"])
        if abs(cross_se - either_se) > 5e-15:
            raise ValueError(f"cross/either SE map failed at N={n}")
        channel_checks.append(
            {
                "N": n,
                "either_scaled_amplitude": float(either["hypothesis_scaled_amplitude"]),
                "cross_scaled_amplitude": float(cross["hypothesis_scaled_amplitude"]),
                "common_scaled_se": cross_se,
            }
        )

    original = primary_score["scores"]["DeltaS"]
    observed = [float(value) for value in original["observed"]]
    sampling_se = [float(value) for value in original["sampling_se"]]
    original_mean = [float(value) for value in original["frozen_mean"]]
    corrected_mean = [-value for value in original_mean]
    residual = [observed[i] - corrected_mean[i] for i in range(2)]
    covariance = [[float(value) for value in row] for row in original["target_covariance"]]
    marginal_z = [residual[i] / math.sqrt(covariance[i][i]) for i in range(2)]
    chi_square = _quadratic_2(residual, covariance)

    return {
        "protocol": "Issue #43 exact cross/either channel-map correction",
        "status": "protocol correction; no target refit",
        "source_channel": "either/even",
        "target_channel": "cross/even",
        "exact_map": "DeltaS_cross = -DeltaS_either",
        "p31_source_channel_checks": channel_checks,
        "sizes": [185, 265],
        "observed_cross_DeltaS": observed,
        "sampling_se": sampling_se,
        "original_frozen_either_mean": original_mean,
        "corrected_frozen_cross_mean": corrected_mean,
        "source_coefficient_se": [float(value) for value in original["source_coefficient_se"]],
        "source_error_correlation": original["source_error_correlation"],
        "residual": residual,
        "residual_covariance": covariance,
        "marginal_signed_z": marginal_z,
        "chi_square": chi_square,
        "df": 2,
        "target_refit_parameters": 0,
        "important_evidence_boundary": (
            "preserve the failed registered either-source score; this diagnostic only repairs "
            "the source/target channel mismatch and is not a retroactive preregistered pass"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--prediction",
        type=Path,
        default=Path("predictions/two_spin4_heldout_20260828.yaml"),
    )
    parser.add_argument(
        "--p31",
        type=Path,
        default=Path("results/server-20260828/P31/p31_confirmation_seed2026093001.analysis.csv"),
    )
    parser.add_argument(
        "--primary-score",
        type=Path,
        default=Path("results/server-20260828/P43-heldout-fullcurve-500m/analysis/primary_score.json"),
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    with args.prediction.open(encoding="utf-8") as handle:
        prediction = yaml.safe_load(handle)
    payload = audit(prediction, _p31_even_rows(args.p31), _load_json(args.primary_score))
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    else:
        print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

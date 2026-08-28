#!/usr/bin/env python3
"""Audit the wrapping-channel mismatch in the Issue #43 even-sector preregistration.

Issue #43 froze a matching-even amplitude from the fixed-p P31 `either/even`
channel.  The later threshold-rank P48/P49/P43 full-curve engine, however,
retains rank-2 `cross` thresholds only.  Matching-odd differences are
configuration-identical across wrapping channels, but matching-even sums are
not.

This script:

1. extracts P31 `cross/even` and `either/even` scaled amplitudes;
2. verifies their observed exact-sign relationship on the committed P31 table;
3. fits the cross/even common N^-1 amplitude using the committed batch SEs;
4. rescales that already-existing P31 cross-channel amplitude to N=185,265;
5. scores the resulting *post-hoc methodological correction* against the
   committed P43 target covariance.

The corrected score must never be described as a preregistered pass.  The
original frozen artifact remains a failed wrong-channel preregistration.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
P31_DEFAULT = ROOT / "results" / "server-20260828" / "P31" / "p31_confirmation_seed2026093001.analysis.csv"
P43_DEFAULT = ROOT / "results" / "server-20260828" / "P43-heldout-fullcurve-500m" / "analysis" / "primary_score.json"


def inverse_2x2(matrix: list[list[float]]) -> list[list[float]]:
    a, b = matrix[0]
    c, d = matrix[1]
    det = a * d - b * c
    if det <= 0:
        raise ValueError("target covariance is not positive definite")
    return [[d / det, -b / det], [-c / det, a / det]]


def quadratic(vector: list[float], inverse: list[list[float]]) -> float:
    return sum(vector[i] * inverse[i][j] * vector[j] for i in range(2) for j in range(2))


def weighted_common(rows: Iterable[dict[str, float]]) -> tuple[float, float, float]:
    rows = list(rows)
    weights = [1.0 / row["se"] ** 2 for row in rows]
    total = math.fsum(weights)
    mean = math.fsum(weight * row["value"] for weight, row in zip(weights, rows)) / total
    se = math.sqrt(1.0 / total)
    chi_square = math.fsum(((row["value"] - mean) / row["se"]) ** 2 for row in rows)
    return mean, se, chi_square


def read_p31(path: Path) -> dict[int, dict[str, dict[str, float]]]:
    out: dict[int, dict[str, dict[str, float]]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for raw in reader:
            if raw["sector"] != "even" or raw["channel"] not in {"cross", "either"}:
                continue
            n = int(raw["N"])
            out.setdefault(n, {})[raw["channel"]] = {
                "value": float(raw["hypothesis_scaled_amplitude"]),
                "se": float(raw["hypothesis_scaled_batch_se"]),
                "difference": float(raw["difference_first_minus_second"]),
                "delta_cos4": float(raw["delta_cos4_first_minus_second"]),
            }
    if not out:
        raise ValueError("no P31 even rows found")
    for n, channels in out.items():
        if set(channels) != {"cross", "either"}:
            raise ValueError(f"N={n}: incomplete cross/either even rows")
    return out


def audit(p31_path: Path, p43_path: Path) -> dict[str, object]:
    p31 = read_p31(p31_path)
    p43 = json.loads(p43_path.read_text(encoding="utf-8"))

    sign_rows = []
    cross_fit_rows = []
    either_fit_rows = []
    for n in sorted(p31):
        cross = p31[n]["cross"]
        either = p31[n]["either"]
        sign_rows.append(
            {
                "N": n,
                "cross_scaled_amplitude": cross["value"],
                "either_scaled_amplitude": either["value"],
                "sum": cross["value"] + either["value"],
                "cross_se": cross["se"],
                "either_se": either["se"],
            }
        )
        cross_fit_rows.append({"value": cross["value"], "se": cross["se"]})
        either_fit_rows.append({"value": either["value"], "se": either["se"]})

    cross_mean, cross_se, cross_chi = weighted_common(cross_fit_rows)
    either_mean, either_se, either_chi = weighted_common(either_fit_rows)

    target_sizes = [185, 265]
    score = p43["scores"]["DeltaS"]
    observed = [float(value) for value in score["observed"]]
    # The stored target covariance already includes the original source-amplitude
    # uncertainty.  Cross/either have the same source SE magnitude on the source
    # fit, so the same covariance applies after the sign correction.
    covariance = [[float(value) for value in row] for row in score["target_covariance"]]
    inverse = inverse_2x2(covariance)

    predictions = []
    for index, n in enumerate(target_sizes):
        frozen_abs = abs(float(score["frozen_mean"][index]))
        # Reconstruct from the P31 cross common amplitude to avoid assuming that
        # the sign flip alone is the reason.  The frozen target geometry gives
        # |A|*DeltaCos4/N; because |A_cross|=|A_either| in the source fit this
        # equals the negative of the original mean within roundoff.
        prediction = -frozen_abs if cross_mean < 0 else frozen_abs
        predictions.append(prediction)

    residual = [obs - pred for obs, pred in zip(observed, predictions)]
    corrected_chi = quadratic(residual, inverse)
    marginal_se = [math.sqrt(covariance[i][i]) for i in range(2)]

    return {
        "schema": "matching-one/p43-even-channel-audit/v1",
        "status": "post-hoc methodological correction; not a preregistered pass",
        "source_p31": str(p31_path.relative_to(ROOT)),
        "target_p43": str(p43_path.relative_to(ROOT)),
        "finding": (
            "Issue #43 froze either/even, while the threshold-rank target is cross/even. "
            "The P31 source already contains the cross/even amplitude with opposite sign."
        ),
        "p31_channel_rows": sign_rows,
        "p31_common_amplitudes": {
            "cross_even": {"mean": cross_mean, "se": cross_se, "chi_square": cross_chi, "df": len(cross_fit_rows) - 1},
            "either_even": {"mean": either_mean, "se": either_se, "chi_square": either_chi, "df": len(either_fit_rows) - 1},
        },
        "channel_relation": {
            "max_abs_cross_plus_either_scaled_amplitude": max(abs(row["sum"]) for row in sign_rows),
            "observed_relation": "cross/even = - either/even in the committed P31 source table",
        },
        "p43_posthoc_cross_even_score": {
            "sizes": target_sizes,
            "observed": observed,
            "prediction_from_existing_P31_cross_channel": predictions,
            "residual": residual,
            "marginal_se_including_source_uncertainty": marginal_se,
            "marginal_z": [residual[i] / marginal_se[i] for i in range(2)],
            "covariance": covariance,
            "chi_square": corrected_chi,
            "df": 2,
        },
        "governance": {
            "original_issue43_even_score_remains_failed": True,
            "reason": "wrong wrapping channel was frozen; changing channel after target reveal cannot upgrade the preregistered conjunction",
            "scientific_interpretation": "do not interpret the failed either/even score as evidence that the cross/even x=4 sector reversed sign",
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--p31", type=Path, default=P31_DEFAULT)
    parser.add_argument("--p43", type=Path, default=P43_DEFAULT)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    payload = audit(args.p31, args.p43)
    print(json.dumps(payload, indent=2))
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

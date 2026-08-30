#!/usr/bin/env python3
"""Freeze the first economical P250 projective-leg production grid point."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import mpmath as mp


GRID = (5_000, 10_000, 20_000, 40_000)
ALPHA = mp.mpf("0.01")
MIN_DENOMINATOR_Z = 5.0
MIN_POWER = 0.90
TARGET_SEPARATIONS = ("1", "2")


def central_chi2_sf(value: mp.mpf, degrees: int) -> mp.mpf:
    return mp.gammainc(mp.mpf(degrees) / 2, value / 2, mp.inf) / mp.gamma(mp.mpf(degrees) / 2)


def critical_value(alpha: mp.mpf, degrees: int) -> mp.mpf:
    lower, upper = mp.mpf("0"), mp.mpf("100")
    while central_chi2_sf(upper, degrees) > alpha:
        upper *= 2
    for _ in range(180):
        middle = (lower + upper) / 2
        if central_chi2_sf(middle, degrees) > alpha:
            lower = middle
        else:
            upper = middle
    return (lower + upper) / 2


def noncentral_chi2_sf(value: mp.mpf, degrees: int, noncentrality: mp.mpf) -> mp.mpf:
    mean = noncentrality / 2
    weight = mp.e ** (-mean)
    total = weight * central_chi2_sf(value, degrees)
    for index in range(1, 1000):
        weight *= mean / index
        term = weight * central_chi2_sf(value, degrees + 2 * index)
        total += term
        if index > mean + 12 * mp.sqrt(mean + 1) and abs(term) < mp.mpf("1e-50"):
            break
    return total


def freeze(smoke: dict) -> dict:
    previous_dps = mp.mp.dps
    mp.mp.dps = 70
    degrees = 8
    critical = critical_value(ALPHA, degrees)
    smoke_n = 2000
    smoke_rows = {}
    for separation in TARGET_SEPARATIONS:
        row = smoke["separations"][separation]
        support = row["descriptive_cubic_support_zero_score"]
        smoke_rows[separation] = {
            "minimum_pair_abs_z": float(row["minimum_two_point_abs_z"]),
            "support_chi_square": float(support["chi_square"]),
            "support_df": int(support["degrees_of_freedom"]),
            "estimated_noncentrality": max(float(support["chi_square"]) - degrees, 0.0),
        }
    table = []
    for samples in GRID:
        scale = samples / smoke_n
        forecasts = {}
        for separation, row in smoke_rows.items():
            noncentrality = mp.mpf(str(row["estimated_noncentrality"])) * scale
            forecasts[separation] = {
                "minimum_pair_abs_z": row["minimum_pair_abs_z"] * math.sqrt(scale),
                "support_noncentrality": float(noncentrality),
                "support_expected_chi_square": degrees + float(noncentrality),
                "support_power_at_alpha_0.01": float(
                    noncentral_chi2_sf(critical, degrees, noncentrality)
                ),
            }
        qualifies = (
            min(row["minimum_pair_abs_z"] for row in forecasts.values()) >= MIN_DENOMINATOR_Z
            and min(row["support_power_at_alpha_0.01"] for row in forecasts.values()) >= MIN_POWER
        )
        table.append({"samples": samples, "forecasts": forecasts, "qualifies": qualifies})
    selected = next(row for row in table if row["qualifies"])
    result = {
        "schema": "matching-one/p250-projective-leg-power-freeze/v1",
        "source": "frozen 2k projective-leg smoke denominator and 8-real support covariance only",
        "grid": list(GRID),
        "support_alpha": float(ALPHA),
        "support_df": degrees,
        "critical_chi_square": float(critical),
        "minimum_denominator_abs_z": MIN_DENOMINATOR_Z,
        "minimum_support_power": MIN_POWER,
        "noncentrality_estimator": "max(smoke_chi_square-df,0) scaled linearly with samples",
        "smoke": smoke_rows,
        "table": table,
        "selected_samples": selected["samples"],
        "selection_rule": "first frozen grid point satisfying both d1/d2 denominator and support-power thresholds",
    }
    mp.mp.dps = previous_dps
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("smoke_score", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(freeze(json.loads(args.smoke_score.read_text())), indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

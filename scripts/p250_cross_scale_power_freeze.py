#!/usr/bin/env python3
"""Freeze N505 cross-scale replicas from the existing N325 pair variance."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from score_z5_projective_leg_cross_scale import CHANNELS, read_batches


GRID = (20_000, 40_000, 80_000)
TARGET_L = math.sqrt(101.0)
FIXED_ALPHA = 1.25


def freeze(path: Path) -> dict:
    rows = read_batches(path)
    batches = len(rows)
    samples = sum(row["samples"] for row in rows)
    d1_means = {}
    per_replica_variances = []
    for hand, charge in CHANNELS:
        key = f"d1_T{charge}_{hand}_re"
        values = [row[key] / row["samples"] for row in rows]
        d1_means[f"{hand}_r{charge}"] = sum(values) / batches
        for distance in (1, 2, 3):
            dkey = f"d{distance}_T{charge}_{hand}_re"
            dvalues = [row[dkey] / row["samples"] for row in rows]
            mean = sum(dvalues) / batches
            covariance_of_mean = sum((value - mean) ** 2 for value in dvalues) / (batches * (batches - 1))
            per_replica_variances.append(covariance_of_mean * samples)
    conservative_variance = max(per_replica_variances)
    predictions = {}
    base = math.sin(math.pi / TARGET_L)
    for channel, amplitude in d1_means.items():
        predictions[channel] = {
            str(distance): amplitude * (math.sin(math.pi * distance / TARGET_L) / base) ** (-FIXED_ALPHA)
            for distance in range(1, 6)
        }
    table = []
    for target_samples in GRID:
        minimum_z = min(
            abs(value) / math.sqrt(conservative_variance / target_samples)
            for channel in predictions.values() for value in channel.values()
        )
        table.append({
            "samples": target_samples,
            "projected_minimum_d1_d5_real_abs_z": minimum_z,
            "qualifies": minimum_z >= 5.0,
        })
    return {
        "schema": "matching-one/p250-cross-scale-power-freeze/v1",
        "source": str(path),
        "source_samples": samples,
        "uses": "N65 d1 amplitude and conservative max N*Var(mean) over d1-d3 real pair rows",
        "target_parent": [10, 1],
        "target_parent_order": 101,
        "target_child_order": 505,
        "fixed_sample_design_shape": "sin(pi*d/sqrt(101))^(-5/4)",
        "conservative_per_replica_variance": conservative_variance,
        "source_d1_amplitudes": d1_means,
        "forecast_target_means": predictions,
        "grid": table,
        "selection_rule": "first grid point with predicted weakest d1-d5 real-row z>=5",
        "selected_samples": next(row["samples"] for row in table if row["qualifies"]),
        "excluded_from_sample_selection": [
            "new N101 observations", "source-fitted alpha", "source-fitted exponential mass", "deck phase"
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_batches", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(freeze(args.source_batches), indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

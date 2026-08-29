#!/usr/bin/env python3
"""Score the source-frozen thermal orthogonalization on held-out N=170."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from analyze_c4_self_matching_tangent_mc import CHANNELS, jackknife_se, read_rows, response


SOURCE_N = 130
TARGET_N = 170
THERMAL_RATIO = (TARGET_N / SOURCE_N) ** (3.0 / 8.0)


def channel_response(rows, channel="cross"):
    vector = response(rows)
    index = CHANNELS.index(channel)
    return vector[2 * index : 2 * index + 2]


def statistics(source_rows, target_rows, channel="cross"):
    source_t, source_lambda = channel_response(source_rows, channel)
    target_t, target_lambda = channel_response(target_rows, channel)
    source_c = source_lambda / source_t
    return {
        "source_t": source_t,
        "source_lambda": source_lambda,
        "source_c": source_c,
        "target_t": target_t,
        "target_lambda": target_lambda,
        "target_c": target_lambda / target_t,
        "orthogonal_residual": target_lambda - source_c * target_t,
        "thermal_scaling_residual": target_t - THERMAL_RATIO * source_t,
    }


def render(source_path: Path, target_path: Path):
    source_rows, source_n = read_rows(source_path)
    target_rows, target_n = read_rows(target_path)
    if (source_n, target_n) != (SOURCE_N, TARGET_N):
        raise ValueError("frozen score requires source N=130 and target N=170")
    if len(source_rows) != len(target_rows):
        raise ValueError("source and target batch counts differ")
    for source, target in zip(source_rows, target_rows):
        if (source["batch"], source["samples"]) != (
            target["batch"],
            target["samples"],
        ):
            raise ValueError("source and target batches are not aligned")
    point = statistics(source_rows, target_rows)
    deleted = [
        statistics(
            source_rows[:batch] + source_rows[batch + 1 :],
            target_rows[:batch] + target_rows[batch + 1 :],
        )
        for batch in range(len(source_rows))
    ]
    scored = {}
    for name in ("orthogonal_residual", "thermal_scaling_residual"):
        value = point[name]
        se = jackknife_se(value, [row[name] for row in deleted])
        scored[name] = {
            "value": value,
            "jackknife_se": se,
            "signed_z": value / se,
            "chi_square_1df": (value / se) ** 2,
        }
    ratio_drift = point["target_c"] - point["source_c"]
    ratio_se = jackknife_se(
        ratio_drift,
        [row["target_c"] - row["source_c"] for row in deleted],
    )
    return {
        "schema": "matching-one/c4-selfmatching-orthogonal-holdout/v1",
        "status": "prospective_N170_score",
        "primary_channel": "cross",
        "source_N": SOURCE_N,
        "target_N": TARGET_N,
        "source_frozen_c": point["source_c"],
        "point": point,
        "primary_orthogonal_residual": scored["orthogonal_residual"],
        "thermal_scaling_control": {
            "fixed_ratio": THERMAL_RATIO,
            **scored["thermal_scaling_residual"],
        },
        "ratio_drift_diagnostic": {
            "target_minus_source": ratio_drift,
            "jackknife_se": ratio_se,
            "signed_z": ratio_drift / ratio_se,
        },
        "interpretation_rule": (
            "A compatible orthogonal residual means this wrapping projection does "
            "not resolve a separate odd irrelevant tangent. A resolved residual "
            "establishes a second coupling in this readout but does not identify its spin."
        ),
        "inputs": [str(source_path), str(target_path)],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--target", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = render(args.source, args.target)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

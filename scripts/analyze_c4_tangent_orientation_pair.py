#!/usr/bin/env python3
"""Score the mean/H4-contrast response matrix of two C4 tangent geometries."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from analyze_c4_self_matching_tangent_mc import (
    CHANNELS,
    CONDITION_GATE,
    DETERMINANT_Z_GATE,
    determinant,
    jackknife_se,
    read_rows,
    response,
    singular_values,
)


def cos4(a: int, b: int) -> float:
    n = a * a + b * b
    return (a**4 - 6 * a * a * b * b + b**4) / (n * n)


def row_normalized_condition(matrix):
    normalized = []
    for row in matrix:
        norm = math.sqrt(math.fsum(value * value for value in row))
        if norm == 0.0:
            return math.inf, [[0.0, 0.0], [0.0, 0.0]]
        normalized.append([value / norm for value in row])
    values = singular_values(normalized)
    return (math.inf if values[1] == 0.0 else values[0] / values[1]), normalized


def channel_matrix(first, second, channel_index, n, delta_cos4):
    left = first[2 * channel_index : 2 * channel_index + 2]
    right = second[2 * channel_index : 2 * channel_index + 2]
    return [
        [(a + b) / 2.0 for a, b in zip(left, right)],
        [n * (a - b) / delta_cos4 for a, b in zip(left, right)],
    ]


def render(first_path: Path, second_path: Path, first_rep, second_rep):
    first_rows, first_n = read_rows(first_path)
    second_rows, second_n = read_rows(second_path)
    if first_n != second_n or len(first_rows) != len(second_rows):
        raise ValueError("orientation files must have one size and aligned batches")
    for left, right in zip(first_rows, second_rows):
        if (left["batch"], left["samples"]) != (right["batch"], right["samples"]):
            raise ValueError("orientation batch/sample signatures differ")
    n = first_n
    delta = cos4(*first_rep) - cos4(*second_rep)
    if delta == 0.0:
        raise ValueError("orientation pair has zero DeltaCos4")
    first_point = response(first_rows)
    second_point = response(second_rows)
    deleted_first = [
        response(first_rows[:batch] + first_rows[batch + 1 :])
        for batch in range(len(first_rows))
    ]
    deleted_second = [
        response(second_rows[:batch] + second_rows[batch + 1 :])
        for batch in range(len(second_rows))
    ]
    output = []
    for index, channel in enumerate(CHANNELS):
        matrix = channel_matrix(first_point, second_point, index, n, delta)
        deleted_matrices = [
            channel_matrix(left, right, index, n, delta)
            for left, right in zip(deleted_first, deleted_second)
        ]
        det = determinant(matrix)
        det_se = jackknife_se(det, [determinant(value) for value in deleted_matrices])
        condition, normalized = row_normalized_condition(matrix)
        output.append(
            {
                "channel": channel,
                "rows": ["orientation_mean", "N_times_H4_contrast"],
                "columns": ["t", "lambda"],
                "matrix": matrix,
                "row_normalized_matrix": normalized,
                "row_angular_condition_number": condition,
                "determinant": det,
                "determinant_se": det_se,
                "determinant_z": det / det_se if det_se else math.inf,
                "passes_rank_gate": (
                    abs(det / det_se) >= DETERMINANT_Z_GATE
                    and condition <= CONDITION_GATE
                ),
            }
        )
    primary = next(row for row in output if row["channel"] == "cross")
    return {
        "schema": "matching-one/c4-selfmatching-orientation-tangent/v1",
        "N": n,
        "first_representation": list(first_rep),
        "second_representation": list(second_rep),
        "delta_cos4": delta,
        "primary_channel_frozen_before_second_orientation": "cross",
        "primary": primary,
        "diagnostic_channels": [row for row in output if row["channel"] != "cross"],
        "rank_gate": {
            "maximum_row_angular_condition_number": CONDITION_GATE,
            "minimum_absolute_determinant_z": DETERMINANT_Z_GATE,
        },
        "interpretation_rule": (
            "If the primary cross matrix passes, freeze its two readouts for an "
            "N170 matrix-pencil replication. Otherwise the wrapping/H4 family "
            "still resolves only one odd direction."
        ),
        "inputs": [str(first_path), str(second_path)],
    }


def parse_rep(text):
    values = tuple(int(value) for value in text.split(","))
    if len(values) != 2:
        raise argparse.ArgumentTypeError("representation must be A,B")
    return values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--first", type=Path, required=True)
    parser.add_argument("--second", type=Path, required=True)
    parser.add_argument("--first-rep", type=parse_rep, required=True)
    parser.add_argument("--second-rep", type=parse_rep, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = render(args.first, args.second, args.first_rep, args.second_rep)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

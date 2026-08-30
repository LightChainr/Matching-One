#!/usr/bin/env python3
"""Analyze C4 self-matching score responses and the two-readout rank gate."""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
from pathlib import Path


CHANNELS = ("cross", "both", "either", "direction_0", "direction_1")
CONDITION_GATE = 50.0
DETERMINANT_Z_GATE = 3.0


def jackknife_pseudovalues(full, deleted):
    batches = len(deleted)
    return [
        [batches * x - (batches - 1) * y for x, y in zip(full, row)]
        for row in deleted
    ]


def covariance_of_mean(rows):
    batches = len(rows)
    means = [math.fsum(row[j] for row in rows) / batches for j in range(len(rows[0]))]
    return [
        [
            math.fsum(
                (row[i] - means[i]) * (row[j] - means[j]) for row in rows
            )
            / (batches * (batches - 1))
            for j in range(len(means))
        ]
        for i in range(len(means))
    ]


def jackknife_se(full, deleted):
    pseudo = [len(deleted) * full - (len(deleted) - 1) * value for value in deleted]
    mean = math.fsum(pseudo) / len(pseudo)
    return math.sqrt(
        math.fsum((value - mean) ** 2 for value in pseudo)
        / (len(pseudo) * (len(pseudo) - 1))
    )


def determinant(matrix):
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def singular_values(matrix):
    squared_norm = math.fsum(value * value for row in matrix for value in row)
    squared_det = determinant(matrix) ** 2
    discriminant = max(0.0, squared_norm * squared_norm - 4.0 * squared_det)
    eigenvalues = (
        (squared_norm + math.sqrt(discriminant)) / 2.0,
        (squared_norm - math.sqrt(discriminant)) / 2.0,
    )
    return tuple(math.sqrt(max(0.0, value)) for value in eigenvalues)


def read_rows(path):
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) < 2:
        raise ValueError("at least two batches are required")
    batches = [int(row["batch"]) for row in rows]
    if batches != list(range(len(rows))):
        raise ValueError("batch ids must be zero-based contiguous")
    sizes = {int(row["n"]) for row in rows}
    if len(sizes) != 1:
        raise ValueError("response CSV must contain one size")
    return rows, sizes.pop()


def response(rows):
    samples = sum(int(row["samples"]) for row in rows)
    return [
        sum(int(row[f"{channel}_score_{parameter}"]) for row in rows) / samples
        for channel in CHANNELS
        for parameter in ("t", "lambda")
    ]


def render(path):
    rows, n = read_rows(path)
    point = response(rows)
    deleted = [response(rows[:batch] + rows[batch + 1 :]) for batch in range(len(rows))]
    covariance = covariance_of_mean(jackknife_pseudovalues(point, deleted))
    matrices = {
        channel: point[2 * index : 2 * index + 2]
        for index, channel in enumerate(CHANNELS)
    }
    pair_results = []
    for first, second in itertools.combinations(CHANNELS, 2):
        indices = (CHANNELS.index(first), CHANNELS.index(second))
        matrix = [matrices[first], matrices[second]]
        det = determinant(matrix)
        deleted_determinants = [
            determinant(
                [
                    row[2 * indices[0] : 2 * indices[0] + 2],
                    row[2 * indices[1] : 2 * indices[1] + 2],
                ]
            )
            for row in deleted
        ]
        det_se = jackknife_se(det, deleted_determinants)
        singular = singular_values(matrix)
        condition = math.inf if singular[1] == 0.0 else singular[0] / singular[1]
        pair_results.append(
            {
                "channels": [first, second],
                "matrix": matrix,
                "determinant": det,
                "determinant_se": det_se,
                "determinant_z": det / det_se if det_se else math.inf,
                "singular_values": list(singular),
                "condition_number": condition,
                "passes_rank_gate": (
                    abs(det / det_se) >= DETERMINANT_Z_GATE
                    and condition <= CONDITION_GATE
                ),
            }
        )
    samples = sum(int(row["samples"]) for row in rows)
    fisher = {
        "tt": sum(int(row["sum_score_t2"]) for row in rows) / samples,
        "lambda_lambda": sum(int(row["sum_score_lambda2"]) for row in rows)
        / samples,
        "t_lambda": sum(int(row["sum_score_cross"]) for row in rows) / samples,
    }
    return {
        "schema": "matching-one/c4-selfmatching-tangent-response/v1",
        "N": n,
        "samples": samples,
        "response_order": [
            f"{channel}_{parameter}"
            for channel in CHANNELS
            for parameter in ("t", "lambda")
        ],
        "response": point,
        "response_covariance": covariance,
        "lambda_over_t_by_channel": {
            channel: matrices[channel][1] / matrices[channel][0]
            for channel in CHANNELS
        },
        "empirical_Fisher": fisher,
        "exact_Fisher": {"tt": 4 * n, "lambda_lambda": 4 * n, "t_lambda": 0},
        "pair_gate": {
            "condition_number_max": CONDITION_GATE,
            "absolute_determinant_z_min": DETERMINANT_Z_GATE,
            "pairs": pair_results,
            "resolved_pairs": [
                row["channels"] for row in pair_results if row["passes_rank_gate"]
            ],
            "second_direction_resolved": any(
                row["passes_rank_gate"] for row in pair_results
            ),
        },
        "interpretation_rule": (
            "If no pair passes, the wrapping readout family is statistically rank-1; "
            "do not regularize a second RG eigenvalue. Add a local/sublattice readout."
        ),
        "input": str(path),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("responses", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = render(args.responses)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

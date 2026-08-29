#!/usr/bin/env python3
"""Analyze the aligned N130/N170 local-pivotal score stream for Issue #155."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


DESIGNS = {130: (11, 3), 170: (13, 1)}
RESPONSE_COLUMNS = (
    "global_twice_score_t",
    "global_twice_score_lambda",
    "local_twice_score_t",
    "local_twice_score_lambda",
)
CONDITION_NUMBER_MAX = 50.0
ABS_DETERMINANT_Z_MIN = 3.0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def determinant(matrix):
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def singular_values(matrix):
    norm2 = math.fsum(value * value for row in matrix for value in row)
    det2 = determinant(matrix) ** 2
    discriminant = max(0.0, norm2 * norm2 - 4.0 * det2)
    values = ((norm2 + math.sqrt(discriminant)) / 2.0,
              (norm2 - math.sqrt(discriminant)) / 2.0)
    return tuple(math.sqrt(max(0.0, value)) for value in values)


def matrix_metrics(matrix):
    singular = singular_values(matrix)
    condition = math.inf if singular[1] == 0.0 else singular[0] / singular[1]
    return {
        "matrix": matrix,
        "determinant": determinant(matrix),
        "singular_values": list(singular),
        "condition_number": condition,
    }


def jackknife_se(deleted):
    mean = math.fsum(deleted) / len(deleted)
    return math.sqrt(
        (len(deleted) - 1) / len(deleted) *
        math.fsum((value - mean) ** 2 for value in deleted)
    )


def jackknife_covariance(deleted_vectors):
    count = len(deleted_vectors)
    means = [math.fsum(row[index] for row in deleted_vectors) / count
             for index in range(len(deleted_vectors[0]))]
    return [[
        (count - 1) / count * math.fsum(
            (row[i] - means[i]) * (row[j] - means[j])
            for row in deleted_vectors
        )
        for j in range(len(means))
    ] for i in range(len(means))]


def _finite(value):
    return value if math.isfinite(value) else None


def _public_metrics(metrics):
    return {
        "matrix": metrics["matrix"],
        "determinant": metrics["determinant"],
        "singular_values": metrics["singular_values"],
        "condition_number": _finite(metrics["condition_number"]),
    }


def read_inputs(batch_path: Path, metadata_path: Path):
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if metadata.get("schema") != "matching-one/c4-local-odd-pivotal-score-stream/v1":
        raise ValueError("unexpected score-stream metadata schema")
    if metadata.get("radius") != 3:
        raise ValueError("Issue #155 production analysis requires frozen radius R=3")
    if metadata.get("samples_per_size") != 200000 or metadata.get("batches") != 100:
        raise ValueError("Issue #155 production analysis requires 200k samples and 100 batches")
    expected_designs = [
        {"N": size, "a": design[0], "b": design[1]}
        for size, design in DESIGNS.items()
    ]
    if metadata.get("designs") != expected_designs:
        raise ValueError("metadata designs differ from frozen N130/N170 contract")

    by_size = {size: [] for size in DESIGNS}
    with batch_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            size = int(row["n"])
            if size not in by_size:
                raise ValueError(f"unexpected size N={size}")
            a, b = DESIGNS[size]
            if (int(row["a"]), int(row["b"])) != (a, b):
                raise ValueError(f"N={size} Gaussian representation changed")
            by_size[size].append(row)

    for size, rows in by_size.items():
        rows.sort(key=lambda row: int(row["batch"]))
        if [int(row["batch"]) for row in rows] != list(range(100)):
            raise ValueError(f"N={size} lacks the 100 contiguous aligned batches")
        if sum(int(row["samples"]) for row in rows) != 200000:
            raise ValueError(f"N={size} sample count changed")
    for batch in range(100):
        alignment = {
            (int(by_size[size][batch]["counter_first"]),
             int(by_size[size][batch]["counter_last_exclusive"]),
             int(by_size[size][batch]["samples"]))
            for size in DESIGNS
        }
        if len(alignment) != 1:
            raise ValueError(f"batch {batch} is not counter-aligned across sizes")
    return metadata, by_size


def response_vector(rows):
    samples = sum(int(row["samples"]) for row in rows)
    return [
        math.fsum(int(row[column]) for row in rows) / (2.0 * samples)
        for column in RESPONSE_COLUMNS
    ]


def as_matrix(vector):
    return [vector[:2], vector[2:]]


def analyze_size(rows):
    point_vector = response_vector(rows)
    point = matrix_metrics(as_matrix(point_vector))
    deleted_vectors = [
        response_vector(rows[:batch] + rows[batch + 1:])
        for batch in range(len(rows))
    ]
    deleted_metrics = [matrix_metrics(as_matrix(vector)) for vector in deleted_vectors]
    covariance = jackknife_covariance(deleted_vectors)
    determinant_se = jackknife_se([row["determinant"] for row in deleted_metrics])
    singular_se = [
        jackknife_se([row["singular_values"][index] for row in deleted_metrics])
        for index in range(2)
    ]
    finite_conditions = [row["condition_number"] for row in deleted_metrics
                         if math.isfinite(row["condition_number"])]
    condition_se = (jackknife_se(finite_conditions)
                    if len(finite_conditions) == len(deleted_metrics) else None)
    determinant_z = (point["determinant"] / determinant_se
                     if determinant_se > 0.0 else math.inf)
    passes = (abs(determinant_z) >= ABS_DETERMINANT_Z_MIN and
              point["condition_number"] <= CONDITION_NUMBER_MAX)
    samples = sum(int(row["samples"]) for row in rows)
    fisher = {
        "tt": math.fsum(int(row["sum_score_t2"]) for row in rows) / samples,
        "lambda_lambda": math.fsum(
            int(row["sum_score_lambda2"]) for row in rows) / samples,
        "t_lambda": math.fsum(int(row["sum_score_cross"]) for row in rows) / samples,
    }
    return {
        "samples": samples,
        "response_order": [
            "global_cross_d_t", "global_cross_d_lambda",
            "local_pivotal_h4_d_t", "local_pivotal_h4_d_lambda",
        ],
        "point": _public_metrics(point),
        "response_covariance": covariance,
        "response_standard_errors": [math.sqrt(max(0.0, covariance[i][i]))
                                      for i in range(4)],
        "metric_standard_errors": {
            "determinant": determinant_se,
            "singular_values": singular_se,
            "condition_number": condition_se,
        },
        "determinant_z": _finite(determinant_z),
        "gate": {
            "absolute_determinant_z_min": ABS_DETERMINANT_Z_MIN,
            "condition_number_max": CONDITION_NUMBER_MAX,
            "passes": passes,
        },
        "empirical_fisher": fisher,
        "exact_fisher": {"tt": 4 * int(rows[0]["n"]),
                         "lambda_lambda": 4 * int(rows[0]["n"]),
                         "t_lambda": 0},
        "delete_one": [
            {"excluded_batch": batch, **_public_metrics(metrics)}
            for batch, metrics in enumerate(deleted_metrics)
        ],
        "_deleted_vectors": deleted_vectors,
    }


def inverse_2x2(matrix):
    det = determinant(matrix)
    if det == 0.0:
        raise ArithmeticError("source response matrix is singular")
    return [[matrix[1][1] / det, -matrix[0][1] / det],
            [-matrix[1][0] / det, matrix[0][0] / det]]


def multiply_2x2(left, right):
    return [[math.fsum(left[i][k] * right[k][j] for k in range(2))
             for j in range(2)] for i in range(2)]


def eigenvalue_rows(matrix, size_ratio=170 / 130):
    trace = matrix[0][0] + matrix[1][1]
    det = determinant(matrix)
    discriminant = trace * trace - 4.0 * det
    if discriminant >= 0.0:
        root = math.sqrt(discriminant)
        values = (complex((trace + root) / 2.0, 0.0),
                  complex((trace - root) / 2.0, 0.0))
    else:
        root = math.sqrt(-discriminant)
        values = (complex(trace / 2.0, root / 2.0),
                  complex(trace / 2.0, -root / 2.0))
    values = sorted(values, key=lambda value: abs(value), reverse=True)
    rows = []
    for value in values:
        modulus = abs(value)
        if modulus <= 0.0:
            raise ArithmeticError("zero generalized eigenvalue has no effective exponent")
        rows.append({
            "real": value.real,
            "imag": value.imag,
            "modulus": modulus,
            "effective_y": 2.0 * math.log(modulus) / math.log(size_ratio),
        })
    return rows


def generalized_diagnostic(size_results):
    if not all(size_results[size]["gate"]["passes"] for size in DESIGNS):
        return {
            "computed": False,
            "resolved_dimension": 1,
            "reason": "both N130 and N170 must pass determinant-z and condition-number gates",
        }
    source = size_results[130]["point"]["matrix"]
    target = size_results[170]["point"]["matrix"]
    transfer = multiply_2x2(target, inverse_2x2(source))
    point = eigenvalue_rows(transfer)
    deleted = []
    for batch in range(100):
        source_deleted = as_matrix(size_results[130]["_deleted_vectors"][batch])
        target_deleted = as_matrix(size_results[170]["_deleted_vectors"][batch])
        matrix = multiply_2x2(target_deleted, inverse_2x2(source_deleted))
        deleted.append(eigenvalue_rows(matrix))
    for branch in range(2):
        for field in ("real", "imag", "modulus", "effective_y"):
            point[branch][f"{field}_jackknife_se"] = jackknife_se(
                [row[branch][field] for row in deleted]
            )
    return {
        "computed": True,
        "resolved_dimension": 2,
        "size_ratio": 170 / 130,
        "transfer_R170_times_inverse_R130": transfer,
        "branches_sorted_by_modulus": point,
        "delete_one": [
            {"excluded_batch": batch, "branches_sorted_by_modulus": rows}
            for batch, rows in enumerate(deleted)
        ],
        "interpretation": "exploratory matrix-pencil diagnostic; no free correction fit",
    }


def render(batch_path: Path, metadata_path: Path):
    metadata, by_size = read_inputs(batch_path, metadata_path)
    size_results = {size: analyze_size(rows) for size, rows in by_size.items()}
    generalized = generalized_diagnostic(size_results)
    for result in size_results.values():
        del result["_deleted_vectors"]
    return {
        "schema": "matching-one/c4-local-odd-pivotal-largeN-gate/v1",
        "issue": 155,
        "design": {
            "sizes": [{"N": size, "a": a, "b": b}
                      for size, (a, b) in DESIGNS.items()],
            "radius": 3,
            "observable_rows": [
                "global_cross_half_difference",
                "local_pivotal_h4_half_difference",
            ],
            "score_columns": ["S_t", "S_lambda"],
            "cross_size_coupling": metadata["cross_size_coupling"],
        },
        "sizes": {str(size): result for size, result in size_results.items()},
        "generalized_eigensystem": generalized,
        "decision": (
            "two_dimensions_resolved_exploratory_pencil_reported"
            if generalized["computed"] else
            "only_one_dimension_resolved_no_second_eigenvalue_claim"
        ),
        "inputs": {
            "batches": str(batch_path),
            "batches_sha256": sha256(batch_path),
            "metadata": str(metadata_path),
            "metadata_sha256": sha256(metadata_path),
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batches", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = render(args.batches, args.metadata)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Verify a supplied rational covariance-subspace discrepancy exactly."""

from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results/model-certificates/framework/confidence-subspace/latest.json"
SCHEMA = "matching-one/exact-confidence-subspace-certificate/v1"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _vector(values: Sequence[Any]) -> list[Fraction]:
    return [Fraction(value) for value in values]


def _matrix(values: Sequence[Sequence[Any]]) -> list[list[Fraction]]:
    matrix = [_vector(row) for row in values]
    _require(matrix and matrix[0] and all(len(row) == len(matrix[0]) for row in matrix), "matrix must be rectangular")
    return matrix


def transpose(matrix: Sequence[Sequence[Fraction]]) -> list[list[Fraction]]:
    return [list(column) for column in zip(*matrix)]


def multiply(left: Sequence[Sequence[Fraction]], right: Sequence[Sequence[Fraction]]) -> list[list[Fraction]]:
    _require(left and right and len(left[0]) == len(right), "matrix product shape mismatch")
    return [
        [sum((left[row][inner] * right[inner][column] for inner in range(len(right))), Fraction()) for column in range(len(right[0]))]
        for row in range(len(left))
    ]


def matvec(matrix: Sequence[Sequence[Fraction]], vector: Sequence[Fraction]) -> list[Fraction]:
    _require(matrix and all(len(row) == len(vector) for row in matrix), "matrix/vector shape mismatch")
    return [sum((entry * value for entry, value in zip(row, vector)), Fraction()) for row in matrix]


def matrix_rank(values: Sequence[Sequence[Fraction]]) -> int:
    work = [list(row) for row in values]
    _require(work and work[0] and all(len(row) == len(work[0]) for row in work), "rank requires a rectangular matrix")
    rank = 0
    for column in range(len(work[0])):
        pivot = next((row for row in range(rank, len(work)) if work[row][column]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        pivot_value = work[rank][column]
        work[rank] = [value / pivot_value for value in work[rank]]
        for row in range(len(work)):
            if row == rank or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [left - factor * right for left, right in zip(work[row], work[rank])]
        rank += 1
    return rank


def determinant(values: Sequence[Sequence[Fraction]]) -> Fraction:
    dimension = len(values)
    _require(dimension and all(len(row) == dimension for row in values), "determinant requires a square matrix")
    work = [list(row) for row in values]
    result = Fraction(1)
    for column in range(dimension):
        pivot = next((row for row in range(column, dimension) if work[row][column]), None)
        if pivot is None:
            return Fraction()
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            result = -result
        pivot_value = work[column][column]
        result *= pivot_value
        for row in range(column + 1, dimension):
            factor = work[row][column] / pivot_value
            for inner in range(column, dimension):
                work[row][inner] -= factor * work[column][inner]
    return result


def _is_positive_semidefinite(matrix: Sequence[Sequence[Fraction]]) -> bool:
    dimension = len(matrix)
    for size in range(1, dimension + 1):
        for indices in combinations(range(dimension), size):
            principal = [[matrix[row][column] for column in indices] for row in indices]
            if determinant(principal) < 0:
                return False
    return True


def verify_discrepancy(descriptor: Mapping[str, Any]) -> Mapping[str, Any]:
    required = {"observed", "predicted", "covariance", "covariance_pseudoinverse", "cutoff", "calibration"}
    _require(set(descriptor) == required, "descriptor fields drift")
    observed = _vector(descriptor["observed"])
    predicted = _vector(descriptor["predicted"])
    covariance = _matrix(descriptor["covariance"])
    pseudoinverse = _matrix(descriptor["covariance_pseudoinverse"])
    dimension = len(observed)
    _require(len(predicted) == dimension, "observed/predicted dimension mismatch")
    _require(len(covariance) == dimension and all(len(row) == dimension for row in covariance), "covariance dimension mismatch")
    _require(len(pseudoinverse) == dimension and all(len(row) == dimension for row in pseudoinverse), "pseudoinverse dimension mismatch")
    _require(covariance == transpose(covariance), "covariance must be symmetric")
    _require(pseudoinverse == transpose(pseudoinverse), "pseudoinverse must be symmetric")
    _require(_is_positive_semidefinite(covariance), "covariance must be positive semidefinite")
    covariance_projection = multiply(covariance, pseudoinverse)
    pseudoinverse_projection = multiply(pseudoinverse, covariance)
    _require(multiply(multiply(covariance, pseudoinverse), covariance) == covariance, "C C+ C identity failed")
    _require(multiply(multiply(pseudoinverse, covariance), pseudoinverse) == pseudoinverse, "C+ C C+ identity failed")
    _require(covariance_projection == transpose(covariance_projection), "C C+ symmetry identity failed")
    _require(pseudoinverse_projection == transpose(pseudoinverse_projection), "C+ C symmetry identity failed")
    residual = [left - right for left, right in zip(observed, predicted)]
    projected = matvec(covariance_projection, residual)
    _require(projected == residual, "residual lies outside covariance range")
    weighted = matvec(pseudoinverse, residual)
    discrepancy = sum((left * right for left, right in zip(residual, weighted)), Fraction())
    cutoff = Fraction(descriptor["cutoff"])
    _require(cutoff >= 0, "cutoff must be nonnegative")
    _require(descriptor["calibration"] == "supplied_external_not_verified", "calibration boundary drift")
    relation = "inside_or_on" if discrepancy <= cutoff else "outside"
    return {
        "ambient_dimension": dimension,
        "covariance_rank": matrix_rank(covariance),
        "covariance_positive_semidefinite": True,
        "moore_penrose_identities": "all_four_exact",
        "residual": [str(value) for value in residual],
        "range_projection_residual": [str(left - right) for left, right in zip(projected, residual)],
        "quadratic_discrepancy": str(discrepancy),
        "cutoff": str(cutoff),
        "relation_to_supplied_outer_set": relation,
        "status": "exact_covariance_subspace_discrepancy_verified",
    }


def frozen_descriptor() -> dict[str, Any]:
    return {
        "observed": ["2", "1", "4"],
        "predicted": ["1", "2", "4"],
        "covariance": [["1", "0", "0"], ["0", "4", "0"], ["0", "0", "0"]],
        "covariance_pseudoinverse": [["1", "0", "0"], ["0", "1/4", "0"], ["0", "0", "0"]],
        "cutoff": "3/2",
        "calibration": "supplied_external_not_verified",
    }


def build_result() -> dict[str, Any]:
    descriptor = frozen_descriptor()
    return {
        "schema": SCHEMA,
        "issue": 370,
        "claim_level": "exact_control",
        "dependency_group": "synthetic-exact-covariance-subspace",
        "descriptor": descriptor,
        "verification": verify_discrepancy(descriptor),
        "claim_boundary": {
            "included": "exact rational verification of a supplied singular covariance, pseudoinverse, range condition, discrepancy, and cutoff comparison",
            "excluded": "confidence-set calibration, covariance estimation, finite-sample coverage, interval enclosure of measured inputs, model acceptance, or physical interpretation",
            "parent_issue": "remain open",
        },
    }


def validate_result(result: Mapping[str, Any]) -> Mapping[str, Any]:
    expected = build_result()
    _require(result == expected, "confidence-subspace certificate does not exactly reproduce")
    verification = result["verification"]
    return {
        "schema": result["schema"],
        "status": "valid_exact_confidence_subspace_certificate",
        "covariance_rank": verification["covariance_rank"],
        "discrepancy": verification["quadratic_discrepancy"],
        "relation": verification["relation_to_supplied_outer_set"],
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--validate", type=Path)
    args = parser.parse_args(argv)
    if args.validate:
        value = json.loads(args.validate.read_text(encoding="utf-8"))
        print(json.dumps(validate_result(value), indent=2, sort_keys=True))
        return 0
    rendered = json.dumps(build_result(), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

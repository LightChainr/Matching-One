#!/usr/bin/env python3
"""Exact rational 2x2 Jordan-geometry and coalescence oracle."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Any, Mapping

from exact_jet_algebra import canonical_fraction


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = ROOT / "analysis" / "jordan_geometry_contract.json"
EXPECTED_SCHEMA = "matching-one/jordan-geometry/v1"
Matrix = tuple[tuple[Fraction, Fraction], tuple[Fraction, Fraction]]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def parse_matrix(raw: Any, label: str) -> Matrix:
    _require(
        isinstance(raw, list)
        and len(raw) == 2
        and all(isinstance(row, list) and len(row) == 2 for row in raw),
        f"{label} must be 2x2",
    )
    return tuple(
        tuple(canonical_fraction(value, f"{label}[{i}][{j}]") for j, value in enumerate(row))
        for i, row in enumerate(raw)
    )  # type: ignore[return-value]


def identity() -> Matrix:
    return ((Fraction(1), Fraction(0)), (Fraction(0), Fraction(1)))


def zero_matrix() -> Matrix:
    return ((Fraction(0), Fraction(0)), (Fraction(0), Fraction(0)))


def matmul(left: Matrix, right: Matrix) -> Matrix:
    return tuple(
        tuple(
            sum((left[i][k] * right[k][j] for k in range(2)), Fraction(0))
            for j in range(2)
        )
        for i in range(2)
    )  # type: ignore[return-value]


def matrix_add(left: Matrix, right: Matrix) -> Matrix:
    return tuple(
        tuple(left[i][j] + right[i][j] for j in range(2)) for i in range(2)
    )  # type: ignore[return-value]


def matrix_scale(value: Fraction, matrix: Matrix) -> Matrix:
    return tuple(
        tuple(value * matrix[i][j] for j in range(2)) for i in range(2)
    )  # type: ignore[return-value]


def trace(matrix: Matrix) -> Fraction:
    return matrix[0][0] + matrix[1][1]


def determinant(matrix: Matrix) -> Fraction:
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def inverse_2x2(matrix: Matrix) -> Matrix:
    det = determinant(matrix)
    _require(det != 0, "similarity matrix is singular")
    a, b = matrix[0]
    c, d = matrix[1]
    return ((d / det, -b / det), (-c / det, a / det))


def subtract_scalar_identity(matrix: Matrix, value: Fraction) -> Matrix:
    return (
        (matrix[0][0] - value, matrix[0][1]),
        (matrix[1][0], matrix[1][1] - value),
    )


def discriminant(matrix: Matrix) -> Fraction:
    return trace(matrix) ** 2 - 4 * determinant(matrix)


def jordan_certificate(matrix: Matrix) -> dict[str, Any]:
    eigenvalue = trace(matrix) / 2
    nilpotent_part = subtract_scalar_identity(matrix, eigenvalue)
    nilpotent_square = matmul(nilpotent_part, nilpotent_part)
    repeated = discriminant(matrix) == 0
    non_scalar = nilpotent_part != zero_matrix()
    square_zero = nilpotent_square == zero_matrix()
    return {
        "eigenvalue": str(eigenvalue),
        "repeated_eigenvalue": repeated,
        "nonzero_nilpotent_part": non_scalar,
        "nilpotent_part_squares_to_zero": square_zero,
        "is_size_two_jordan_block": repeated and non_scalar and square_zero,
    }


def similarity_transform(matrix: Matrix, similarity: Matrix) -> Matrix:
    return matmul(matmul(inverse_2x2(similarity), matrix), similarity)


def cocycle_matrix(generator: Matrix, parameter: Fraction) -> Matrix:
    return matrix_add(identity(), matrix_scale(parameter, generator))


def validate_nilpotent_cocycle(generator: Matrix, first: Fraction, second: Fraction) -> bool:
    _require(generator != zero_matrix(), "nilpotent generator must be nonzero")
    _require(matmul(generator, generator) == zero_matrix(), "generator must square to zero")
    left = matmul(cocycle_matrix(generator, first), cocycle_matrix(generator, second))
    right = cocycle_matrix(generator, first + second)
    return left == right


def coalescence_record(eigenvalue: Fraction, epsilon: Fraction) -> dict[str, str]:
    _require(epsilon > 0, "coalescing epsilon must be positive")
    matrix: Matrix = ((eigenvalue, Fraction(1)), (Fraction(0), eigenvalue + epsilon))
    angle_squared = epsilon**2 / (1 + epsilon**2)
    return {
        "epsilon": str(epsilon),
        "discriminant": str(discriminant(matrix)),
        "eigenbasis_determinant": str(epsilon),
        "squared_sine_angle": str(angle_squared),
    }


def validate_contract(contract: Mapping[str, Any]) -> dict[str, Any]:
    _require(contract.get("schema") == EXPECTED_SCHEMA, "schema drifted")
    _require(contract.get("issue") == 218, "issue must be 218")
    _require(contract.get("status") == "exact_synthetic_2x2_only", "status drifted")

    eigenvalue = canonical_fraction(contract.get("eigenvalue"), "eigenvalue")
    jordan = parse_matrix(contract.get("jordan_matrix"), "Jordan matrix")
    expected_jordan: Matrix = (
        (eigenvalue, Fraction(1)),
        (Fraction(0), eigenvalue),
    )
    _require(jordan == expected_jordan, "declared Jordan normal form drifted")
    certificate = jordan_certificate(jordan)
    _require(certificate["is_size_two_jordan_block"], "matrix is not a size-two Jordan block")

    similarity = parse_matrix(contract.get("similarity"), "similarity")
    transformed = similarity_transform(jordan, similarity)
    transformed_certificate = jordan_certificate(transformed)
    _require(transformed_certificate["is_size_two_jordan_block"], "Jordan structure was not similarity invariant")
    _require(trace(transformed) == trace(jordan), "trace was not similarity invariant")
    _require(determinant(transformed) == determinant(jordan), "determinant was not similarity invariant")
    _require(discriminant(transformed) == discriminant(jordan), "discriminant was not similarity invariant")

    generator = parse_matrix(contract.get("nilpotent_generator"), "nilpotent generator")
    raw_parameters = contract.get("cocycle_parameters")
    _require(isinstance(raw_parameters, list) and len(raw_parameters) == 2, "two cocycle parameters are required")
    parameters = tuple(
        canonical_fraction(value, f"cocycle parameter {index}")
        for index, value in enumerate(raw_parameters)
    )
    _require(validate_nilpotent_cocycle(generator, *parameters), "nilpotent cocycle law failed")

    raw_epsilons = contract.get("coalescing_epsilons")
    _require(isinstance(raw_epsilons, list) and len(raw_epsilons) >= 2, "at least two epsilons are required")
    epsilons = [
        canonical_fraction(value, f"coalescing epsilon {index}")
        for index, value in enumerate(raw_epsilons)
    ]
    _require(all(left > right for left, right in zip(epsilons, epsilons[1:])), "epsilons must strictly decrease")
    records = [coalescence_record(eigenvalue, epsilon) for epsilon in epsilons]
    angles = [canonical_fraction(record["squared_sine_angle"], "squared sine angle") for record in records]
    _require(all(left > right for left, right in zip(angles, angles[1:])), "eigenvectors did not coalesce monotonically")

    return {
        "schema": EXPECTED_SCHEMA,
        "status": "valid_exact_jordan_geometry",
        "jordan_certificate": certificate,
        "similarity_invariants_preserved": True,
        "nilpotent_cocycle_composes": True,
        "coalescence_family": records,
        "coalescence_monotone": True,
        "contains_empirical_transfer_matrix_claim": False,
        "parent_issue": "remain open",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    args = parser.parse_args()
    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    print(json.dumps(validate_contract(contract), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

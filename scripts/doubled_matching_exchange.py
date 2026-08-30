#!/usr/bin/env python3
"""Exact doubled matching-pair exchange and RG-intertwiner oracle."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

from exact_jet_algebra import canonical_fraction


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = ROOT / "analysis" / "doubled_matching_exchange_contract.json"
EXPECTED_SCHEMA = "matching-one/doubled-matching-exchange/v1"
Matrix = tuple[tuple[Fraction, ...], ...]
Vector = tuple[Fraction, ...]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def parse_matrix(raw: Any, label: str) -> Matrix:
    _require(isinstance(raw, list) and bool(raw), f"{label} must be a nonempty matrix")
    width = len(raw[0]) if isinstance(raw[0], list) else 0
    _require(width > 0 and all(isinstance(row, list) and len(row) == width for row in raw), f"{label} is ragged")
    return tuple(
        tuple(canonical_fraction(value, f"{label}[{i}][{j}]") for j, value in enumerate(row))
        for i, row in enumerate(raw)
    )


def parse_vector(raw: Any, length: int, label: str) -> Vector:
    _require(isinstance(raw, list) and len(raw) == length, f"{label} must have length {length}")
    return tuple(canonical_fraction(value, f"{label}[{index}]") for index, value in enumerate(raw))


def identity(size: int) -> Matrix:
    return tuple(tuple(Fraction(i == j) for j in range(size)) for i in range(size))


def matmul(left: Matrix, right: Matrix) -> Matrix:
    _require(bool(left) and bool(right), "matrices must not be empty")
    _require(len(left[0]) == len(right), "matrix dimensions do not align")
    _require(all(len(row) == len(left[0]) for row in left), "left matrix is ragged")
    _require(all(len(row) == len(right[0]) for row in right), "right matrix is ragged")
    return tuple(
        tuple(
            sum((left[i][k] * right[k][j] for k in range(len(right))), Fraction(0))
            for j in range(len(right[0]))
        )
        for i in range(len(left))
    )


def matvec(matrix: Matrix, vector: Vector) -> Vector:
    _require(bool(matrix) and len(matrix[0]) == len(vector), "matrix/vector dimensions do not align")
    return tuple(sum((entry * value for entry, value in zip(row, vector)), Fraction(0)) for row in matrix)


def inverse_2x2(matrix: Matrix) -> Matrix:
    _require(len(matrix) == 2 and all(len(row) == 2 for row in matrix), "matrix must be 2x2")
    a, b = matrix[0]
    c, d = matrix[1]
    determinant = a * d - b * c
    _require(determinant != 0, "matrix is singular")
    return ((d / determinant, -b / determinant), (-c / determinant, a / determinant))


def block_diagonal(left: Matrix, right: Matrix) -> Matrix:
    _require(len(left) == len(left[0]) and len(right) == len(right[0]), "blocks must be square")
    zero_left = (Fraction(0),) * len(right)
    zero_right = (Fraction(0),) * len(left)
    return tuple(tuple(row) + zero_left for row in left) + tuple(zero_right + tuple(row) for row in right)


def exchange_matrix(a: Matrix, a_inverse: Matrix) -> Matrix:
    size = len(a)
    _require(size == len(a[0]) == len(a_inverse) == len(a_inverse[0]), "matching maps must be equally sized square matrices")
    zero = (Fraction(0),) * size
    return tuple(zero + tuple(row) for row in a_inverse) + tuple(tuple(row) + zero for row in a)


def lift_exchange_vector(vector: Vector, a: Matrix, parity: int) -> Vector:
    _require(parity in (-1, 1), "exchange parity must be +1 or -1")
    image = matvec(a, vector)
    return vector + tuple(parity * value for value in image)


def validate_contract(contract: Mapping[str, Any]) -> dict[str, Any]:
    _require(contract.get("schema") == EXPECTED_SCHEMA, "schema drifted")
    _require(contract.get("issue") == 61, "issue must be 61")
    _require(
        contract.get("status") == "exact_finite_dimensional_linearization_only",
        "status drifted",
    )
    a = parse_matrix(contract.get("matching_map_A"), "A")
    declared_inverse = parse_matrix(contract.get("matching_map_A_inverse"), "A inverse")
    rg = parse_matrix(contract.get("RG_G"), "RG_G")
    rg_star = parse_matrix(contract.get("RG_G_star"), "RG_G_star")
    size = len(a)
    _require(size == 2 and all(len(matrix) == 2 and len(matrix[0]) == 2 for matrix in (a, rg, rg_star)), "oracle matrices must remain 2x2")
    a_inverse = inverse_2x2(a)
    _require(a_inverse == declared_inverse, "declared inverse drifted")
    _require(matmul(a_inverse, a) == identity(size), "left inverse failed")
    _require(matmul(a, a_inverse) == identity(size), "right inverse failed")

    intertwines = matmul(a, rg) == matmul(rg_star, a)
    _require(intertwines, "A does not intertwine the RG maps")
    exchange = exchange_matrix(a, a_inverse)
    doubled_rg = block_diagonal(rg, rg_star)
    _require(matmul(exchange, exchange) == identity(2 * size), "doubled exchange is not involutive")
    _require(matmul(exchange, doubled_rg) == matmul(doubled_rg, exchange), "doubled exchange does not commute with RG")

    probe = parse_vector(contract.get("probe_vector"), size, "probe vector")
    _require(any(probe), "probe vector must be nonzero")
    eigenspaces = {}
    for parity in (-1, 1):
        lifted = lift_exchange_vector(probe, a, parity)
        _require(matvec(exchange, lifted) == tuple(parity * value for value in lifted), "exchange eigenspace lift failed")
        eigenspaces[str(parity)] = [str(value) for value in lifted]

    degenerate = parse_matrix(contract.get("degenerate_RG"), "degenerate RG")
    b = parse_matrix(contract.get("alternative_identification_B"), "B")
    _require(b != a, "alternative identification must differ from A")
    _require(inverse_2x2(b), "alternative identification must be invertible")
    _require(matmul(a, degenerate) == matmul(degenerate, a), "A must intertwine the degenerate block")
    _require(matmul(b, degenerate) == matmul(degenerate, b), "B must intertwine the degenerate block")
    _require(lift_exchange_vector(probe, a, 1) != lift_exchange_vector(probe, b, 1), "degenerate identifications did not change the parity basis")

    return {
        "schema": EXPECTED_SCHEMA,
        "status": "valid_exact_doubled_exchange",
        "dimension_per_tangent_space": size,
        "matching_map_invertible": True,
        "matching_map_intertwines_rg": True,
        "doubled_exchange_involutive": True,
        "doubled_exchange_commutes_with_rg": True,
        "exchange_eigenvectors": eigenspaces,
        "degenerate_rg_does_not_fix_identification": True,
        "contains_physical_parity_assignment": False,
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

#!/usr/bin/env python3
"""Exact normalized odd-jet invariants and covariance propagation."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence, Tuple


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = ROOT / "analysis" / "odd_jet_invariant_contract.json"
FORBIDDEN_KEYS = frozenset(
    {
        "estimate",
        "fit",
        "p_value",
        "root",
        "samples",
        "target_data",
        "universality_claim",
    }
)
Matrix = Tuple[Tuple[Fraction, ...], ...]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _walk_forbidden(value: Any, path: str = "$") -> None:
    if isinstance(value, Mapping):
        bad = sorted(FORBIDDEN_KEYS.intersection(value))
        _require(not bad, "%s contains empirical fields: %s" % (path, ",".join(bad)))
        for key, child in value.items():
            _walk_forbidden(child, "%s.%s" % (path, key))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            _walk_forbidden(child, "%s[%d]" % (path, index))


def canonical_fraction(value: Any, label: str) -> Fraction:
    _require(type(value) is str, "%s must be an exact rational string" % label)
    try:
        parsed = Fraction(value)
    except (ValueError, ZeroDivisionError) as error:
        raise ValueError("%s is not rational" % label) from error
    _require(str(parsed) == value, "%s is not canonically encoded" % label)
    return parsed


def validate_orders(orders: Sequence[int]) -> Tuple[int, ...]:
    normalized = tuple(orders)
    _require(bool(normalized) and normalized[0] == 1, "odd jet must start at derivative order 1")
    _require(all(type(order) is int and order > 0 and order % 2 == 1 for order in normalized), "derivative orders must be positive odd integers")
    _require(tuple(sorted(set(normalized))) == normalized, "derivative orders must be unique and increasing")
    return normalized


def validate_jet(orders: Sequence[int], derivatives: Sequence[Fraction]) -> Tuple[Tuple[int, ...], Tuple[Fraction, ...]]:
    valid_orders = validate_orders(orders)
    values = tuple(derivatives)
    _require(len(values) == len(valid_orders), "derivative count does not match orders")
    _require(values[0] != 0, "first derivative must be nonzero")
    return valid_orders, values


def normalized_invariants(orders: Sequence[int], derivatives: Sequence[Fraction]) -> Mapping[int, Fraction]:
    valid_orders, values = validate_jet(orders, derivatives)
    first = values[0]
    return {order: value / first**order for order, value in zip(valid_orders[1:], values[1:])}


def invariant_jacobian(orders: Sequence[int], derivatives: Sequence[Fraction]) -> Matrix:
    valid_orders, values = validate_jet(orders, derivatives)
    first = values[0]
    rows = []
    for row_index, order in enumerate(valid_orders[1:], start=1):
        row = [Fraction(0) for _ in valid_orders]
        row[0] = -order * values[row_index] / first ** (order + 1)
        row[row_index] = Fraction(1) / first**order
        rows.append(tuple(row))
    return tuple(rows)


def validate_covariance(covariance: Sequence[Sequence[Fraction]], dimension: int) -> Matrix:
    matrix = tuple(tuple(row) for row in covariance)
    _require(len(matrix) == dimension and all(len(row) == dimension for row in matrix), "covariance dimension mismatch")
    _require(all(matrix[i][j] == matrix[j][i] for i in range(dimension) for j in range(dimension)), "covariance must be symmetric")
    _require(all(matrix[i][i] >= 0 for i in range(dimension)), "covariance diagonal must be nonnegative")
    return matrix


def propagate_covariance(
    orders: Sequence[int],
    derivatives: Sequence[Fraction],
    covariance: Sequence[Sequence[Fraction]],
) -> Matrix:
    valid_orders, _values = validate_jet(orders, derivatives)
    source = validate_covariance(covariance, len(valid_orders))
    jacobian = invariant_jacobian(valid_orders, derivatives)
    return tuple(
        tuple(
            sum(
                jacobian[i][a] * source[a][b] * jacobian[j][b]
                for a in range(len(valid_orders))
                for b in range(len(valid_orders))
            )
            for j in range(len(jacobian))
        )
        for i in range(len(jacobian))
    )


def rescale_coordinate(
    orders: Sequence[int], derivatives: Sequence[Fraction], scale: Fraction
) -> Tuple[Fraction, ...]:
    valid_orders, values = validate_jet(orders, derivatives)
    _require(scale != 0, "coordinate scale must be nonzero")
    return tuple(scale**order * value for order, value in zip(valid_orders, values))


def rescale_covariance(
    orders: Sequence[int], covariance: Sequence[Sequence[Fraction]], scale: Fraction
) -> Matrix:
    valid_orders = validate_orders(orders)
    source = validate_covariance(covariance, len(valid_orders))
    _require(scale != 0, "coordinate scale must be nonzero")
    return tuple(
        tuple(scale ** (valid_orders[i] + valid_orders[j]) * source[i][j] for j in range(len(valid_orders)))
        for i in range(len(valid_orders))
    )


def rescale_observable(derivatives: Sequence[Fraction], scale: Fraction) -> Tuple[Fraction, ...]:
    _require(scale != 0, "observable scale must be nonzero")
    return tuple(scale * value for value in derivatives)


def _parse_matrix(value: Any, label: str) -> Matrix:
    _require(isinstance(value, list), "%s must be a matrix" % label)
    return tuple(
        tuple(canonical_fraction(entry, "%s[%d][%d]" % (label, i, j)) for j, entry in enumerate(row))
        for i, row in enumerate(value)
    )


def _matrix_text(matrix: Matrix) -> list[list[str]]:
    return [[str(value) for value in row] for row in matrix]


def validate_contract(contract: Mapping[str, Any]) -> Mapping[str, Any]:
    _walk_forbidden(contract)
    _require(contract.get("schema") == "matching-one/odd-jet-invariants/v1", "unknown schema")
    _require(contract.get("issue") == 16, "wrong issue")
    _require(contract.get("status") == "exact_algebra_only_no_estimate", "scope status drift")
    _require(contract.get("coordinate_convention") == "M_tilde(u)=M(lambda*u)", "coordinate convention drift")
    orders = validate_orders(contract.get("derivative_orders", []))
    _require(orders == (1, 3, 5), "fixture derivative orders drift")

    fixture = contract.get("synthetic_fixture", {})
    derivative_values = fixture.get("derivatives", [])
    _require(isinstance(derivative_values, list), "derivatives must be a list")
    derivatives = tuple(canonical_fraction(value, "derivative[%d]" % index) for index, value in enumerate(derivative_values))
    covariance = _parse_matrix(fixture.get("covariance"), "covariance")
    validate_jet(orders, derivatives)
    validate_covariance(covariance, len(orders))

    invariants = normalized_invariants(orders, derivatives)
    invariant_covariance = propagate_covariance(orders, derivatives, covariance)
    expected = contract.get("expected", {})
    expected_invariants = expected.get("invariants", {})
    _require(
        {str(order): str(value) for order, value in invariants.items()} == expected_invariants,
        "stored normalized invariants drift",
    )
    _require(_matrix_text(invariant_covariance) == expected.get("invariant_covariance"), "stored invariant covariance drift")

    scale_values = fixture.get("coordinate_rescales")
    _require(isinstance(scale_values, list) and bool(scale_values), "coordinate rescale fixtures missing")
    audited_scales = []
    for index, value in enumerate(scale_values):
        scale = canonical_fraction(value, "coordinate_rescales[%d]" % index)
        _require(scale != 0, "coordinate rescale fixture must be nonzero")
        changed_derivatives = rescale_coordinate(orders, derivatives, scale)
        changed_covariance = rescale_covariance(orders, covariance, scale)
        _require(normalized_invariants(orders, changed_derivatives) == invariants, "invariants changed under coordinate rescaling")
        _require(
            propagate_covariance(orders, changed_derivatives, changed_covariance) == invariant_covariance,
            "invariant covariance changed under coordinate rescaling",
        )
        audited_scales.append(str(scale))

    boundary = contract.get("claim_boundary", {})
    _require(boundary.get("parent_issue") == "remain open", "parent boundary changed")
    return {
        "schema": contract["schema"],
        "status": "valid_exact_algebra_only",
        "invariants": {str(order): str(value) for order, value in invariants.items()},
        "invariant_covariance": _matrix_text(invariant_covariance),
        "coordinate_rescales_verified": audited_scales,
        "contains_estimate": False,
        "parent_issue": "remain open",
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    args = parser.parse_args(argv)
    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    print(json.dumps(validate_contract(contract), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Exact finite derivative chains for the Bernoulli log-odds coordinate."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

from exact_jet_algebra import (
    canonical_fraction,
    compose_derivatives,
    derivatives_to_series,
    inverse_series,
    parse_derivatives,
    series_to_derivatives,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = ROOT / "analysis" / "logodds_derivative_contract.json"
EXPECTED_SCHEMA = "matching-one/logodds-derivative-chain/v1"
FORBIDDEN_KEYS = {
    "basis_rotation",
    "covariance",
    "estimate",
    "histograms",
    "model_score",
    "samples",
    "transfer_eigenvalue",
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _walk_forbidden(value: Any, path: str = "$") -> None:
    if isinstance(value, dict):
        bad = sorted(set(value) & FORBIDDEN_KEYS)
        _require(not bad, f"{path} contains forbidden empirical keys: {','.join(bad)}")
        for key, child in value.items():
            _walk_forbidden(child, f"{path}.{key}")
    elif isinstance(value, list):
        for index, child in enumerate(value):
            _walk_forbidden(child, f"{path}[{index}]")


def _poly_multiply(left: Sequence[Fraction], right: Sequence[Fraction]) -> list[Fraction]:
    output = [Fraction(0) for _ in range(len(left) + len(right) - 1)]
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            output[i + j] += left_value * right_value
    return output


def _poly_derivative(coefficients: Sequence[Fraction]) -> list[Fraction]:
    return [degree * coefficients[degree] for degree in range(1, len(coefficients))]


def _poly_evaluate(coefficients: Sequence[Fraction], value: Fraction) -> Fraction:
    result = Fraction(0)
    for coefficient in reversed(coefficients):
        result = result * value + coefficient
    return result


def logistic_delta_derivatives(center: Fraction, order: int) -> list[Fraction]:
    """Return derivatives of ``p(eta)-p(eta0)`` at a rational center."""

    _require(0 < center < 1, "center probability must lie strictly between zero and one")
    _require(order >= 1, "order must be positive")
    derivative_polynomial = [Fraction(0), Fraction(1)]
    values: list[Fraction] = []
    logistic_vector_field = [Fraction(0), Fraction(1), Fraction(-1)]
    for _ in range(order + 1):
        values.append(_poly_evaluate(derivative_polynomial, center))
        derivative_polynomial = _poly_multiply(
            _poly_derivative(derivative_polynomial), logistic_vector_field
        )
    values[0] = Fraction(0)
    return values


def response_power_coefficients(conditional_values: Sequence[Fraction]) -> list[Fraction]:
    """Convert ``E[q_K]`` under ``Bin(N,p)`` to ascending power coefficients."""

    _require(bool(conditional_values), "conditional response must not be empty")
    n = len(conditional_values) - 1
    output = [Fraction(0) for _ in range(n + 1)]
    for k, value in enumerate(conditional_values):
        bernstein = math.comb(n, k) * value
        for degree in range(k, n + 1):
            output[degree] += bernstein * (-1) ** (degree - k) * math.comb(
                n - k, degree - k
            )
    return output


def polynomial_derivatives(
    coefficients: Sequence[Fraction], center: Fraction, order: int
) -> list[Fraction]:
    output = []
    for derivative_order in range(order + 1):
        value = Fraction(0)
        for degree in range(derivative_order, len(coefficients)):
            falling = math.factorial(degree) // math.factorial(degree - derivative_order)
            value += coefficients[degree] * falling * center ** (degree - derivative_order)
        output.append(value)
    return output


def p_to_eta_derivatives(
    response_p_derivatives: Sequence[Fraction], center: Fraction, order: int
) -> list[Fraction]:
    _require(len(response_p_derivatives) >= order + 1, "p derivative jet is too short")
    logistic = logistic_delta_derivatives(center, order)
    return compose_derivatives(response_p_derivatives, logistic, order)


def eta_to_p_derivatives(
    response_eta_derivatives: Sequence[Fraction], center: Fraction, order: int
) -> list[Fraction]:
    _require(len(response_eta_derivatives) >= order + 1, "eta derivative jet is too short")
    logistic = logistic_delta_derivatives(center, order)
    eta_delta_series = inverse_series(derivatives_to_series(logistic), order)
    eta_delta_derivatives = series_to_derivatives(eta_delta_series)
    return compose_derivatives(response_eta_derivatives, eta_delta_derivatives, order)


def response_chain(
    conditional_values: Sequence[Fraction], center: Fraction, order: int
) -> tuple[list[Fraction], list[Fraction]]:
    coefficients = response_power_coefficients(conditional_values)
    p_derivatives = polynomial_derivatives(coefficients, center, order)
    return p_derivatives, p_to_eta_derivatives(p_derivatives, center, order)


def complement_values(
    conditional_values: Sequence[Fraction], parity: int
) -> list[Fraction]:
    _require(parity in (-1, 1), "complement parity must be +1 or -1")
    return [parity * value for value in reversed(conditional_values)]


def validate_contract(contract: Mapping[str, Any]) -> dict[str, Any]:
    _walk_forbidden(contract)
    _require(contract.get("schema") == EXPECTED_SCHEMA, "schema drifted")
    _require(contract.get("issue") == 182, "issue must be 182")
    _require(contract.get("status") == "exact_finite_check_only", "status drifted")
    order = contract.get("maximum_order")
    _require(order == 6, "maximum order must remain six")
    center = canonical_fraction(contract.get("center_probability"), "center probability")
    _require(0 < center < 1, "center probability must lie strictly between zero and one")

    response = contract.get("synthetic_bernstein_response", {})
    n = response.get("N")
    _require(n == 6, "synthetic response size must remain six")
    conditional = parse_derivatives(
        response.get("conditional_values", []), n, "conditional response"
    )
    parity = response.get("complement_parity")
    _require(parity in (-1, 1), "complement parity must be +1 or -1")

    p_derivatives, eta_derivatives = response_chain(conditional, center, order)
    expected_p = parse_derivatives(contract.get("expected_p_derivatives", []), order, "expected p jet")
    expected_eta = parse_derivatives(
        contract.get("expected_eta_derivatives", []), order, "expected eta jet"
    )
    _require(p_derivatives == expected_p, "p derivative chain drifted")
    _require(eta_derivatives == expected_eta, "eta derivative chain drifted")
    _require(
        eta_to_p_derivatives(eta_derivatives, center, order) == p_derivatives,
        "p/eta derivative round trip failed",
    )

    complement = complement_values(conditional, parity)
    complement_center = 1 - center
    complement_p, complement_eta = response_chain(complement, complement_center, order)
    expected_complement = [parity * (-1) ** r * eta_derivatives[r] for r in range(order + 1)]
    declared_complement = parse_derivatives(
        contract.get("expected_complement_eta_derivatives", []),
        order,
        "expected complement eta jet",
    )
    _require(complement_eta == expected_complement, "general complement parity failed")
    _require(declared_complement == expected_complement, "declared complement eta jet drifted")
    _require(
        eta_to_p_derivatives(complement_eta, complement_center, order) == complement_p,
        "complement p/eta derivative round trip failed",
    )

    return {
        "schema": EXPECTED_SCHEMA,
        "status": "valid_exact_logodds_chain",
        "maximum_order": order,
        "center_probability": str(center),
        "complement_center_probability": str(complement_center),
        "p_derivatives": [str(value) for value in p_derivatives],
        "eta_derivatives": [str(value) for value in eta_derivatives],
        "complement_eta_derivatives": [str(value) for value in complement_eta],
        "round_trip_exact": True,
        "general_complement_parity_exact": True,
        "contains_empirical_result": False,
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

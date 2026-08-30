#!/usr/bin/env python3
"""Exact centered Taylor-parity oracle for the square-bond duality control."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

from exact_jet_algebra import canonical_fraction
from square_bond_duality_exact import (
    CHANNELS,
    channel_indicators,
    primal_dual_wrapping,
    square_bond_pairs,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = ROOT / "analysis" / "square_bond_centered_parity_contract.json"
EXPECTED_SCHEMA = "matching-one/square-bond-centered-parity/v1"
Polynomial = tuple[Fraction, ...]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def polynomial_multiply(left: Polynomial, right: Polynomial) -> Polynomial:
    result = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, left_value in enumerate(left):
        for j, right_value in enumerate(right):
            result[i + j] += left_value * right_value
    return tuple(result)


def polynomial_power(base: Polynomial, exponent: int) -> Polynomial:
    _require(exponent >= 0, "polynomial exponent must be nonnegative")
    result: Polynomial = (Fraction(1),)
    for _ in range(exponent):
        result = polynomial_multiply(result, base)
    return result


def centered_from_bernstein_sums(sums: Sequence[Fraction]) -> Polynomial:
    """Expand sum_k a_k p^k(1-p)^(n-k) at p=1/2+t."""

    _require(bool(sums), "Bernstein sums must not be empty")
    degree = len(sums) - 1
    result = [Fraction(0)] * (degree + 1)
    plus: Polynomial = (Fraction(1, 2), Fraction(1))
    minus: Polynomial = (Fraction(1, 2), Fraction(-1))
    for occupied, aggregate in enumerate(sums):
        term = polynomial_multiply(
            polynomial_power(plus, occupied),
            polynomial_power(minus, degree - occupied),
        )
        for power, value in enumerate(term):
            result[power] += aggregate * value
    return tuple(result)


def enumerate_bernstein_sums(length: int) -> dict[str, dict[str, tuple[Fraction, ...]]]:
    _require(length == 2, "this bounded oracle is frozen to L=2")
    pairs = square_bond_pairs(length)
    bond_count = len(pairs)
    sums = {
        channel: {
            "S": [Fraction(0)] * (bond_count + 1),
            "D": [Fraction(0)] * (bond_count + 1),
        }
        for channel in CHANNELS
    }
    for mask in range(1 << bond_count):
        occupied = bin(mask).count("1")
        primal, dual = primal_dual_wrapping(length, mask, pairs)
        primal_bits = channel_indicators(primal)
        dual_bits = channel_indicators(dual)
        for channel in CHANNELS:
            sums[channel]["S"][occupied] += Fraction(
                primal_bits[channel] + dual_bits[channel], 2
            )
            sums[channel]["D"][occupied] += Fraction(
                primal_bits[channel] - dual_bits[channel]
            )
    return {
        channel: {name: tuple(values) for name, values in sectors.items()}
        for channel, sectors in sums.items()
    }


def enumerate_centered_parity(length: int = 2) -> dict[str, Any]:
    sums = enumerate_bernstein_sums(length)
    bond_count = 2 * length * length
    channels: dict[str, Any] = {}
    d_polynomials: set[Polynomial] = set()
    for channel in CHANNELS:
        s_sums = sums[channel]["S"]
        d_sums = sums[channel]["D"]
        s_complement = all(s_sums[k] == s_sums[bond_count - k] for k in range(bond_count + 1))
        d_complement = all(d_sums[k] == -d_sums[bond_count - k] for k in range(bond_count + 1))
        s_centered = centered_from_bernstein_sums(s_sums)
        d_centered = centered_from_bernstein_sums(d_sums)
        s_odd_zero = all(value == 0 for order, value in enumerate(s_centered) if order % 2 == 1)
        d_even_zero = all(value == 0 for order, value in enumerate(d_centered) if order % 2 == 0)
        d_polynomials.add(d_centered)
        channels[channel] = {
            "S_bernstein_complement_even": s_complement,
            "D_bernstein_complement_odd": d_complement,
            "S_centered_coefficients": [str(value) for value in s_centered],
            "D_centered_coefficients": [str(value) for value in d_centered],
            "S_all_odd_derivatives_vanish_at_half": s_odd_zero,
            "D_all_even_derivatives_vanish_at_half": d_even_zero,
            "D_at_half": str(d_centered[0]),
            "D_first_derivative_at_half": str(d_centered[1]),
        }
    passed = all(
        row["S_bernstein_complement_even"]
        and row["D_bernstein_complement_odd"]
        and row["S_all_odd_derivatives_vanish_at_half"]
        and row["D_all_even_derivatives_vanish_at_half"]
        for row in channels.values()
    )
    return {
        "length": length,
        "bond_count": bond_count,
        "configurations": 1 << bond_count,
        "channels": channels,
        "distinct_D_centered_polynomials": len(d_polynomials),
        "passed": passed,
    }


def validate_contract(contract: Mapping[str, Any]) -> dict[str, Any]:
    _require(contract.get("schema") == EXPECTED_SCHEMA, "schema drifted")
    _require(contract.get("issue") == 42, "issue must be 42")
    _require(contract.get("status") == "exact_l2_centered_polynomial_only", "status drifted")
    _require(contract.get("channels") == list(CHANNELS), "channel registry drifted")
    length = contract.get("length")
    _require(length == 2, "length must remain two")
    result = enumerate_centered_parity(length)
    _require(result["passed"], "centered parity checks failed")

    declared_d = contract.get("expected_centered_D_coefficients")
    _require(isinstance(declared_d, list), "expected D polynomial must be a list")
    expected_d = [
        str(canonical_fraction(value, f"expected D coefficient {index}"))
        for index, value in enumerate(declared_d)
    ]
    expected_derivative = str(
        canonical_fraction(
            contract.get("expected_D_first_derivative_at_half"),
            "expected D first derivative",
        )
    )
    for channel, row in result["channels"].items():
        _require(row["D_centered_coefficients"] == expected_d, f"{channel} D polynomial drifted")
        _require(row["D_first_derivative_at_half"] == expected_derivative, f"{channel} D derivative drifted")
        _require(row["D_at_half"] == "0", f"{channel} D center did not vanish")
    _require(expected_derivative != "0", "D first derivative must remain nonzero")

    return {
        "schema": EXPECTED_SCHEMA,
        "status": "valid_exact_square_bond_centered_parity",
        "length": length,
        "configurations": result["configurations"],
        "channels_checked": len(result["channels"]),
        "S_is_centered_even_in_every_channel": True,
        "D_is_centered_odd_in_every_channel": True,
        "common_D_centered_polynomial": expected_d,
        "common_D_first_derivative_at_half": expected_derivative,
        "D_center_zero_does_not_force_D_derivative_zero": True,
        "contains_orientation_production_result": False,
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

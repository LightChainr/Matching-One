#!/usr/bin/env python3
"""Exact formal asymptotics of the adjacent-size matching annihilator.

Assume at the critical point

    M_L = a L^(-13/4) (1 + c L^(-q) + ...),
    M'_L = b L^(3/4) (1 + ...).

The equation L^(13/4) M_L = (L-1)^(13/4) M_(L-1) has linearized root
shift

    delta p = (a c / b) L^(-(4+q)) [q/4 + O(1/L)].

This oracle derives the bracket as an exact Fraction power series using the
generalized binomial theorem.  It does not read matching-root data.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Any, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = ROOT / "analysis" / "adjacent_annihilator_asymptotics_contract.json"
EXPECTED_SCHEMA = "matching-one/adjacent-annihilator-asymptotics/v1"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def parse_fraction(value: Any, label: str) -> Fraction:
    _require(isinstance(value, (str, int)), f"{label} must be an exact rational")
    try:
        return Fraction(value)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError(f"{label} must be an exact rational") from exc


def generalized_binomial(exponent: Fraction, order: int) -> Fraction:
    _require(order >= 0, "binomial order must be nonnegative")
    result = Fraction(1)
    for index in range(order):
        result *= exponent - index
        result /= index + 1
    return result


def one_minus_power(exponent: Fraction, max_order: int) -> tuple[Fraction, ...]:
    """Coefficients of 1-(1-x)^exponent through x^max_order."""

    _require(max_order >= 1, "series order must be positive")
    return tuple(
        Fraction(0)
        if order == 0
        else -generalized_binomial(exponent, order) * (-1) ** order
        for order in range(max_order + 1)
    )


def series_divide(
    numerator: Sequence[Fraction], denominator: Sequence[Fraction], terms: int
) -> tuple[Fraction, ...]:
    _require(terms >= 1, "quotient must contain at least one term")
    _require(bool(denominator) and denominator[0] != 0, "series denominator must have nonzero constant")
    result = [Fraction(0)] * terms
    for order in range(terms):
        value = numerator[order] if order < len(numerator) else Fraction(0)
        for lower in range(order):
            if order - lower < len(denominator):
                value -= result[lower] * denominator[order - lower]
        result[order] = value / denominator[0]
    return tuple(result)


def annihilator_profile(q: Fraction, terms: int = 5) -> dict[str, Any]:
    _require(q > 0, "relative correction q must be positive")
    _require(terms >= 2, "at least two correction terms are required")
    # After x=1/L, strip the common leading x from
    # A=1-(1-x)^(-q) and B=1-(1-x)^4.  The root shift is -x^(q+4) A/B.
    numerator_full = one_minus_power(-q, terms)
    denominator_full = one_minus_power(Fraction(4), terms)
    numerator = numerator_full[1:]
    denominator = denominator_full[1:]
    quotient = tuple(-value for value in series_divide(numerator, denominator, terms - 1))
    _require(quotient[0] == q / 4, "leading annihilator coefficient drifted")
    return {
        "q": str(q),
        "root_power_w": str(q + 4),
        "normalized_shift_series": [str(value) for value in quotient],
        "leading_coefficient": str(quotient[0]),
    }


def validate_contract(contract: Mapping[str, Any]) -> dict[str, Any]:
    _require(contract.get("schema") == EXPECTED_SCHEMA, "schema drifted")
    _require(contract.get("issue") == 47, "issue must be 47")
    _require(contract.get("status") == "formal_series_only", "status drifted")
    raw_candidates = contract.get("candidates")
    _require(isinstance(raw_candidates, list) and raw_candidates, "candidates must be nonempty")
    candidates = [parse_fraction(value, f"candidates[{index}]") for index, value in enumerate(raw_candidates)]
    _require(candidates == sorted(candidates) and len(set(candidates)) == len(candidates), "candidates must be sorted and distinct")
    expected = contract.get("expected")
    _require(isinstance(expected, Mapping), "expected profiles must be a mapping")

    profiles = []
    for q in candidates:
        profile = annihilator_profile(q)
        frozen = expected.get(str(q))
        _require(isinstance(frozen, Mapping), f"missing expected profile for q={q}")
        _require(profile["root_power_w"] == str(parse_fraction(frozen.get("w"), f"q={q} w")), f"q={q} root power drifted")
        _require(profile["leading_coefficient"] == str(parse_fraction(frozen.get("leading"), f"q={q} leading")), f"q={q} leading coefficient drifted")
        profiles.append(profile)

    return {
        "schema": EXPECTED_SCHEMA,
        "status": "valid_adjacent_annihilator_formal_asymptotics",
        "leading_matching_power": "13/4",
        "slope_power": "3/4",
        "profiles": profiles,
        "maps_q_to_w_as_4_plus_q": True,
        "uses_root_or_production_data": False,
        "selects_winning_exponent": False,
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

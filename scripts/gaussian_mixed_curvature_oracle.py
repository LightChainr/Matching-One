#!/usr/bin/env python3
"""Formal norm-2/norm-5 mixed-curvature identities in exact arithmetic."""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Dict, Mapping


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis" / "gaussian_mixed_curvature_contract.json"

TERM_KEYS = (
    "constant",
    "logN",
    "logN_sq",
    "log2",
    "log5",
    "logN_log2",
    "logN_log5",
    "log2_sq",
    "log2_log5",
    "log5_sq",
)


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def quadratic_log_at(
    constant: Fraction,
    linear: Fraction,
    quadratic: Fraction,
    power2: int,
    power5: int,
) -> Dict[str, Fraction]:
    """Expand A+B log(sN)+C log(sN)^2 in a formal logarithm basis."""

    a = Fraction(power2)
    b = Fraction(power5)
    return {
        "constant": constant,
        "logN": linear,
        "logN_sq": quadratic,
        "log2": linear * a,
        "log5": linear * b,
        "logN_log2": 2 * quadratic * a,
        "logN_log5": 2 * quadratic * b,
        "log2_sq": quadratic * a * a,
        "log2_log5": 2 * quadratic * a * b,
        "log5_sq": quadratic * b * b,
    }


def add_scaled(
    target: Dict[str, Fraction], source: Mapping[str, Fraction], scale: int
) -> None:
    for key in TERM_KEYS:
        target[key] += scale * source[key]


def mixed_log_curvature(
    constant: Fraction, linear: Fraction, quadratic: Fraction
) -> Dict[str, Fraction]:
    """Apply F(10N)-F(5N)-F(2N)+F(N) formally."""

    result = {key: Fraction(0) for key in TERM_KEYS}
    for power2, power5, sign in ((1, 1, 1), (0, 1, -1), (1, 0, -1), (0, 0, 1)):
        add_scaled(
            result,
            quadratic_log_at(constant, linear, quadratic, power2, power5),
            sign,
        )
    return {key: value for key, value in result.items() if value}


def ordinary_power_mixed_response(beta: int) -> Fraction:
    """Return the normalized response for Y(N)=N^(-beta)."""

    if beta <= 0:
        raise ValueError("beta must be a positive integer")
    direct = (
        Fraction(1, 10**beta)
        - Fraction(1, 5**beta)
        - Fraction(1, 2**beta)
        + 1
    )
    factorized = (1 - Fraction(1, 2**beta)) * (1 - Fraction(1, 5**beta))
    if direct != factorized:
        raise AssertionError("mixed-power factorization failed")
    return direct


def serialize_terms(terms: Mapping[str, Fraction]) -> Dict[str, str]:
    return {key: fraction_text(value) for key, value in terms.items()}


def build_contract() -> Dict[str, object]:
    return {
        "schema": "matching-one/gaussian-mixed-curvature/v1",
        "status": "valid_exact_norm2_norm5_mixed_curvature",
        "parent_issue": "remain open",
        "operator": "F(10N)-F(5N)-F(2N)+F(N)",
        "constant_response": serialize_terms(
            mixed_log_curvature(Fraction(7), Fraction(0), Fraction(0))
        ),
        "rank2_affine_log_response": serialize_terms(
            mixed_log_curvature(Fraction(7), Fraction(11), Fraction(0))
        ),
        "unit_quadratic_log_response": serialize_terms(
            mixed_log_curvature(Fraction(0), Fraction(0), Fraction(1))
        ),
        "quadratic_rule": "2*C*log(2)*log(5)",
        "ordinary_power_responses": {
            str(beta): fraction_text(ordinary_power_mixed_response(beta))
            for beta in (1, 2, 3)
        },
        "ordinary_power_rule": "(1-2^(-beta))*(1-5^(-beta))",
        "uses_data": False,
        "identifies_a_jordan_or_physical_mechanism": False,
    }


def validate_contract(path: Path = CONTRACT) -> Dict[str, object]:
    with path.open(encoding="utf-8") as handle:
        frozen = json.load(handle)
    actual = build_contract()
    if frozen != actual:
        raise AssertionError("checked-in Gaussian mixed-curvature contract drifted")
    return actual


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=CONTRACT)
    args = parser.parse_args()
    print(json.dumps(validate_contract(args.contract), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

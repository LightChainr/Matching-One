#!/usr/bin/env python3
"""Exact reflection, median, and mode certificate for a threshold profile."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
from typing import Sequence

from threshold_histogram_profile import (
    density_coefficients,
    evaluate_polynomial,
    integrate_density,
    mixture_weights,
    parse_histogram,
    validate_contract,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = ROOT / "analysis" / "threshold_histogram_profile_contract.json"


def reflection_coefficients(coefficients: Sequence[Fraction]) -> list[Fraction]:
    """Return coefficients of p(1-x) in the same power basis."""
    values = [Fraction(value) for value in coefficients]
    reflected = [Fraction(0) for _ in values]
    for degree, coefficient in enumerate(values):
        for power in range(degree + 1):
            reflected[power] += (
                coefficient
                * math.comb(degree, power)
                * ((-1) ** power)
            )
    return reflected


def derivative_coefficients(coefficients: Sequence[Fraction]) -> list[Fraction]:
    derivative = [
        Fraction(degree) * Fraction(coefficient)
        for degree, coefficient in enumerate(coefficients)
        if degree
    ]
    while derivative and derivative[-1] == 0:
        derivative.pop()
    return derivative


def weights_are_reflection_symmetric(weights: Sequence[Fraction]) -> bool:
    values = tuple(Fraction(value) for value in weights)
    return values == tuple(reversed(values))


def density_is_reflection_symmetric(coefficients: Sequence[Fraction]) -> bool:
    values = [Fraction(value) for value in coefficients]
    return reflection_coefficients(values) == values


def cdf_has_complement_reflection(coefficients: Sequence[Fraction]) -> bool:
    values = [Fraction(value) for value in coefficients]
    reflected = reflection_coefficients(values)
    total = [left + right for left, right in zip(values, reflected)]
    return total == [Fraction(1)] + [Fraction(0) for _ in values[1:]]


def certify_unique_midpoint_mode(density: Sequence[Fraction]) -> dict[str, str | bool]:
    """Certify the strict mode using the exact linear derivative gate.

    This deliberately handles only a nonconstant quadratic density. Higher
    degree profiles require a separate exact root-isolation certificate.
    """
    values = [Fraction(value) for value in density]
    derivative = derivative_coefficients(values)
    if len(derivative) != 2:
        raise ValueError("unique-mode gate requires an exactly linear derivative")
    intercept, slope = derivative
    if not intercept > 0 or not slope < 0:
        raise ValueError("linear derivative does not change from positive to negative")
    root = -intercept / slope
    if root != Fraction(1, 2):
        raise ValueError("unique mode is not at the reflection midpoint")
    if evaluate_polynomial(values, Fraction(0)) <= 0 or evaluate_polynomial(values, Fraction(1)) <= 0:
        raise ValueError("density endpoints must be positive for a unique median certificate")
    return {
        "unique_mode_certified": True,
        "mode": "1/2",
        "derivative_power_coefficients": [str(value) for value in derivative],
        "derivative_sign": "positive on [0,1/2), negative on (1/2,1]",
    }


def build_artifact(contract_path: Path = DEFAULT_CONTRACT) -> dict:
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    validate_contract(contract)
    n = contract["N"]
    minus = parse_histogram(contract["K_minus_counts"], n, "K_minus")
    plus = parse_histogram(contract["K_plus_counts"], n, "K_plus")
    weights = mixture_weights(minus, plus, n)
    density = density_coefficients(weights)
    cdf = integrate_density(density)

    gates = {
        "rank_weights_reflection_symmetric": weights_are_reflection_symmetric(weights),
        "density_reflection_symmetric": density_is_reflection_symmetric(density),
        "cdf_complement_reflection": cdf_has_complement_reflection(cdf),
        "midpoint_cdf_is_half": evaluate_polynomial(cdf, Fraction(1, 2)) == Fraction(1, 2),
    }
    if not all(gates.values()):
        failed = ", ".join(name for name, passed in gates.items() if not passed)
        raise ValueError("reflection certificate failed: " + failed)
    mode = certify_unique_midpoint_mode(density)
    return {
        "schema": "matching-one/threshold-profile-symmetry-certificate/v1",
        "issue": 28,
        "data_class": "exact synthetic histogram only",
        "N": n,
        "mixture_weights_by_rank": [str(value) for value in weights],
        "density_power_coefficients": [str(value) for value in density],
        "cdf_power_coefficients": [str(value) for value in cdf],
        "reflection_gates": gates,
        "unique_median_certified": True,
        "median": "1/2",
        **mode,
        "boundary": (
            "Exact synthetic low-degree certificate only: no production histogram, empirical "
            "median/mode, bootstrap distance, tail fit, or universality claim."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rendered = json.dumps(build_artifact(args.contract), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


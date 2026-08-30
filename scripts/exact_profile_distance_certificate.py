#!/usr/bin/env python3
"""Exact L2 density and CDF distances for rational polynomial profiles."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Sequence

from exact_polynomial_root_certificate import fraction_text
from threshold_histogram_profile import (
    density_coefficients,
    integrate_density,
    mixture_weights,
    parse_histogram,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = ROOT / "analysis" / "threshold_histogram_profile_contract.json"


def pad(polynomial: Sequence[Fraction], size: int) -> list[Fraction]:
    values = [Fraction(value) for value in polynomial]
    return values + [Fraction(0) for _ in range(size - len(values))]


def subtract(first: Sequence[Fraction], second: Sequence[Fraction]) -> list[Fraction]:
    size = max(len(first), len(second))
    left, right = pad(first, size), pad(second, size)
    return [a - b for a, b in zip(left, right)]


def polynomial_inner_product(
    first: Sequence[Fraction], second: Sequence[Fraction]
) -> Fraction:
    return sum(
        (
            Fraction(left) * Fraction(right) / Fraction(i + j + 1)
            for i, left in enumerate(first)
            for j, right in enumerate(second)
        ),
        Fraction(0),
    )


def integral(polynomial: Sequence[Fraction]) -> Fraction:
    return sum(
        (Fraction(value) / Fraction(index + 1) for index, value in enumerate(polynomial)),
        Fraction(0),
    )


def validate_density(density: Sequence[Fraction], label: str) -> list[Fraction]:
    values = [Fraction(value) for value in density]
    if not values:
        raise ValueError(f"{label} density must not be empty")
    if integral(values) != 1:
        raise ValueError(f"{label} density must integrate exactly to one")
    return values


def exact_profile_distance(
    first_density: Sequence[Fraction], second_density: Sequence[Fraction]
) -> dict[str, Fraction]:
    first = validate_density(first_density, "first")
    second = validate_density(second_density, "second")
    density_difference = subtract(first, second)
    cdf_difference = subtract(integrate_density(first), integrate_density(second))
    density_l2 = polynomial_inner_product(density_difference, density_difference)
    cdf_l2 = polynomial_inner_product(cdf_difference, cdf_difference)
    if density_l2 < 0 or cdf_l2 < 0:
        raise ArithmeticError("squared profile distance must be nonnegative")
    return {
        "density_L2_squared": density_l2,
        "cdf_Cramer_von_Mises_squared": cdf_l2,
    }


def serialize_distance(distance: dict[str, Fraction]) -> dict[str, str]:
    return {name: fraction_text(value) for name, value in distance.items()}


def build_artifact(contract_path: Path = DEFAULT_CONTRACT) -> dict:
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    n = contract["N"]
    minus = parse_histogram(contract["K_minus_counts"], n, "K_minus")
    plus = parse_histogram(contract["K_plus_counts"], n, "K_plus")
    frozen = density_coefficients(mixture_weights(minus, plus, n))
    uniform = [Fraction(1)]
    beta22 = [Fraction(0), Fraction(6), Fraction(-6)]
    distance = exact_profile_distance(frozen, uniform)
    reverse = exact_profile_distance(uniform, frozen)
    identity = exact_profile_distance(frozen, frozen)
    if distance != reverse or any(identity.values()):
        raise ValueError("profile distance symmetry/identity gate failed")

    frozen_difference = subtract(frozen, uniform)
    beta_difference = subtract(beta22, uniform)
    mixture_size = max(len(frozen), len(uniform), len(beta22))
    mixture_identity = pad(frozen, mixture_size) == [
        (left + right) / 2
        for left, right in zip(
            pad(uniform, mixture_size), pad(beta22, mixture_size)
        )
    ]
    if not mixture_identity:
        raise ValueError("frozen uniform/Beta(2,2) mixture identity drifted")
    gram = [
        [
            polynomial_inner_product(frozen_difference, frozen_difference),
            polynomial_inner_product(frozen_difference, beta_difference),
        ],
        [
            polynomial_inner_product(beta_difference, frozen_difference),
            polynomial_inner_product(beta_difference, beta_difference),
        ],
    ]
    gram_determinant = gram[0][0] * gram[1][1] - gram[0][1] * gram[1][0]
    if gram[0][1] != gram[1][0] or gram_determinant < 0:
        raise ValueError("profile-difference Gram gate failed")
    return {
        "schema": "matching-one/exact-profile-distance-certificate/v1",
        "issue": 28,
        "data_class": "exact synthetic rational polynomial profiles",
        "frozen_density_power_coefficients": [fraction_text(value) for value in frozen],
        "reference_density": "uniform_on_[0,1]",
        "reference_density_power_coefficients": ["1"],
        "distance": serialize_distance(distance),
        "symmetry_certified": True,
        "identity_zero_certified": True,
        "exact_affine_mixture_identity": "frozen = (uniform + Beta(2,2))/2",
        "control_gram": {
            "basis": ["frozen-minus-uniform", "beta22-minus-uniform"],
            "matrix": [[fraction_text(value) for value in row] for row in gram],
            "determinant": fraction_text(gram_determinant),
            "positive_semidefinite_certified": (
                gram[0][0] >= 0 and gram[1][1] >= 0 and gram_determinant >= 0
            ),
        },
        "boundary": (
            "Exact synthetic polynomial distances only: no production histogram, bootstrap "
            "uncertainty, cross-model universality, or tail claim."
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

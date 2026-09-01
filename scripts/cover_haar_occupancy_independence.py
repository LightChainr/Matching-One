#!/usr/bin/env python3
"""Exact occupancy independence for additive Haar cover couplings."""

from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import combinations
import json
from math import factorial
from pathlib import Path
from typing import Dict, List, Sequence


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis" / "cover_haar_occupancy_independence_contract.json"


def validate_inputs(degree: int, probability: Fraction) -> None:
    if degree < 2:
        raise ValueError("Haar independence requires cover degree at least two")
    if not 0 < probability < 1:
        raise ValueError("occupancy probability must lie strictly between zero and one")


def bounded_box_cdf(threshold: Fraction, caps: Sequence[Fraction]) -> Fraction:
    """Lebesgue volume of 0<=x_i<=caps_i and sum x_i<threshold."""

    dimension = len(caps)
    total = Fraction()
    for size in range(dimension + 1):
        for subset in combinations(range(dimension), size):
            remainder = threshold - sum((caps[index] for index in subset), Fraction())
            if remainder > 0:
                total += (-1) ** size * remainder**dimension
    return total / factorial(dimension)


def fractional_window_volume(caps: Sequence[Fraction], width: Fraction) -> Fraction:
    """Volume where the fractional part of the coordinate sum lies in [0,width)."""

    if not 0 <= width <= 1:
        raise ValueError("fractional window width must lie in [0,1]")
    dimension = len(caps)
    return sum(
        (
            bounded_box_cdf(Fraction(integer) + width, caps)
            - bounded_box_cdf(Fraction(integer), caps)
            for integer in range(dimension)
        ),
        Fraction(),
    )


def joint_parent_child(
    degree: int, probability: Fraction, *, antithetic: bool = False
) -> Fraction:
    """P(parent occupied and child coordinate zero occupied)."""

    validate_inputs(degree, probability)
    caps = [probability] + [Fraction(1)] * (degree - 1)
    if antithetic:
        return probability - fractional_window_volume(caps, 1 - probability)
    return fractional_window_volume(caps, probability)


def occupancy_covariance(
    degree: int, probability: Fraction, *, antithetic: bool = False
) -> Fraction:
    """Cov(parent occupancy, the child fiber occupancy mean)."""

    joint = joint_parent_child(degree, probability, antithetic=antithetic)
    return joint - probability * probability


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def frozen_case(degree: int, probability: Fraction) -> Dict[str, object]:
    additive_joint = joint_parent_child(degree, probability)
    antithetic_joint = joint_parent_child(degree, probability, antithetic=True)
    return {
        "cover_degree": degree,
        "occupancy_probability": fraction_text(probability),
        "product_of_marginals": fraction_text(probability * probability),
        "h0_joint_parent_child": fraction_text(additive_joint),
        "h0_covariance_with_fiber_mean": fraction_text(
            occupancy_covariance(degree, probability)
        ),
        "h1_joint_parent_child": fraction_text(antithetic_joint),
        "h1_covariance_with_fiber_mean": fraction_text(
            occupancy_covariance(degree, probability, antithetic=True)
        ),
    }


def build_contract() -> Dict[str, object]:
    return {
        "schema": "matching-one/cover-haar-occupancy-independence/v1",
        "status": "valid_exact_haar_independence_certificate",
        "parent_issue": "remain open",
        "degree_domain": "Q>=2",
        "h0_parent_uniform": "V=frac(sum_m U_m)",
        "h1_parent_uniform": "1-V",
        "theorem": "V is independent of every U_i because the modulo-one sum of the other coordinates is Haar uniform",
        "occupancy_conclusion": "both H0 and H1 parent occupancies have zero covariance with the child fiber occupancy mean",
        "frozen_cases": [
            frozen_case(degree, probability)
            for degree in (2, 5)
            for probability in (Fraction(2, 5), Fraction(3, 5))
        ],
        "exact_rational_volume_oracle": True,
        "claim_boundary": {
            "proves_nonlinear_topological_covariance": False,
            "covers_threshold_rank_permutations": False,
            "includes_wall_time_model": False,
            "proves_residual_variance_gain": False,
            "makes_production_recommendation": False,
        },
    }


def validate_contract(path: Path = CONTRACT) -> Dict[str, object]:
    with path.open(encoding="utf-8") as handle:
        frozen = json.load(handle)
    actual = build_contract()
    if frozen != actual:
        raise AssertionError("checked-in Haar occupancy contract drifted")
    return actual


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=CONTRACT)
    args = parser.parse_args()
    print(json.dumps(validate_contract(args.contract), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Exact occupancy covariance for section and antithetic cover couplings."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Dict


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis" / "cover_section_occupancy_covariance_contract.json"


def validate_inputs(degree: int, probability: Fraction) -> None:
    if degree < 1:
        raise ValueError("cover degree must be positive")
    if not 0 < probability < 1:
        raise ValueError("occupancy probability must lie strictly between zero and one")


def section_covariance(degree: int, probability: Fraction) -> Fraction:
    """Cov(1[U0<p], Q^-1 sum_m 1[Um<p])."""

    validate_inputs(degree, probability)
    return probability * (1 - probability) / degree


def antithetic_covariance(degree: int, probability: Fraction) -> Fraction:
    """Cov(1[1-U0<p], Q^-1 sum_m 1[Um<p])."""

    validate_inputs(degree, probability)
    overlap = max(Fraction(), 2 * probability - 1)
    return (overlap - probability * probability) / degree


def squared_correlation(covariance: Fraction, degree: int, probability: Fraction) -> Fraction:
    validate_inputs(degree, probability)
    bernoulli_variance = probability * (1 - probability)
    fiber_mean_variance = bernoulli_variance / degree
    return covariance * covariance / (bernoulli_variance * fiber_mean_variance)


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def coupling_case(degree: int, probability: Fraction) -> Dict[str, object]:
    direct = section_covariance(degree, probability)
    antithetic = antithetic_covariance(degree, probability)
    direct_r2 = squared_correlation(direct, degree, probability)
    antithetic_r2 = squared_correlation(antithetic, degree, probability)
    return {
        "cover_degree": degree,
        "occupancy_probability": fraction_text(probability),
        "section_covariance": fraction_text(direct),
        "section_squared_correlation": fraction_text(direct_r2),
        "antithetic_covariance": fraction_text(antithetic),
        "antithetic_squared_correlation": fraction_text(antithetic_r2),
        "antithetic_is_strictly_smaller_in_magnitude": abs(antithetic) < direct,
    }


def build_contract() -> Dict[str, object]:
    probabilities = (Fraction(2, 5), Fraction(3, 5))
    degrees = (2, 5)
    return {
        "schema": "matching-one/cover-section-occupancy-covariance/v1",
        "status": "valid_exact_bernoulli_fiber_certificate",
        "parent_issue": "remain open",
        "child_statistic": "fiber occupancy mean Q^-1 sum_m 1[U_m<p]",
        "section_parent": "1[U_0<p]",
        "antithetic_parent": "1[1-U_0<p]",
        "exact_formulas": {
            "section_covariance": "p(1-p)/Q",
            "section_squared_correlation": "1/Q",
            "antithetic_covariance": "-min(p,1-p)^2/Q",
            "antithetic_squared_correlation": "min(p,1-p)^4/(Q p^2(1-p)^2)",
        },
        "comparison_scope": "the declared section and antithetic-section pair only",
        "section_squared_correlation_ceiling_within_pair": "1/Q",
        "frozen_cases": [
            coupling_case(degree, probability)
            for degree in degrees
            for probability in probabilities
        ],
        "claim_boundary": {
            "covers_h0_h1_additive_couplings": False,
            "proves_topological_observable_covariance": False,
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
        raise AssertionError("checked-in cover-section occupancy contract drifted")
    return actual


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=CONTRACT)
    args = parser.parse_args()
    print(json.dumps(validate_contract(args.contract), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

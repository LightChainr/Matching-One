#!/usr/bin/env python3
"""Reconstruct an exact finite threshold profile from integer rank histograms."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

from exact_jet_algebra import canonical_fraction


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = ROOT / "analysis" / "threshold_histogram_profile_contract.json"
EXPECTED_SCHEMA = "matching-one/threshold-histogram-profile/v1"
FORBIDDEN_KEYS = {
    "bootstrap",
    "empirical_kappa",
    "heldout_score",
    "production_histogram",
    "source_data",
    "tail_fit",
    "wasserstein_distance",
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


def parse_histogram(raw: Mapping[str, Any], n: int, label: str) -> dict[int, int]:
    _require(isinstance(raw, dict), f"{label} must be an object")
    output: dict[int, int] = {}
    for raw_rank, count in raw.items():
        _require(isinstance(raw_rank, str) and raw_rank.isdigit(), f"{label} ranks must be decimal strings")
        rank = int(raw_rank)
        _require(str(rank) == raw_rank, f"{label} rank keys must be canonical")
        _require(1 <= rank <= n, f"{label} rank is out of range")
        _require(isinstance(count, int) and not isinstance(count, bool), f"{label} counts must be integers")
        _require(count >= 0, f"{label} counts must be nonnegative")
        output[rank] = count
    _require(sum(output.values()) > 0, f"{label} must contain samples")
    return output


def mixture_weights(
    minus_counts: Mapping[int, int], plus_counts: Mapping[int, int], n: int
) -> tuple[Fraction, ...]:
    _require(n >= 1, "N must be positive")
    for label, counts in (("K_minus", minus_counts), ("K_plus", plus_counts)):
        _require(bool(counts), f"{label} must contain samples")
        for rank, count in counts.items():
            _require(isinstance(rank, int) and 1 <= rank <= n, f"{label} rank is out of range")
            _require(
                isinstance(count, int) and not isinstance(count, bool) and count >= 0,
                f"{label} counts must be nonnegative integers",
            )
    minus_samples = sum(minus_counts.values())
    plus_samples = sum(plus_counts.values())
    _require(minus_samples > 0, "histograms must contain samples")
    _require(minus_samples == plus_samples, "K_minus and K_plus sample counts must match")
    denominator = 2 * minus_samples
    return tuple(
        Fraction(minus_counts.get(rank, 0) + plus_counts.get(rank, 0), denominator)
        for rank in range(1, n + 1)
    )


def beta_density_coefficients(n: int, rank: int) -> list[Fraction]:
    _require(n >= 1, "N must be positive")
    _require(1 <= rank <= n, "rank is out of range")
    coefficients = [Fraction(0) for _ in range(n)]
    normalization = n * math.comb(n - 1, rank - 1)
    for offset in range(n - rank + 1):
        degree = rank - 1 + offset
        coefficients[degree] += (
            normalization * math.comb(n - rank, offset) * ((-1) ** offset)
        )
    return coefficients


def density_coefficients(weights: Sequence[Fraction]) -> list[Fraction]:
    n = len(weights)
    _require(n >= 1, "mixture must contain at least one rank")
    _require(all(weight >= 0 for weight in weights), "mixture weights must be nonnegative")
    _require(sum(weights) == 1, "mixture weights must sum to one")
    output = [Fraction(0) for _ in range(n)]
    for rank, weight in enumerate(weights, start=1):
        beta = beta_density_coefficients(n, rank)
        output = [left + weight * right for left, right in zip(output, beta)]
    return output


def integrate_density(coefficients: Sequence[Fraction]) -> list[Fraction]:
    return [Fraction(0)] + [coefficient / (degree + 1) for degree, coefficient in enumerate(coefficients)]


def evaluate_polynomial(coefficients: Sequence[Fraction], value: Fraction) -> Fraction:
    output = Fraction(0)
    for coefficient in reversed(coefficients):
        output = output * value + coefficient
    return output


def beta_raw_moment(n: int, rank: int, order: int) -> Fraction:
    _require(order >= 0, "moment order must be nonnegative")
    output = Fraction(1)
    for offset in range(order):
        output *= Fraction(rank + offset, n + 1 + offset)
    return output


def raw_moments(weights: Sequence[Fraction], maximum_order: int) -> list[Fraction]:
    n = len(weights)
    _require(maximum_order >= 0, "maximum moment order must be nonnegative")
    return [
        sum(
            (weight * beta_raw_moment(n, rank, order) for rank, weight in enumerate(weights, start=1)),
            Fraction(0),
        )
        for order in range(maximum_order + 1)
    ]


def parse_fraction_vector(values: Any, length: int, label: str) -> list[Fraction]:
    _require(isinstance(values, list) and len(values) == length, f"{label} must have length {length}")
    return [canonical_fraction(value, f"{label}[{index}]") for index, value in enumerate(values)]


def validate_contract(contract: Mapping[str, Any]) -> dict[str, Any]:
    _walk_forbidden(contract)
    _require(contract.get("schema") == EXPECTED_SCHEMA, "schema drifted")
    _require(contract.get("issue") == 28, "issue must be 28")
    _require(contract.get("status") == "exact_synthetic_histogram_only", "status drifted")
    n = contract.get("N")
    _require(isinstance(n, int) and n >= 1, "N must be positive")
    minus = parse_histogram(contract.get("K_minus_counts", {}), n, "K_minus")
    plus = parse_histogram(contract.get("K_plus_counts", {}), n, "K_plus")
    weights = mixture_weights(minus, plus, n)
    density = density_coefficients(weights)
    cdf = integrate_density(density)
    moments = raw_moments(weights, 6)

    expected_density = parse_fraction_vector(
        contract.get("expected_density_power_coefficients"), n, "expected density"
    )
    expected_cdf = parse_fraction_vector(
        contract.get("expected_cdf_power_coefficients"), n + 1, "expected CDF"
    )
    expected_moments = parse_fraction_vector(
        contract.get("expected_raw_moments_0_through_6"), 7, "expected moments"
    )
    expected_mean = canonical_fraction(contract.get("expected_mean_from_ranks"), "expected mean")

    _require(density == expected_density, "density coefficients drifted")
    _require(cdf == expected_cdf, "CDF coefficients drifted")
    _require(moments == expected_moments, "raw moments drifted")
    _require(evaluate_polynomial(cdf, Fraction(0)) == 0, "CDF does not start at zero")
    _require(evaluate_polynomial(cdf, Fraction(1)) == 1, "CDF does not end at one")
    _require(evaluate_polynomial(cdf, Fraction(1, 2)) == Fraction(1, 2), "synthetic median drifted")
    _require(evaluate_polynomial(density, Fraction(1, 2)) == Fraction(5, 4), "synthetic peak drifted")

    mean_from_ranks = sum(
        (weight * rank for rank, weight in enumerate(weights, start=1)), Fraction(0)
    ) / (n + 1)
    _require(mean_from_ranks == moments[1] == expected_mean, "rank mean identity failed")

    return {
        "schema": EXPECTED_SCHEMA,
        "status": "valid_exact_threshold_profile",
        "N": n,
        "samples_per_histogram": sum(minus.values()),
        "mixture_weights_by_rank": [str(weight) for weight in weights],
        "density_power_coefficients": [str(value) for value in density],
        "cdf_power_coefficients": [str(value) for value in cdf],
        "raw_moments_0_through_6": [str(value) for value in moments],
        "mean_from_ranks": str(mean_from_ranks),
        "cdf_normalized_exactly": True,
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

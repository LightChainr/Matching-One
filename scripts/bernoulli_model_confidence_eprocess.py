#!/usr/bin/env python3
"""Exact finite-horizon e-process confidence set for two Bernoulli models.

Each declared model is tested by the likelihood ratio of the other model
against it.  Inverting the two e-values at threshold 1/alpha gives an
anytime-valid confidence sequence over the declared model family.  The oracle
enumerates all paths, verifies fixed-time and bounded-stopping expectations,
and computes the exact probability that the true declared model is ever
excluded.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Any, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = ROOT / "analysis" / "bernoulli_model_confidence_contract.json"
EXPECTED_SCHEMA = "matching-one/bernoulli-model-confidence/v1"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def parse_fraction(value: Any, label: str) -> Fraction:
    _require(isinstance(value, (str, int)), f"{label} must be an exact rational")
    try:
        return Fraction(value)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError(f"{label} must be an exact rational") from exc


def path_probability(path: Sequence[int], p: Fraction) -> Fraction:
    _require(all(bit in (0, 1) for bit in path), "path bits must be zero or one")
    occupied = sum(path)
    return p**occupied * (1 - p) ** (len(path) - occupied)


def likelihood_ratio(path: Sequence[int], null: Fraction, alternative: Fraction) -> Fraction:
    _require(0 < null < 1 and 0 < alternative < 1, "model probabilities must be interior")
    _require(null != alternative, "null and alternative models must differ")
    occupied = sum(path)
    vacant = len(path) - occupied
    return (
        (alternative / null) ** occupied
        * ((1 - alternative) / (1 - null)) ** vacant
    )


def all_paths(length: int) -> tuple[tuple[int, ...], ...]:
    _require(length >= 0, "path length must be nonnegative")
    return tuple(
        tuple((mask >> index) & 1 for index in range(length))
        for mask in range(1 << length)
    )


def confidence_sets(
    path: Sequence[int], models: Sequence[Fraction], threshold: Fraction
) -> tuple[tuple[int, ...], ...]:
    _require(len(models) == 2, "this bounded oracle requires exactly two models")
    _require(threshold > 1, "e-value threshold must exceed one")
    result = []
    for time in range(len(path) + 1):
        prefix = path[:time]
        retained = tuple(
            null_index
            for null_index in range(2)
            if likelihood_ratio(
                prefix, models[null_index], models[1 - null_index]
            )
            < threshold
        )
        result.append(retained)
    return tuple(result)


def model_certificate(
    models: Sequence[Fraction], null_index: int, horizon: int, alpha: Fraction
) -> dict[str, Any]:
    _require(len(models) == 2, "this bounded oracle requires exactly two models")
    _require(null_index in (0, 1), "null index must be zero or one")
    _require(horizon >= 1, "horizon must be positive")
    _require(0 < alpha < 1, "alpha must lie between zero and one")
    null = models[null_index]
    alternative = models[1 - null_index]
    threshold = 1 / alpha

    fixed_time_expectations: list[Fraction] = []
    for time in range(horizon + 1):
        expectation = sum(
            path_probability(path, null)
            * likelihood_ratio(path, null, alternative)
            for path in all_paths(time)
        )
        fixed_time_expectations.append(expectation)
    _require(all(value == 1 for value in fixed_time_expectations), "fixed-time e-mean drifted")

    crossing_probability = Fraction(0)
    stopped_expectation = Fraction(0)
    crossing_time_mass = [Fraction(0)] * (horizon + 1)
    for path in all_paths(horizon):
        probability = path_probability(path, null)
        stop = horizon
        crossed = False
        for time in range(1, horizon + 1):
            if likelihood_ratio(path[:time], null, alternative) >= threshold:
                stop = time
                crossed = True
                break
        stopped_value = likelihood_ratio(path[:stop], null, alternative)
        stopped_expectation += probability * stopped_value
        if crossed:
            crossing_probability += probability
            crossing_time_mass[stop] += probability

        sets = confidence_sets(path, models, threshold)
        ever_excluded = any(null_index not in retained for retained in sets)
        _require(ever_excluded == crossed, "confidence-set inversion disagreed with e-value crossing")

    _require(stopped_expectation == 1, "bounded optional-stopping expectation drifted")
    _require(crossing_probability <= alpha, "Ville bound failed")
    return {
        "null_index": null_index,
        "null_p": str(null),
        "alternative_p": str(alternative),
        "fixed_time_expectations": [str(value) for value in fixed_time_expectations],
        "first_crossing_mass_by_time": [str(value) for value in crossing_time_mass],
        "ever_excluded_probability": str(crossing_probability),
        "alpha": str(alpha),
        "ville_bound_holds": True,
        "stopped_e_expectation": str(stopped_expectation),
        "bounded_optional_stopping_holds": True,
    }


def validate_contract(contract: Mapping[str, Any]) -> dict[str, Any]:
    _require(contract.get("schema") == EXPECTED_SCHEMA, "schema drifted")
    _require(contract.get("issue") == 126, "issue must be 126")
    _require(contract.get("status") == "finite_declared_model_family_only", "status drifted")
    raw_models = contract.get("models")
    _require(isinstance(raw_models, list) and len(raw_models) == 2, "models must have length two")
    models = [parse_fraction(value, f"models[{index}]") for index, value in enumerate(raw_models)]
    _require(models[0] < models[1], "models must be strictly increasing")
    horizon = contract.get("horizon")
    _require(isinstance(horizon, int), "horizon must be an integer")
    alpha = parse_fraction(contract.get("alpha"), "alpha")
    expected_crossing = parse_fraction(
        contract.get("expected_ever_excluded_probability"),
        "expected ever-excluded probability",
    )

    certificates = [model_certificate(models, index, horizon, alpha) for index in range(2)]
    for row in certificates:
        _require(
            parse_fraction(row["ever_excluded_probability"], "computed crossing probability")
            == expected_crossing,
            "ever-excluded probability drifted",
        )
    return {
        "schema": EXPECTED_SCHEMA,
        "status": "valid_finite_bernoulli_model_confidence_sequence",
        "models": [str(value) for value in models],
        "horizon": horizon,
        "threshold": str(1 / alpha),
        "certificates": certificates,
        "simultaneous_true_model_coverage_at_least": str(1 - alpha),
        "fixed_time_e_means_equal_one": True,
        "bounded_optional_stopping_exact": True,
        "off_model_coverage_defined": False,
        "permits_adaptive_model_generation_from_same_data": False,
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

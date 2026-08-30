#!/usr/bin/env python3
"""Exact logit-parameter and Russo identities for a finite monotone event.

For independent Bernoulli bits, let K be the occupation count and A a monotone
event with probability F(p).  With eta=log(p/(1-p)), finite differentiation
gives

    dF/deta = Cov(1_A, K) = p(1-p) F'(p)

and therefore

    d logit(F)/deta = E[K | A] - E[K | not A].

The bounded oracle uses the threshold event K >= r and verifies the identities,
including the Russo pivotal sum, with Fraction arithmetic.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Any, Mapping


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = ROOT / "analysis" / "bernoulli_logit_hazard_contract.json"
EXPECTED_SCHEMA = "matching-one/bernoulli-logit-hazard/v1"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def parse_fraction(value: Any, label: str) -> Fraction:
    _require(isinstance(value, (str, int)), f"{label} must be an exact rational")
    try:
        result = Fraction(value)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError(f"{label} must be an exact rational") from exc
    return result


def bernoulli_weight(occupied: int, total: int, p: Fraction) -> Fraction:
    return p**occupied * (1 - p) ** (total - occupied)


def threshold_oracle(total: int, threshold: int, p: Fraction) -> dict[str, Fraction]:
    _require(total >= 1, "total must be positive")
    _require(1 <= threshold <= total, "threshold must lie between one and total")
    _require(0 < p < 1, "p must lie strictly between zero and one")

    event_probability = Fraction(0)
    occupied_event_moment = Fraction(0)
    occupied_moment = Fraction(0)
    derivative_p = Fraction(0)
    for mask in range(1 << total):
        occupied = bin(mask).count("1")
        weight = bernoulli_weight(occupied, total, p)
        occupied_moment += occupied * weight
        if occupied >= threshold:
            event_probability += weight
            occupied_event_moment += occupied * weight
            derivative_p += weight * (
                Fraction(occupied, 1) / p
                - Fraction(total - occupied, 1) / (1 - p)
            )

    complement_probability = 1 - event_probability
    _require(0 < event_probability < 1, "event must be nontrivial")
    conditional_event = occupied_event_moment / event_probability
    conditional_complement = (
        occupied_moment - occupied_event_moment
    ) / complement_probability
    score_covariance = occupied_event_moment - event_probability * occupied_moment
    derivative_eta = p * (1 - p) * derivative_p
    logit_derivative = derivative_eta / (
        event_probability * complement_probability
    )

    pivotal_sum = Fraction(0)
    for _edge in range(total):
        for other_mask in range(1 << (total - 1)):
            other_occupied = bin(other_mask).count("1")
            if other_occupied == threshold - 1:
                pivotal_sum += bernoulli_weight(other_occupied, total - 1, p)

    _require(occupied_moment == total * p, "Bernoulli score mean drifted")
    _require(derivative_p == pivotal_sum, "Russo pivotal identity failed")
    _require(derivative_eta == score_covariance, "natural-parameter score identity failed")
    _require(
        logit_derivative == conditional_event - conditional_complement,
        "logit derivative did not equal the conditional occupation gap",
    )
    return {
        "F": event_probability,
        "dF_dp": derivative_p,
        "russo_pivotal_sum": pivotal_sum,
        "dF_deta": derivative_eta,
        "score_covariance": score_covariance,
        "d_logitF_deta": logit_derivative,
        "E_K_given_A": conditional_event,
        "E_K_given_not_A": conditional_complement,
        "conditional_occupation_gap": conditional_event - conditional_complement,
    }


def string_payload(row: Mapping[str, Fraction]) -> dict[str, str]:
    return {name: str(value) for name, value in row.items()}


def validate_contract(contract: Mapping[str, Any]) -> dict[str, Any]:
    _require(contract.get("schema") == EXPECTED_SCHEMA, "schema drifted")
    _require(contract.get("issue") == 100, "issue must be 100")
    _require(contract.get("status") == "exact_finite_event_only", "status drifted")
    total = contract.get("bits")
    threshold = contract.get("threshold")
    _require(isinstance(total, int), "bits must be an integer")
    _require(isinstance(threshold, int), "threshold must be an integer")
    raw_points = contract.get("p_points")
    _require(isinstance(raw_points, list) and raw_points, "p_points must be nonempty")
    points = [parse_fraction(value, f"p_points[{index}]") for index, value in enumerate(raw_points)]
    _require(points == sorted(points) and len(set(points)) == len(points), "p_points must be sorted and distinct")

    expected = contract.get("expected")
    _require(isinstance(expected, Mapping), "expected values must be a mapping")
    rows: dict[str, Any] = {}
    for p in points:
        row = threshold_oracle(total, threshold, p)
        key = str(p)
        frozen = expected.get(key)
        _require(isinstance(frozen, Mapping), f"missing expected row for p={key}")
        for field in ("F", "dF_dp", "dF_deta", "d_logitF_deta"):
            declared = parse_fraction(frozen.get(field), f"expected {key} {field}")
            _require(row[field] == declared, f"p={key} {field} drifted")
        rows[key] = string_payload(row)

    return {
        "schema": EXPECTED_SCHEMA,
        "status": "valid_exact_bernoulli_logit_hazard_oracle",
        "bits": total,
        "threshold": threshold,
        "points": rows,
        "russo_identity_exact": True,
        "natural_parameter_score_identity_exact": True,
        "logit_derivative_is_conditional_occupation_gap": True,
        "contains_monte_carlo": False,
        "identifies_four_arm_exponent": False,
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

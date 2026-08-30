#!/usr/bin/env python3
"""Exact finite Bernoulli oracle for predictable arm-selection e-processes."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

from exact_jet_algebra import canonical_fraction


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = ROOT / "analysis" / "predictable_arm_eprocess_contract.json"
EXPECTED_SCHEMA = "matching-one/predictable-arm-eprocess/v1"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


@dataclass(frozen=True)
class Arm:
    name: str
    null_success: Fraction
    alternative_success: Fraction

    def probability(self, outcome: int, *, alternative: bool = False) -> Fraction:
        _require(outcome in (0, 1), "outcome must be zero or one")
        success = self.alternative_success if alternative else self.null_success
        return success if outcome else 1 - success

    def likelihood_ratio(self, outcome: int) -> Fraction:
        return self.probability(outcome, alternative=True) / self.probability(outcome)


def parse_arms(raw: Any) -> dict[str, Arm]:
    _require(isinstance(raw, Mapping) and set(raw) == {"A", "B"}, "arms must be exactly A and B")
    arms: dict[str, Arm] = {}
    for name in ("A", "B"):
        row = raw[name]
        _require(isinstance(row, Mapping), f"arm {name} must be a mapping")
        null_success = canonical_fraction(row.get("null_success"), f"arm {name} null success")
        alternative_success = canonical_fraction(
            row.get("alternative_success"), f"arm {name} alternative success"
        )
        _require(0 < null_success < 1, f"arm {name} null probability must be interior")
        _require(0 < alternative_success < 1, f"arm {name} alternative probability must be interior")
        _require(null_success != alternative_success, f"arm {name} hypotheses must differ")
        arms[name] = Arm(name, null_success, alternative_success)
    return arms


def arm_for_history(history: Sequence[int], policy: Mapping[str, Any]) -> str:
    if not history:
        arm = policy.get("initial_arm")
    elif history[-1] == 0:
        arm = policy.get("after_failure")
    else:
        arm = policy.get("after_success")
    _require(arm in {"A", "B"}, "policy selected an unknown arm")
    return str(arm)


def path_values(
    outcomes: Sequence[int], arms: Mapping[str, Arm], policy: Mapping[str, Any]
) -> tuple[Fraction, Fraction, tuple[str, ...]]:
    null_probability = Fraction(1)
    e_value = Fraction(1)
    selected: list[str] = []
    history: list[int] = []
    for outcome in outcomes:
        arm_name = arm_for_history(history, policy)
        arm = arms[arm_name]
        null_probability *= arm.probability(outcome)
        e_value *= arm.likelihood_ratio(outcome)
        selected.append(arm_name)
        history.append(outcome)
    return null_probability, e_value, tuple(selected)


def all_binary_sequences(length: int) -> list[tuple[int, ...]]:
    _require(length >= 0, "sequence length must be nonnegative")
    return [tuple((mask >> index) & 1 for index in range(length)) for mask in range(1 << length)]


def terminal_enumeration(
    horizon: int, arms: Mapping[str, Arm], policy: Mapping[str, Any]
) -> dict[str, Any]:
    rows = []
    total_probability = Fraction(0)
    mean_e_value = Fraction(0)
    for outcomes in all_binary_sequences(horizon):
        probability, e_value, selected = path_values(outcomes, arms, policy)
        total_probability += probability
        mean_e_value += probability * e_value
        rows.append(
            {
                "outcomes": "".join(str(value) for value in outcomes),
                "arms": "".join(selected),
                "null_probability": str(probability),
                "e_value": str(e_value),
            }
        )
    return {
        "paths": rows,
        "path_count": len(rows),
        "null_probability_sum": str(total_probability),
        "null_mean_e_value": str(mean_e_value),
    }


def stopped_enumeration(
    horizon: int,
    threshold: Fraction,
    arms: Mapping[str, Arm],
    policy: Mapping[str, Any],
) -> dict[str, Any]:
    _require(threshold > 1, "stopping threshold must exceed one")
    leaves: list[dict[str, str | int | bool]] = []

    def visit(
        history: tuple[int, ...], probability: Fraction, e_value: Fraction
    ) -> None:
        stopped = e_value >= threshold
        if stopped or len(history) == horizon:
            _, _, selected = path_values(history, arms, policy)
            leaves.append(
                {
                    "outcomes": "".join(str(value) for value in history),
                    "arms": "".join(selected),
                    "time": len(history),
                    "threshold_hit": stopped,
                    "null_probability": str(probability),
                    "e_value": str(e_value),
                }
            )
            return
        arm = arms[arm_for_history(history, policy)]
        for outcome in (0, 1):
            visit(
                history + (outcome,),
                probability * arm.probability(outcome),
                e_value * arm.likelihood_ratio(outcome),
            )

    visit((), Fraction(1), Fraction(1))
    probability_sum = sum((Fraction(str(row["null_probability"])) for row in leaves), Fraction(0))
    stopped_mean = sum(
        (
            Fraction(str(row["null_probability"])) * Fraction(str(row["e_value"]))
            for row in leaves
        ),
        Fraction(0),
    )
    return {
        "leaves": leaves,
        "leaf_count": len(leaves),
        "null_probability_sum": str(probability_sum),
        "null_mean_stopped_e_value": str(stopped_mean),
    }


def conditional_mean_checks(arms: Mapping[str, Arm]) -> dict[str, str]:
    checks = {}
    for name, arm in arms.items():
        mean = sum(
            (arm.probability(outcome) * arm.likelihood_ratio(outcome) for outcome in (0, 1)),
            Fraction(0),
        )
        checks[name] = str(mean)
    return checks


def peek_both_expected_factor(arms: Mapping[str, Arm]) -> Fraction:
    """Invalid control: observe both current outcomes, then report the larger LR."""

    expected = Fraction(0)
    for outcome_a in (0, 1):
        for outcome_b in (0, 1):
            probability = arms["A"].probability(outcome_a) * arms["B"].probability(outcome_b)
            reported = max(
                arms["A"].likelihood_ratio(outcome_a),
                arms["B"].likelihood_ratio(outcome_b),
            )
            expected += probability * reported
    return expected


def validate_contract(contract: Mapping[str, Any]) -> dict[str, Any]:
    _require(contract.get("schema") == EXPECTED_SCHEMA, "schema drifted")
    _require(contract.get("issue") == 126, "issue must be 126")
    _require(contract.get("status") == "exact_finite_bernoulli_only", "status drifted")
    horizon = contract.get("horizon")
    _require(isinstance(horizon, int) and 1 <= horizon <= 12, "horizon must be an integer in [1,12]")
    threshold = canonical_fraction(contract.get("stopping_threshold"), "stopping threshold")
    arms = parse_arms(contract.get("arms"))
    policy = contract.get("predictable_policy")
    _require(isinstance(policy, Mapping), "predictable policy must be a mapping")

    conditional = conditional_mean_checks(arms)
    _require(set(conditional.values()) == {"1"}, "conditional LR mean is not one")
    terminal = terminal_enumeration(horizon, arms, policy)
    _require(terminal["null_probability_sum"] == "1", "terminal path probabilities do not sum to one")
    _require(terminal["null_mean_e_value"] == "1", "terminal e-value mean is not one")
    stopped = stopped_enumeration(horizon, threshold, arms, policy)
    _require(stopped["null_probability_sum"] == "1", "stopped path probabilities do not sum to one")
    _require(stopped["null_mean_stopped_e_value"] == "1", "bounded stopping changed the e-value mean")

    peek_factor = peek_both_expected_factor(arms)
    declared_peek = canonical_fraction(
        contract.get("peek_both_negative_control_expected_factor"), "peek-both expected factor"
    )
    _require(peek_factor == declared_peek, "peek-both negative control drifted")
    _require(peek_factor > 1, "peek-both control did not break calibration")

    return {
        "schema": EXPECTED_SCHEMA,
        "status": "valid_exact_predictable_arm_eprocess",
        "conditional_null_lr_means": conditional,
        "terminal": {
            "path_count": terminal["path_count"],
            "null_probability_sum": terminal["null_probability_sum"],
            "null_mean_e_value": terminal["null_mean_e_value"],
        },
        "bounded_stopping": {
            "leaf_count": stopped["leaf_count"],
            "null_probability_sum": stopped["null_probability_sum"],
            "null_mean_stopped_e_value": stopped["null_mean_stopped_e_value"],
        },
        "peek_both_expected_factor": str(peek_factor),
        "nonpredictable_selection_breaks_calibration": True,
        "contains_production_stopping_rule": False,
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

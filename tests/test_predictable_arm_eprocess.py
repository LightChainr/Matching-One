from __future__ import annotations

from copy import deepcopy
from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from predictable_arm_eprocess import (  # noqa: E402
    Arm,
    path_values,
    peek_both_expected_factor,
    stopped_enumeration,
    validate_contract,
)


class PredictableArmEprocessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = json.loads(
            (ROOT / "analysis" / "predictable_arm_eprocess_contract.json").read_text(
                encoding="utf-8"
            )
        )

    def test_checked_in_contract_closes_exactly(self) -> None:
        result = validate_contract(self.contract)
        self.assertEqual(result["conditional_null_lr_means"], {"A": "1", "B": "1"})
        self.assertEqual(result["terminal"]["path_count"], 16)
        self.assertEqual(result["terminal"]["null_mean_e_value"], "1")
        self.assertEqual(result["bounded_stopping"]["null_mean_stopped_e_value"], "1")
        self.assertEqual(result["peek_both_expected_factor"], "4/3")
        self.assertFalse(result["contains_production_stopping_rule"])

    def test_arm_choice_depends_only_on_past_outcomes(self) -> None:
        arms = {
            "A": Arm("A", Fraction(1, 3), Fraction(2, 3)),
            "B": Arm("B", Fraction(1, 4), Fraction(1, 2)),
        }
        policy = self.contract["predictable_policy"]
        _, _, selected = path_values((1, 0, 1, 1), arms, policy)
        self.assertEqual(selected, ("A", "B", "A", "B"))

    def test_bounded_stopping_partitions_null_probability(self) -> None:
        arms = {
            "A": Arm("A", Fraction(1, 3), Fraction(2, 3)),
            "B": Arm("B", Fraction(1, 4), Fraction(1, 2)),
        }
        result = stopped_enumeration(4, Fraction(4), arms, self.contract["predictable_policy"])
        self.assertEqual(result["null_probability_sum"], "1")
        self.assertEqual(result["null_mean_stopped_e_value"], "1")
        self.assertTrue(any(row["threshold_hit"] for row in result["leaves"]))

    def test_current_outcome_peek_is_an_invalid_negative_control(self) -> None:
        arms = {
            "A": Arm("A", Fraction(1, 3), Fraction(2, 3)),
            "B": Arm("B", Fraction(1, 4), Fraction(1, 2)),
        }
        self.assertEqual(peek_both_expected_factor(arms), Fraction(4, 3))

    def test_contract_drift_fails_closed(self) -> None:
        changed = deepcopy(self.contract)
        changed["peek_both_negative_control_expected_factor"] = "1"
        with self.assertRaisesRegex(ValueError, "negative control drifted"):
            validate_contract(changed)

        changed = deepcopy(self.contract)
        changed["predictable_policy"]["after_success"] = "C"
        with self.assertRaisesRegex(ValueError, "unknown arm"):
            validate_contract(changed)


if __name__ == "__main__":
    unittest.main()

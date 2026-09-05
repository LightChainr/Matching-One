
from __future__ import annotations
from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from bernoulli_logit_hazard_oracle import threshold_oracle, validate_contract  # noqa: E402


class BernoulliLogitHazardOracleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = json.loads(
            (ROOT / "analysis" / "bernoulli_logit_hazard_contract.json").read_text(
                encoding="utf-8"
            )
        )

    def test_checked_in_contract_closes_exactly(self) -> None:
        result = validate_contract(self.contract)
        self.assertTrue(result["russo_identity_exact"])
        self.assertTrue(result["natural_parameter_score_identity_exact"])
        self.assertTrue(result["logit_derivative_is_conditional_occupation_gap"])
        self.assertFalse(result["contains_monte_carlo"])
        self.assertFalse(result["identifies_four_arm_exponent"])

    def test_three_identity_routes_agree(self) -> None:
        for p in (Fraction(1, 5), Fraction(1, 2), Fraction(4, 5)):
            row = threshold_oracle(5, 3, p)
            self.assertEqual(row["dF_dp"], row["russo_pivotal_sum"])
            self.assertEqual(row["dF_deta"], row["score_covariance"])
            self.assertEqual(
                row["d_logitF_deta"], row["conditional_occupation_gap"]
            )

    def test_threshold_edge_cases_remain_exact(self) -> None:
        p = Fraction(2, 5)
        any_open = threshold_oracle(3, 1, p)
        all_open = threshold_oracle(3, 3, p)
        self.assertEqual(any_open["F"], 1 - (1 - p) ** 3)
        self.assertEqual(all_open["F"], p**3)


if __name__ == "__main__":
    unittest.main()

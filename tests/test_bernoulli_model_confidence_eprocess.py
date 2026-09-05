
from __future__ import annotations
from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from bernoulli_model_confidence_eprocess import (  # noqa: E402
    confidence_sets,
    likelihood_ratio,
    model_certificate,
    validate_contract,
)


class BernoulliModelConfidenceEprocessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = json.loads(
            (ROOT / "analysis" / "bernoulli_model_confidence_contract.json").read_text(
                encoding="utf-8"
            )
        )
        cls.models = (Fraction(1, 3), Fraction(2, 3))

    def test_checked_in_contract_closes_exactly(self) -> None:
        result = validate_contract(self.contract)
        self.assertEqual(result["threshold"], "4")
        self.assertEqual(result["simultaneous_true_model_coverage_at_least"], "3/4")
        self.assertTrue(result["fixed_time_e_means_equal_one"])
        self.assertTrue(result["bounded_optional_stopping_exact"])
        self.assertFalse(result["off_model_coverage_defined"])
        for row in result["certificates"]:
            self.assertEqual(row["ever_excluded_probability"], "13/81")
            self.assertEqual(row["stopped_e_expectation"], "1")

    def test_likelihood_ratios_are_reciprocal(self) -> None:
        for path in ((1, 1), (0, 1, 0), (1, 0, 1, 1)):
            forward = likelihood_ratio(path, self.models[0], self.models[1])
            reverse = likelihood_ratio(path, self.models[1], self.models[0])
            self.assertEqual(forward * reverse, 1)

    def test_confidence_set_can_exclude_and_reinclude_a_model(self) -> None:
        sets = confidence_sets((1, 1, 0, 0), self.models, Fraction(4))
        self.assertNotIn(0, sets[2])
        self.assertIn(0, sets[4])

    def test_ville_bound_holds_for_each_declared_truth(self) -> None:
        for null_index in (0, 1):
            row = model_certificate(self.models, null_index, 4, Fraction(1, 4))
            self.assertTrue(row["ville_bound_holds"])
            self.assertEqual(row["fixed_time_expectations"], ["1"] * 5)


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

from copy import deepcopy
from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from threshold_histogram_profile import (  # noqa: E402
    beta_density_coefficients,
    beta_raw_moment,
    density_coefficients,
    evaluate_polynomial,
    integrate_density,
    mixture_weights,
    parse_histogram,
    raw_moments,
    validate_contract,
)


class ThresholdHistogramProfileTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contract = json.loads(
            (ROOT / "analysis" / "threshold_histogram_profile_contract.json").read_text(
                encoding="utf-8"
            )
        )

    def test_checked_in_synthetic_profile_is_exact(self) -> None:
        result = validate_contract(self.contract)
        self.assertEqual(result["mixture_weights_by_rank"], ["1/8", "3/8", "3/8", "1/8"])
        self.assertEqual(result["mean_from_ranks"], "1/2")
        self.assertTrue(result["cdf_normalized_exactly"])
        self.assertFalse(result["contains_empirical_result"])

    def test_integer_order_statistic_density_normalizes(self) -> None:
        for n in range(1, 8):
            for rank in range(1, n + 1):
                density = beta_density_coefficients(n, rank)
                cdf = integrate_density(density)
                self.assertEqual(evaluate_polynomial(cdf, Fraction(0)), 0)
                self.assertEqual(evaluate_polynomial(cdf, Fraction(1)), 1)
                self.assertEqual(beta_raw_moment(n, rank, 0), 1)
                self.assertEqual(beta_raw_moment(n, rank, 1), Fraction(rank, n + 1))

    def test_arbitrary_histogram_mean_matches_rank_identity(self) -> None:
        minus = {1: 2, 3: 1, 5: 2}
        plus = {2: 1, 4: 3, 5: 1}
        weights = mixture_weights(minus, plus, 5)
        moments = raw_moments(weights, 3)
        rank_mean = sum(
            (weight * rank for rank, weight in enumerate(weights, start=1)), Fraction(0)
        ) / 6
        self.assertEqual(moments[1], rank_mean)
        cdf = integrate_density(density_coefficients(weights))
        self.assertEqual(evaluate_polynomial(cdf, Fraction(1)), 1)

    def test_histogram_validation_fails_closed(self) -> None:
        for raw, message in [({"0": 1}, "out of range"), ({"01": 1}, "canonical"), ({"1": -1}, "nonnegative")]:
            with self.subTest(raw=raw):
                with self.assertRaisesRegex(ValueError, message):
                    parse_histogram(raw, 4, "test")
        with self.assertRaisesRegex(ValueError, "sample counts must match"):
            mixture_weights({1: 1}, {2: 2}, 4)
        with self.assertRaisesRegex(ValueError, "rank is out of range"):
            mixture_weights({5: 1}, {2: 1}, 4)

    def test_contract_coefficient_and_moment_drift_is_rejected(self) -> None:
        changed = deepcopy(self.contract)
        changed["expected_density_power_coefficients"][1] = "4"
        with self.assertRaisesRegex(ValueError, "density coefficients drifted"):
            validate_contract(changed)
        changed = deepcopy(self.contract)
        changed["expected_raw_moments_0_through_6"][4] = "0"
        with self.assertRaisesRegex(ValueError, "raw moments drifted"):
            validate_contract(changed)

    def test_empirical_result_fields_are_rejected(self) -> None:
        changed = deepcopy(self.contract)
        changed["tail_fit"] = {"exponent": "4/3"}
        with self.assertRaisesRegex(ValueError, "forbidden empirical keys"):
            validate_contract(changed)


if __name__ == "__main__":
    unittest.main()

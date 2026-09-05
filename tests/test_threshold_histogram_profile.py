
from __future__ import annotations
from fractions import Fraction
import json
import math
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from threshold_histogram_profile import (  # noqa: E402
    beta_density_coefficients,
    beta_raw_moment,
    central_moments,
    density_coefficients,
    exact_shape_invariants,
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
        self.assertEqual(
            result["central_moments_0_through_6"],
            ["1", "0", "1/15", "0", "1/112", "0", "1/672"],
        )
        self.assertEqual(
            result["exact_shape_invariants"],
            {
                "variance": "1/15",
                "signed_skewness_squared": "0",
                "kurtosis": "225/112",
                "excess_kurtosis": "-111/112",
                "signed_standardized_fifth_squared": "0",
                "standardized_sixth": "1125/224",
            },
        )
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

    def test_shape_invariants_are_exact_under_positive_affine_change(self) -> None:
        weights = mixture_weights({1: 2, 3: 1, 5: 2}, {2: 1, 4: 3, 5: 1}, 5)
        raw = raw_moments(weights, 6)
        scale, offset = Fraction(7, 3), Fraction(-2, 5)
        transformed = [
            sum(
                (
                    Fraction(math.comb(order, index))
                    * (scale ** index)
                    * (offset ** (order - index))
                    * raw[index]
                    for index in range(order + 1)
                ),
                Fraction(0),
            )
            for order in range(7)
        ]
        original_shape = exact_shape_invariants(raw)
        transformed_shape = exact_shape_invariants(transformed)
        self.assertEqual(
            transformed_shape["variance"], scale * scale * original_shape["variance"]
        )
        for key in set(original_shape) - {"variance"}:
            self.assertEqual(transformed_shape[key], original_shape[key])
        central = central_moments(raw)
        transformed_central = central_moments(transformed)
        self.assertEqual(
            transformed_central,
            [central[order] * (scale ** order) for order in range(7)],
        )

    def test_shape_invariants_reject_degenerate_or_inconsistent_moments(self) -> None:
        with self.assertRaisesRegex(ValueError, "positive variance"):
            exact_shape_invariants([Fraction(1), *([Fraction(0)] * 6)])
        with self.assertRaisesRegex(ValueError, "positive variance"):
            exact_shape_invariants(
                [Fraction(1), Fraction(0), Fraction(-1), *([Fraction(0)] * 4)]
            )
        with self.assertRaisesRegex(ValueError, "fourth moment"):
            exact_shape_invariants(
                [
                    Fraction(1),
                    Fraction(0),
                    Fraction(1),
                    Fraction(0),
                    Fraction(0),
                    Fraction(0),
                    Fraction(1),
                ]
            )


if __name__ == "__main__":
    unittest.main()

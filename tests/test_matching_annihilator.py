
from __future__ import annotations
from pathlib import Path
import sys
import unittest

import mpmath as mp


ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT))
from scripts import matching_annihilator as module  # noqa: E402


class MatchingAnnihilatorTests(unittest.TestCase):
    def setUp(self) -> None:
        mp.mp.dps = 80

    def test_two_size_weights_match_known_ratio(self) -> None:
        exponent = mp.mpf(13) / 4
        matrix = module.constraint_matrix([15, 16], [exponent])
        weights = module.solve_weights(matrix)
        expected_ratio = -(mp.mpf(15) / 16) ** exponent
        self.assertLess(abs(weights[0] / weights[1] - expected_ratio), mp.mpf("1e-70"))
        self.assertLess(abs(weights[0] + weights[1] - 1), mp.mpf("1e-70"))

    def test_three_size_annihilates_two_terms(self) -> None:
        result = module.compute(
            [14, 15, 16], [mp.mpf(13) / 4, mp.mpf(25) / 4], None
        )
        self.assertTrue(
            all(abs(mp.mpf(value)) < mp.mpf("1e-60") for value in result.residuals)
        )
        self.assertGreater(mp.mpf(result.l1_norm), 1)

    def test_weighted_value_removes_declared_powers(self) -> None:
        sizes = [12, 13, 14, 15, 16]
        exponents = [mp.mpf(13) / 4, mp.mpf(25) / 4]
        target = mp.mpf("0.59274605079")
        values = [
            target
            + mp.mpf(7) * mp.power(size, -exponents[0])
            - mp.mpf(3) * mp.power(size, -exponents[1])
            for size in sizes
        ]
        result = module.compute(sizes, exponents, values)
        self.assertIsNotNone(result.weighted_value)
        assert result.weighted_value is not None
        self.assertLess(abs(mp.mpf(result.weighted_value) - target), mp.mpf("1e-35"))


if __name__ == "__main__":
    unittest.main()

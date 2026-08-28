from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import unittest

import mpmath as mp
import yaml


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "predictions" / "gaussian_norm5_radial_competitors_20260828.yaml"


class Norm5RadialPredictionArtifactTests(unittest.TestCase):
    def test_deltaM_and_root_ratios_follow_frozen_exponents(self) -> None:
        mp.mp.dps = 80
        with ARTIFACT.open(encoding="utf-8") as handle:
            frozen = yaml.safe_load(handle)

        self.assertEqual(frozen["claim_level"], "C0")
        angular = -mp.mpf(14) / 25
        slope_ratio = mp.power(5, mp.mpf(3) / 8)
        expected = {
            "x21_over_4_thermal_level4": Fraction(13, 8),
            "x14_over_3_V13_parity_failure": Fraction(4, 3),
            "x17_over_4_W22_log_leakage": Fraction(9, 8),
        }

        for name, alpha_fraction in expected.items():
            model = frozen["models"][name]
            alpha = mp.mpf(alpha_fraction.numerator) / alpha_fraction.denominator
            delta_expected = angular * mp.power(5, -alpha)
            delta_reported = mp.mpf(model["deltaM_child_over_parent"])
            self.assertLess(abs(delta_reported - delta_expected), mp.mpf("1e-48"))

            root_expected = delta_expected / slope_ratio
            root_reported = mp.mpf(
                model["full_curve_if_run"]["root_gap_child_over_parent"]
            )
            self.assertLess(abs(root_reported - root_expected), mp.mpf("1e-48"))

            root_alpha_reported = Fraction(
                model["full_curve_if_run"]["root_gap_exponent_in_N"]
            )
            self.assertEqual(root_alpha_reported, alpha_fraction + Fraction(3, 8))


if __name__ == "__main__":
    unittest.main()

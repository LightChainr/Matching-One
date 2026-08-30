from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import unittest

import mpmath as mp
import yaml


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "predictions" / "gaussian_normalized_p4_multiplier_spectrum_20260828.yaml"


class GaussianNormalizedP4MultiplierSpectrumTests(unittest.TestCase):
    def setUp(self) -> None:
        mp.mp.dps = 80
        with ARTIFACT.open(encoding="utf-8") as handle:
            self.frozen = yaml.safe_load(handle)

    def test_pure_h4_normalized_projectors_have_positive_radial_ratio(self) -> None:
        exponents = {
            "P4_S": Fraction(1, 1),
            "P4_Dprime": Fraction(5, 8),
            "P4_D": Fraction(13, 8),
            "P4_Sprime": Fraction(5, 4),
        }
        for section, q in (("norm2", 2), ("norm5", 5)):
            reported = self.frozen[section]["normalized_pure_H4_child_over_parent"]
            for name, alpha_fraction in exponents.items():
                alpha = mp.mpf(alpha_fraction.numerator) / alpha_fraction.denominator
                expected = mp.power(q, -alpha)
                actual = mp.mpf(str(reported[name]["decimal"]))
                self.assertGreater(actual, 0)
                self.assertLess(abs(actual - expected), mp.mpf("1e-48"))

    def test_raw_and_normalized_norm5_harmonic_factors_are_distinct(self) -> None:
        r4 = Fraction(-14, 25)
        r12 = Fraction(23506, 15625)
        normalized = r12 / r4
        self.assertEqual(normalized, Fraction(-1679, 625))
        stored = Fraction(
            self.frozen["norm5"]["normalized_H12_over_H4_angular_ratio"]["exact"]
        )
        self.assertEqual(stored, normalized)

    def test_norm5_h12_normalized_targets_use_r12_over_r4(self) -> None:
        angular = mp.mpf(-1679) / 625
        for name, alpha_fraction in (
            ("P4_D", Fraction(13, 8)),
            ("P4_Sprime", Fraction(5, 4)),
        ):
            alpha = mp.mpf(alpha_fraction.numerator) / alpha_fraction.denominator
            expected = angular * mp.power(5, -alpha)
            actual = mp.mpf(
                str(self.frozen["norm5"]["normalized_same_radial_H12_adversary"][name]["decimal"])
            )
            self.assertLess(abs(actual - expected), mp.mpf("1e-48"))

    def test_derivative_exponent_propagation_for_radial_adversaries(self) -> None:
        models = self.frozen["thermal_H4_radial_adversaries"]
        thermal_step = Fraction(3, 8)
        for name in (
            "x21_over_4",
            "x14_over_3_V13_parity_failure",
            "x17_over_4_W22_log_leakage",
        ):
            model = models[name]
            self.assertEqual(
                Fraction(model["D_alpha"]) - Fraction(model["Sprime_alpha"]),
                thermal_step,
            )


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import unittest

import mpmath as mp
import yaml


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "predictions" / "gaussian_multiplier_derivative_spectrum_20260828.yaml"


class GaussianMultiplierDerivativeSpectrumTests(unittest.TestCase):
    def load(self) -> dict:
        with ARTIFACT.open(encoding="utf-8") as handle:
            return yaml.safe_load(handle)

    def test_normalized_h4_targets_cancel_angular_ratio(self) -> None:
        mp.mp.dps = 80
        frozen = self.load()
        self.assertEqual(frozen["version"], 3)
        self.assertEqual(frozen["claim_level"], "C0")

        channel_exponents = {
            name: Fraction(payload["exponent_alpha_in_N"])
            for name, payload in frozen["channel_models"].items()
        }
        for section, q in (("norm2", 2), ("norm5", 5)):
            for name, alpha_fraction in channel_exponents.items():
                alpha = mp.mpf(alpha_fraction.numerator) / alpha_fraction.denominator
                expected = mp.power(q, -alpha)
                reported = mp.mpf(
                    frozen[section]["normalized_P4_child_over_parent"][name]["decimal"]
                )
                self.assertGreater(reported, 0)
                self.assertLess(abs(reported - expected), mp.mpf("1e-48"))

    def test_raw_h4_contrast_retains_angular_sign(self) -> None:
        mp.mp.dps = 80
        frozen = self.load()
        norm2_raw = -mp.power(2, -mp.mpf(13) / 8)
        norm5_raw = -mp.mpf(14) / 25 * mp.power(5, -mp.mpf(13) / 8)
        self.assertEqual(
            frozen["norm2"]["raw_H4_contrast_child_over_parent"]["thermal_odd_D"],
            "-2^(-13/8)",
        )
        reported5 = mp.mpf(
            frozen["norm5"]["raw_H4_contrast_child_over_parent"]["thermal_odd_D"][
                "decimal"
            ]
        )
        self.assertLess(abs(reported5 - norm5_raw), mp.mpf("1e-48"))
        self.assertLess(norm2_raw, 0)
        self.assertLess(reported5, 0)

    def test_derivative_exponents_are_one_thermal_step_above_center_channels(self) -> None:
        frozen = self.load()
        exponents = {
            name: Fraction(payload["exponent_alpha_in_N"])
            for name, payload in frozen["channel_models"].items()
        }
        thermal_step = Fraction(3, 8)
        self.assertEqual(
            exponents["identity_even_S"] - exponents["identity_even_Dprime"],
            thermal_step,
        )
        self.assertEqual(
            exponents["thermal_odd_D"] - exponents["thermal_odd_Sprime"],
            thermal_step,
        )

    def test_fixed_thermal_radial_adversaries_propagate_to_normalized_Sprime(self) -> None:
        mp.mp.dps = 80
        frozen = self.load()
        thermal_step = Fraction(3, 8)
        competitors = frozen["thermal_radial_competitor_rule"]["candidates"]
        for name, payload in competitors.items():
            alpha_d = Fraction(payload["central_D_alpha_N"])
            alpha_sprime = Fraction(payload["Sprime_alpha_N"])
            self.assertEqual(alpha_d - alpha_sprime, thermal_step, name)

        for section, q in (("norm2", 2), ("norm5", 5)):
            adversaries = frozen[section]["normalized_P4_thermal_radial_adversaries"]
            for name, payload in adversaries.items():
                alpha_sprime_fraction = Fraction(competitors[name]["Sprime_alpha_N"])
                alpha_sprime = (
                    mp.mpf(alpha_sprime_fraction.numerator)
                    / alpha_sprime_fraction.denominator
                )
                expected = mp.power(q, -alpha_sprime)
                reported = mp.mpf(payload["Sprime_decimal"])
                self.assertGreater(reported, 0)
                self.assertLess(abs(reported - expected), mp.mpf("1e-48"))

    def test_norm5_h12_normalized_adversary_uses_ratio_of_angular_levers(self) -> None:
        mp.mp.dps = 80
        frozen = self.load()
        adversary = frozen["norm5_same_radial_H12_adversary_for_normalized_thermal_P4"]
        angular = mp.mpf(-1679) / 625
        self.assertEqual(adversary["normalized_angular_ratio_H12_over_H4"], "-1679/625")
        for key, alpha in (
            ("thermal_odd_D", mp.mpf(13) / 8),
            ("thermal_odd_Sprime", mp.mpf(5) / 4),
        ):
            expected = angular * mp.power(5, -alpha)
            reported = mp.mpf(adversary[key]["decimal"])
            self.assertLess(abs(reported - expected), mp.mpf("1e-48"))
            self.assertLess(reported, 0)


if __name__ == "__main__":
    unittest.main()

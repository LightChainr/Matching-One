from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import unittest

import mpmath as mp
import yaml


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "predictions" / "gaussian_multiplier_derivative_spectrum_20260828.yaml"


class GaussianMultiplierDerivativeSpectrumTests(unittest.TestCase):
    def test_h4_multiplier_targets_follow_channel_exponents(self) -> None:
        mp.mp.dps = 80
        with ARTIFACT.open(encoding="utf-8") as handle:
            frozen = yaml.safe_load(handle)

        self.assertEqual(frozen["claim_level"], "C0")
        channel_exponents = {
            name: Fraction(payload["exponent_alpha_in_N"])
            for name, payload in frozen["channel_models"].items()
        }
        for section, q, angular in (
            ("norm2", 2, -mp.mpf(1)),
            ("norm5", 5, -mp.mpf(14) / 25),
        ):
            for name, alpha_fraction in channel_exponents.items():
                alpha = mp.mpf(alpha_fraction.numerator) / alpha_fraction.denominator
                expected = angular * mp.power(q, -alpha)
                reported = mp.mpf(
                    frozen[section]["fixed_child_over_parent"][name]["decimal"]
                )
                self.assertLess(abs(reported - expected), mp.mpf("1e-48"))

    def test_derivative_exponents_are_one_thermal_step_above_center_channels(self) -> None:
        with ARTIFACT.open(encoding="utf-8") as handle:
            frozen = yaml.safe_load(handle)

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

    def test_fixed_thermal_radial_adversaries_propagate_to_Sprime(self) -> None:
        mp.mp.dps = 80
        with ARTIFACT.open(encoding="utf-8") as handle:
            frozen = yaml.safe_load(handle)

        thermal_step = Fraction(3, 8)
        competitors = frozen["thermal_radial_competitor_rule"]["candidates"]
        for name, payload in competitors.items():
            alpha_d = Fraction(payload["central_D_alpha_N"])
            alpha_sprime = Fraction(payload["Sprime_alpha_N"])
            self.assertEqual(alpha_d - alpha_sprime, thermal_step, name)

        for section, q, angular in (
            ("norm2", 2, -mp.mpf(1)),
            ("norm5", 5, -mp.mpf(14) / 25),
        ):
            adversaries = frozen[section]["thermal_radial_adversaries"]
            for name, payload in adversaries.items():
                source_name = (
                    "x14_over_3_V13" if name == "x14_over_3_V13" else "x17_over_4_W22"
                )
                alpha_sprime_fraction = Fraction(
                    competitors[source_name]["Sprime_alpha_N"]
                )
                alpha_sprime = (
                    mp.mpf(alpha_sprime_fraction.numerator)
                    / alpha_sprime_fraction.denominator
                )
                expected = angular * mp.power(q, -alpha_sprime)
                reported = mp.mpf(payload["Sprime_decimal"])
                self.assertLess(abs(reported - expected), mp.mpf("1e-48"))


if __name__ == "__main__":
    unittest.main()

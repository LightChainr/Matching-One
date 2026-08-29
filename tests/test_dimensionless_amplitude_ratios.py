from __future__ import annotations

from fractions import Fraction
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from dimensionless_amplitude_ratios import (  # noqa: E402
    descriptive_p48_ratios,
    metric_free_exponents,
    ratio_exponent,
    scaled_ratio,
)


class DimensionlessAmplitudeRatioTests(unittest.TestCase):
    def test_metric_free_ratios_have_vanishing_n_power(self) -> None:
        exponents = metric_free_exponents()
        self.assertEqual(exponents["R_I"], 0)
        self.assertEqual(exponents["R_T"], 0)

    def test_raw_odd_even_ratio_does_not_cancel_n(self) -> None:
        self.assertEqual(ratio_exponent("P4_D", ("P4_S",)), Fraction(-5, 8))

    def test_scaled_formula_matches_unscaled_cancellation(self) -> None:
        self.assertAlmostEqual(scaled_ratio(-0.0237, -0.00955, 1.746), 1.4215, places=3)

    def test_descriptive_p48_reconstruction_is_labeled(self) -> None:
        payload = descriptive_p48_ratios()
        self.assertTrue(payload["not_a_numerical_target_for_issue_57"])
        self.assertEqual(payload["role"], "development_descriptive_only")
        rows = payload["rows"]
        self.assertEqual([row["N"] for row in rows], [185, 265])
        self.assertAlmostEqual(rows[0]["R_I"], 1.42, places=1)
        self.assertAlmostEqual(rows[0]["R_T"], 5.27, places=1)
        self.assertAlmostEqual(rows[1]["R_I"], 1.20, places=1)
        self.assertAlmostEqual(rows[1]["R_T"], 5.12, places=1)

    def test_R_T_is_more_stable_than_raw_scaled_Sprime_in_descriptive_view(self) -> None:
        rows = descriptive_p48_ratios()["rows"]
        sprime = [abs(row["A_S_prime"]) for row in rows]
        thermal = [row["R_T"] for row in rows]
        sprime_rel = abs(sprime[1] - sprime[0]) / abs(sprime[0])
        thermal_rel = abs(thermal[1] - thermal[0]) / abs(thermal[0])
        self.assertLess(thermal_rel, sprime_rel)


if __name__ == "__main__":
    unittest.main()

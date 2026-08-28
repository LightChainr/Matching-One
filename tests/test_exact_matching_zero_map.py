#!/usr/bin/env python3
from pathlib import Path
import sys
import unittest

import mpmath as mp


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import exact_matching_zero_map as zero_map  # noqa: E402


class ExactMatchingZeroMapTests(unittest.TestCase):
    def test_exact_matching_partner_transform(self):
        coefficients = [-1, 0, 0, 6, 0, 0, 0, -18, 18, -4]
        partner = zero_map.matching_partner_coefficients(coefficients)
        with mp.workdps(80):
            for value in (mp.mpc("0.37", "0.2"), mp.mpc("-1.1", "0.7")):
                self.assertAlmostEqual(
                    abs(zero_map.evaluate(partner, value) + zero_map.evaluate(coefficients, 1 - value)),
                    0.0,
                    places=60,
                )

    def test_high_precision_roots_pass_pairing_audits(self):
        coefficients = [-1, 0, 0, 6, 0, 0, 0, -18, 18, -4]
        summary, roots = zero_map.analyze_polynomial("axis", 3, coefficients, 70)
        self.assertEqual(summary["status"], "OK")
        self.assertEqual(len(roots), 9)
        self.assertEqual(summary["metrics"]["real_root_count"], 3)
        self.assertEqual(summary["metrics"]["nonreal_root_count"], 6)
        self.assertLess(summary["audit"]["max_conjugate_pair_distance"], mp.mpf("1e-50"))
        self.assertLess(summary["audit"]["max_matching_partner_pair_distance"], mp.mpf("1e-50"))
        self.assertLess(summary["audit"]["max_normalized_polynomial_residual"], mp.mpf("1e-50"))

    def test_inverse_n_prediction_scores_a_heldout_value(self):
        rows = [
            {"L": 1, "N": 1, "metrics": {"physical_root_0_1": mp.mpf("0.5")}},
            {"L": 2, "N": 4, "metrics": {"physical_root_0_1": mp.mpf("0.54")}},
        ]
        result = zero_map.linear_inverse_n_prediction(rows[0], rows[1], 9, "physical_root_0_1")
        self.assertEqual(result["training_L"], [1, 2])
        self.assertEqual(result["target_N"], 9)
        self.assertTrue(mp.isfinite(result["prediction"]))


if __name__ == "__main__":
    unittest.main()

import json
import math
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import e4_balanced_pell_estimator as estimator  # noqa: E402


class E4BalancedPellEstimatorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = json.loads((ROOT / "analysis" / "e4_balanced_pell_manifest.yaml").read_text())
        cls.result = estimator.analyze(cls.config)

    def test_pell_families_and_e4_signs(self):
        for row in self.result["rows"]:
            self.assertEqual(row["negative"]["pell_residual"], -2)
            self.assertEqual(row["positive"]["pell_residual"], 1)
            self.assertLess(row["negative"]["E4"], 0)
            self.assertGreater(row["positive"]["E4"], 0)

    def test_frozen_weights_match_issue_values(self):
        expected = [
            (0.0714582925, 0.9285417075),
            (0.0714531199, 0.9285468801),
            (0.0714531180, 0.9285468820),
        ]
        for row, (negative, positive) in zip(self.result["rows"], expected):
            self.assertAlmostEqual(row["weights"]["negative"], negative, places=9)
            self.assertAlmostEqual(row["weights"]["positive"], positive, places=9)

    def test_weights_are_positive_normalized_and_cancel_h4(self):
        for row in self.result["rows"]:
            self.assertGreater(row["weights"]["negative"], 0)
            self.assertGreater(row["weights"]["positive"], 0)
            self.assertAlmostEqual(row["weight_sum"], 1.0, places=15)
            self.assertLess(row["h4_cancellation_relative"], 1e-14)

    def test_scalar_ratio_converges_to_fundamental_unit_prediction(self):
        target = (2 + math.sqrt(3)) ** -7
        ratios = [row["next_scalar_coefficient_ratio"] for row in self.result["rows"][:-1]]
        self.assertLess(abs(ratios[-1] - target), abs(ratios[0] - target))
        self.assertAlmostEqual(self.result["asymptotic_scalar_ratio"], target, places=15)

    def test_selection_boundary_excludes_target_outcomes(self):
        boundary = self.config["selection_boundary"]
        self.assertIn("Root estimates", boundary)
        self.assertIn("forbidden", boundary)


if __name__ == "__main__":
    unittest.main()

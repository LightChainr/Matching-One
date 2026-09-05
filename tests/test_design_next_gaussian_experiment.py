from fractions import Fraction
from pathlib import Path
import json
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import design_next_gaussian_experiment as design  # noqa: E402


class GaussianExperimentDesignTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = json.loads((ROOT / "analysis" / "gaussian_experiment_design_manifest.yaml").read_text())
        cls.result = design.design(cls.config)

    def test_norm5_exact_harmonic_ratios(self):
        first, second = (8, 1), (7, 4)
        multiplier = (2, -1)
        self.assertEqual(design.angular_ratio(first, second, multiplier, 4), Fraction(-14, 25))
        self.assertEqual(design.angular_ratio(first, second, multiplier, 8), Fraction(-1054, 625))
        self.assertEqual(design.angular_ratio(first, second, multiplier, 12), Fraction(23506, 15625))

    def test_norm5_is_near_optimal_for_h4_h12(self):
        rows = self.result["norm5_reference_rows"]
        self.assertEqual(len(rows), 2)
        self.assertTrue(all(row["h4_vs_h12_rank"] <= 10 for row in rows))

    def test_candidates_obey_caps_and_have_exact_geometry(self):
        self.assertEqual(self.result["candidate_count"], len(self.result["candidates"]))
        for row in self.result["candidates"][:self.result["output_top_n"]]:
            self.assertLessEqual(row["target_N"], self.config["caps"]["target_N_max"])
            self.assertLessEqual(row["multiplier_norm"], self.config["caps"]["multiplier_norm_max"])
            self.assertIn("numerator", row["harmonic_angular_ratios"]["H4"])
            self.assertGreater(row["estimated_cpu_seconds"], 0)


if __name__ == "__main__":
    unittest.main()

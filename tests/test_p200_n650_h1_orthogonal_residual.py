import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from analyze_p200_n650_h1_orthogonal_residual import analyze


class P200N650H1OrthogonalResidualTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.payload = analyze(
            ROOT / "results/local-20260829/P200-n650-mixed-join-smoke/n650_20k.batches.csv",
            ROOT / "results/local-20260829/P200-n650-mixed-join-smoke/n650_20k.metadata.json",
            ROOT / "predictions/p200_n650_mixed_join_phaseB_20260829.json",
        )

    def test_joint_covariance_is_complete_and_symmetric(self):
        block = self.payload["joint_same_stream_state"]
        self.assertEqual(len(block["state_order"]), 8)
        self.assertEqual(len(block["delete_one_covariance"]), 8)
        for first in range(8):
            for second in range(8):
                self.assertAlmostEqual(
                    block["delete_one_covariance"][first][second],
                    block["delete_one_covariance"][second][first],
                    places=14,
                )

    def test_h1_orthogonal_mode_survives_but_geometry_does_not(self):
        projection = self.payload["one_gain_H1_projection"]
        self.assertLess(projection["orthogonal_common_mode"]["p_value"], 1e-100)
        self.assertGreater(
            projection["geometry_difference_after_same_gain"]["p_value"], 0.9
        )
        self.assertAlmostEqual(
            projection["derived_mean"][1], -34.495921874765514, places=10
        )

    def test_identifiability_boundary_is_explicit(self):
        self.assertIn(
            "not identified",
            self.payload["observable_boundary"]["path_or_state_memory"],
        )
        self.assertAlmostEqual(
            self.payload["one_gain_H1_projection"]["gain_free_collinearity_certificate"]["z"],
            139.37330658266495,
            places=10,
        )


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_threshold_rank_root_doubling import score  # noqa: E402


class ThresholdRankRootDoublingTests(unittest.TestCase):
    def test_exact_lineage_relation_has_zero_residual(self) -> None:
        sizes = [65, 85, 130, 145, 170]
        values = [-4.0, -8.0, -1.0, 0.5, -2.0]
        covariance = [
            [0.01 if i == j else 0.0 for j in range(len(sizes))]
            for i in range(len(sizes))
        ]
        summary = {
            "format_version": 2,
            "sizes": sizes,
            "nonlinear_estimator": {
                "root_gap_method": "delete_one_jackknife_pseudovalues"
            },
            "metrics": {
                "root_gap": {
                    "means": dict(zip(map(str, sizes), values)),
                    "covariance_of_means": covariance,
                }
            },
        }
        result = score(summary)
        full = result["full_cross_size_covariance"]
        self.assertEqual(full["target_ratio"], -0.25)
        self.assertAlmostEqual(full["lineages"][0]["observed_ratio"], -0.25)
        self.assertAlmostEqual(full["lineages"][1]["observed_ratio"], -0.25)
        self.assertAlmostEqual(full["joint_residual_chi_square"], 0.0)

    def test_full_and_diagonal_scores_are_both_retained(self) -> None:
        sizes = [65, 85, 130, 145, 170]
        values = [-4.0, -8.0, -1.1, 0.5, -1.8]
        covariance = [
            [0.04 if i == j else 0.01 for j in range(len(sizes))]
            for i in range(len(sizes))
        ]
        summary = {
            "format_version": 2,
            "sizes": sizes,
            "metrics": {
                "root_gap": {
                    "means": dict(zip(map(str, sizes), values)),
                    "covariance_of_means": covariance,
                }
            },
        }
        result = score(summary)
        self.assertIn("full_cross_size_covariance", result)
        self.assertIn("diagonal_cross_size_covariance", result)
        self.assertNotEqual(
            result["full_cross_size_covariance"]["joint_residual_chi_square"],
            result["diagonal_cross_size_covariance"]["joint_residual_chi_square"],
        )


if __name__ == "__main__":
    unittest.main()

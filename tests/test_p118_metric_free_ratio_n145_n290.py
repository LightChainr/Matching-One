from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_p118_metric_free_ratio_n145_n290 import (  # noqa: E402
    FEATURE_ORDER,
    add_covariances,
    marginal_score,
    ratio_coordinates,
    validate_manifest,
)


class P118MetricFreeRatioTests(unittest.TestCase):
    def test_ratio_coordinates_use_the_frozen_monomials(self) -> None:
        stat = {
            "mean_slope": 2.0,
            "P4_S": 4.0,
            "P4_D": 5.0,
            "P4_D_prime": 16.0,
            "P4_S_prime": 30.0,
        }
        self.assertEqual(ratio_coordinates(stat), {"R_I": 2.0, "R_T": 3.0})
        self.assertEqual(FEATURE_ORDER, ("R_I", "R_T"))

    def test_zero_nonlinear_denominator_fails_closed(self) -> None:
        stat = {
            "mean_slope": 2.0,
            "P4_S": 0.0,
            "P4_D": 5.0,
            "P4_D_prime": 16.0,
            "P4_S_prime": 30.0,
        }
        with self.assertRaisesRegex(ValueError, "denominator"):
            ratio_coordinates(stat)

    def test_independent_residual_covariance_is_a_sum(self) -> None:
        self.assertEqual(
            add_covariances([[1.0, 2.0], [2.0, 4.0]], [[3.0, 4.0], [4.0, 5.0]]),
            [[4.0, 6.0], [6.0, 9.0]],
        )

    def test_marginal_score_uses_two_sided_frozen_alpha(self) -> None:
        score = marginal_score(2.0, 4.0, 0.01)
        self.assertAlmostEqual(score["signed_z"], 1.0)
        self.assertAlmostEqual(score["p_value"], 0.31731050786291415)
        self.assertIn("survives", score["decision"])

    def test_manifest_hashes_and_chronology_validate(self) -> None:
        manifest_path = ROOT / "analysis/p118_metric_free_ratio_n145_n290_manifest.json"
        manifest = validate_manifest(manifest_path, ROOT)
        self.assertEqual(manifest["inputs"]["parent_N"], 145)
        self.assertEqual(manifest["inputs"]["child_N"], 290)
        self.assertTrue(manifest["inputs"]["independent_random_streams"])
        self.assertEqual(
            manifest["chronology"]["definition_freeze_commit"],
            "920c6393f6db7927887df0905dfaed81838ca062",
        )
        self.assertEqual(
            manifest["chronology"]["first_target_result_commit"],
            "9675bce5b406247e15c03bca20abef954f26a3a2",
        )

    def test_checked_in_result_retains_frozen_order_when_present(self) -> None:
        result = (
            ROOT
            / "results/server-20260829/P118-metric-free-ratio-N145-N290/score.json"
        )
        if not result.exists():
            self.skipTest("result has not been revealed yet")
        payload = json.loads(result.read_text(encoding="utf-8"))
        self.assertEqual(payload["feature_order"], ["R_I", "R_T"])
        self.assertEqual(payload["constant_response_null"]["joint"]["degrees_of_freedom"], 2)


if __name__ == "__main__":
    unittest.main()

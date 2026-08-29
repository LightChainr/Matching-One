from __future__ import annotations

import json
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
RESULT = (
    ROOT
    / "results"
    / "server-20260829"
    / "mean-jd4-q2-chain"
    / "analysis"
    / "score.json"
)


class MeanJD4Q2ChainScoreTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.score = json.loads(RESULT.read_text(encoding="utf-8"))

    def test_frozen_stream_and_clean_audits(self) -> None:
        self.assertEqual(
            self.score["preregistration_commit"],
            "0b3cebf9cd1b536c859d3fad1591345725a068f7",
        )
        for n in (65, 130, 260):
            record = self.score["provenance"][str(n)]
            self.assertEqual(record["metadata"]["samples_per_pair"], 5_000_000)
            self.assertEqual(record["metadata"]["batches"], 100)
            self.assertEqual(record["metadata"]["replica_counter_first"], 9_300_000_000)
            self.assertFalse(any(record["complement_audit"].values()))
            self.assertEqual(
                record["binary_sha256"],
                "46d8a2690b9a3b1899b3fe61e9a2c16019cb39487d493998c477ca302eaa1223",
            )

    def test_full_covariance_is_symmetric(self) -> None:
        covariance = self.score["full_cross_size_delete_one_covariance"]
        self.assertEqual(len(covariance), 6)
        self.assertTrue(all(len(row) == 6 for row in covariance))
        for i in range(6):
            self.assertGreater(float(covariance[i][i]), 0.0)
            for j in range(6):
                self.assertAlmostEqual(float(covariance[i][j]), float(covariance[j][i]))

    def test_both_preregistered_one_step_phases_are_rejected(self) -> None:
        scores = self.score["primary_phase_first"]
        self.assertGreater(abs(float(scores["N65_to_N130"]["phase_residual_z"])), 20.0)
        self.assertGreater(abs(float(scores["N130_to_N260"]["phase_residual_z"])), 5.0)
        self.assertLess(float(scores["N65_to_N130"]["target_chi2_p_2dof"]), 1e-20)
        self.assertLess(float(scores["N130_to_N260"]["target_chi2_p_2dof"]), 1e-3)

    def test_failed_zy_setup_is_not_the_primary_path(self) -> None:
        failed = self.score["failed_setup_excluded"]
        self.assertEqual(failed["path"], "/workspace/Matching-One-mean-jd4-q2-chain")
        self.assertTrue(self.score["remote_raw"]["65"].startswith("/workspace/mean-jd4-N65-5m/"))

    def test_post_reveal_gate_is_labelled_diagnostic(self) -> None:
        gate = self.score["post_reveal_gate_diagnostic"]
        self.assertIn("not a rescue", gate["status"])
        for ratio in gate["magnitude_ratios"].values():
            self.assertGreater(float(ratio), 1.2)
            self.assertLess(float(ratio), 1.4)


if __name__ == "__main__":
    unittest.main()

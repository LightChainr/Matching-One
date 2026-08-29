#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from analyze_annulus_channel_recurrence import analyze, synthetic_oracles  # noqa: E402


class AnnulusChannelRecurrenceTests(unittest.TestCase):
    def test_synthetic_rank_two_classes_and_fractional_prediction(self) -> None:
        payload = synthetic_oracles()
        self.assertAlmostEqual(payload["J2"]["recurrence"]["Delta"], 0.0, places=12)
        self.assertGreater(payload["R2"]["recurrence"]["Delta"], 0.0)
        self.assertLess(payload["C2"]["recurrence"]["Delta"], 0.0)
        for record in payload.values():
            self.assertLess(record["max_abs_error"], 1e-12)

    def test_committed_block_identifiability_is_not_overclaimed(self) -> None:
        source = ROOT / "results/server-20260829/P225-norm5-multiradius/analysis.json"
        if not source.exists():
            self.skipTest("PR247 production artifact is not present")
        result = analyze(source)
        self.assertIn("saturated", result["identifiability"]["R2"])
        self.assertIn("saturated", result["identifiability"]["C2"])
        self.assertEqual(
            result["model_profiles"]["minus"]["R2_gap1"]["degrees_of_freedom"], 1)
        plus = result["channel_recurrences"]["A_plus"]
        minus = result["channel_recurrences"]["A_minus"]
        self.assertEqual(plus["point_class"], "R2_positive_roots")
        self.assertEqual(minus["point_class"], "C2_principal_branch")
        self.assertLess(abs(plus["Delta_over_SE"]), 1.0)
        self.assertLess(abs(minus["Delta_over_SE"]), 1.0)
        self.assertEqual(result["next_acquisition"]["geometry"]["N"], 365)

    def test_result_is_machine_readable(self) -> None:
        path = ROOT / "results/annulus-radial-design/latest.json"
        if not path.exists():
            self.skipTest("generated result not present")
        payload = json.loads(path.read_text(encoding="utf-8"))
        self.assertEqual(payload["schema"], "matching-one/annulus-channel-recurrence/v1")
        self.assertEqual(payload["issue"], 253)


if __name__ == "__main__":
    unittest.main()

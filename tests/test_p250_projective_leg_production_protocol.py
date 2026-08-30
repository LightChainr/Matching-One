from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p250_projective_leg_power_freeze import freeze  # noqa: E402
from score_z5_charged_multiseparation import read_batches  # noqa: E402
from score_z5_projective_leg_production import support_then_phase  # noqa: E402


SMOKE = ROOT / "results/local-20260830/P250-z5-projective-leg-smoke"


class P250ProjectiveLegProductionProtocolTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = json.loads((SMOKE / "response_2k.json").read_text())
        cls.batches = read_batches(SMOKE / "response_2k.batches.csv")
        cls.smoke_score = json.loads((SMOKE / "score_2k.json").read_text())

    def test_first_qualifying_grid_point_is_10k(self) -> None:
        result = freeze(self.smoke_score)
        self.assertEqual(result["selected_samples"], 10_000)
        self.assertFalse(result["table"][0]["qualifies"])
        self.assertTrue(result["table"][1]["qualifies"])
        self.assertGreater(result["table"][1]["forecasts"]["1"]["support_power_at_alpha_0.01"], 0.98)

    def test_phase_code_path_stays_locked_when_support_gate_fails(self) -> None:
        result = support_then_phase(
            self.payload, self.batches, minimum_pair_z=100.0, support_alpha=0.01
        )
        self.assertFalse(result["support_gate_passed"])
        self.assertFalse(result["phase_closure"]["computed"])
        self.assertEqual(result["phase_closure"]["status"], "locked_support_gate_failed")

    def test_phase_code_path_can_only_unlock_after_both_support_rows(self) -> None:
        result = support_then_phase(
            self.payload, self.batches, minimum_pair_z=0.0, support_alpha=1.0
        )
        self.assertTrue(result["support_gate_passed"])
        self.assertTrue(result["phase_closure"]["computed"])
        self.assertEqual(result["phase_closure"]["joint_zero_score"]["degrees_of_freedom"], 4)

    def test_committed_manifest_freezes_fresh_seed_and_phase_lock(self) -> None:
        manifest = json.loads(
            (ROOT / "analysis/p250_projective_leg_fresh_production_manifest.json").read_text()
        )
        power = json.loads(
            (ROOT / "analysis/p250_projective_leg_power_freeze.json").read_text()
        )
        self.assertEqual(manifest["run"]["samples"], power["selected_samples"])
        self.assertEqual(manifest["run"]["seed"], 25033433720260930)
        self.assertEqual(manifest["run"]["replica_offset"], 0)
        self.assertEqual(manifest["phase_policy"], "locked_until_both_support_stages_pass")


if __name__ == "__main__":
    unittest.main()

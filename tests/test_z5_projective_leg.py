from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_z5_charged_multiseparation import read_batches  # noqa: E402
from score_z5_projective_leg import score  # noqa: E402
from z5_charged_multiseparation_mc import FIELD_ORDER  # noqa: E402
from z5_projective_leg_multiseparation_mc import SMOKE_CAP, exact_gate, run  # noqa: E402


RESULT = ROOT / "results/local-20260830/P250-z5-projective-leg-smoke"


class Z5ProjectiveLegTests(unittest.TestCase):
    def test_exact_gate_has_black_and_white_propagation(self) -> None:
        gate = exact_gate()
        self.assertTrue(gate["passed"])
        for hand in gate["constructed_same_component_propagation"]:
            self.assertEqual(hand["component_rank"], 1)
            self.assertEqual(hand["white_component_rank"], 1)
            self.assertEqual(hand["line_gcd"], 1)
            self.assertEqual(hand["black_separated_scalar_values"], [1, 1, 1, 1])
            self.assertEqual(hand["white_separated_scalar_values"], [-1, -1, -1, -1])
            self.assertEqual(hand["basis_change_physical_residual"], [0, 0])

    def test_deck_character_is_exact_to_roundoff(self) -> None:
        gate = exact_gate()
        self.assertLess(max(gate["deck_character_residuals"].values()), 3e-16)

    def test_tiny_stream_retains_joint_interface(self) -> None:
        rows, analysis = run(4, 2, 1, 0.59274605079, 25033433720260830, 0)
        self.assertEqual(len(rows), 2)
        self.assertEqual(len(FIELD_ORDER), 72)
        self.assertEqual(len(analysis["covariance_of_mean"]), 72)
        self.assertEqual(SMOKE_CAP, 2000)
        self.assertLess(analysis["pair_imaginary_max"], 1e-12)
        self.assertLess(analysis["DFT_conjugacy_max_abs"], 1e-12)

    def test_checked_score_is_reproducible(self) -> None:
        response = json.loads((RESULT / "response_2k.json").read_text())
        batches = read_batches(RESULT / "response_2k.batches.csv")
        self.assertEqual(score(response, batches), json.loads((RESULT / "score_2k.json").read_text()))

    def test_two_separations_pass_without_phase_selection(self) -> None:
        checked = json.loads((RESULT / "score_2k.json").read_text())
        self.assertEqual(checked["usable_separation_count"], 2)
        self.assertTrue(checked["promotion_gate_at_least_two"])
        self.assertTrue(checked["separations"]["1"]["usable_for_promotion"])
        self.assertTrue(checked["separations"]["2"]["usable_for_promotion"])
        self.assertFalse(checked["separations"]["3"]["usable_for_promotion"])
        for row in checked["separations"].values():
            self.assertFalse(row["phase_score_computed"])


if __name__ == "__main__":
    unittest.main()

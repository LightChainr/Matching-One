
from __future__ import annotations
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_p159_pell_hex_filter import build_score  # noqa: E402


RESULTS = ROOT / "results" / "local-20260829" / "P156-square-bond-primitive-pilot"


class P159PellHexFilterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.score = build_score(
            RESULTS / "result.batches.csv", RESULTS / "result.json"
        )

    def test_exact_oracle_and_basis_transport_pass(self):
        self.assertTrue(self.score["gates"]["exact_oracle"]["passed"])
        transport = self.score["basis_transport"]
        self.assertTrue(transport["passed"])
        self.assertEqual(transport["transport_Dminus2_to_Dplus1"], [[1, 0], [0, 1]])
        self.assertEqual(
            transport["ordered_unoriented_line_cycle"],
            ["l0->l1", "l1->l2", "l2->l0"],
        )

    def test_character_detects_signal_but_fails_h4_phase_gate(self):
        gates = self.score["gates"]
        self.assertTrue(gates["nontrivial_character_detection"]["passed"])
        self.assertTrue(gates["reflection_null"]["passed"])
        self.assertTrue(gates["ordinary_H4_simple_zero_phase"]["observed_C_has_same_sign"])
        self.assertTrue(gates["ordinary_H4_simple_zero_phase"]["phase_coordinates_have_opposite_sign"])
        self.assertFalse(gates["ordinary_H4_simple_zero_phase"]["passed"])
        self.assertFalse(gates["H4_specific_filter"]["passed"])

    def test_full_covariance_and_historical_minus2_residual(self):
        covariance = self.score["joint_contrast_covariance_of_mean"]
        self.assertEqual((len(covariance), len(covariance[0])), (6, 6))
        self.assertTrue(all(covariance[i][j] == 0.0
                            for i in range(3) for j in range(3, 6)))
        score = self.score["exploratory_historical_minus2_score"]
        self.assertAlmostEqual(score["observed_ratio_Dminus2_over_Dplus1"], 0.5721444746428873)
        self.assertAlmostEqual(score["residual_z"], 4.960540996712738)
        self.assertFalse(score["passes_abs_z_le_2"])
        self.assertFalse(
            self.score["governance"]["minus2_score_is_preregistered_for_this_C_observable"]
        )


if __name__ == "__main__":
    unittest.main()

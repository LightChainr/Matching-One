from __future__ import annotations

from fractions import Fraction
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from p337_n13_r_odd_charged_source import build_charged_certificate  # noqa: E402


class N13ROddChargedSourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = build_charged_certificate()

    def test_projective_quarter_turn_charges_and_selection_rules(self) -> None:
        representation = self.payload["representation"]
        self.assertIn("charge 2", representation["C4_charge"])
        self.assertIn("B1", representation["D4_refinement"])
        self.assertIn("B2", representation["D4_refinement"])
        self.assertIn("linear H and D responses remain zero",
                      self.payload["selection_rules"]["A_source"])
        self.assertIn("linear H and A responses remain zero",
                      self.payload["selection_rules"]["D_source"])
        matrix = self.payload["linear_response_matrix_at_p_ref"]["exact"]
        self.assertEqual(matrix[0], ["0", "0"])
        self.assertEqual(matrix[1][1], "0")
        self.assertEqual(matrix[2][0], "0")
        self.assertGreater(Fraction(matrix[1][0]), 0)
        self.assertGreater(Fraction(matrix[2][1]), 0)

    def test_f3_phase_defect_activates_both_odd_channels(self) -> None:
        for channel in self.payload["channels"].values():
            self.assertGreater(
                channel["reference_evaluation"]["decimal"]["susceptibility"], 0.0
            )
            for row in channel["state_response_rows"]:
                self.assertEqual(Fraction(row["unweighted_charge_coefficient"]), 0)
                defect = row["F3_phase_defect"]
                self.assertEqual(
                    Fraction(defect["Z_u_minus_1"])
                    + Fraction(defect["Z_u_0"])
                    + Fraction(defect["Z_u_plus_1"]),
                    1,
                )
                self.assertEqual(
                    Fraction(defect["O_u_minus_1"])
                    + Fraction(defect["O_u_0"])
                    + Fraction(defect["O_u_plus_1"]),
                    0,
                )

    def test_charged_current_continuity_is_coefficientwise(self) -> None:
        for channel in self.payload["channels"].values():
            for row in channel["flux_response_rows"]:
                self.assertEqual(Fraction(row["unweighted_birth_charge"]), 0)
                self.assertEqual(Fraction(row["unweighted_exit_charge"]), 0)
                self.assertEqual(
                    Fraction(row["charged_derivative_response"]),
                    Fraction(row["charged_birth_response"])
                    - Fraction(row["charged_exit_response"]),
                )
                self.assertTrue(row["continuity_pass"])

    def test_reference_responses_match_exact_n13_values(self) -> None:
        a = self.payload["channels"]["A_axis_odd"]["reference_evaluation"]["decimal"]
        d = self.payload["channels"]["D_diagonal_odd"]["reference_evaluation"]["decimal"]
        self.assertAlmostEqual(a["susceptibility"], 0.3198939547688)
        self.assertAlmostEqual(d["susceptibility"], 0.020121410529484324)
        self.assertGreater(a["derivative_response"], 0.0)
        self.assertLess(d["derivative_response"], 0.0)
        self.assertTrue(self.payload["gates"]["all_pass"])


if __name__ == "__main__":
    unittest.main()

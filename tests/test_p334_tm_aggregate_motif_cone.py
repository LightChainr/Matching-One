import json
from pathlib import Path
import unittest

from p334_tm_aggregate_motif_cone import build_result


ROOT = Path(__file__).resolve().parents[1]


class P334TMAggregateMotifConeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = build_result()

    def test_four_motif_polynomial_is_exact_TM(self):
        audit = self.result["bounded_audit"]
        self.assertEqual(audit["line_layer_rows"], 984)
        self.assertEqual(audit["four_motif_identity_pass"], 984)
        self.assertEqual(audit["corrected_Rayleigh_pass"], 984)

    def test_both_cover_mechanisms_are_independently_needed(self):
        audit = self.result["bounded_audit"]
        self.assertEqual(audit["ordinary_Rayleigh_fail"], 16)
        self.assertEqual(audit["synergy_only_fail"], 68)
        cone = self.result["mechanism_independence"]
        self.assertEqual(len(cone["Pareto_frontier"]), 9)
        self.assertEqual(len(cone["exact_convex_lower_hull"]), 4)
        self.assertEqual(cone["minimum_corrected_cover_over_hard"], "32/17")

    def test_relative_displacement_and_order_locality_fail(self):
        audit = self.result["bounded_audit"]
        self.assertEqual(audit["relative_pair_tables"], 51_912)
        self.assertEqual(audit["relative_pair_corrected_fail"], 3_900)
        self.assertEqual(audit["displacement_order_tables"], 3_330)
        self.assertEqual(audit["displacement_order_corrected_fail"], 220)
        first = self.result["failed_stronger_routes"][
            "relative_displacement_locality"
        ]["first_failure"]
        self.assertEqual(first["N"], 6)
        self.assertEqual(first["polynomial"], -16)

    def test_spectral_and_complement_shortcuts_fail(self):
        audit = self.result["bounded_audit"]
        self.assertEqual(
            audit["signed_kernel_simple_contrast_negative_rows"], 802
        )
        self.assertEqual(audit["Alexander_reflected_row_pairs"], 492)
        self.assertEqual(audit["Alexander_reflected_margin_equal"], 18)
        self.assertEqual(audit["Alexander_reflected_margin_unequal"], 474)

    def test_checked_artifact_reproduces(self):
        checked = json.loads(
            (ROOT / "results/p334-tm-aggregate-motif-cone/latest.json").read_text()
        )
        self.assertEqual(checked, self.result)


if __name__ == "__main__":
    unittest.main()

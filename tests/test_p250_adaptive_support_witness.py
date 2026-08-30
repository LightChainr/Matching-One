#!/usr/bin/env python3
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import p250_adaptive_support_witness as oracle  # noqa: E402


class AdaptiveSupportWitnessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = oracle.build_result()

    def test_support_involution_is_exact(self):
        gate = self.result["support_involution_gate"]
        self.assertTrue(gate["passed"])
        self.assertEqual(gate["defined_and_undefined_support_pairs_checked"], 746_496)

    def test_every_defined_l3_rectangle_retains_order(self):
        result = self.result["exhaustive_L3"]
        self.assertEqual(result["defined_rectangles"], 8_136)
        self.assertEqual(result["noncommuting_final_fields"], 8_136)
        self.assertEqual(result["nonzero_R_minus"], 8_136)
        self.assertEqual(result["R_minus_histogram"], {"1": 4_320, "2": 3_816})

    def test_first_witness_changes_second_support(self):
        witness = self.result["exhaustive_L3"]["first_exact_witness"]
        supports = witness["supports"]
        self.assertEqual(supports["D0"]["site_coordinate"], [0, 1])
        self.assertEqual(supports["J0"]["site_coordinate"], [2, 1])
        self.assertEqual(supports["J_after_D"]["site_coordinate"], [0, 1])
        self.assertEqual(witness["R_minus"], 1)

    def test_witness_survives_l3_and_l4(self):
        rows = self.result["isometric_witness_orbits"]
        self.assertEqual([row["L"] for row in rows], [3, 4])
        self.assertEqual(
            [row["translation_dihedral_covariance_checks"] for row in rows],
            [72, 128],
        )
        self.assertTrue(all(row["rectangle"]["R_minus"] == 1 for row in rows))

    def test_pilot_is_frozen_without_production(self):
        interface = self.result["frozen_minimal_pilot_interface"]
        self.assertIn("J_after_D_site", interface["branch_specific_supports"])
        self.assertIn("D_after_J_site", interface["branch_specific_supports"])
        self.assertIn("no stochastic production", interface["production_status"])


if __name__ == "__main__":
    unittest.main()

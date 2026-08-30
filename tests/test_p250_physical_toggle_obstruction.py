#!/usr/bin/env python3
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import p250_physical_toggle_obstruction as oracle  # noqa: E402


class PhysicalToggleObstructionTests(unittest.TestCase):
    def test_overwrite_algebra(self):
        result = oracle.overwrite_algebra()
        self.assertGreater(result["same_support_state_checks"], 0)
        self.assertGreater(result["disjoint_support_state_checks"], 0)

    def test_typed_involution(self):
        result = oracle.typed_involution_gate(oracle.axis_integer_torus(3))
        self.assertTrue(result["passed"])
        self.assertEqual(result["leg_checks"], 512 * 9)

    def test_l3_has_symmetric_but_no_order_signal(self):
        result = oracle.exhaustive_marked_triples_l3()
        self.assertEqual(result["eligible_rectangles_checked"], 64_512)
        self.assertTrue(result["R_plus_nonzero_exists"])
        self.assertFalse(result["R_minus_nonzero_exists"])
        self.assertEqual(result["R_minus_histogram"], {"0": 64_512})

    def test_minimal_witness(self):
        witness = oracle.exhaustive_marked_triples_l3()["minimal_R_plus_witness"]
        self.assertEqual(witness["base_occupied_vertices"], [0, 1])
        self.assertEqual(witness["responses"]["L_J"], 1)
        self.assertEqual(witness["R_plus"], 1)
        self.assertEqual(witness["R_minus"], 0)

    def test_decision_stops_before_inventing_semantics(self):
        decision = oracle.build_result()["decision"]
        self.assertFalse(decision["pilot_runner_added"])
        self.assertIn("state", decision["minimum_escape_from_no_go"][0])


if __name__ == "__main__":
    unittest.main()

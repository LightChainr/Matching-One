import json
from pathlib import Path
import tempfile
import unittest

from p334_lorentzian_support_gate import (
    build_result,
    homogeneous_exponent,
    main,
    m_convex_exchange_witness,
    real_root_audit,
)
from fractions import Fraction


ROOT = Path(__file__).resolve().parents[1]


class P334LorentzianSupportGateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = build_result()

    def test_minimal_support_fails_m_convex_exchange(self):
        gate = self.result["lorentzian_support_gate"]
        self.assertEqual(gate["N"], 4)
        self.assertEqual(gate["matrix"], [[2, 0], [0, 2]])
        self.assertEqual(gate["family_masks"], [3, 12])
        witness = gate["witness"]
        self.assertEqual(witness["left_exponent"], [2, 1, 1, 0, 0])
        self.assertEqual(witness["right_exponent"], [2, 0, 0, 1, 1])
        self.assertFalse(
            any(
                row["left_in_support"] and row["right_in_support"]
                for row in witness["candidate_exchanges"]
            )
        )

    def test_exchange_helper_is_exact_on_a_uniform_layer(self):
        family = {mask for mask in range(1 << 4) if mask.bit_count() == 2}
        self.assertIsNone(m_convex_exchange_witness(family, 4))
        self.assertEqual(homogeneous_exponent(3, 4), (2, 1, 1, 0, 0))

    def test_natural_real_rooted_strengthenings_fail_exactly(self):
        normalized = self.result["normalized_q_real_root_gate"]
        self.assertEqual(normalized["N"], 6)
        self.assertEqual(
            normalized["audit"]["coefficients_low_to_high"],
            ["1/5", "3/5", "3/5"],
        )
        self.assertEqual(normalized["audit"]["distinct_real_roots"], 0)
        raw = self.result["raw_count_real_root_gate"]
        self.assertEqual(raw["N"], 8)
        self.assertEqual(raw["audit"]["degree"], 4)
        self.assertEqual(raw["audit"]["distinct_real_roots"], 2)
        self.assertFalse(real_root_audit([Fraction(1), Fraction(3), Fraction(3)])["all_roots_real"])

    def test_normalized_matching_has_an_exact_dead_end(self):
        gate = self.result["normalized_matching_gate"]
        self.assertEqual(gate["N"], 11)
        self.assertEqual(gate["matrix"], [[11, 3], [0, 1]])
        self.assertEqual(gate["lower_layer"], 7)
        self.assertEqual(gate["maximum_flow"], 605)
        self.assertEqual(gate["required_flow"], 726)
        self.assertEqual(gate["dead_end_mask"], 471)
        self.assertEqual(gate["neighbor_masks"], [])
        self.assertTrue(
            all(row["new_mark"][0] == 2 for row in gate["dead_end_additions"])
        )

    def test_checked_in_artifact_reproduces(self):
        checked = json.loads(
            (ROOT / "results/p334-lorentzian-support-gate/latest.json").read_text()
        )
        self.assertEqual(checked, self.result)


if __name__ == "__main__":
    unittest.main()

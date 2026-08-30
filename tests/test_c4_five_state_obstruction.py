from __future__ import annotations

import copy
import json
import sys
import unittest
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import c4_five_state_obstruction as oracle  # noqa: E402


class C4FiveStateObstructionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.artifact = oracle.build_certificate()

    def test_all_exact_model_checks(self) -> None:
        self.assertTrue(all(self.artifact["exact_checks"].values()))
        self.assertEqual(self.artifact["passed_exact_checks"], 20)

    def test_orbit_arithmetic_for_reduced_five_points(self) -> None:
        self.assertEqual(self.artifact["reduced_support_orbit_counts_1_2_4"], [(1, 0, 1)])

    def test_reduced_character(self) -> None:
        self.assertEqual(self.artifact["reduced_control"]["character_R_powers"], [5, 1, 1, 1])
        self.assertEqual(self.artifact["reduced_control"]["multiplicities_1_i_minus1_minusi"], [2, 1, 1, 1])

    def test_nonreduced_character_is_separated(self) -> None:
        self.assertEqual(self.artifact["nonreduced_control"]["character_R_powers"], [5, -1, 1, -1])
        self.assertEqual(self.artifact["nonreduced_control"]["multiplicities_1_i_minus1_minusi"], [1, 1, 2, 1])

    def test_trace_character_multiplicities_sum_to_dimension(self) -> None:
        for name in ("reduced_control", "nonreduced_control"):
            m0, m1, m2, m3 = self.artifact[name]["multiplicities_1_i_minus1_minusi"]
            self.assertEqual(m0 + m1 + m2 + m3, 5)
            self.assertEqual(m1, m3)
            self.assertEqual([5, m0 - m2, m0 - m1 + m2 - m3, m0 - m2],
                             self.artifact[name]["character_R_powers"])

    def test_nonreduced_generator_really_has_nilpotency_order_three(self) -> None:
        x = oracle.mat(self.artifact["nonreduced_control"]["log_U"])
        zero = oracle.add((F(0), oracle.eye()))
        self.assertNotEqual(oracle.power(x, 2), zero)
        self.assertEqual(oracle.power(x, 3), zero)

    def test_input_contract(self) -> None:
        with self.assertRaises(ValueError):
            oracle.mat(((1, 0), (0, 1)))
        bad = [[0] * 5 for _ in range(5)]
        bad[0][0] = 1.0
        with self.assertRaises(TypeError):
            oracle.mat(bad)

    def test_separate_verifier(self) -> None:
        self.assertEqual(oracle.verify_certificate(self.artifact)["status"], "verified_exactly")

    def test_corrupted_rotation_is_rejected(self) -> None:
        bad = copy.deepcopy(self.artifact)
        bad["nonreduced_control"]["R"][0][0] = "2"
        with self.assertRaises(ValueError):
            oracle.verify_certificate(bad)

    def test_corrupted_character_is_rejected(self) -> None:
        bad = copy.deepcopy(self.artifact)
        bad["nonreduced_control"]["character_R_powers"] = [5, 1, 1, 1]
        with self.assertRaises(ValueError):
            oracle.verify_certificate(bad)

    def test_corrupted_translation_is_rejected(self) -> None:
        bad = copy.deepcopy(self.artifact)
        bad["nonreduced_control"]["U"][0][0] = "2"
        with self.assertRaises(ValueError):
            oracle.verify_certificate(bad)

    def test_missing_reduced_point_is_rejected(self) -> None:
        bad = copy.deepcopy(self.artifact)
        bad["reduced_control"]["joint_points"].pop()
        with self.assertRaises(ValueError):
            oracle.verify_certificate(bad)

    def test_checked_in_certificate_reproduces(self) -> None:
        stored = json.loads((ROOT / "results/jordan-nonseparation/c4-five-state.json").read_text())
        self.assertEqual(stored, json.loads(json.dumps(self.artifact)))


if __name__ == "__main__":
    unittest.main()

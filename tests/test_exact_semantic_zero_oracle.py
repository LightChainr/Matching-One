import copy
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from exact_semantic_zero_oracle import DEFAULT_OUTPUT, build_result, validate_result  # noqa: E402


class ExactSemanticZeroOracleTests(unittest.TestCase):
    def test_primitive_bezout_witness_is_exactly_one(self) -> None:
        witness = build_result()["primitive_bezout_infeasibility_witness"]
        self.assertEqual(witness["multipliers"], ["7", "-1"])
        self.assertEqual(witness["result_coefficients_constant_first"], ["1", "0"])
        self.assertEqual(witness["gcd_of_integer_multipliers"], 1)
        self.assertEqual(witness["status"], "exactly_infeasible")

    def test_semantic_zero_is_a_preoptimization_hard_gate(self) -> None:
        gate = build_result()["gate_result"]
        self.assertIs(gate["semantic_zero_enforced_as_hard_equality"], True)
        self.assertIs(gate["contradictory_row_caught_before_optimization"], True)
        self.assertIs(gate["solver_invoked"], False)

    def test_claim_boundary_does_not_promote_the_synthetic_zero(self) -> None:
        excluded = build_result()["claim_boundary"]["excluded"]
        self.assertIn("physical selection rule", excluded)
        self.assertIn("model validation", excluded)

    def test_checked_in_certificate_reproduces_and_tampering_fails_closed(self) -> None:
        checked = json.loads(DEFAULT_OUTPUT.read_text(encoding="utf-8"))
        self.assertEqual(checked, build_result())
        self.assertIs(validate_result(checked)["solver_invoked"], False)
        tampered = copy.deepcopy(checked)
        tampered["primitive_bezout_infeasibility_witness"]["result_coefficients_constant_first"] = ["0", "0"]
        with self.assertRaisesRegex(ValueError, "does not exactly reproduce"):
            validate_result(tampered)


if __name__ == "__main__":
    unittest.main()

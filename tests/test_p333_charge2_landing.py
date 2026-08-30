import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from p333_charge2_landing import (  # noqa: E402
    build_result,
    charge2_character,
    landing_reference_state,
)


RESULT = ROOT / "results/p333-charge2-landing/latest.json"


class TestP333Charge2Landing(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = json.loads(RESULT.read_text(encoding="utf-8"))

    def test_frozen_charge_two_orbit(self):
        character = [charge2_character(site) for site in range(4)]
        self.assertEqual(character, [1, -1, 1, -1])
        self.assertEqual(character[1:] + character[:1], [-value for value in character])

    def test_landing_reference_states(self):
        self.assertEqual(
            [landing_reference_state(site) for site in range(4)],
            [(0, 1, 1, 1), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)],
        )

    def test_non_scalar_and_translation_gates(self):
        self.assertTrue(self.result["non_scalar_gate"]["passes"])
        self.assertEqual(self.result["non_scalar_gate"]["new_non_scalar_rank"], 1)
        self.assertEqual(self.result["translation_covariance"]["G0_residual_rank"], 0)
        self.assertEqual(self.result["translation_covariance"]["G1_residual_rank"], 0)

    def test_exact_dimension_ladder_stops_at_gram(self):
        stages = self.result["stages"]
        self.assertEqual(
            [
                stages[name]["affine_tangent_dimension"]
                for name in (
                    "affine_q_jet",
                    "endpoint_radical_normalized",
                    "gram_self_adjoint",
                    "source_landing_normalized",
                )
            ],
            [3, 1, None, None],
        )
        self.assertEqual(self.result["canonical_restricted_gram_skew_rank"], 4)

    def test_first_empty_restriction_has_exact_left_null_witness(self):
        gate = self.result["first_empty_restriction"]
        self.assertEqual(gate["from_stage"], "endpoint_radical_normalized")
        self.assertEqual(gate["added_stage"], "gram_self_adjoint")
        self.assertEqual(gate["coefficient_rank"], 0)
        self.assertEqual(gate["augmented_rank"], 1)
        witness = gate["inconsistency_witness"]
        self.assertTrue(all(value == 0 for value in witness["left_times_parameter_matrix"]))
        self.assertEqual(witness["left_times_rhs"], 1)

    def test_failure_exhausts_terminal_landing_irreps_only(self):
        self.assertEqual(self.result["decision"], "charge2_fails")
        exhaustion = self.result["irrep_exhaustion_if_failure"]
        self.assertIn("e7e6c80", exhaustion["scalar_endpoint_irreps"])
        self.assertIn("7b40ec7", exhaustion["C4_charge1_landing_irrep"])
        self.assertIn("fails here", exhaustion["C4_charge2_landing_irrep"])
        self.assertIn("Retain rooted/landing connectivity", exhaustion["next_semantic_change"])

    def test_full_certificate_recomputes(self):
        self.assertEqual(build_result(), self.result)


if __name__ == "__main__":
    unittest.main()

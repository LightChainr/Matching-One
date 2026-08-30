from fractions import Fraction
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from noncrossing_connectivity_codec import noncrossing_states  # noqa: E402
from p321_homology_trace_certificate import (  # noqa: E402
    action_matrix,
    matrix_multiply,
    rotate_state,
)
from p333_generic_q_detach_intertwiner import (  # noqa: E402
    build_result,
    detach_jet,
    detach_state,
)


RESULT = ROOT / "results/p333-generic-q-detach-intertwiner/latest.json"


class TestP333GenericQDetachIntertwiner(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = json.loads(RESULT.read_text(encoding="utf-8"))

    def test_detach_state(self):
        self.assertEqual(detach_state((0, 0, 1), 0), (0, 1, 2))
        self.assertEqual(detach_state((0, 0, 1), 2), (0, 0, 1))
        self.assertEqual(detach_state((0, 1, 0, 1), 1), (0, 1, 0, 2))

    def test_detach_polynomial_has_frozen_loop_weight(self):
        states = noncrossing_states(3)
        constant, velocity = detach_jet(3, 0)
        for column, state in enumerate(states):
            is_singleton = state.count(state[0]) == 1
            self.assertEqual(sum(row[column] for row in constant), 1)
            self.assertEqual(sum(row[column] for row in velocity), int(is_singleton))
            if is_singleton:
                self.assertEqual(
                    [constant[row][column] for row in range(len(states))],
                    [velocity[row][column] for row in range(len(states))],
                )

    def test_translation_covariance_for_both_q_coefficients(self):
        for width in (2, 3, 4):
            translation = action_matrix(width, lambda state: rotate_state(state, 1))
            for site in range(width):
                constant, velocity = detach_jet(width, site)
                shifted_constant, shifted_velocity = detach_jet(width, (site + 1) % width)
                self.assertEqual(
                    matrix_multiply(translation, constant),
                    matrix_multiply(shifted_constant, translation),
                )
                self.assertEqual(
                    matrix_multiply(translation, velocity),
                    matrix_multiply(shifted_velocity, translation),
                )

    def test_exact_dimension_ladder(self):
        rows = self.result["widths"]
        self.assertEqual([row["width"] for row in rows], [2, 3, 4])
        self.assertEqual(
            [row["stages"]["generic_q_affine_hom_jet"]["affine_tangent_dimension"] for row in rows],
            [2, 2, 2],
        )
        self.assertEqual(
            [row["stages"]["endpoint_normalized"]["affine_tangent_dimension"] for row in rows],
            [0, 0, 0],
        )
        self.assertEqual(
            [row["stages"]["gram_radical_self_adjoint"]["consistent"] for row in rows],
            [True, False, False],
        )

    def test_endpoint_selects_canonical_full_q_translation(self):
        self.assertTrue(
            self.result["endpoint_uniquely_selects_canonical_all_widths"]
        )
        for row in self.result["widths"]:
            endpoint = row["stages"]["endpoint_normalized"]
            self.assertEqual(endpoint["affine_tangent_dimension"], 0)
            self.assertEqual(endpoint["canonical_X0_T_V0_nonzero_equations"], 0)
            self.assertTrue(row["endpoint_uniquely_selects_X0_T_V0"])

    def test_nondegenerate_gram_no_go_is_zeroth_order(self):
        for row in self.result["widths"][1:]:
            self.assertTrue(row["zeroth_order_projection"]["inherited_no_go"])
            witness = row["gram_restriction_on_endpoint_moduli"]["inconsistency_witness"]
            self.assertEqual(witness["left_times_rhs"], 1)
            self.assertEqual(witness["left_times_parameter_matrix"], [])

    def test_width_two_has_zero_velocity_overlap(self):
        row = self.result["widths"][0]
        self.assertEqual(row["decision"], "reopened_unique")
        overlap = row["velocity_radical_overlap"]
        self.assertEqual(overlap["particular_rank"], 0)
        self.assertEqual(overlap["particular_induced_radical_velocity"], [[0]])
        self.assertTrue(overlap["all_final_tangents_have_zero_radical_velocity"])

    def test_secondary_full_gram_cannot_reopen(self):
        self.assertEqual(
            [row["secondary_full_gram_source"]["consistent"] for row in self.result["widths"]],
            [True, False, False],
        )

    def test_full_certificate_recomputes(self):
        self.assertEqual(build_result(), self.result)


if __name__ == "__main__":
    unittest.main()

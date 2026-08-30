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
from p333_one_mark_endpoint_jet import (  # noqa: E402
    block_diagonal_with_mark,
    build_result,
    extended_gram_jet,
    marked_detach_jet,
    nullspace_basis,
    singleton_covector,
)


RESULT = ROOT / "results/p333-one-mark-endpoint-jet/latest.json"


class TestP333OneMarkEndpointJet(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = json.loads(RESULT.read_text(encoding="utf-8"))

    def test_singleton_covector(self):
        states = noncrossing_states(3)
        values = singleton_covector(states, 0)
        self.assertEqual(len(values), 5)
        self.assertEqual(values[states.index((0, 1, 2))], 1)
        self.assertEqual(values[states.index((0, 0, 0))], 0)

    def test_marked_detach_has_terminal_emission(self):
        states = noncrossing_states(3)
        zero, velocity = marked_detach_jet(3, 0)
        size = len(states) + 1
        self.assertEqual(zero[-1][-1], 1)
        self.assertEqual(velocity[-1][-1], 0)
        for column, state in enumerate(states):
            expected = int(state.count(state[0]) == 1)
            self.assertEqual(velocity[-1][column], expected)
        self.assertTrue(all(zero[row][-1] == 0 for row in range(size - 1)))

    def test_marked_detach_translation_covariance(self):
        for width in (3, 4):
            translation = block_diagonal_with_mark(
                action_matrix(width, lambda state: rotate_state(state, 1))
            )
            for site in range(width):
                zero, velocity = marked_detach_jet(width, site)
                next_zero, next_velocity = marked_detach_jet(
                    width, (site + 1) % width
                )
                self.assertEqual(
                    matrix_multiply(translation, zero),
                    matrix_multiply(next_zero, translation),
                )
                self.assertEqual(
                    matrix_multiply(translation, velocity),
                    matrix_multiply(next_velocity, translation),
                )

    def test_extended_gram_rank_and_radical(self):
        expected = {3: (3, 3), 4: (3, 12)}
        for width, (rank, radical_dimension) in expected.items():
            states = noncrossing_states(width)
            g0, _ = extended_gram_jet(states)
            basis = nullspace_basis(g0)
            self.assertEqual(len(basis[0]), radical_dimension)
            for row in range(len(g0)):
                for column in range(radical_dimension):
                    self.assertEqual(
                        sum(g0[row][pivot] * basis[pivot][column] for pivot in range(len(g0))),
                        Fraction(0),
                    )
            self.assertEqual(len(g0) - radical_dimension, rank)

    def test_dimension_table_and_first_empty_gate(self):
        rows = self.result["widths"]
        self.assertEqual([row["width"] for row in rows], [3, 4])
        for row in rows:
            stages = row["stages"]
            self.assertEqual(stages["marked_affine_hom_jet"]["affine_tangent_dimension"], 5)
            self.assertEqual(stages["endpoint_radical_normalized"]["affine_tangent_dimension"], 3)
            self.assertFalse(stages["gram_self_adjoint"]["consistent"])
            witness = row["first_empty_restriction"]
            self.assertEqual(witness["added_stage"], "gram_self_adjoint")
            self.assertEqual(witness["coefficient_rank"], 0)
            self.assertEqual(witness["augmented_rank"], 1)
            self.assertEqual(witness["inconsistency_witness"]["left_times_rhs"], 1)

    def test_q_velocity_really_reaches_mark_from_radical(self):
        self.assertTrue(self.result["genuine_q_velocity_gate_passed"])
        for row in self.result["widths"]:
            for site in row["mark_velocity_gate"]["sites"]:
                self.assertEqual(site["mark_injection_rank_on_radical"], 1)

    def test_smallest_counterexample_and_lower_bound(self):
        width3, width4 = self.result["widths"]
        self.assertEqual(width3["extended_dimension"], 6)
        self.assertEqual(width3["decision"], "one_mark_insufficient")
        self.assertEqual(width3["canonical_translation"]["restricted_gram_skew_rank"], 2)
        self.assertEqual(width3["scalar_mark_lower_bound_if_canonical"]["total_marks_including_current"], 2)
        self.assertEqual(width4["canonical_translation"]["restricted_gram_skew_rank"], 4)
        self.assertEqual(width4["scalar_mark_lower_bound_if_canonical"]["total_marks_including_current"], 3)

    def test_full_certificate_recomputes(self):
        self.assertEqual(build_result(), self.result)


if __name__ == "__main__":
    unittest.main()

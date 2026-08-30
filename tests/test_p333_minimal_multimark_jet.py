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
from p333_minimal_multimark_jet import (  # noqa: E402
    block_diagonal_with_marks,
    build_result,
    extended_multimark_gram_jet,
    falling_factorial,
    multimark_detach_jet,
)
from p333_one_mark_endpoint_jet import nullspace_basis  # noqa: E402


RESULT = ROOT / "results/p333-minimal-multimark-jet/latest.json"


class TestP333MinimalMultimarkJet(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = json.loads(RESULT.read_text(encoding="utf-8"))

    def test_falling_factorials(self):
        self.assertEqual([falling_factorial(4, r) for r in range(5)], [1, 4, 12, 24, 24])
        self.assertEqual(falling_factorial(2, 3), 0)

    def test_frozen_case_gram_ranks(self):
        for width, marks, expected_rank, expected_radical in (
            (3, 2, 5, 2),
            (4, 3, 7, 10),
        ):
            states = noncrossing_states(width)
            g0, _ = extended_multimark_gram_jet(states, marks)
            basis = nullspace_basis(g0)
            self.assertEqual(len(g0) - len(basis[0]), expected_rank)
            self.assertEqual(len(basis[0]), expected_radical)

    def test_detach_emission_is_frozen_falling_factorial(self):
        states = noncrossing_states(4)
        _, velocity = multimark_detach_jet(4, 0, 3)
        ordinary = len(states)
        singleton = states.index((0, 1, 2, 3))
        connected = states.index((0, 0, 0, 0))
        self.assertEqual(
            [velocity[ordinary + mark][singleton] for mark in range(3)],
            [1, 3, 6],
        )
        self.assertEqual(
            [velocity[ordinary + mark][connected] for mark in range(3)],
            [0, 0, 0],
        )

    def test_multimark_translation_covariance(self):
        for width, marks in ((3, 2), (4, 3)):
            translation = block_diagonal_with_marks(
                action_matrix(width, lambda state: rotate_state(state, 1)), marks
            )
            for site in range(width):
                zero, velocity = multimark_detach_jet(width, site, marks)
                next_zero, next_velocity = multimark_detach_jet(
                    width, (site + 1) % width, marks
                )
                self.assertEqual(
                    matrix_multiply(translation, zero),
                    matrix_multiply(next_zero, translation),
                )
                self.assertEqual(
                    matrix_multiply(translation, velocity),
                    matrix_multiply(next_velocity, translation),
                )

    def test_exact_dimension_ladders_and_witnesses(self):
        cases = self.result["cases"]
        expected = {(3, 2): (10, 8), (4, 3): (17, 15)}
        for row in cases:
            affine, endpoint = expected[(row["width"], row["marks"])]
            self.assertEqual(row["stages"]["affine_q_jet"]["affine_tangent_dimension"], affine)
            self.assertEqual(row["stages"]["endpoint_radical_normalized"]["affine_tangent_dimension"], endpoint)
            self.assertFalse(row["stages"]["gram_self_adjoint"]["consistent"])
            witness = row["first_empty_restriction"]
            self.assertEqual(witness["added_stage"], "gram_self_adjoint")
            self.assertEqual(witness["coefficient_rank"], 0)
            self.assertEqual(witness["augmented_rank"], 1)
            self.assertEqual(witness["inconsistency_witness"]["left_times_rhs"], 1)

    def test_raised_bounds_and_width_four_exhaustion(self):
        width3, width4 = self.result["cases"]
        self.assertEqual(width3["raised_bound_if_empty"]["total_scalar_marks"], 3)
        self.assertFalse(width3["raised_bound_if_empty"]["exceeds_available_independent_marks"])
        self.assertEqual(width4["raised_bound_if_empty"]["total_scalar_marks"], 5)
        self.assertEqual(width4["raised_bound_if_empty"]["maximum_nonzero_falling_factorial_marks"], 4)
        self.assertTrue(width4["raised_bound_if_empty"]["exceeds_available_independent_marks"])

    def test_full_certificate_recomputes(self):
        self.assertEqual(build_result(), self.result)


if __name__ == "__main__":
    unittest.main()

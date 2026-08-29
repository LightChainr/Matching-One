#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from exact_vjs_collision_tangent import render  # noqa: E402


RESULT = render()


class ExactVJSCollisionTangentTests(unittest.TestCase):
    def test_tiny_probabilities_and_partition_are_exact(self) -> None:
        probabilities = RESULT["vjs_geometric_probabilities"]
        self.assertEqual(probabilities["pair0_distinct"]["Q1"], "37/256")
        self.assertEqual(probabilities["P0_four_clusters"]["Q1"], "1/256")
        self.assertEqual(probabilities["P1_one_propagating_cluster"]["Q1"], "3/128")
        self.assertEqual(probabilities["P2_two_propagating_clusters"]["Q1"], "9/256")
        self.assertTrue(RESULT["exact_partition_check"]["equal"])

    def test_measure_projector_and_explicit_terms_close_exactly(self) -> None:
        tensor = RESULT["vjs_centered_two_point_tensor"]
        decomposition = tensor["derivative_decomposition"]
        self.assertTrue(tensor["all_three_contributions_nonzero"])
        self.assertTrue(decomposition["sum_equals_direct"])
        self.assertEqual(
            decomposition["finite_confluent_projector_derivative"],
            {"I_pair": "0/1", "X_shared_colour": "9/128", "J_all_ones": "0/1"},
        )

    def test_issue262_regular_block_is_used_not_singular_projectors(self) -> None:
        gate = RESULT["issue262_reuse_gate"]
        self.assertEqual(gate["P_regular_derivative_at_Q1"], "d_Q(P_singlet+P_[2])=X")
        self.assertFalse(gate["individual_projector_derivatives_used"])

    def test_scaling_collision_is_not_inferred_from_one_distance(self) -> None:
        boundary = RESULT["collision_boundary"]
        self.assertIn(
            "the logarithmic coefficient 2sqrt(3)/pi from a one-distance L2 oracle",
            boundary["not_determined_here"],
        )


if __name__ == "__main__":
    unittest.main()

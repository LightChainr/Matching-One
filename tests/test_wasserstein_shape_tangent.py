#!/usr/bin/env python3
"""Lock the affine-null and rank-recovery controls of the Wasserstein shape tangent."""

from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import wasserstein_shape_tangent as shape  # noqa: E402


class AffineTangent(unittest.TestCase):
    def test_the_tangent_is_the_two_directions_an_affine_map_can_move_along(self) -> None:
        """Stops us believing a tangent that is not span{1, Q}.

        If the basis drifts -- to span{1, u}, say, or to Q alone -- then a pure
        width change stops being projected out and every shape norm reported
        afterwards silently contains it.
        """
        base = [3.0, 5.0, 11.0, 17.0]
        basis = shape.affine_tangent(base)
        self.assertEqual(basis, [[1.0, 1.0, 1.0, 1.0], base])

    def test_the_base_quantile_is_not_itself_affine_in_u(self) -> None:
        """Stops the controls passing for a reason that has nothing to do with the projection.

        A base proportional to u would make span{1, Q} equal span{1, u}, and a
        deformation written in u would land inside the tangent by accident.  The
        wrong result this catches is an affine-null control that is null because
        the base was degenerate.
        """
        levels = shape.quantile_levels()
        base = shape._base_quantile(levels)
        size = len(levels)
        residuals = shape.subspace_residual(
            base, [[1.0 if i == j else 0.0 for j in range(size)] for i in range(size)],
            [[1.0] * size, list(levels)])
        self.assertGreater(residuals["statistic"], 1e-3)

    def test_every_affine_displacement_has_exactly_zero_shape_flow(self) -> None:
        """Stops us believing a shape flow that is a center or width change.

        This is the control the whole method rests on.  A nonzero statistic on a
        pure translation or a pure rescaling would mean every later shape norm
        is contaminated by the drift it was supposed to remove.
        """
        for row in shape.affine_null_controls(shape.quantile_levels()):
            self.assertLess(row["statistic"], 1e-20, row["case"])
            self.assertLess(row["unweighted_shape_norm"], 1e-12, row["case"])

    def test_a_declared_deformation_is_seen_and_its_generator_closes_it(self) -> None:
        """Stops us believing a projector that annihilates everything.

        A statistic that is zero on affine displacements is only half of the
        test; it must be nonzero on a deformation genuinely outside the tangent,
        and must return to zero once that deformation joins the basis.  The
        wrong number here is a shape norm of zero on data that has shape.
        """
        rows = {row["case"]: row for row in
                shape.rank_recovery_control(shape.quantile_levels())}
        one = rows["one_shape_direction"]
        self.assertGreater(one["statistic_affine_basis_only"], 1e-3)
        self.assertLess(one["statistic_after_one_generator"], 1e-20)

    def test_one_generator_does_not_close_a_rank_two_deformation(self) -> None:
        """Stops us believing a rank-1 shape flow where the flow is rank 2.

        The decision table of #582 turns on exactly this distinction, so the
        control has to show that adding the first generator leaves a residual
        and adding the second removes it.  A method that closed the rank-2 case
        with one direction would report every flow as rank 1.
        """
        rows = {row["case"]: row for row in
                shape.rank_recovery_control(shape.quantile_levels())}
        two = rows["two_shape_directions"]
        self.assertGreater(two["statistic_after_one_generator"], 1e-3)
        self.assertLess(two["statistic_after_all_generators"], 1e-20)
        self.assertLess(two["statistic_after_one_generator"],
                        two["statistic_affine_basis_only"])

    def test_the_weighted_statistic_is_not_carried_by_the_noisiest_quantile(self) -> None:
        """Stops us reading a shape rank off a norm that a thin tail supplied.

        #582 proposes an unweighted L2 shape norm and then has to freeze a
        central window u_min to keep the tails from dominating it.  This control
        exhibits the case: an affine displacement with a one-sigma excursion in
        the noisiest coordinate.  The wrong number is the unweighted norm's
        verdict -- a large shape flow -- taken as evidence of shape.
        """
        row = shape.tail_dominance_control(shape.quantile_levels())
        self.assertAlmostEqual(row["tail_excursion_in_units_of_its_own_sigma"], 1.0, places=9)
        self.assertGreater(row["unweighted_shape_norm"], 100.0)
        self.assertLess(row["weighted_statistic"], row["degrees_of_freedom"])
        self.assertLess(row["weighted_equivalent_sigma"], 1.0)

    def test_the_degrees_of_freedom_are_the_levels_minus_the_affine_tangent(self) -> None:
        """Stops us believing a chi-square scored against the wrong df.

        Two columns of tangent are removed, so nine levels leave seven.  A df
        that silently grew would make every shape flow look more significant
        than it is.
        """
        for count in (5, 9, 13):
            rows = shape.affine_null_controls(shape.quantile_levels(count))
            for row in rows:
                self.assertEqual(row["degrees_of_freedom"], count - 2, (count, row["case"]))

    def test_three_levels_is_the_floor(self) -> None:
        """Stops us building a shape statistic with no room for a shape.

        Two levels are exactly spanned by the affine tangent, so the residual is
        identically zero and no deformation could ever be detected.
        """
        with self.assertRaises(ValueError):
            shape.quantile_levels(2)
        self.assertEqual(len(shape.quantile_levels(3)), 3)

    def test_the_assembled_report_refuses_the_claim_it_could_be_read_as_making(self) -> None:
        """Stops the synthetic controls being cited as a result about a threshold law."""
        payload = shape.assemble()
        joined = " ".join(payload["not_established"])
        self.assertIn("synthetic", joined)
        self.assertIn("percolation", joined)


if __name__ == "__main__":
    unittest.main()

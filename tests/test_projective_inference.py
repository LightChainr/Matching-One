#!/usr/bin/env python3
"""Tests for denominator-free model comparison.

The wrong numbers these exist to stop us believing are the ones the N=580
ladder produced: verdicts computed by dividing a well-measured response by a
poorly-measured one.  The anchor test below is the one that matters -- it pins
this module to Fieller's test, which is already known to be right in the two
entry case, so the generalisation cannot quietly be a different statistic.
"""

from __future__ import annotations

import math
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import projective_inference as projective  # noqa: E402


class ProjectiveInferenceTests(unittest.TestCase):
    # The committed N=580 run: A4 at r = 1, 2, 4 with the measured covariance
    # of #575.  Real numbers, because a synthetic covariance would not have the
    # feature that caused the problem -- a first entry 3.6 sigma from zero.
    Y = [9.01643304e-04, 2.91097703e-03, 4.13180763e-03]
    VAR = [6.207021e-08, 6.713835e-08, 3.807515e-08]
    RHO_14 = -0.1648

    def _two_entry_covariance(self):
        s0 = math.sqrt(self.VAR[0])
        s2 = math.sqrt(self.VAR[2])
        c = self.RHO_14 * s0 * s2
        return [[self.VAR[0], c], [c, self.VAR[2]]]

    def test_one_ray_in_two_entries_is_exactly_fieller_squared(self) -> None:
        """The anchor. Fieller is the 2x1 case, so this must reproduce it.

        The wrong number this stops us believing is a generalisation that has
        drifted into being a different test.  Fieller's z on the ratio of two
        entries is already established; if the covariance-weighted distance to
        the model ray is not its square, one of the two is wrong, and this test
        does not say which -- it says stop.
        """
        covariance = self._two_entry_covariance()
        observed = [self.Y[0], self.Y[2]]
        for predicted in (1.0, 4.0, 10.9908008589, 16.0, 120.79770352, 2080.30719731):
            fieller = (
                (self.Y[2] - predicted * self.Y[0])
                / math.sqrt(self.VAR[2] + predicted ** 2 * self.VAR[0]
                            - 2 * predicted * covariance[0][1])
            )
            result = projective.ray_residual(observed, covariance, [1.0, predicted])
            self.assertEqual(result["degrees_of_freedom"], 1)
            self.assertAlmostEqual(
                result["statistic"], fieller ** 2,
                delta=abs(fieller ** 2) * 1e-9 + 1e-12,
            )

    def test_no_denominator_is_ever_chosen(self) -> None:
        """Scaling the model direction must not move the statistic.

        The wrong number here is one that depends on which entry was written
        first.  A ray is a ray: (1, R) and (c, cR) are the same model, and a
        test that distinguished them would be testing our bookkeeping.
        """
        covariance = self._two_entry_covariance()
        observed = [self.Y[0], self.Y[2]]
        base = projective.ray_residual(observed, covariance, [1.0, 4.0])["statistic"]
        for scale in (1e-6, 0.5, 7.0, 1e6, -3.0):
            scaled = projective.ray_residual(
                observed, covariance, [scale, 4.0 * scale]
            )["statistic"]
            self.assertAlmostEqual(scaled, base, delta=abs(base) * 1e-9 + 1e-12)

    def test_a_model_through_the_data_has_no_residual(self) -> None:
        """A ray that passes through the observation must score zero.

        The wrong number this catches is a statistic that has picked up an
        offset -- from a mis-signed cross term, say -- which would inflate every
        model equally and be invisible in a comparison between models.
        """
        covariance = self._two_entry_covariance()
        observed = [self.Y[0], self.Y[2]]
        exact = projective.ray_residual(
            observed, covariance, [self.Y[0], self.Y[2]]
        )
        self.assertLess(exact["statistic"], 1e-20)
        self.assertGreater(exact["p_value"], 0.999)

    def test_a_nuisance_direction_costs_one_degree_of_freedom(self) -> None:
        """Carrying a known systematic must be cheaper than discarding an entry.

        The frozen ladder dropped the r=2 rung from its decision because that
        rung carries spin-8 leakage. The wrong belief this test guards is that
        dropping was the only option: a systematic with a known direction is a
        second basis vector, which costs one degree of freedom rather than a
        whole measurement.
        """
        s = [math.sqrt(v) for v in self.VAR]
        c14 = self.RHO_14 * s[0] * s[2]
        c12 = 0.0209 * s[0] * s[1]
        covariance = [[self.VAR[0], c12, c14],
                      [c12, self.VAR[1], 0.0],
                      [c14, 0.0, self.VAR[2]]]
        ray = projective.subspace_residual(self.Y, covariance, [[1.0, 2.0, 4.0]])
        with_nuisance = projective.subspace_residual(
            self.Y, covariance, [[1.0, 2.0, 4.0], [-1.0, 1.0, -1.0]]
        )
        self.assertEqual(ray["degrees_of_freedom"], 2)
        self.assertEqual(with_nuisance["degrees_of_freedom"], 1)
        # Absorbing a direction can only reduce the distance to the model set.
        self.assertLessEqual(with_nuisance["statistic"], ray["statistic"] + 1e-15)

    def test_dependent_model_directions_are_refused_rather_than_fitted(self) -> None:
        """Stops a model set that is secretly smaller than it looks.

        Two directions that span one line give a singular normal-equation
        matrix.  Solving it anyway returns an arbitrary split of one amplitude
        between two coefficients, and a degrees-of-freedom count that is one too
        low -- which makes every p-value too small.
        """
        covariance = self._two_entry_covariance()
        with self.assertRaisesRegex(ValueError, "linearly dependent"):
            projective.subspace_residual(
                [self.Y[0], self.Y[2]], covariance, [[1.0, 4.0], [2.0, 8.0]]
            )

    def test_a_null_direction_of_the_covariance_is_dropped_not_inverted(self) -> None:
        """Stops manufactured certainty along a direction with no information.

        A jackknife over 100 batches cannot support more than 99 directions, so
        a measured covariance can be singular by construction.  Inverting a zero
        eigenvalue would make the residual in that direction infinitely
        significant. The rank must fall instead, and the degrees of freedom with
        it.
        """
        covariance = [[1.0, 1.0, 0.0], [1.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        result = projective.ray_residual([1.0, 1.0, 2.0], covariance, [1.0, 1.0, 1.0])
        self.assertEqual(result["covariance_rank"], 2)
        self.assertEqual(result["degrees_of_freedom"], 1)
        self.assertTrue(math.isfinite(result["statistic"]))

    def test_the_chi_square_tail_matches_values_that_can_be_checked_by_hand(self) -> None:
        """The p-values are computed, not tabulated, so check them once.

        At one degree of freedom the tail is the two-sided normal tail, which is
        known in closed form; at two it is a plain exponential.
        """
        for sigma in (1.0, 2.0, 3.0, 5.0):
            expected = math.erfc(sigma / math.sqrt(2.0))
            self.assertAlmostEqual(
                projective.chi_square_upper_tail(sigma ** 2, 1) / expected, 1.0, places=9
            )
        for statistic in (0.5, 2.0, 9.0):
            self.assertAlmostEqual(
                projective.chi_square_upper_tail(statistic, 2),
                math.exp(-statistic / 2.0), places=12,
            )

    def test_equivalent_sigma_is_a_tail_not_a_square_root(self) -> None:
        """Stops sqrt(chi-square) being quoted as sigma at more than one df.

        The wrong number here is the one everybody writes: sqrt(D) reported as
        a sigma when D has two or more degrees of freedom, which overstates the
        significance because a larger D is expected under the null.
        """
        self.assertAlmostEqual(projective._equivalent_sigma(9.0, 1), 3.0, places=9)
        two = projective._equivalent_sigma(9.0, 2)
        self.assertLess(two, 3.0)
        self.assertGreater(two, 2.0)


if __name__ == "__main__":
    unittest.main()

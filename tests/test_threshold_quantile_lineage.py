#!/usr/bin/env python3
"""Lock the fast threshold-quantile path against the exact one, and the spin-4 weights."""

from __future__ import annotations

import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import threshold_quantile_lineage as lineage  # noqa: E402


def _synthetic(n: int) -> tuple[dict[int, int], dict[int, int]]:
    """A histogram with support spread over many ranks, and matched sample counts."""
    minus = {k: (k % 5) + 1 for k in range(1, n // 2 + 1)}
    plus = {k: (k % 3) + 1 for k in range(n // 2 + 1, n + 1)}
    gap = sum(minus.values()) - sum(plus.values())
    if gap > 0:
        plus[n] = plus.get(n, 0) + gap
    elif gap < 0:
        minus[1] = minus.get(1, 0) - gap
    return minus, plus


class FastPathAgainstExact(unittest.TestCase):
    def test_the_contract_histogram_reproduces_its_exact_cdf(self) -> None:
        """Stops us believing a profile CDF that is not the committed channel.

        The wrong number is a threshold law reconstructed by a different mixture
        convention from ``threshold_histogram_profile``'s, which would make every
        quantile below a different observable from the one #28 froze.
        """
        contract = json.loads(
            (ROOT / "analysis" / "threshold_histogram_profile_contract.json")
            .read_text(encoding="utf-8"))
        minus = {int(k): v for k, v in contract["K_minus_counts"].items()}
        plus = {int(k): v for k, v in contract["K_plus_counts"].items()}
        report = lineage.exactness_control(minus, plus, contract["N"])
        self.assertLess(report["largest_cdf_difference"], 1e-14)
        self.assertEqual(report["largest_quantile_excursion_outside_the_exact_bracket"], 0.0)

    def test_the_binomial_identity_holds_on_a_spread_out_histogram(self) -> None:
        """Stops us believing a fast path validated only where the support is tiny.

        The production identity is ``F(p) = sum_j B(N,j,p) G(j)``, and the wrong
        number it could produce is a CDF whose binomial window silently clipped
        mass.  A four-rank contract cannot exercise that; twenty ranks can.
        """
        minus, plus = _synthetic(20)
        report = lineage.exactness_control(minus, plus, 20)
        self.assertLess(report["largest_cdf_difference"], 1e-13)
        self.assertEqual(report["largest_quantile_excursion_outside_the_exact_bracket"], 0.0)

    def test_the_cdf_is_monotone_and_hits_both_ends(self) -> None:
        """Stops us bisecting a function that is not a CDF."""
        minus, plus = _synthetic(20)
        cumulative = lineage.rank_cdf(minus, plus, 20)
        values = [lineage.profile_cdf(cumulative, 20, p / 50.0) for p in range(1, 50)]
        for left, right in zip(values, values[1:]):
            self.assertLessEqual(left, right + 1e-15)
        self.assertEqual(lineage.profile_cdf(cumulative, 20, 0.0), 0.0)
        self.assertEqual(lineage.profile_cdf(cumulative, 20, 1.0), 1.0)

    def test_the_binomial_window_does_not_underflow_at_large_n(self) -> None:
        """Stops the N=1300 underflow from coming back.

        Anchoring at j = 0 makes ``(1-p)**N`` exactly zero near N = 790 at the
        threshold, and every later term with it.  The wrong number that produced
        was a channel that returned 0.0 with nothing raised.
        """
        for size in (800, 1500, 4000):
            low, weights = lineage._binomial_weights(size, 0.5927)
            self.assertGreater(max(weights), 0.0)
            self.assertLess(low, int(size * 0.5927))
            self.assertGreater(low + len(weights), int(size * 0.5927))


class SpinFourWeights(unittest.TestCase):
    def test_cos_four_theta_is_the_exact_rational_it_should_be(self) -> None:
        """Stops a sign or a factor entering the coefficient that sets the correction.

        ``cos 4theta = Re(w^4)/|w|^4``.  A wrong sign here would move the spin-4
        correction the wrong way and would look like a shape flow.
        """
        self.assertAlmostEqual(lineage.cos_four_theta(1, 0), 1.0, places=15)
        self.assertAlmostEqual(lineage.cos_four_theta(1, 1), -1.0, places=15)
        self.assertAlmostEqual(lineage.cos_four_theta(8, 1), 3713 / 4225, places=15)
        self.assertAlmostEqual(lineage.cos_four_theta(7, 4), -2047 / 4225, places=15)

    def test_the_weights_cancel_cos_four_theta_and_sum_to_one(self) -> None:
        """Stops us believing a spin-0 law that still carries spin 4.

        The equal-weight average leaves +0.197 of it at N=65, and that residue
        alternates in sign along a lineage.  The wrong result is a correction
        that does not actually zero the net coefficient.
        """
        for first, second in ((0.8788, -0.4845), (0.2120, 0.9755), (-0.9162, -0.0233)):
            cos4 = {"first": first, "second": second}
            weights = lineage.spin_zero_weights(cos4)
            self.assertAlmostEqual(sum(weights.values()), 1.0, places=12)
            self.assertAlmostEqual(
                weights["first"] * first + weights["second"] * second, 0.0, places=12)

    def test_same_sign_orientations_are_flagged_as_an_extrapolation(self) -> None:
        """Stops an extrapolated correction being reported as if it were a mixture.

        N=325 and N=425 carry the same sign of cos4theta on both orientations, so
        their spin-0 weights leave [0,1].  The correction is still linear and
        still valid, but it is a weaker object and must not pass silently.
        """
        interior = lineage.spin_zero_weights({"first": 0.8788, "second": -0.4845})
        self.assertTrue(lineage.is_interpolation(interior))
        outside = lineage.spin_zero_weights({"first": 0.2120, "second": 0.9755})
        self.assertFalse(lineage.is_interpolation(outside))

    def test_equal_orientations_refuse_rather_than_divide_by_a_small_gap(self) -> None:
        """Stops a spin-4 correction being manufactured from two identical geometries."""
        with self.assertRaises(ValueError):
            lineage.spin_zero_weights({"first": 0.5, "second": 0.5})


class Jackknife(unittest.TestCase):
    def test_the_covariance_is_of_one_object_across_the_whole_grid(self) -> None:
        """Stops us building a diagonal covariance for nine points of one curve.

        Nine quantiles of one distribution, estimated from one sample, are
        strongly correlated.  A diagonal S would make almost any displacement
        look significant, which is the failure mode this whole line of work is
        about.  The wrong number is an off-diagonal that came out zero.
        """
        full = [0.4, 0.5, 0.6]
        deleted = {b: [0.4 + 0.001 * b, 0.5 + 0.0009 * b, 0.6 + 0.0011 * b]
                   for b in range(20)}
        covariance = lineage.jackknife_covariance(full, deleted)
        for i in range(3):
            self.assertGreater(covariance[i][i], 0.0)
            for j in range(3):
                self.assertAlmostEqual(covariance[i][j], covariance[j][i], places=18)
        correlation = covariance[0][2] / math.sqrt(covariance[0][0] * covariance[2][2])
        self.assertGreater(correlation, 0.99)

    def test_deleting_a_batch_removes_exactly_that_batch(self) -> None:
        """Stops a delete-one that leaves counts behind or removes too many."""
        total = {1: 5, 2: 7, 3: 2}
        self.assertEqual(lineage._subtract(total, {2: 7}), {1: 5, 3: 2})
        with self.assertRaises(ValueError):
            lineage._subtract(total, {1: 6})


if __name__ == "__main__":
    unittest.main()

import math
import sys
import unittest
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from analyze_p48_retrospective import (  # noqa: E402
    Histogram,
    covariance_of_mean,
    gls_score,
    pseudovalues,
    tail_and_derivative,
    validate_alignment,
)


class P48RetrospectiveTests(unittest.TestCase):
    def test_threshold_tail_and_derivative(self):
        # A deterministic K=2 threshold at N=3 has tail 3p^2(1-p)+p^3.
        p = 0.37
        value, derivative = tail_and_derivative([0, 0, 1, 0], 1, p)
        self.assertAlmostEqual(value, 3 * p * p * (1 - p) + p ** 3, places=14)
        self.assertAlmostEqual(derivative, 6 * p * (1 - p), places=14)

    def test_the_binomial_tail_survives_the_sizes_we_want_to_run(self):
        """The wrong number here is a silent zero, not an exception.

        The obvious recurrence starts at ``(1-p)**n`` and walks up. At the
        percolation threshold ``1-p`` is about 0.407, so that first term
        underflows to exactly zero somewhere near 790 sites, and the
        multiplicative recurrence then stays at zero for every later term: the
        tail, its derivative, and every channel built on them come back as 0.0
        with nothing raised. A 1M pilot at N=1300 is where we found this.

        A deterministic threshold at K makes the tail an exact binomial
        survival function, so it can be checked against a value nothing in this
        module computes.
        """
        for n in (290, 800, 1300, 2600):
            hist = [0] * (n + 1)
            hist[n // 2] = 1
            p = 0.5927460508
            value, derivative = tail_and_derivative(hist, 1, p)
            with self.subTest(n=n):
                # the median-rank tail sits strictly inside (0,1) for these sizes
                self.assertGreater(value, 1e-6, "tail collapsed to zero")
                self.assertLess(value, 1.0)
                self.assertGreater(derivative, 0.0)
                # d/dp P(Bin >= K) = n * C(n-1, K-1) p^(K-1) q^(n-K), checked in logs
                rank = n // 2
                expected = (
                    math.lgamma(n + 1)
                    - math.lgamma(rank)
                    - math.lgamma(n - rank + 1)
                    + (rank - 1) * math.log(p)
                    + (n - rank) * math.log(1.0 - p)
                )
                self.assertAlmostEqual(math.log(derivative), expected, places=8)

    def test_the_tail_is_a_probability_at_every_rank(self):
        """A monotone survival function is what every downstream root-find assumes.

        `project_size` bisects on the sign of a difference of these tails; a
        non-monotone or out-of-range tail would send it off the end of (0,1),
        which is the symptom the N=1300 pilot actually showed.

        The tolerance is not slack. Summing n+1 recurrence terms leaves a
        rounding residue of order n*eps -- about 1e-12 at n=1300 -- so the tail
        can sit a few parts in 1e12 above 1. The quantities built on it are
        differences of order 1e-3, which is nine orders of magnitude away.
        """
        n = 1300
        p = 0.4
        tolerance = 1e-9
        previous = 1.0
        for rank in (1, 100, 400, 520, 700, 1000, n):
            hist = [0] * (n + 1)
            hist[rank] = 1
            value, _ = tail_and_derivative(hist, 1, p)
            with self.subTest(rank=rank):
                self.assertGreaterEqual(value, 0.0)
                self.assertLessEqual(value, 1.0 + tolerance)
                self.assertLessEqual(value, previous + tolerance)
            previous = value
        self.assertLess(previous, 1e-100, "the top-rank tail should be vanishingly small")

    def test_pseudovalue_covariance_normalization(self):
        self.assertEqual(pseudovalues(2.0, [1.5, 2.0, 2.5]), [3.0, 2.0, 1.0])
        covariance = covariance_of_mean([[3.0, 6.0], [2.0, 4.0], [1.0, 2.0]])
        self.assertAlmostEqual(covariance[0][0], 1.0 / 3.0)
        self.assertAlmostEqual(covariance[0][1], 2.0 / 3.0)
        self.assertAlmostEqual(covariance[1][1], 4.0 / 3.0)

    def test_gls_retains_train_heldout_cross_covariance(self):
        values = [1.0, 1.2, 1.1]
        covariance = [
            [0.04, 0.00, 0.01],
            [0.00, 0.09, 0.02],
            [0.01, 0.02, 0.16],
        ]
        score = gls_score(values, covariance, [0, 1], [2])
        expected_amplitude = (1.0 / 0.04 * 1.0 + 1.0 / 0.09 * 1.2) / (1.0 / 0.04 + 1.0 / 0.09)
        self.assertAlmostEqual(score["amplitude"], expected_amplitude)
        influence = score["influence"]
        cross = covariance[0][2] * influence[0] + covariance[1][2] * influence[1]
        expected_variance = 0.16 + score["amplitude_se"] ** 2 - 2 * cross
        self.assertAlmostEqual(score["residual_covariance"][0][0], expected_variance)

    def test_alignment_rejects_unsynchronized_batches(self):
        records = {}
        for orientation in ("first", "second"):
            for batch in (0, 1):
                n = 2
                records[(n, orientation, batch)] = Histogram(
                    n, 1, 1, orientation, batch, 1, [0, 1, 0], [0, 1, 0]
                )
        # A second N is missing batch 1, so its counter blocks cannot be
        # deleted synchronously with the first N.
        for orientation in ("first", "second"):
            records[(5, orientation, 0)] = Histogram(
                5, 2, 1, orientation, 0, 1, [0, 1, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0]
            )
        with self.assertRaisesRegex(ValueError, "same contiguous batch ids"):
            validate_alignment(records)


if __name__ == "__main__":
    unittest.main()

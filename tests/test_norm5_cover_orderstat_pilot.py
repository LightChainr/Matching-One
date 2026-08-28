from __future__ import annotations

import itertools
import math
from pathlib import Path
import random
import statistics
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from gaussian_covering_map import GaussianPair, canonical_cover  # noqa: E402
from norm5_cover_orderstat_pilot import (  # noqa: E402
    beta_order_cdf,
    child_priorities,
    orderstat_parent_fields,
    priorities_to_permutation,
)
from integer_period_torus import gaussian_integer_torus  # noqa: E402
from threshold_rank_nz import ThresholdRankEngine, enumerate_exact  # noqa: E402


class OrderStatisticCouplingTests(unittest.TestCase):
    def test_integer_beta_cdf_endpoints_and_median_symmetry(self) -> None:
        for degree in (2, 3, 5):
            for order in range(1, degree + 1):
                self.assertEqual(beta_order_cdf(order, degree, 0.0), 0.0)
                self.assertEqual(beta_order_cdf(order, degree, 1.0), 1.0)
                self.assertAlmostEqual(
                    beta_order_cdf(order, degree, 0.37),
                    1.0 - beta_order_cdf(degree + 1 - order, degree, 0.63),
                    places=14,
                )

    def test_each_order_field_has_uniform_moments(self) -> None:
        generator = random.Random(17)
        degree = 5
        fields = orderstat_parent_fields(
            [tuple(generator.random() for _ in range(degree)) for _ in range(30000)]
        )
        for field in fields:
            self.assertLess(abs(statistics.mean(field) - 0.5), 0.006)
            self.assertLess(abs(statistics.variance(field) - 1.0 / 12.0), 0.003)

    def test_each_parent_permutation_is_uniform_on_tiny_system(self) -> None:
        generator = random.Random(29)
        degree = 5
        permutations = tuple(itertools.permutations(range(3)))
        counts = [{permutation: 0 for permutation in permutations} for _ in range(degree)]
        samples = 18000
        for _ in range(samples):
            fibers = [tuple(generator.random() for _ in range(degree)) for _ in range(3)]
            for order, field in enumerate(orderstat_parent_fields(fibers)):
                counts[order][priorities_to_permutation(field)] += 1
        expected = samples / math.factorial(3)
        # Six cells and a fixed seed: this broad Pearson bound guards gross
        # marginal defects without turning the regression into a flaky test.
        for order_counts in counts:
            chi2 = sum((value - expected) ** 2 / expected for value in order_counts.values())
            self.assertLess(chi2, 20.0)

    def test_child_assignment_is_a_bijective_relabeling(self) -> None:
        cover = canonical_cover(GaussianPair(2, 1), GaussianPair(4, 3))
        generator = random.Random(41)
        fibers = [
            tuple(generator.random() for _ in range(cover.degree))
            for _ in range(cover.parent.n)
        ]
        child = child_priorities(cover, fibers)
        self.assertEqual(sorted(child), sorted(value for row in fibers for value in row))
        self.assertEqual(len(priorities_to_permutation(child)), cover.child.n)

    def test_each_field_is_unbiased_for_tiny_threshold_rank_observable(self) -> None:
        geometry = gaussian_integer_torus(2, 1)
        engine = ThresholdRankEngine(geometry)
        exact = enumerate_exact(geometry)
        exact_mean = exact.sum_kplus / exact.sample_count
        generator = random.Random(97)
        degree = 5
        sums = [0.0] * degree
        samples = 6000
        for _ in range(samples):
            fibers = [
                tuple(generator.random() for _ in range(degree))
                for _ in range(geometry.n)
            ]
            for order, field in enumerate(orderstat_parent_fields(fibers)):
                permutation = priorities_to_permutation(field)
                sums[order] += engine.threshold_ranks(permutation)[1]
        for total in sums:
            self.assertLess(abs(total / samples - exact_mean), 0.04)


if __name__ == "__main__":
    unittest.main()

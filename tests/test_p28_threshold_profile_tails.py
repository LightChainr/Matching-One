import math
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_p28_threshold_profile_tails import (  # noqa: E402
    beta_density_basis,
    binomial_tail_basis,
    mixture_location_scale,
    orthogonal_contrasts,
)


class P28ThresholdProfileTailTests(unittest.TestCase):
    def test_uniform_rank_mixture_is_uniform(self):
        n = 5
        counts = {rank: 2 for rank in range(1, n + 1)}
        center, scale = mixture_location_scale(n, counts, sum(counts.values()))
        self.assertAlmostEqual(center, 0.5)
        self.assertAlmostEqual(scale, math.sqrt(1 / 12))
        basis = beta_density_basis(n, 0.37)
        density = sum(counts[rank] * basis[rank - 1] for rank in counts) / sum(counts.values())
        self.assertAlmostEqual(density, 1.0, places=12)

    def test_binomial_tail_basis(self):
        tails = binomial_tail_basis(2, 0.5)
        self.assertAlmostEqual(tails[0], 0.75)
        self.assertAlmostEqual(tails[1], 0.25)

    def test_contrasts_annihilate_intercept_and_slope(self):
        x = [value ** (4 / 3) for value in (2.5, 2.75, 3.0, 3.25, 3.5)]
        contrasts = orthogonal_contrasts(x)
        self.assertEqual(len(contrasts), 3)
        for row in contrasts:
            self.assertAlmostEqual(sum(row), 0.0, places=12)
            self.assertAlmostEqual(sum(a * b for a, b in zip(row, x)), 0.0, places=12)


if __name__ == "__main__":
    unittest.main()

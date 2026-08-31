from fractions import Fraction as F
from math import comb, lcm
from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"scripts"))
from p267_thermal_warp_invariants import compose, mul, scale
from p267_scalar_clock_transport import exact_empirical_root_intervals, feature_vector


def to_bernstein(power, degree):
    return [sum(power[k]*F(comb(i, k), comb(degree, k))
                for k in range(min(i+1, len(power)))) for i in range(degree+1)]


class ScalarClockTests(unittest.TestCase):
    def test_one_monotone_scalar_map_preserves_both_branches_and_valley(self):
        d = mul([F(0), F(1), F(-1)], [F(7, 20), F(-1), F(1)])
        phi = [F(0), F(6, 5), F(-1, 5)]
        u = scale(compose(d, phi), F(-3, 7))
        n = len(u)-1
        co = np.zeros((2, 2, n+1))
        co[0, 0], co[1, 0] = to_bernstein(d, n), to_bernstein(u, n)
        out, curvatures = feature_vector(co, [(.1, .4), (.4, .6), (.6, .9)])
        np.testing.assert_allclose(out[-2:], 0, atol=1e-12)
        self.assertTrue(np.all(np.sign(curvatures[0]) == -np.sign(curvatures[1])))

    def test_exact_root_count_with_ordered_rational_intervals(self):
        roots = [F(1, 3), F(2, 5), F(4, 5)]
        p = [F(1)]
        for r in roots:
            p = mul(p, [-r, F(1)])
        b = to_bernstein(p, 3)
        denominator = lcm(*(x.denominator for x in b))
        intervals = exact_empirical_root_intervals([int(x*denominator) for x in b])
        self.assertEqual(len(intervals), 3)
        for root, (left, right) in zip(roots, intervals):
            self.assertLess(F(left), root)
            self.assertLess(root, F(right))


if __name__ == "__main__":
    unittest.main()

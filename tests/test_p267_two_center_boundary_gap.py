from fractions import Fraction as F
from itertools import combinations
from math import prod
import unittest


def exact_gap(centers, weights):
    mean = sum(w*x for w, x in zip(weights, centers))
    y = [x-mean for x in centers]
    v, m3, m4 = [sum(w*x**r for w, x in zip(weights, y)) for r in (2, 3, 4)]
    hankel = (v*m4-m3*m3-v**3)/v**3
    vandermonde = prod(weights)*prod((b-a)**2 for a, b in combinations(centers, 2))/v**3
    energy = sum(w*(x*x-m3/v*x-v)**2 for w, x in zip(weights, y))/v**2
    return hankel, vandermonde, energy


class ExactRankTwoGap(unittest.TestCase):
    def test_identities_affine_boundary_and_broader_two_lobe_example(self):
        centers, weights = [F(-2), F(0), F(3)], [F(1, 5), F(1, 2), F(3, 10)]
        result = exact_gap(centers, weights)
        self.assertEqual(len(set(result)), 1)
        moved = [F(5, 3)*x-F(7, 5) for x in centers]
        self.assertEqual(exact_gap(moved, weights), result)
        self.assertEqual(exact_gap(centers, [F(0), F(2, 5), F(3, 5)]), (0, 0, 0))
        # B=+-1, Z=+-1 have a three-point sum but are two translates of one
        # symmetric kernel. Adding the same Gaussian preserves this example.
        self.assertEqual(exact_gap([F(-2), F(0), F(2)], [F(1, 4), F(1, 2), F(1, 4)]), (1, 1, 1))


if __name__ == "__main__":
    unittest.main()

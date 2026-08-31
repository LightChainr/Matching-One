import sys
from pathlib import Path
import unittest
from fractions import Fraction
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from p437_high_pass_mc import COEFFICIENTS, MASK, MATRICES, coupled_masks, degree5, random_mask, non_low_degree_certificate
from integer_period_torus import integer_torus_geometry


class PilotTests(unittest.TestCase):
    def test_exact_operator_targets(self):
        for degree in range(5):
            self.assertEqual(sum(Fraction(c, 2 ** (level * degree)) for level, c in enumerate(COEFFICIENTS)), 0)
        self.assertEqual(sum(Fraction(c, 2 ** (level * 5)) for level, c in enumerate(COEFFICIENTS)), Fraction(9765, 32768))

    def test_nested_common_noise(self):
        for replica in range(10):
            levels, replacement = coupled_masks(437, replica)
            self.assertEqual(levels[0], random_mask(437, replica, 0))
            changed = 0
            for value in levels:
                now = levels[0] ^ value
                self.assertEqual(changed & now, changed)
                self.assertEqual((value ^ levels[0]) & (value ^ replacement), 0)
                changed = now
            self.assertEqual(replacement & MASK, replacement)

    def test_geometry_and_normalized_control(self):
        for matrix in MATRICES:
            geometry = integer_torus_geometry(matrix)
            self.assertEqual((geometry.n, len(geometry.primal_edges)), (112, 224))
        self.assertEqual(sum(degree5(mask) for mask in range(32)), 0)
        self.assertEqual(sum(degree5(mask) ** 2 for mask in range(32)), 32)

    def test_topology_is_not_degree_at_most_four(self):
        self.assertEqual(non_low_degree_certificate()["mixed_difference_re"], "-1/3")


if __name__ == "__main__":
    unittest.main()

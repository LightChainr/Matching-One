import sys
from fractions import Fraction
from pathlib import Path
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from p437_positive_difference_bridge import certificate, divided_difference_coefficients, mixed_derivative


class PositiveBridgeTests(unittest.TestCase):
    def test_exact_multiplier_bridge(self):
        self.assertEqual(divided_difference_coefficients(), [1, -31, 310, -1240, 1984, -1024])
        self.assertEqual(len(certificate()["multiplier_identity"]), 21)

    def test_pointwise_degree_cutoff_and_signed_pair(self):
        def chi(mask, n):
            return Fraction((-1) ** (n - bin(mask & ((1 << n) - 1)).count("1")))
        sites = range(5)
        self.assertEqual(mixed_derivative(0, sites, lambda m: chi(m, 4)), 0)
        self.assertEqual(mixed_derivative(0, sites, lambda m: chi(m, 5)), 1)
        left = mixed_derivative(0, sites, lambda m: chi(m, 6))
        right = mixed_derivative(32, sites, lambda m: chi(m, 6))
        self.assertEqual(left * right, -1)

    def test_uniform_subset_variance_is_exact(self):
        row = certificate()["sparse_degree5_uniform_subset_obstruction"]
        q = Fraction(row["nonzero_probability"])
        z = Fraction(row["nonzero_response"])
        mean = Fraction(row["mean"])
        self.assertEqual(q * z, mean)
        self.assertEqual(q * z * z - mean * mean, Fraction(row["variance"]))


if __name__ == "__main__":
    unittest.main()

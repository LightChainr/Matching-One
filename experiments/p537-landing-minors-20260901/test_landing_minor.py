import unittest
from fractions import Fraction as F

import landing_minor as a


class LandingMinorTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = a.exact_landing_minor()

    def test_physical_branch_free_c4_orbit(self):
        self.assertEqual(len(self.result["orbit"]), 4)
        for row in self.result["orbit"]:
            self.assertIn(row["neighbor_occupancy_NESW"], ([1, 0, 1, 0], [0, 1, 0, 1]))
            self.assertEqual(row["rank_before_after"], [0, 1])
            self.assertEqual(row["kernel_before_after"], ["0", "1/4"])

    def test_nonzero_minor_at_a_rational_control(self):
        self.assertEqual(F(self.result["determinant_at_half"]), F(1001, 536870912))

    def test_finite_root_is_not_a_minor_root(self):
        self.assertEqual(self.result["gcd_root_and_determinant"], ["1"])
        self.assertEqual(self.result["decision"], "nonzero_exact_root_conditioned_2x2_minor")


if __name__ == "__main__":
    unittest.main()

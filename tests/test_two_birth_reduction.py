"""Targeted controls; these do not determine the universal sharp-threshold constant."""
from fractions import Fraction
from pathlib import Path
import sys
import unittest
import mpmath as mp
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import two_birth_reduction as t


class TwoBirthReductionTests(unittest.TestCase):
    def test_lifted_parallel_edges(self):
        self.assertEqual(t.lifted_rank(2, 2, 0), 0)
        self.assertEqual(t.lifted_rank(2, 2, 1), 0)
        self.assertEqual(t.lifted_rank(2, 2, 3), 1)
        self.assertEqual(t.lifted_rank(2, 2, 7), 2)
        self.assertEqual(t.lifted_rank(2, 2, 15), 2)

    def test_exact_census_and_permutation_integral(self):
        with mp.workdps(70):
            data = t.exact_small_controls()
        self.assertEqual(data["total_configurations"], 336)
        self.assertEqual(data["permutation_control"]["permutations"], 720)
        self.assertEqual(data["permutation_control"]["birth_time_covariance"], "1123/58800")

    def test_beta_integral(self):
        self.assertEqual(t.beta_integral([1, 4, 6, 4, 1]), Fraction(1))
        self.assertEqual(t.beta_integral([0, 0, 6, 0, 0]), Fraction(1, 5))

    def test_median_and_quartile_bracketing(self):
        with mp.workdps(65):
            for m in (2, 4, 16):
                a = t.quantile(m, "first", ".5")
                b = t.quantile(m, "second", ".5")
                q = t.quantile(m, "mixture", ".5")
                self.assertLessEqual(t.quantile(m, "mixture", ".25"), a)
                self.assertLessEqual(a, q)
                self.assertLessEqual(q, b)
                self.assertLessEqual(b, t.quantile(m, "mixture", ".75"))

    def test_probability_range(self):
        with mp.workdps(65):
            for m in (2, 8, 128):
                for p in (mp.mpf(".01"), mp.mpf(".5"), mp.mpf(".99")):
                    probs = t.sector_probabilities(m, p)
                    self.assertLess(abs(sum(probs) - 1), mp.mpf("1e-60"))
                    self.assertTrue(all(-mp.mpf("1e-60") <= x <= 1 + mp.mpf("1e-60") for x in probs))

    def test_fixed_length_coupling_and_gap(self):
        with mp.workdps(65):
            data = t.finite_record(4)
        self.assertGreater(float(data["integral_P1"]), 0)
        self.assertLessEqual(float(data["mixture_W1_to_two_median_atoms"]), float(data["component_coupling_W1_upper"]))

    def test_illegal_inputs(self):
        with self.assertRaises(ValueError):
            t.lifted_rank(1, 4, 0)
        with self.assertRaises(ValueError):
            t.quantile(2, "first", 1)
        with self.assertRaises(ValueError):
            t.sector_probabilities(1, mp.mpf(".5"))


if __name__ == "__main__":
    unittest.main()

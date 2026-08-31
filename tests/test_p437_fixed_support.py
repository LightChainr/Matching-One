import sys
from fractions import Fraction
from pathlib import Path
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from integer_period_torus import integer_torus_geometry
from p437_high_pass_mc import MATRICES
from p437_fixed_support_mc import DENOMINATOR, child_differences, energy_numerator
from p437_positive_difference_bridge import H5, spectral_multiplier


class FixedSupportTests(unittest.TestCase):
    def test_previous_exact_witness_normalization(self):
        geometry = [integer_torus_geometry(matrix) for matrix in MATRICES]
        outside = sum(1 << edge for edge in (140,168,196))
        differences = child_differences(outside, geometry)
        self.assertEqual(differences, (-1,0,0))
        self.assertEqual(Fraction(energy_numerator(differences),DENOMINATOR),Fraction(1,9216))

    def test_population_bound_mode_by_mode(self):
        self.assertEqual(spectral_multiplier(5),H5)
        for j in range(5,225):
            self.assertGreaterEqual(spectral_multiplier(j),H5)


if __name__ == "__main__":
    unittest.main()

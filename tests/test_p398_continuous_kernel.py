from fractions import Fraction as F
from math import sqrt
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"scripts"))
from p398_continuous_kernel import derive_generator, sample, stationary_limit


class ContinuousKernelTests(unittest.TestCase):
    def test_exact_generator_and_critical_amplitude(self):
        derive_generator()
        row=stationary_limit(F(1))
        self.assertEqual(row["C0_complex_re_im"],
            [[(F(6,7),F(0)),(F(-4,7),F(4,7))],[(F(-4,7),F(-4,7)),(F(6,7),F(0))]])
        self.assertEqual(row["C0_conjugate_M_offdiagonal_skew"],(0,0))

    def test_amplitude_metric_free_fingerprints(self):
        for k in (F(1,4),F(1),F(4)):
            row=sample(k)
            self.assertAlmostEqual(25*row["metric_free_mass_split"]**2*(2-row["diagonal_gauge_free_mixing"]),2,places=12)
            for point in row["kernel_samples"]:
                self.assertAlmostEqual(point["metric_free_mass_split"],row["metric_free_mass_split"],places=10)
                self.assertAlmostEqual(point["diagonal_gauge_free_mixing"],row["diagonal_gauge_free_mixing"],places=10)
                self.assertAlmostEqual(point["signed_channel_axis"],float(k-1)/sqrt(float(k*k+6*k+1)),places=10)


if __name__=="__main__":
    unittest.main()

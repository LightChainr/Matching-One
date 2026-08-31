import sys
from pathlib import Path
import unittest
from fractions import Fraction
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"scripts"))
from score_p437_coherent_decomposition import coherent_weight,symmetry_certificate
from p437_positive_difference_bridge import H5,spectral_multiplier


class CoherentTests(unittest.TestCase):
    def test_u_statistic_and_variance_correction(self):
        z=np.array([[1.,2.],[3.,-1.],[-2.,4.],[.5,-3.]])
        pair=sum(float(z[i]@z[j]) for i in range(4) for j in range(4) if i!=j)/12
        self.assertAlmostEqual(coherent_weight(z),pair)
        self.assertAlmostEqual(pair,float(z.mean(0)@z.mean(0)-np.trace(np.cov(z,rowvar=False)/4)))

    def test_unresolved_can_be_negative(self):
        self.assertLess(coherent_weight(np.array([[1.,0.],[-1.,0.]])),0)

    def test_exact_mirror_phase_not_total_zero(self):
        c=symmetry_certificate()
        self.assertEqual(c["independent_constraint_on_real_child_coefficients"],[[0,1,-1]])
        self.assertEqual(len(c["allowed_child_coefficient_basis"]),2)

    def test_refined_multiplier(self):
        self.assertEqual(spectral_multiplier(6),Fraction(615195,1048576))
        self.assertEqual(spectral_multiplier(6)/H5,Fraction(63,32))


if __name__=="__main__":
    unittest.main()

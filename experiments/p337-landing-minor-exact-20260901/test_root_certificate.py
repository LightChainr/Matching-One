from __future__ import annotations
import unittest
from fractions import Fraction as F
import root_certificate as rc

class RootCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls): cls.result=rc.build()
    def test_matching_root_is_unique(self): self.assertEqual(self.result['matching_root_count'],1)
    def test_no_root_in_decision_denominators(self):
        self.assertEqual(self.result['determinant_root_count'],0)
        self.assertEqual(self.result['thermal_sum_root_count'],0)
        self.assertEqual(self.result['mixed_wronskian_root_count'],0)
    def test_endpoint_signs(self): self.assertTrue(all(F(x)<0 for x in self.result['determinant_signs']+self.result['thermal_sum_signs']+self.result['mixed_wronskian_signs']))
    def test_root_schur_residual(self):
        self.assertEqual(F(self.result['half_root_schur_even_residual']),F(533831111,1539745775616))
        self.assertGreater(float(self.result['matrix_at_root_midpoint']['root_schur_even_residual']),0.0)
    def test_root_conditioned_mixed_hessian(self):
        self.assertEqual(F(self.result['half_root_conditioned_mixed_hessian']),F(1397902943671,32208917889024))
        self.assertGreater(float(self.result['matrix_at_root_midpoint']['root_conditioned_mixed_hessian']),0.0)

if __name__=='__main__':unittest.main()

from __future__ import annotations
import unittest
from fractions import Fraction as F
import root_certificate as rc

class RootCertificateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls): cls.result=rc.build()
    def test_matching_root_is_unique(self): self.assertEqual(self.result['matching_root_count'],1)
    def test_landing_determinant_has_no_root(self): self.assertEqual(self.result['determinant_root_count'],0)
    def test_determinant_is_negative_on_interval_endpoints(self): self.assertTrue(all(F(x)<0 for x in self.result['determinant_signs']))
    def test_root_midpoint_matrix_is_nonzero(self): self.assertLess(float(self.result['matrix_at_root_midpoint']['determinant']),0.0)

if __name__=='__main__':unittest.main()

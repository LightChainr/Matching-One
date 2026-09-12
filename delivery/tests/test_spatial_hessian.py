import json
from fractions import Fraction as F
from pathlib import Path
import sys
import unittest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import torus_source_hessian as h

class SpatialHessianTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result=h.analyze()
    def test_known_rank_controls(self):
        self.assertEqual([h.physical_rank(x) for x in (0,15,(1<<16)-1,3+(3<<4))],[0,1,2,0])
    def test_rank_totals(self):
        self.assertEqual(self.result['rank_totals'],[36559,19932,9045])
    def test_uniform_thermal_sum_rule(self):
        a=self.result['origin_M_hessian_coefficients']
        self.assertEqual([16*sum(a[d][i] for d in range(16)) for i in range(15)],self.result['Msecond_coefficients'])
    def test_root_bracket(self):
        lo,hi=map(F,self.result['root_bracket'])
        self.assertLess(h.bernstein(self.result['M_coefficients'],lo),0)
        self.assertGreater(h.bernstein(self.result['M_coefficients'],hi),0)
        self.assertGreater(F(self.result['Mprime_root_interval'][0]),0)
    def test_fourier_exact_half_values(self):
        expected={(0,0):7429,(1,0):2665,(2,0):3733,(1,1):-2483,(2,1):-2567,(2,2):-5355}
        for row in self.result['fourier_modes']:
            k=tuple(row['wavevector'])
            if k in expected:self.assertEqual(F(row['eigenvalue_at_half']),F(expected[k],8192))
    def test_full_root_hessian_inertia(self):
        signs=[]
        for row in self.result['fourier_modes']:
            lo,hi=map(F,row['additive_root_curvature_interval'])
            signs.append(1 if lo>0 else -1 if hi<0 else 0)
        self.assertEqual([signs.count(1),signs.count(-1),signs.count(0)],[9,6,1])
    def test_logit_sign_certificates(self):
        signs=[]
        for row in self.result['fourier_modes']:
            if row['wavevector']==[0,0]:continue
            lo,hi=map(F,row['logit_root_curvature_interval'])
            signs.append(1 if lo>0 else -1 if hi<0 else 0)
        self.assertEqual([signs.count(1),signs.count(-1),signs.count(0)],[9,6,0])
    def test_root_uniform_gauge(self):
        row=self.result['fourier_modes'][0]
        self.assertEqual(row['additive_root_curvature_interval'],['0','0'])
        self.assertEqual(row['logit_source_logit_root_second_derivative_unit_RMS'],0)
        self.assertLess(row['logit_source_p_root_second_derivative_unit_RMS'],0)
    def test_hessian_trace_and_disorder_bias(self):
        total=sum(x['additive_p_root_second_derivative_unit_RMS']/16 for x in self.result['fourier_modes'])
        self.assertAlmostEqual(total,2*self.result['independent_additive_disorder_bias_per_variance'],places=14)
    def test_invalid_masks(self):
        for mask in (-1,1<<16):
            with self.assertRaises(ValueError):h.physical_rank(mask)

if __name__=='__main__':unittest.main()

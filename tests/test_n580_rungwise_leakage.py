"""Cheap mathematical controls of the cone projection, not a data-rescoring suite."""
from pathlib import Path
import sys
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import n580_rungwise_leakage as m
import projective_inference as p

class ConeTests(unittest.TestCase):
    def test_exact_identity_covariance_toy(self):
        y=[4,1,1];s=np.eye(3);v=[1,1,1];leak=[1,1,1]
        got=m.cone_projection(y,s,v,.5,leak)
        checked=m.high_precision_check(y,s,v,.5,got,leak)
        self.assertAlmostEqual(checked['distance'],2/11,places=13)
        np.testing.assert_allclose(checked['mean'],[42/11,14/11,14/11],rtol=1e-13)

    def test_bound_zero_equals_signed_gls_line(self):
        y=[.3,.7,.9];s=[[1,.1,0],[.1,2,.1],[0,.1,1]];v=[1,2,4]
        got=m.cone_projection(y,s,v,0)
        expected=p.ray_residual(y,s,v)['statistic']
        self.assertAlmostEqual(got['distance'],expected,places=13)

    def test_sign_reversal_and_nested_bounds(self):
        y=np.array([4,1,1]);s=np.eye(3);v=[1,1,1];leak=[1,1,1]
        previous=float('inf')
        for b in (0,.1,.3,.5,.7):
            value=m.cone_projection(y,s,v,b,leak)['distance']
            self.assertLessEqual(value,previous+1e-12)
            self.assertAlmostEqual(value,m.cone_projection(-y,s,v,b,leak)['distance'],places=13)
            previous=value
        self.assertLess(previous,1e-20)

if __name__=='__main__':
    unittest.main()

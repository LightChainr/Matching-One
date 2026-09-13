import sys
from pathlib import Path
import unittest
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from shape_reflection_order import summarize

class ReflectionOrder(unittest.TestCase):
    def test_common_logit_reflections_keep_order(self):
        f=lambda p,k:k*(1-p)/(p+k*(1-p))
        for p in (.3,.5,.7):
            self.assertLess(f(p,1),f(p,2));self.assertLess(f(p,2),f(p,3))
    def test_crossing_valid_involution_is_possible_without_common_chart(self):
        h=lambda x:x+.2*x*(1-x)*((x-.5)**2-.02)
        r=lambda x:brentq(lambda y:h(y)-(1-h(x)),0,1)
        self.assertLess(r(.1)-.9,0);self.assertGreater(r(.5)-.5,0)
        for p in (.1,.5,.9):self.assertAlmostEqual(r(r(p)),p,places=11)
    def test_bands_do_not_call_unresolved_signs_a_crossing(self):
        y=np.array([1,1,1,-1]*3,dtype=float)
        self.assertTrue(summarize(y,np.eye(12)*1e-4)['any_resolved_crossing'])
        self.assertFalse(summarize(y,np.eye(12))['any_resolved_crossing'])

if __name__=='__main__':unittest.main()

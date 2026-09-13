from pathlib import Path
import sys
import unittest
from fractions import Fraction
from mpmath import mp
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from shape_reflection_commutator import Law,Component

class InverseMixture(unittest.TestCase):
    def law(self):
        pool={'first':{'minus':{2:1},'plus':{2:1}},'second':{'minus':{4:1},'plus':{4:1}}}
        return Law(pool,5,{'first':Fraction(2,5),'second':Fraction(3,5)})
    def test_one_level_cdf_inverts_quantile_mixture(self):
        law=self.law()
        for u in (.1,.3,.5,.8,.9):
            p=law.quantile(u)
            self.assertAlmostEqual(law.cdf(p),u,places=12)
            self.assertAlmostEqual(law.reflect(law.reflect(p)),p,places=12)
    def test_bernstein_cdf_matches_beta_exact_formula(self):
        with mp.workdps(55):
            c=Component({2:1},{2:1},5);p=mp.mpf('.3')
            expected=mp.betainc(2,4,0,p,regularized=True)
            self.assertLess(abs(c.mp_cdf(p)-expected),mp.mpf('1e-50'))
    def test_high_precision_mixture_inverse(self):
        with mp.workdps(55):
            law=self.law();u=mp.mpf('.3');p=law.mp_quantile(u)
            self.assertLess(abs(law.mp_cdf(p)-u),mp.mpf('1e-23'))

if __name__=='__main__':unittest.main()

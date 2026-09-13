from fractions import Fraction as F
from pathlib import Path
import sys
import unittest
from mpmath import mp
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from shape_common_chart import reflection_commutator

class CommonChart(unittest.TestCase):
    def test_shared_logit_chart_exactly_commutes(self):
        maps=[lambda p,k=k:k*(1-p)/(p+k*(1-p)) for k in (F(1),F(2),F(3))]
        for p in (F(3,10),F(1,2),F(7,10)):
            self.assertEqual(reflection_commutator(*maps,p)['difference'],0)
    def test_individually_valid_reflections_need_not_share_chart(self):
        with mp.workdps(40):
            p=mp.mpf('.5')
            r1=lambda x:1-x
            r2=lambda x:2*(1-x)/(x+2*(1-x))
            r3=lambda x:mp.root(1-x**3,3)
            self.assertLess(abs(r3(r3(p))-p),mp.mpf('1e-35'))
            self.assertGreater(reflection_commutator(r1,r2,r3,p)['difference'],mp.mpf('.018'))
    def test_degenerate_zero_does_not_identify_a_chart(self):
        # R1=R2 makes the test zero, regardless of R3: power may be absent.
        with mp.workdps(40):
            r=lambda x:1-x
            r3=lambda x:mp.root(1-x**3,3)
            self.assertLess(abs(reflection_commutator(r,r,r3,mp.mpf('.5'))['difference']),mp.mpf('1e-35'))

if __name__=='__main__':unittest.main()

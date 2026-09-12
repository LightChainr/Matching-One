import sys
from pathlib import Path
import unittest
from fractions import Fraction as F
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import rectangular_rank_odds as r

class TestRankOdds(unittest.TestCase):
    def test_local_box_guard(self):
        with self.assertRaises(ValueError):r.local_arm(0,3,8,0,1,False)
        with self.assertRaises(ValueError):r.sign_certificate(3,1,F(1,100))
    def test_path_bound(self):
        self.assertEqual(r.path_arm_bound(4,F(1,10),2),F(3,250))
        self.assertEqual(r.path_arm_bound(8,F(1,20),2),F(7,1000))
    def test_exact_positive_margin(self):
        c=r.sign_certificate(16,2,F(3,250))
        self.assertTrue(c['uniform_sign_certified'])
        self.assertGreater(c['gamma_diagnostic'],0)
    def test_failed_margin_not_a_result(self):
        self.assertFalse(r.sign_certificate(16,2,F(1,10))['uniform_sign_certified'])
    def test_two_cross_cycles(self):
        w=m=8;mask=sum(1<<x for x in range(w))|sum(1<<(w*y) for y in range(m))
        self.assertEqual(r.lifted_rank(mask,w,m),2)
        self.assertTrue(r.local_arm(mask,w,m,0,2,False))
        self.assertTrue(all(r.slab_crossing(mask,w,m,y,2,False) for y in (0,3)))
    def test_slabs_disjoint_sites(self):
        w,m,R=8,13,2
        sets=[{w*y+x for y in range(a,a+R+1) for x in range(w)}
              for a in range(0,m-R,R+1)]
        self.assertEqual(sum(map(len,sets)),len(set.union(*sets)))
    def test_information_identity(self):
        v=r.fisher_decomposition(F(1,200),F(1,200),F(-1,200),F(3,200))
        self.assertTrue(v['exact_decomposition'])
        self.assertEqual(F(v['conditional_H']),F(1,2))
    def test_counts_derivative(self):
        value,der=r.evaluate_counts([0,0,1],F(1,3))
        self.assertEqual((value,der),(F(1,9),F(2,3)))

if __name__=='__main__':unittest.main()

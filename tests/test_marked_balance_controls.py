import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import marked_balance_controls as m
from fractions import Fraction as F

class TestMarkedBalance(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.r=m.report()
    def test_three_root_signs(self):
        self.assertEqual([x['sign'] for x in self.r['rational_sign_checks']],[-1,1,-1,1])
        self.assertTrue(self.r['at_least_three_distinct_roots'])
    def test_positive_normalization(self):
        for x in self.r['rational_sign_checks']:
            self.assertGreater(F(x['normalizer']),0)
            self.assertLessEqual(abs(F(x['normalized_balance'])),1)
    def test_odds_oscillation_bound(self):
        self.assertTrue(all(x['multiplicative_odds_bound_checked'] for x in self.r['rational_sign_checks']))

if __name__=='__main__':unittest.main()

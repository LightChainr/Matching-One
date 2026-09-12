"""The crossing comparison and all-length sign use rational bounds."""
from pathlib import Path
from fractions import Fraction
import sys
import unittest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import width4_cylinder_certificate as c

class CylinderCertificateTests(unittest.TestCase):
    def test_exact_intervals_and_all_m_sign(self):
        d=c.certificate()
        self.assertTrue(d['one_positive_root_Descartes'])
        self.assertGreater(Fraction(d['hprime_p_interval'][0]),0)
        self.assertLess(16*Fraction(d['relative_remainder_ratio_upper'])**6,2)
        self.assertEqual([r['m'] for r in d['all_m_root_sign']['m2_to_5']],[2,3,4,5])
        self.assertTrue(d['matches_Jacobsen_2015_table2_n4_within_1e_minus40'])

if __name__=='__main__':unittest.main()

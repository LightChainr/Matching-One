import sys
import unittest
from fractions import Fraction as F
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import coalescent_bulk_clock as c

class CoalescentBulkClockTests(unittest.TestCase):
    def test_01_cycle_histories(self):
        for n in range(2,7):
            result=c.enumerate_histories(n)
            self.assertEqual(result['multiplicity_per_history'],2**(n-1))
    def test_02_uniform_full_history_transitions(self):
        self.assertEqual(c.enumerate_histories(6)['ranked_histories'],2700)
    def test_03_eppf(self):
        self.assertEqual(c.enumerate_histories(6)['eppf_equalities'],203)
    def test_04_two_parallel_edges(self):
        self.assertEqual(c.cycle_history((0,1),(0,1)),c.cycle_history((0,1),(1,0)))
    def test_05_fixed_order_is_not_uniform_pairs(self):
        first=set()
        from itertools import permutations
        for perm in permutations(range(4)):
            first.add(c.cycle_history((0,1,2,3),perm)[1])
        self.assertEqual(len(first),4) # not choose(4,2)=6
    def test_06_pair_survival_two_methods(self):
        for n in range(2,13):
            for p in (F(0),F(1,7),F(1,3),F(1,2),F(4,5),F(1)):
                self.assertEqual(c.pair_separated_by_edges(n,p),c.pair_separated_by_levels(n,p))
    def test_07_jet_inverse(self):
        x={(0,0):F(3),(1,0):F(2),(0,1):F(-1),(1,1):F(5)}
        self.assertEqual(c.mul(x,c.inv(x)),{(0,0):F(1)})
    def test_08_score_fourth_moments(self):
        for k,a in ((F(1),F(1)),(F(16),F(1,8)),(F(256),F(1,32)),(F(65536),F(1,512))):
            self.assertEqual(c.score_moments(k,a)['fourth'],18)
    def test_09_transform_same_time(self):
        for r in (F(0),F(4,19),F(8)):
            v=c.six_variable_transform(1,1,r,F(1,3),F(1,5),F(1,2),F(3,4),F(2,3),F(-1,4))
            expected=(1+F(1,3)+F(1,5)+r*(1-F(3,8))+(F(2,3)-F(1,4))**2/2)**-2
            self.assertEqual(v,expected)
    def test_10_transform_marginals(self):
        for k,a in ((F(2),F(1,2)),(F(4),F(1,3))):
            v=c.six_variable_transform(k,a,8,0,0,1,1,F(2,3),0)
            self.assertEqual(v,(1+F(2,3)**2/2)**-2)
            v=c.six_variable_transform(k,a,8,0,0,1,1,0,F(2,3))
            self.assertEqual(v,(1+F(2,3)**2/2)**-2)
    def test_11_transform_negative_cross(self):
        self.assertGreater(c.six_variable_transform(2,F(7,10),0,0,0,1,1,10,-7),0)
    def test_12_bad_covariance_rejected(self):
        with self.assertRaises(ValueError): c.score_moments(4,1)
    def test_13_nested_labels(self):
        self.assertEqual(len(c.label_covariance_checks()),9)
    def test_14_retention_NB_Beta(self):
        for k in (2,3,4,8):
            self.assertAlmostEqual(c.retention_moment_mixture(k,1),c.retention_mean(k),places=12)
    def test_15_retention_has_atom(self):
        self.assertEqual(c.retention_mean(1),1)
        self.assertTrue(0<c.retention_mean(4)<1)
    def test_16_four_coordinate_closure(self):
        self.assertEqual(c.closure_stationarity_checks(),768)
    def test_17_complete_report(self):
        r=c.report()
        self.assertEqual(sum(x['realizations'] for x in r['genealogy']),89438)
        self.assertEqual(r['pair_survival_equalities'],66)
    def test_18_invalid_enumeration(self):
        with self.assertRaises(ValueError): c.enumerate_histories(7)

if __name__=='__main__':
    unittest.main()

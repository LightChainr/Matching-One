import unittest
from fractions import Fraction as F
import thermal_gate as a

class ThermalGateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.join=a.join_audit()
        cls.case=a.enumerate_case(3,0)
    def test_bell_counts(self):
        self.assertEqual(len(a.partitions(8)),4140)
        self.assertEqual(sum(a.nshared(p)==2 for p in a.partitions(8)),1922)
    def test_canonical_known_witnesses(self):
        self.assertEqual(a.g16(tuple(map(int,'00110011'))),1)
        self.assertEqual(a.g16(tuple(map(int,'00110012'))),-2)
    def test_zero_one_shared(self):
        self.assertTrue(all(a.g16(p)==0 for p in a.partitions(8) if a.nshared(p)<=1))
    def test_join_count_bound(self):
        self.assertEqual(self.join['checks'],64954)
        self.assertEqual(self.join['max_abs_delta16'],68)
        self.assertEqual(self.join['support_failures'],0)
    def test_rewiring_sign_flip(self):
        p=(0,0,1,1,0,0,1,2);q=a.merge_pattern(p,(1,2))
        self.assertEqual((a.nshared(p),a.nshared(q)),(2,2))
        self.assertEqual((a.g16(p),a.g16(q)),(-2,1))
    def test_dual_topology(self):
        self.assertEqual(self.case['topology_failures'],0)
        self.assertEqual(self.case['states'],512)
    def test_local_merge_equals_recomputation(self):
        self.assertEqual(self.case['state_edge_pair_checks'],25344)
        self.assertEqual(self.case['pair_update_failures'],0)
    def test_complete_covariance_identity(self):
        for r in self.case['covariance_controls']:
            x=F(r['kernel_rearrangement']['exact'])+F(r['observable_pivotal']['exact'])
            self.assertEqual(x,F(r['covariance_derivative']['exact']))
            self.assertEqual(r['score_residual'],'0')
    def test_endpoint_dilution(self):
        for r in self.case['covariance_controls']:
            self.assertEqual(F(r['raw_derivative_endpoint_dilution']['exact']),
                -2*F(r['mean_source']['exact'])/(1-F(r['p'])))
    def test_readout_term_can_reverse_sign(self):
        r=next(v for v in self.case['covariance_controls'] if v['p']=='3/5' and v['observable']=='E')
        self.assertGreater(F(r['kernel_rearrangement']['exact']),0)
        self.assertLess(F(r['covariance_derivative']['exact']),0)
    def test_mean_derivative_zero_does_not_control_covariance_derivative(self):
        # n=2, O=X, source=1[X!=Y], p=1/2.
        n=2;obs=[m&1 for m in range(4)];src=[64*int(m in (1,2)) for m in range(4)]
        r=a.conditional_covariance_decomposition(n,obs,src,F(1,2))
        self.assertEqual(F(r['source_derivative']['exact']),0)
        self.assertEqual(F(r['covariance_derivative']['exact']),-F(1,2))
    def test_geometry_rejects_aliased_ports(self):
        with self.assertRaises(ValueError):a.Torus(2,0)

if __name__=='__main__':unittest.main()

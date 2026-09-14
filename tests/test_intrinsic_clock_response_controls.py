import sys
import unittest
from fractions import Fraction as F
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import intrinsic_clock_response_controls as c

class ExactInterfaces(unittest.TestCase):
    def test_single_ring(self):
        self.assertEqual(c.winding_count(7),1)
        self.assertEqual(c.winding_count(3),0)
    def test_two_rings_and_late_merge(self):
        early=7|(7<<6)
        self.assertEqual(c.winding_count(early),2)
        self.assertEqual(c.winding_count(early|(1<<3)),1)
    def test_clock_endpoints(self):
        for r in (F(1,3),F(1,2),F(3,4)):
            self.assertEqual(c.beta(F(0),r),0)
            self.assertEqual(c.beta(r,r),c.nu(r))
    def test_positive_bernstein_derivative(self):
        for r in (F(1,3),F(1,2),F(3,4)):
            a=c.clock_bernstein(r)
            self.assertTrue(all(a[j+1]>=a[j] for j in range(c.N)))
            for j in range(1,8):
                self.assertGreater(c.beta_prime(r*F(j,8),r),0)
                self.assertLessEqual(c.beta_prime(r*F(j,8),r),c.N)
    def test_common_label_mass_transport(self):
        for p,r in ((F(1,6),F(1,3)),(F(1,4),F(1,2))):
            beta,loss,pairs,n=c.direct_pair_values(p,r)
            self.assertEqual(n,3**c.N)
            self.assertEqual(beta,c.beta(p,r))
            self.assertEqual(beta+loss,c.nu(p))
            self.assertLessEqual(loss,pairs)
    def test_inverse_clock(self):
        r,p0=F(3,4),F(1,4)
        for s in (F(1,2),F(1),F(2)):
            target=s*c.beta(p0,r)
            lo,hi=c.invert_beta(target,r)
            self.assertLessEqual(c.beta(lo,r),target)
            self.assertLessEqual(target,c.beta(hi,r))
            self.assertLessEqual(hi-lo,r/F(2**48))
    def test_nonlinear_tangent_and_normal(self):
        result=c.alignment_controls()
        self.assertEqual(result['exact_controls'],20)
    def test_higher_spin_alias(self):
        result=c.harmonic_controls()
        self.assertEqual(result['hidden_H8_into_scalar'],F(429,625))
        self.assertEqual(result['hidden_H8_into_spin4'],F(196,625))
    def test_two_harmonic_recovery_and_alias(self):
        h=F(-527,625)
        for b0,b4,b8 in ((F(2),F(3),F(0)),(F(0),F(0),F(1)),(F(2),F(-3),F(5))):
            ya=b0+b4+b8
            yb=b0+b4*h+b8*(2*h*h-1)
            fit4=(ya-yb)/(1-h)
            fit0=(yb-h*ya)/(1-h)
            self.assertEqual(fit0,b0+F(429,625)*b8)
            self.assertEqual(fit4,b4+F(196,625)*b8)
    def test_invalid_parameters(self):
        with self.assertRaises(ValueError): c.beta(F(3,4),F(1,2))
        with self.assertRaises(ValueError): c.clock_bernstein(F(1))
        with self.assertRaises(ValueError): c.invert_beta(F(-1),F(1,2))
class PhysicalSourceChecks(unittest.TestCase):
    def test_torus_empty_full(self):
        self.assertEqual(c.torus_rank(0),0)
        self.assertEqual(c.torus_rank(511),2)
    def test_pair_source_equivalence(self):
        self.assertEqual(c.rank_response_controls()['normalized_measure_equalities'],4608)
    def test_normalization_repair(self):
        r=c.repaired_table()
        for row in r['rows']:
            self.assertGreater(row['physical_black_N'],0)
            self.assertLess(abs(row['pair_normal_difference']),F(2,10**15))
    def test_spinful_zero_mode(self):
        self.assertEqual(c.spin_zero_mode_control()['chiral_spin4_I3_eigenvalue'],F(-55,9216))

if __name__=='__main__': unittest.main()

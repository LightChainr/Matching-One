"""Concrete false-verdict regressions and full-rank Fieller controls."""
import math
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import projective_inference as p


class ProjectiveInferenceTests(unittest.TestCase):
    Y=[9.01643304e-4,2.91097703e-3,4.13180763e-3]
    VAR=[6.207021e-8,6.713835e-8,3.807515e-8]
    S=[[1,0,0],[0,1,0],[0,0,0]]

    def pair(self):
        c=-.1648*math.sqrt(self.VAR[0]*self.VAR[2])
        return [self.Y[0],self.Y[2]],[[self.VAR[0],c],[c,self.VAR[2]]]

    def test_fieller_identity(self):
        y,s=self.pair()
        for r in (1,4,10.9908008589,16,120.79770352,2080.30719731):
            z=(y[1]-r*y[0])/math.sqrt(s[1][1]-2*r*s[0][1]+r*r*s[0][0])
            got=p.ray_residual(y,s,[1,r])
            self.assertEqual(got['degrees_of_freedom'],1)
            self.assertAlmostEqual(got['statistic'],z*z,delta=1e-9*z*z+1e-12)

    def test_signed_line_scaling(self):
        y,s=self.pair()
        reference=p.ray_residual(y,s,[1,4])['statistic']
        for c in (1e-6,.5,7,1e6,-3):
            self.assertAlmostEqual(p.ray_residual(y,s,[c,4*c])['statistic'],reference,
                                   delta=1e-9*reference+1e-12)

    def test_exact_fit(self):
        y,s=self.pair()
        self.assertLess(p.ray_residual(y,s,y)['statistic'],1e-20)

    def test_nuisance_lowers_df_and_distance(self):
        y=[1,3,4];s=[[1,0,0],[0,1,0],[0,0,1]]
        one=p.subspace_residual(y,s,[[1,2,4]])
        two=p.subspace_residual(y,s,[[1,2,4],[-1,1,-1]])
        self.assertEqual((one['degrees_of_freedom'],two['degrees_of_freedom']),(2,1))
        self.assertLessEqual(two['statistic'],one['statistic'])

    def test_dependent_full_rank_design_refused(self):
        y,s=self.pair()
        with self.assertRaisesRegex(ValueError,'linearly dependent'):
            p.subspace_residual(y,s,[[1,4],[2,8]])

    def test_singular_covariance_requires_semantics(self):
        with self.assertRaisesRegex(ValueError,'must be declared'):
            p.ray_residual([0,0,1],self.S,[1,0,0])

    def test_exact_support_contradiction_is_not_p_one(self):
        with self.assertRaisesRegex(ValueError,'incompatible'):
            p.ray_residual([0,0,1],self.S,[1,0,0],covariance_mode='exact_support')

    def test_exact_support_constrains_amplitude_and_df(self):
        got=p.ray_residual([0,0,1],self.S,[1,0,1],covariance_mode='exact_support')
        self.assertEqual(got['amplitudes'],[1])
        self.assertEqual(got['statistic'],1)
        self.assertEqual(got['degrees_of_freedom'],2)

    def test_known_duplicate_coordinate_support(self):
        s=[[1,1,0],[1,1,0],[0,0,1]]
        got=p.ray_residual([1,1,2],s,[1,1,1],covariance_mode='exact_support')
        self.assertEqual(got['covariance_rank'],2)
        self.assertEqual(got['degrees_of_freedom'],1)
        self.assertAlmostEqual(got['statistic'],.5)

    def test_negative_covariance_refused(self):
        with self.assertRaisesRegex(ValueError,'positive semidefinite'):
            p.ray_residual([0,0,1],[[1,0,0],[0,1,0],[0,0,-.1]],[1,0,0])

    def test_cutoff_is_not_structural_null(self):
        with self.assertRaisesRegex(ValueError,'numerical cutoff'):
            p.ray_residual([0,0,1],[[1,0,0],[0,1,0],[0,0,1e-15]],
                           [1,0,1],covariance_mode='exact_support')

    def test_all_deterministic_compatible(self):
        got=p.ray_residual([2,4],[[0,0],[0,0]],[1,2],covariance_mode='exact_support')
        self.assertEqual(got['degrees_of_freedom'],0)
        self.assertEqual(got['statistic'],0)

    def test_chi_square_tail(self):
        for z in (1,2,3,5):
            self.assertAlmostEqual(p.chi_square_upper_tail(z*z,1)/math.erfc(z/math.sqrt(2)),1,places=9)
        for x in (.5,2,9):
            self.assertAlmostEqual(p.chi_square_upper_tail(x,2),math.exp(-x/2),places=12)

    def test_sigma_is_not_sqrt_for_multiple_df(self):
        self.assertAlmostEqual(p._equivalent_sigma(9,1),3,places=9)
        self.assertTrue(2<p._equivalent_sigma(9,2)<3)


if __name__ == '__main__':
    unittest.main()

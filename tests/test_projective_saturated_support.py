"""Prevent saturated p=0 artefacts without clipping positive-df residuals."""
import math
from pathlib import Path
import sys
import unittest
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import projective_inference as p


class SaturatedSupportTests(unittest.TestCase):
    def test_full_rank_saturation_is_an_algebraic_zero(self):
        for scale in (1e-20, 1., 1e20):
            with self.subTest(scale=scale):
                result = p.subspace_residual(
                    [.3*scale, .7*scale], [[1, .23], [.23, 2]],
                    [[1, .123], [.47, .9]])
                self.assertEqual(result['degrees_of_freedom'], 0)
                self.assertEqual(result['statistic'], 0.)
                self.assertEqual(result['p_value'], 1.)
                self.assertEqual(result['equivalent_sigma'], 0.)

    def test_exact_support_saturation(self):
        result = p.subspace_residual([.3, .3, .7],
                    [[1, 1, 0], [1, 1, 0], [0, 0, 2]],
                    [[1, 1, .123], [.47, .47, .9]], covariance_mode='exact_support')
        self.assertEqual(result['degrees_of_freedom'], 0)
        self.assertEqual(result['statistic'], 0.)
        self.assertEqual(result['p_value'], 1.)

    def test_zero_df_does_not_hide_support_incompatibility(self):
        with self.assertRaisesRegex(ValueError, 'incompatible'):
            p.subspace_residual([.3, .4, .7],
                [[1, 1, 0], [1, 1, 0], [0, 0, 2]],
                [[1, 1, .123], [.47, .47, .9]], covariance_mode='exact_support')

    def test_positive_df_residual_is_not_clipped(self):
        result = p.ray_residual([1, 1e-15], [[1, 0], [0, 1]], [1, 0])
        self.assertEqual(result['degrees_of_freedom'], 1)
        self.assertGreater(result['statistic'], 0.)
        self.assertAlmostEqual(result['statistic']/1e-30, 1.)

    def test_covariance_support_and_psd_regressions(self):
        s = [[1,0,0],[0,1,0],[0,0,0]]
        with self.assertRaisesRegex(ValueError, 'incompatible'):
            p.ray_residual([0,0,1],s,[1,0,0],covariance_mode='exact_support')
        result=p.ray_residual([0,0,1],s,[1,0,1],covariance_mode='exact_support')
        self.assertEqual((result['statistic'],result['degrees_of_freedom']), (1.,2))
        with self.assertRaisesRegex(ValueError, 'positive semidefinite'):
            p.ray_residual([0,0,1],[[1,0,0],[0,1,0],[0,0,-.1]],[1,0,0])

    def test_fieller_identity_unchanged(self):
        y = [.3, .7]; s = [[1, .23], [.23, 2]]
        for r in (1., 4., 16., -2.):
            expected = (y[1]-r*y[0])**2/(s[1][1]-2*r*s[0][1]+r*r*s[0][0])
            self.assertAlmostEqual(p.ray_residual(y,s,[1,r])['statistic'], expected, places=14)

if __name__ == '__main__':
    unittest.main()

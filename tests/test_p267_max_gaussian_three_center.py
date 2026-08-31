import sys
from pathlib import Path
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"scripts"))
from p267_max_gaussian_three_center import heat_moments, rank_step_moments, realization


class ThreeCenterMomentTests(unittest.TestCase):
    def test_known_common_gaussian_mixture(self):
        centers = np.array([-1., 0., 1.])/np.sqrt(.7)
        weights, t = np.array([.25, .5, .25]), .2/.7
        discrete = np.array([weights@centers**r for r in range(9)])
        moments = heat_moments(discrete, t)
        found_t, found_centers, found_weights, prediction, _ = realization(moments)
        self.assertAlmostEqual(found_t, t, places=11)
        np.testing.assert_allclose(found_centers, centers, atol=1e-11)
        np.testing.assert_allclose(found_weights, weights, atol=1e-11)
        np.testing.assert_allclose(prediction, moments, atol=1e-10)

    def test_integrated_uniform_rank_bins(self):
        f = np.r_[np.ones(17), 0.]
        moments, mean, variance, _ = rank_step_moments(f)
        self.assertAlmostEqual(mean, .5)
        self.assertAlmostEqual(variance, 1/12)
        np.testing.assert_allclose(moments, [1, 0, 1, 0, 9/5, 0, 27/7, 0, 9], atol=1e-12)


if __name__ == "__main__":
    unittest.main()

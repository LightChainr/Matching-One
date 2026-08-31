from pathlib import Path
import sys
import unittest
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from p267_clock_shape_quotient import contrasts, parameters, jackknife_parameters


class ClockQuotientTests(unittest.TestCase):
    def test_common_offset_and_clock_profile_cancel(self):
        base = np.array([2., -1., 4., 3.])
        d = np.array([3., 7., -2., 5.])
        data = base+np.outer([0., 1., -.3], d)
        values = parameters(data.ravel())
        self.assertAlmostEqual(values[0], -.3)
        np.testing.assert_allclose(values[4:7], 0, atol=1e-14)

    def test_shear_is_exact_and_odd_area_gauge_is_zero(self):
        d = np.array([3., 7., -2., 5.])
        residual = np.array([.1, -.2, 0., .4])
        y = np.outer([0., 1., -.3], d)+np.outer([0., 0., 1.], residual)
        values = parameters(y.ravel())
        np.testing.assert_allclose(values[4:7], residual[[0, 1, 3]])
        dd, uu = contrasts(y.ravel())
        self.assertAlmostEqual(-2*(uu[2]-values[0]*dd[2]), 0)

    def test_jackknife_keeps_ratio_covariance(self):
        mean = np.array([[0., 0., 0., 0.], [3., 7., -2., 5.],
                         [-.8, -2.3, .6, -1.1]]).ravel()
        perturbation = np.arange(12)/10000
        batches = np.array([mean-perturbation, mean+perturbation,
                            mean-2*perturbation, mean+2*perturbation])
        point, cov = jackknife_parameters(mean, batches)
        self.assertEqual(cov.shape, (10, 10))
        self.assertGreater(cov[0, 0], 0)
        self.assertTrue(np.isfinite(point).all() and np.isfinite(cov).all())


if __name__ == "__main__":
    unittest.main()

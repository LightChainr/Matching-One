from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"scripts"))
from p334_conditional_thermal_averaging import weighted_covariance


class ConditionalThermalTests(unittest.TestCase):
    def test_tail_variance_is_not_Bernoulli_variance(self):
        # Equal probability of two already-averaged readouts .2 and .8.
        values = np.array([[.2, 1., .4], [.8, 1., .6]])
        mean, cov = weighted_covariance(values, np.array([.5, .5]))
        np.testing.assert_allclose(mean, [.5, 1., .5])
        np.testing.assert_allclose(cov, [[.09, 0, .03], [0, 0, 0], [.03, 0, .01]], atol=1e-15)
        self.assertNotEqual(cov[0, 0], mean[0]*(1-mean[0]))


if __name__ == "__main__":
    unittest.main()

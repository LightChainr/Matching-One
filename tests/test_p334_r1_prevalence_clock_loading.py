from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"scripts"))
from p334_r1_prevalence_clock_loading import decompose, score_batch_means


class R1SymmetricDecomposition(unittest.TestCase):
    def test_identity_orientation_swap_and_shared_LOO(self):
        mean = np.array([.4, .3, .12, .06, .08, .045])
        value = decompose(mean, .8)
        swapped = decompose(mean[[1, 0, 3, 2, 5, 4]], -.8)
        for start in (6, 13):
            self.assertAlmostEqual(value[start]+value[start+1], value[start+2])
            np.testing.assert_allclose(value[start:start+3], swapped[start:start+3])
        perturbation = np.arange(-9.5, 10, 1)[:, None]*np.array([.001, -.0005, .0003, .0001, .0002, -.0001])
        result = score_batch_means(mean+perturbation, .8)
        self.assertLess(result["max_LOO_additive_identity_residual"], 1e-15)
        cov = np.array(result["full_covariance"])
        for start in (6, 13):
            self.assertAlmostEqual(cov[start+2, start+2], cov[start, start]+cov[start+1, start+1]+2*cov[start, start+1])


if __name__ == "__main__":
    unittest.main()

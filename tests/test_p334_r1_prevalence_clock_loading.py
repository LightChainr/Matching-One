from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"scripts"))
from p334_r1_prevalence_clock_loading import (decompose, four_state_variance,
                                             risk_pair_sums, score_batch_means)


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

    def test_four_state_total_covariance_identity(self):
        risks = np.tile([[0, 0], [0, 1], [1, 0], [1, 1]], (5, 1))
        perturbation = np.arange(20)*.001
        y = np.column_stack((risks[:, 0]*(.1+perturbation), risks[:, 1]*(.11-perturbation),
                             risks[:, 0]*.37, risks[:, 1]*.37))
        point, details = four_state_variance(risk_pair_sums(risks, y, -.8))
        contrast = np.column_stack((y[:, 0]-y[:, 1], y[:, 2]-y[:, 3]))/-.8
        np.testing.assert_allclose(details["total_individual_covariance"], np.cov(contrast, rowvar=False, ddof=0))
        self.assertAlmostEqual(point[7], 1.)
        self.assertGreater(point[1], 0.)


if __name__ == "__main__":
    unittest.main()

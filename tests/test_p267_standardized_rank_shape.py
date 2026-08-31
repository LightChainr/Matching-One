from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"scripts"))
from p267_standardized_rank_shape import beta_raw_moments, shape_moments


class StandardizedShapeTests(unittest.TestCase):
    def test_uniform_cell_and_symmetric_beta(self):
        f = np.array([0., 0., 1., 0., 0.])
        step, canonical, _, _ = shape_moments(f, beta_raw_moments(4))
        np.testing.assert_allclose(step, [0, 9/5, 0, 27/7], atol=1e-12)
        np.testing.assert_allclose(canonical, [0, 7/3, 0, 245/33], atol=2e-11)

    def test_rank_step_invariance_under_free_location_width_and_amplitude(self):
        f = np.array([0., 2., 1., 5., 3., 0., 0., 0., 0.])
        transformed = np.zeros(17)
        # Same ordered step pattern shifted by 1/4 and narrowed by factor 2.
        transformed[4:12] = -3*f[:8]
        old = shape_moments(f, beta_raw_moments(8))[0]
        new = shape_moments(transformed, beta_raw_moments(16))[0]
        np.testing.assert_allclose(old, new, atol=2e-12)


if __name__ == "__main__":
    unittest.main()

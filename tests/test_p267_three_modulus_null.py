from fractions import Fraction as F
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from p267_three_modulus_null import (Q2, SHAPES, WEIGHTS, dot, certificate,
                                    normalized_weights, project_joint, slope)
from p267_third_geometry_certificate import shape_coordinate


class ThreeModulusTests(unittest.TestCase):
    def test_exact_offset_and_amplitude_cancellation(self):
        for model, shape in SHAPES.items():
            response = [Q2(F(17, 13), F(-2, 7))+Q2(F(3, 8), F(5, 9))*g for g in shape]
            self.assertEqual(dot(WEIGHTS[model], response), Q2())
        self.assertEqual(len(certificate()["pairwise_discriminators"]), 3)

    def test_competing_models_intersect_only_on_flat_response(self):
        for truth, shape in SHAPES.items():
            for tested in SHAPES:
                residual = dot(normalized_weights(tested), shape)
                self.assertEqual(residual == Q2(), truth == tested)
                self.assertEqual(residual, (slope(truth)-slope(tested))*(shape[1]-shape[0]))

    def test_exact_cm_values_agree_with_independent_qseries(self):
        for expected, tau in zip(SHAPES["affine_E4"], [(F(0), F(2)), (F(0), F(4)), (F(1, 2), F(1))]):
            self.assertAlmostEqual(expected.as_float(), shape_coordinate(tau), places=13)

    def test_full_covariance_projection_keeps_common_noise_cancellation(self):
        g = SHAPES["affine_E4"]
        mean = [3*j+(j+1)*value.as_float() for value in g for j in range(4)]
        covariance = [[1.0+(0.2 if i == j else 0.0) for j in range(12)] for i in range(12)]
        residual, cov = project_joint(mean, covariance, "affine_E4")
        self.assertTrue(all(abs(x) < 1e-12 for x in residual))
        variance = 0.2*sum(w.as_float()**2 for w in normalized_weights("affine_E4"))
        for i in range(4):
            for j in range(4):
                self.assertAlmostEqual(cov[i][j], variance if i == j else 0.0, places=12)


if __name__ == "__main__":
    unittest.main()

from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path


SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

from control_variate_estimator import (  # noqa: E402
    FrozenEstimator,
    FrozenZeroMeanControls,
    minimum_variance_weights,
    sample_covariance,
)


class MinimumVarianceTests(unittest.TestCase):
    def test_known_diagonal_solution(self) -> None:
        weights, ridge = minimum_variance_weights([[4.0, 0.0], [0.0, 1.0]])
        self.assertEqual(ridge, 0.0)
        self.assertAlmostEqual(weights[0], 0.2)
        self.assertAlmostEqual(weights[1], 0.8)
        self.assertAlmostEqual(math.fsum(weights), 1.0)

    def test_singular_covariance_is_regularized(self) -> None:
        weights, ridge = minimum_variance_weights([[1.0, 1.0], [1.0, 1.0]])
        self.assertGreater(ridge, 0.0)
        self.assertAlmostEqual(weights[0], 0.5, places=6)
        self.assertAlmostEqual(weights[1], 0.5, places=6)

    def test_tiny_covariance_uses_relative_pivot_tolerance(self) -> None:
        weights, ridge = minimum_variance_weights([[4e-30, 0.0], [0.0, 1e-30]])
        self.assertEqual(ridge, 0.0)
        self.assertAlmostEqual(weights[0], 0.2)
        self.assertAlmostEqual(weights[1], 0.8)

    def test_large_covariance_does_not_use_absolute_normalizer_cutoff(self) -> None:
        weights, ridge = minimum_variance_weights([[4e30, 0.0], [0.0, 1e30]])
        self.assertEqual(ridge, 0.0)
        self.assertAlmostEqual(weights[0], 0.2)
        self.assertAlmostEqual(weights[1], 0.8)

    def test_covariance_and_pilot_frozen_evaluation(self) -> None:
        pilot = [[-2.0, -1.0], [-2.0, 1.0], [2.0, -1.0], [2.0, 1.0]]
        covariance = sample_covariance(pilot)
        self.assertAlmostEqual(covariance[0][1], 0.0)
        self.assertAlmostEqual(covariance[0][0] / covariance[1][1], 4.0)

        estimator = FrozenEstimator.fit(("e", "c"), pilot)
        self.assertAlmostEqual(estimator.weights[0], 0.2)
        self.assertAlmostEqual(estimator.weights[1], 0.8)
        result = estimator.evaluate([[1.0, 3.0], [3.0, 1.0]])
        self.assertEqual(result["samples"], 2)
        self.assertAlmostEqual(result["mean"], 2.0)
        self.assertAlmostEqual(result["sample_variance"], 0.72)
        self.assertAlmostEqual(result["standard_error"], 0.6)

    def test_requires_independent_evaluation_sample_size(self) -> None:
        estimator = FrozenEstimator.fit(("a", "b"), [[0.0, 1.0], [1.0, 0.0]])
        with self.assertRaisesRegex(ValueError, "at least two rows"):
            estimator.evaluate([[1.0, 1.0]])

    def test_frozen_zero_mean_control_regression(self) -> None:
        controls = []
        targets = []
        for first in (-1.0, 1.0):
            for second in (-1.0, 1.0):
                for noise in (-1.0, 1.0):
                    controls.append([first, second])
                    targets.append(3.0 * first - 2.0 * second + 0.5 * noise)
        estimator = FrozenZeroMeanControls.fit(
            ("first", "second"), targets, controls
        )
        self.assertAlmostEqual(estimator.coefficients[0], -3.0)
        self.assertAlmostEqual(estimator.coefficients[1], 2.0)
        result = estimator.evaluate(targets, controls)
        self.assertAlmostEqual(result["mean"], 0.0)
        self.assertLess(result["sample_variance"], 0.3)


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
from __future__ import annotations

import math
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from analyze_norm4_variance_pilot import (  # noqa: E402
    covariance_of_mean,
    jackknife_se,
    quadratic_2,
    transform_covariance,
    u_value,
)


class Norm4VariancePilotTests(unittest.TestCase):
    def test_u_definition(self) -> None:
        state = {"P4_S_prime": 2.0, "mean_slope": 4.0}
        self.assertAlmostEqual(u_value(16, state), 45.254833995939045)

    def test_jackknife_se_uses_pseudovalue_normalization(self) -> None:
        self.assertAlmostEqual(jackknife_se(1.0, [0.9, 1.1]), 0.1)

    def test_covariance_and_two_dimensional_quadratic(self) -> None:
        covariance = covariance_of_mean([[1.0, 2.0], [3.0, 0.0]])
        self.assertEqual(covariance, [[1.0, -1.0], [-1.0, 1.0]])
        self.assertAlmostEqual(quadratic_2([2.0, 3.0], [[4.0, 0.0], [0.0, 9.0]]), 2.0)

    def test_covariance_transform(self) -> None:
        matrix = [[4.0, 1.0], [1.0, 9.0]]
        transformed = transform_covariance(matrix, [[1.0, 2.0]])
        self.assertEqual(transformed, [[44.0]])


if __name__ == "__main__":
    unittest.main()

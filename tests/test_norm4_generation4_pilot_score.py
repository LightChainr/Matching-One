#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_norm4_generation4_pilot import (  # noqa: E402
    recurrence_matrix,
    transformed_score,
)


class Norm4Generation4PilotScoreTests(unittest.TestCase):
    def test_half_recurrence_coefficients(self) -> None:
        matrix = recurrence_matrix(1, 0.5)
        self.assertEqual(matrix[0][:4], [-0.5, 2.0, -2.5, 1.0])
        self.assertEqual(matrix[1][4:], [-0.5, 2.0, -2.5, 1.0])

    def test_exact_half_mode_scores_zero(self) -> None:
        vector = [1.0, 3.0, 6.0, 9.5, 2.0, 4.0, 7.0, 10.5]
        covariance = [[float(i == j) for j in range(8)] for i in range(8)]
        score = transformed_score(vector, covariance, 1, 0.5, ["a", "b"])
        self.assertTrue(all(abs(value) < 1e-12 for value in score["residual"]))


if __name__ == "__main__":
    unittest.main()

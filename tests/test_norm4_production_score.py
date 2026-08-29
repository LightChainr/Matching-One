#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_norm4_production import (  # noqa: E402
    METRICS,
    SIZE_ORDER,
    row,
    solve_root,
    transform_covariance,
)


class Norm4ProductionScoreTests(unittest.TestCase):
    def test_root_solver(self) -> None:
        self.assertAlmostEqual(solve_root(lambda x: x - 0.6), 0.6)

    def test_metric_row_places_coefficients(self) -> None:
        result = row("U", (1, -3, 2, 0, 0, 0))
        self.assertEqual(len(result), len(SIZE_ORDER) * len(METRICS))
        self.assertEqual(result[0::3], [1, -3, 2, 0, 0, 0])

    def test_covariance_transform(self) -> None:
        covariance = [[4.0, 1.0], [1.0, 9.0]]
        self.assertEqual(transform_covariance(covariance, [[1.0, 2.0]]), [[44.0]])


if __name__ == "__main__":
    unittest.main()

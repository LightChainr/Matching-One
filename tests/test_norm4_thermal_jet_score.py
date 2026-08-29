#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_norm4_thermal_jet import ORDERS, cocycle_transform  # noqa: E402


class Norm4ThermalJetScoreTests(unittest.TestCase):
    def test_q2_transform_is_one_common_generator(self) -> None:
        matrix = cocycle_transform(1.5)
        self.assertEqual(len(matrix), 2 * len(ORDERS))
        self.assertEqual(matrix[0][0], 0.5)
        self.assertEqual(matrix[0][5], -1.5)
        self.assertEqual(matrix[0][10], 1.0)
        self.assertEqual(sum(value != 0.0 for value in matrix[0]), 3)

    def test_jordan_transform(self) -> None:
        row = cocycle_transform(2.0)[-1]
        self.assertEqual(row[19], 1.0)
        self.assertEqual(row[24], -2.0)
        self.assertEqual(row[29], 1.0)


if __name__ == "__main__":
    unittest.main()

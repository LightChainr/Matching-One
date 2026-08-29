#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from analyze_c4_self_matching_tangent_mc import (  # noqa: E402
    determinant,
    singular_values,
)
from analyze_c4_tangent_orientation_pair import (  # noqa: E402
    channel_matrix,
    row_normalized_condition,
)


class C4SelfMatchingTangentMCTests(unittest.TestCase):
    def test_determinant(self):
        self.assertEqual(determinant([[2.0, 1.0], [3.0, 4.0]]), 5.0)

    def test_singular_values_detect_rank_one(self):
        values = singular_values([[1.0, 2.0], [2.0, 4.0]])
        self.assertAlmostEqual(values[0], 5.0)
        self.assertAlmostEqual(values[1], 0.0)

    def test_orientation_matrix_and_scale_free_condition(self):
        matrix = channel_matrix([2.0, 1.0], [1.0, 2.0], 0, 10, 2.0)
        self.assertEqual(matrix, [[1.5, 1.5], [5.0, -5.0]])
        condition, normalized = row_normalized_condition(matrix)
        self.assertAlmostEqual(condition, 1.0)
        self.assertAlmostEqual(sum(value * value for value in normalized[0]), 1.0)


if __name__ == "__main__":
    unittest.main()

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
from score_c4_tangent_orthogonal_holdout import (  # noqa: E402
    THERMAL_RATIO,
    statistics,
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

    def test_orthogonalization_uses_source_ratio(self):
        source = [
            {"batch": "0", "samples": "1", "cross_score_t": "2", "cross_score_lambda": "1", **{f"{c}_score_{p}": "0" for c in ("both", "either", "direction_0", "direction_1") for p in ("t", "lambda")}},
            {"batch": "1", "samples": "1", "cross_score_t": "2", "cross_score_lambda": "1", **{f"{c}_score_{p}": "0" for c in ("both", "either", "direction_0", "direction_1") for p in ("t", "lambda")}},
        ]
        target = [dict(row, cross_score_t="4", cross_score_lambda="2") for row in source]
        result = statistics(source, target)
        self.assertEqual(result["source_c"], 0.5)
        self.assertEqual(result["orthogonal_residual"], 0.0)
        self.assertGreater(THERMAL_RATIO, 1.0)


if __name__ == "__main__":
    unittest.main()

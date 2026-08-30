#!/usr/bin/env python3
"""Exact geometry gates for the P321 equal-area rectangle design."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "p321_equal_area_rectangle_design.py"
SPEC = importlib.util.spec_from_file_location("p321_design", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class P321EqualAreaRectangleDesignTests(unittest.TestCase):
    def test_each_scale_has_one_common_area_and_five_shapes(self) -> None:
        rows = MODULE.build_design()["rows"]
        for scale in (1, 2, 3):
            selected = [row for row in rows if row["scale"] == scale]
            self.assertEqual(len(selected), 5)
            self.assertEqual({row["N"] for row in selected}, {144 * scale * scale})
            self.assertEqual(
                {row["aspect_ratio"] for row in selected},
                {"1/1", "16/9", "9/4", "4/1", "9/1"},
            )

    def test_matrices_have_declared_positive_determinant(self) -> None:
        for row in MODULE.build_design()["rows"]:
            matrix = row["period_matrix_row_major"]
            determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
            self.assertEqual(determinant, row["N"])

    def test_width_amplitude_conversion_is_frozen(self) -> None:
        contract = MODULE.build_design()["scaling_contract"]
        self.assertIn("C_width(rho)=C_N(rho)/rho^2", contract["width_amplitude_conversion"])


if __name__ == "__main__":
    unittest.main()

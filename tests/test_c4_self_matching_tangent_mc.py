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


class C4SelfMatchingTangentMCTests(unittest.TestCase):
    def test_determinant(self):
        self.assertEqual(determinant([[2.0, 1.0], [3.0, 4.0]]), 5.0)

    def test_singular_values_detect_rank_one(self):
        values = singular_values([[1.0, 2.0], [2.0, 4.0]])
        self.assertAlmostEqual(values[0], 5.0)
        self.assertAlmostEqual(values[1], 0.0)


if __name__ == "__main__":
    unittest.main()

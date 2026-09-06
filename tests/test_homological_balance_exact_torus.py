"""Wrong number this would stop us believing: dual_fail != 0 on the L=3 honest torus,
or p_L^H leaving (0.5865, 0.5866), or the two winding algorithms disagreeing."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "homological_balance"))

import exact_torus_enum as enum  # noqa: E402
import verify_independent as v2  # noqa: E402


class HomologicalBalanceExactTorusTests(unittest.TestCase):
    def test_l3_alexander_census_and_root(self) -> None:
        coeff, dual_fail, rank_pairs = enum.enumerate_ML(3)
        self.assertEqual(dual_fail, 0)
        self.assertEqual(rank_pairs[(0, 2)], 259)
        self.assertEqual(rank_pairs[(1, 1)], 162)
        self.assertEqual(rank_pairs[(2, 0)], 91)
        self.assertEqual(enum.eval_ML(coeff, 0), -1)
        self.assertEqual(enum.eval_ML(coeff, 1), 1)
        half = enum.eval_ML(coeff, enum.Fraction(1, 2))
        self.assertEqual(half, enum.Fraction(-21, 64))
        root = enum.find_root(coeff)
        self.assertGreater(root, 0.5865)
        self.assertLess(root, 0.5866)

    def test_two_algorithms_agree_on_l3(self) -> None:
        self.assertEqual(v2.cross_check(3), 0)

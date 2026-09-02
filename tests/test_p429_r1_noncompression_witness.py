#!/usr/bin/env python3
"""Independent exact-rational check of the n=7 r=1 non-compression witness."""
from __future__ import annotations

import sys
import unittest
from fractions import Fraction
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SYS_DIR = ROOT / "research" / "summary_search"
if str(SYS_DIR) not in sys.path:
    sys.path.insert(0, str(SYS_DIR))

from bounded_summary_search import (  # noqa: E402
    S_coeffs,
    connected_carrier,
    experiments,
    is_plane_two_terminal,
    neighborhood_key,
    safe_table,
)
from verify_witness import A, B, connecting_ksets  # noqa: E402


class TestR1NoncompressionWitness(unittest.TestCase):
    def test_carrier_and_planarity(self):
        for g in (A, B):
            self.assertEqual(g.n, 7)
            self.assertTrue(connected_carrier(g))
            self.assertTrue(is_plane_two_terminal(g))
            self.assertEqual(g.adj_L.bit_count(), 1)
            self.assertEqual(g.adj_R.bit_count(), 1)
            self.assertFalse(g.lr_edge)

    def test_identical_safe_polynomial_and_r1(self):
        tA, tB = safe_table(A), safe_table(B)
        sA, sB = S_coeffs(tA, 7), S_coeffs(tB, 7)
        self.assertEqual(sA, (1, 7, 21, 35, 33, 15, 2, 0))
        self.assertEqual(sB, sA)
        self.assertEqual(sA[1], 7)
        self.assertEqual(sA[2], comb(7, 2))
        self.assertEqual(neighborhood_key(A, 1), neighborhood_key(B, 1))
        self.assertNotEqual(neighborhood_key(A, 2), neighborhood_key(B, 2))

    def test_e2c2_split_and_other_experiments_agree(self):
        eA = experiments(safe_table(A), 7)
        eB = experiments(safe_table(B), 7)
        for k in (
            "E0_c1",
            "E0_c2",
            "E0_mix",
            "E1_c1",
            "E1_c2",
            "E1_mix",
            "E2_c1",
            "E2_mix",
        ):
            self.assertEqual(eA[k], eB[k], k)
        self.assertEqual(eA["E2_c2"], Fraction(937, 1050))
        self.assertEqual(eB["E2_c2"], Fraction(313, 350))
        self.assertEqual(abs(eA["E2_c2"] - eB["E2_c2"]), Fraction(1, 525))
        self.assertEqual(eA["E1_c1"], 1)
        self.assertEqual(eB["E1_c1"], 1)

    def test_mincut_intersection_pattern(self):
        cutsA = connecting_ksets(safe_table(A), 7, 4)
        cutsB = connecting_ksets(safe_table(B), 7, 4)
        self.assertEqual(len(cutsA), 2)
        self.assertEqual(len(cutsB), 2)
        self.assertEqual(set(cutsA[0]) & set(cutsA[1]), {5, 6})
        self.assertEqual(set(cutsB[0]) & set(cutsB[1]), {0, 5, 6})


if __name__ == "__main__":
    unittest.main()

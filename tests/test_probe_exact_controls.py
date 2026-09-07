"""C10 — exact-controls test suite (L=3 fast; L=4 slow-marked)."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path
from fractions import Fraction

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "homological_balance"))
sys.path.insert(0, str(ROOT / "scripts" / "probe"))

import exact_torus_enum as enum  # noqa: E402
import verify_independent as v2  # noqa: E402


class ProbeExactControlsTests(unittest.TestCase):
    # ---- C1 ----
    def test_c1_l3_census(self) -> None:
        coeff, dual_fail, rank_pairs = enum.enumerate_ML(3)
        self.assertEqual(dual_fail, 0)
        self.assertEqual(rank_pairs[(0, 2)], 259)
        self.assertEqual(rank_pairs[(1, 1)], 162)
        self.assertEqual(rank_pairs[(2, 0)], 91)
        self.assertEqual(enum.eval_ML(coeff, 0), -1)
        self.assertEqual(enum.eval_ML(coeff, 1), 1)
        self.assertEqual(enum.eval_ML(coeff, Fraction(1, 2)), Fraction(-21, 64))
        self.assertAlmostEqual(enum.find_root(coeff), 0.586511455113, places=9)

    # ---- C2 ----
    def test_c2_l3_M_polynomial(self) -> None:
        from exact_controls_q_and_m import enumerate_joint, bernstein_from_joint, bernstein_to_power, poly_eval
        N = 9
        jb, _ = enumerate_joint(3)
        coeff = bernstein_from_joint(jb, N)
        a = bernstein_to_power(coeff)
        self.assertEqual(poly_eval(a, 0), -1)
        self.assertEqual(poly_eval(a, 1), 1)
        self.assertEqual(poly_eval(a, Fraction(1, 2)), Fraction(-21, 64))

    # ---- C3 ----
    def test_c3_quantiles_and_self_symmetry_kill(self) -> None:
        from exact_controls_q_and_m import run_L
        r = run_L(3)
        # Q(0.5) == p_L^H within 1e-12
        self.assertAlmostEqual(r["quantiles"]["0.5"], r["p_L_H"], places=12)
        # self-symmetry kill: max |Q(u)+Q(1-u)-1| must NOT be ~0
        mx = max(abs(v) for v in r["Q_plus_Q1minus_1"].values())
        self.assertGreater(mx, 1e-3)

    # ---- C5 ----
    def test_c5_two_algorithms_l3_mismatch_zero(self) -> None:
        self.assertEqual(v2.cross_check(3), 0)

    # ---- C6 ----
    def test_c6_wrap_identity_failure_zero(self) -> None:
        from square_bond_wrap_identity import bond_edges, bond_ambient_rank
        L = 3
        edges = bond_edges(L)
        E = len(edges)
        fail = 0
        for mask in range(1 << E):
            oe = [edges[k] for k in range(E) if (mask >> k) & 1]
            r = bond_ambient_rank(oe, L)
            X = r - 1
            O = 1 if r > 0 else 0
            if O != 1 + (X - X * X) // 2:
                fail += 1
        self.assertEqual(fail, 0)

    # ---- C8 ----
    def test_c8_matching_involution_l3(self) -> None:
        from matching_involution_poly import enumerate_joint, bernstein_from_joint, white_bernstein_from_joint
        N = 9
        jb, jw = enumerate_joint(3)
        cb = bernstein_from_joint(jb, N)
        cw = white_bernstein_from_joint(jw, N)
        for k in range(N + 1):
            self.assertEqual(cw[k], -cb[N - k])

    # ---- L=4 slow ----
    @unittest.skipUnless(False, "slow: L=4 census (optional, run manually)")
    def test_c1_l4_census_slow(self) -> None:
        coeff, dual_fail, rank_pairs = enum.enumerate_ML(4)
        self.assertEqual(dual_fail, 0)
        self.assertEqual(rank_pairs[(0, 2)], 36559)
        self.assertEqual(rank_pairs[(1, 1)], 19932)
        self.assertEqual(rank_pairs[(2, 0)], 9045)
        self.assertEqual(enum.eval_ML(coeff, Fraction(1, 2)), Fraction(-13757, 32768))
        self.assertAlmostEqual(enum.find_root(coeff), 0.590672112331, places=9)


if __name__ == "__main__":
    unittest.main(verbosity=2)

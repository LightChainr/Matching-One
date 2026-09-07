"""probe #625 test suite (L=3 fast; L=4 slow-marked).

Checks the numbers the probe reports, not PR #606's census (cited, and
already re-checked by #627's own tests).
"""
from __future__ import annotations

import sys
import unittest
from pathlib import Path
from fractions import Fraction

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "probe"))
sys.path.insert(0, str(ROOT / "scripts" / "homological_balance"))

from mhalf_common import (  # noqa: E402
    stream_joint, joint_to_bernstein, eval_bernstein, invert_F_bisect,
)


class ProbeMhalfVsShapeTests(unittest.TestCase):
    def test_d1_alexander_precondition_l3(self) -> None:
        counts, N = stream_joint(3)
        self.assertEqual(N, 9)
        self.assertEqual(sum(counts.values()), 512)
        # rank-pair totals must match the cited census
        rank_mass = {r: sum(c for (n, rr), c in counts.items() if rr == r)
                     for r in (0, 1, 2)}
        self.assertEqual(rank_mass[0], 259)
        self.assertEqual(rank_mass[1], 162)
        self.assertEqual(rank_mass[2], 91)

    def test_d2_cited_inputs_reproduced_l3(self) -> None:
        counts, N = stream_joint(3)
        coeff = joint_to_bernstein(counts, N)
        Mh = eval_bernstein(coeff, Fraction(1, 2))
        self.assertEqual(Mh, Fraction(-21, 64))
        F = lambda p: (1 + eval_bernstein(coeff, Fraction(p))) / 2
        self.assertEqual(F(Fraction(1, 2)), Fraction(43, 128))
        pLH = invert_F_bisect(F, Fraction(1, 2))
        self.assertLess(abs(float(pLH) - 0.5865114551126757), 1e-12)

    def test_d2_Z3_monotone_and_anchors(self) -> None:
        counts, N = stream_joint(3)
        coeff = joint_to_bernstein(counts, N)
        F = lambda p: (1 + eval_bernstein(coeff, Fraction(p))) / 2
        UG = [Fraction(k, 10) for k in range(1, 10)]
        Q = {u: invert_F_bisect(F, u) for u in UG}
        Z = {u: (Q[u] - Q[Fraction(2, 10)]) / (Q[Fraction(8, 10)] - Q[Fraction(2, 10)])
             for u in UG}
        self.assertEqual(Z[Fraction(2, 10)], 0)
        self.assertEqual(Z[Fraction(8, 10)], 1)
        # monotone in u
        prev = None
        for u in UG:
            if prev is not None:
                self.assertLess(prev, Z[u])
            prev = Z[u]
        # spot value from the note table
        self.assertLess(abs(float(Z[Fraction(5, 10)]) - 0.528759509), 1e-8)

    def test_d3_gameB_m_pinned_and_Z_moves(self) -> None:
        counts, N = stream_joint(3)
        rows = [(n, r, Fraction(c)) for (n, r), c in counts.items()]

        def F_factory(lam):
            def F(p):
                p = Fraction(p); q = 1 - p
                num = Fraction(0); den = Fraction(0)
                for (n, r, c) in rows:
                    w = c * (p ** n) * (q ** (N - n)) * (lam ** r)
                    num += w * (r - 1); den += w
                return (1 + num / den) / 2
            return F

        F1 = F_factory(Fraction(1))
        F2 = F_factory(Fraction(2))
        # m pinned at 0 on the p_mid line for both lambdas
        for F in (F1, F2):
            p_mid = invert_F_bisect(F, Fraction(1, 2))
            self.assertLess(abs(2 * F(p_mid) - 1), Fraction(1, 10 ** 30))
        # ... but the shape at interior u differs far beyond 1e-8
        for a, b in ((Fraction(1, 10), Fraction(9, 10)),
                     (Fraction(2, 10), Fraction(8, 10))):
            Qs = {}
            for F in (F1, F2):
                Qa = invert_F_bisect(F, a); Qb = invert_F_bisect(F, b)
                Qm = invert_F_bisect(F, Fraction(1, 2))
                Qs = Qm  # placeholder to keep names
                Zm = (Qm - Qa) / (Qb - Qa)
            # compare Z(1/2) between lambda=1 and lambda=2
            def Z_half(F):
                Qa = invert_F_bisect(F, a); Qb = invert_F_bisect(F, b)
                Qm = invert_F_bisect(F, Fraction(1, 2))
                return (Qm - Qa) / (Qb - Qa)
            dz = abs(Z_half(F1) - Z_half(F2))
            self.assertGreater(dz, Fraction(1, 100))

    def test_d4_bond_m_half_zero_and_3_atom(self) -> None:
        from bond_L3_Z import bond_edges, bond_ambient_rank
        L = 3
        edges = bond_edges(L)
        E = len(edges)
        counts = {}
        for mask in range(1 << E):
            n = bin(mask).count("1")
            oe = [edges[k] for k in range(E) if (mask >> k) & 1]
            key = (n, bond_ambient_rank(oe, L))
            counts[key] = counts.get(key, 0) + 1
        M_half = sum(c * (r - 1) * (Fraction(1, 2) ** n) * (Fraction(1, 2) ** (E - n))
                     for (n, r), c in counts.items())
        self.assertEqual(M_half, 0)
        rank_tot = [0, 0, 0]
        for (n, r), c in counts.items():
            rank_tot[r] += c
        self.assertEqual(rank_tot, [75460, 111224, 75460])  # == #627 C6
        # 3-atom law => staircase => Z undefined
        self.assertTrue(all(rank_tot[r] > 0 for r in range(3)))

    # ---- L=4 slow ----
    @unittest.skipUnless(False, "slow: L=4 joint (optional, run manually)")
    def test_d2_l4_cited_inputs_slow(self) -> None:
        counts, N = stream_joint(4)
        coeff = joint_to_bernstein(counts, N)
        self.assertEqual(eval_bernstein(coeff, Fraction(1, 2)),
                         Fraction(-13757, 32768))


if __name__ == "__main__":
    unittest.main(verbosity=2)

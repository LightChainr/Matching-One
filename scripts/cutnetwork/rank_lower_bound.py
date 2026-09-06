#!/usr/bin/env python3
"""Branching-Hankel rank lower bound for the #549 parallel-gadget family.

Protocol semantics (the direct generalisation of #549's published fork):

  * A bank of k future-vertex-disjoint gadgets, each of type A or B; the
    hidden predictive coordinate is a = #{A-type gadgets} in {0..k}.
  * A single "#549 fork" on a group of `g` gadgets is: one shared update
    vertex (uniform among the 8g future vertices of the group), then two
    independently continued clones; success = both clones safe.  Its exact
    success probability given the group's A-count a_g is

        F(g, a_g) = [343 g^3 - 182 g^2 + 25 g + 4 a_g] / [8g (8g-1)^2],

    which is AFFINE in a_g (slope 4/[8g(8g-1)^2] > 0).

  * A "grouped-fork" experiment E_{(g_1,...,g_j)} partitions the k gadgets
    into j disjoint groups and runs one independent fork per group, requiring
    all j forks to succeed.  The groups are future-disjoint, so conditional on
    the per-group A-counts the events are independent.

Theorem (verified exactly below).  For the experiment E_j with group sizes
(1,1,...,1,k-j+1) (j groups), the success probability P_j(a) is a polynomial
of exact degree j in a with nonzero leading coefficient.  Therefore the
(k+1) x (k+1) response matrix with columns E_0 (trivial), E_1, ..., E_k and
rows a = 0..k is Vandermonde-equivalent and has exact rank k+1.

Consequences (the two named major outcomes of the probe):
  A. branching-Hankel rank >= k+1 on the k+1 #549 predictive classes;
  B. for every fixed depth d (at most d independent forks), the response
     space is contained in span{1, a, ..., a^d}, so rank <= d+1 = C_d
     independent of k, while the class count k+1 grows without bound.
     Depth is thus a genuine complexity resource, and r_d(k) = d+1 for d <= k.

Deterministic exact rational arithmetic; no sampling.
"""

from fractions import Fraction
from itertools import product
from math import comb


def fork(k: int, a: int) -> Fraction:
    """#549 single-fork probability for a group of k gadgets with a A-type."""
    return Fraction(343 * k ** 3 - 182 * k ** 2 + 25 * k + 4 * a,
                    8 * k * (8 * k - 1) ** 2)


def group_fork_prob(gs, a: int) -> Fraction:
    """Success probability of grouped-fork experiment with group sizes gs,
    conditioned on total A-count a (per-group counts are hypergeometric)."""
    k = sum(gs)
    denom = comb(k, a)
    j = len(gs)
    total = Fraction(0)

    def enumerate_assignments(l, rem):
        if l == j - 1:
            if 0 <= rem <= gs[l]:
                yield (rem,)
            return
        for x in range(0, min(gs[l], rem) + 1):
            for rest in enumerate_assignments(l + 1, rem - x):
                yield (x,) + rest

    for assign in enumerate_assignments(0, a):
        w = Fraction(1)
        for gl, al in zip(gs, assign):
            w *= Fraction(comb(gl, al), 1)
        w /= Fraction(denom, 1)
        for gl, al in zip(gs, assign):
            w *= fork(gl, al)
        total += w
    return total


def poly_degree(vals) -> int:
    d = list(vals)
    for order in range(1, len(vals) + 1):
        d = [d[i + 1] - d[i] for i in range(len(d) - 1)]
        if all(x == 0 for x in d):
            return order - 1
    return len(vals) - 1


def rank_rows(rows) -> int:
    rows = [[Fraction(x) for x in r] for r in rows]
    m, n = len(rows), len(rows[0])
    rank = 0
    for col in range(n):
        piv = next((r for r in range(rank, m) if rows[r][col] != 0), None)
        if piv is None:
            continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        pv = rows[rank][col]
        rows[rank] = [x / pv for x in rows[rank]]
        for r in range(m):
            if r != rank and rows[r][col] != 0:
                f = rows[r][col]
                rows[r] = [rows[r][c] - f * rows[rank][c] for c in range(n)]
        rank += 1
    return rank


def main() -> None:
    import json
    from pathlib import Path

    out = {}
    for K in range(2, 9):
        cols = []
        degrees = []
        for j in range(K + 1):
            if j == 0:
                probs = [Fraction(1)] * (K + 1)
            else:
                gs = [1] * (j - 1) + [K - (j - 1)]
                probs = [group_fork_prob(gs, a) for a in range(K + 1)]
            cols.append(probs)
            if j > 0:
                degrees.append(poly_degree(probs))
        M = [[cols[j][a] for j in range(K + 1)] for a in range(K + 1)]
        rk = rank_rows(M)
        out[K] = {"degrees_E1_to_Ek": degrees, "rank": rk, "classes": K + 1}
        assert rk == K + 1, f"rank mismatch at K={K}"
        assert degrees == list(range(1, K + 1)), f"degree mismatch at K={K}"
        print(f"K={K}: degrees(E_1..E_K)={degrees}  rank(M)={rk} = classes")

    dest = Path(__file__).resolve().parents[2] / "results" / "cutnetwork-rank-vs-depth" / "latest.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out, indent=1))
    print("wrote", dest)


if __name__ == "__main__":
    main()

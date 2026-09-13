#!/usr/bin/env python3
"""#673 verification: brute-force re-derivation of the five-cell onsets at
COMMITTED sizes only (axis L=3,4; diamond L=2,3), using the repo's own
classifier (torus_homology). No new census: this re-verifies identities at
sizes whose tables are already committed (PRs #653/#657; diamond L=2 is
recomputable in seconds at 2^8 and is used only to test the *L=2 point* of
the closed forms, not as a new census rung).

Checks:
  I1: n*b = C(N,k) below onset            (axis k<L, diamond k<2L)
  I2: b*n = C(N,k) for k >= N-L+1, strict at k=N-L   (both geometries)
  axis I3a: d0(L) = d1(L) = L, and every such set is a full row/column
  axis I3b: b*n(2L-1) = L^2, b*n = 0 below, and every such set is a cross
  axis I4:  C - n*b = d0 + d1 on L <= k < 2L-1
  diamond d0 black side: #black-dir0-only sets at k=2L == L*C(2L,L)
  diamond b*b(2L) = 2L, and every such set is a straight diagonal line
"""
import itertools
import os
import sys
from collections import Counter
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, os.pardir, "scripts"))

from matched_torus_reference import axis_geometry, diamond_geometry
from torus_homology import classify_configuration

failures = []


def check(name, ok):
    print(("PASS" if ok else "FAIL"), name)
    if not ok:
        failures.append(name)


def type_pair(g, combo):
    active = [False] * g.n
    for i in combo:
        active[i] = True
    black, _ = classify_configuration(g, active)
    white, _ = classify_configuration(g, [not a for a in active], matching=True)
    bt = ("both" if black.both else "dir0" if black.direction_0 else
          "dir1" if black.direction_1 else "neither")
    wt = ("both" if white.both else "dir0" if white.direction_0 else
          "dir1" if white.direction_1 else "neither")
    return bt, wt


def cell_counts(g, k):
    counts = {}
    for combo in itertools.combinations(range(g.n), k):
        t = type_pair(g, combo)
        counts[t] = counts.get(t, 0) + 1
    return counts


def is_full_line_axis(g, combo, orientation):
    """orientation 'row' = fixed y, 'col' = fixed x; full = L sites."""
    coords = [g.coordinates[i] for i in combo]
    fixed = {c[1] for c in coords} if orientation == "row" else {c[0] for c in coords}
    return len(coords) == g.L and len(fixed) == 1


def run_axis(L):
    g = axis_geometry(L)
    N = L * L
    print(f"== axis L={L} (N={N})")
    ok = all(all(t == ("neither", "both") for t in
                 (type_pair(g, c) for c in itertools.combinations(range(g.n), k)))
             for k in range(L))
    check(f"I1: n*b=C(N,k) for k<{L}", ok)
    ok = all(all(t == ("both", "neither") for t in
                 (type_pair(g, c) for c in itertools.combinations(range(g.n), k)))
             for k in range(N - L + 1, N + 1))
    ok = ok and any(t != ("both", "neither") for t in cells_at(g, N - L))
    check(f"I2: b*n=C(N,k) for k>={N-L+1}, strict at {N-L}", ok)
    c = cell_counts(g, L)
    check(f"I3a: d0({L}) == {L}", c.get(("dir0", "dir0"), 0) == L)
    check(f"I3a: d1({L}) == {L}", c.get(("dir1", "dir1"), 0) == L)
    ok = True
    for combo in itertools.combinations(range(g.n), L):
        bt, wt = type_pair(g, combo)
        if bt == "dir0" and wt == "dir0":
            ok &= is_full_line_axis(g, combo, "row")
        if bt == "dir1" and wt == "dir1":
            ok &= is_full_line_axis(g, combo, "col")
    check("I3a: every d0 set is a full row and every d1 set a full column", ok)
    if 2 * L - 1 <= N:
        c = cell_counts(g, 2 * L - 1)
        check(f"I3b: b*n({2*L-1}) == L^2 == {L*L}",
              c.get(("both", "neither"), 0) == L * L)
        below = all(cell_counts(g, k).get(("both", "neither"), 0) == 0
                    for k in range(2 * L - 1))
        check(f"I3b: b*n = 0 for k < {2*L-1}", below)
        ok = True
        for combo in itertools.combinations(range(g.n), 2 * L - 1):
            bt, wt = type_pair(g, combo)
            if bt == "both" and wt == "neither":
                coords = [g.coordinates[i] for i in combo]
                cx = Counter(x for x, _ in coords)
                cy = Counter(y for _, y in coords)
                ok &= (max(cx.values()) == L and max(cy.values()) == L)
        check("I3b: every b*n set at k=2L-1 contains a full row and a full column",
              ok)
    ok = True
    for k in range(L, min(2 * L - 1, N + 1)):
        c = cell_counts(g, k)
        ok &= (c.get(("both", "both"), 0) == 0 and c.get(("both", "neither"), 0) == 0)
        ok &= (comb(N, k) - c.get(("neither", "both"), 0)
               == c.get(("dir0", "dir0"), 0) + c.get(("dir1", "dir1"), 0))
    check(f"I4: C-n*b = d0+d1 for {L} <= k < {2*L-1} (b*n = b*b = 0 there)", ok)


def is_diag_line_diamond(g, combo, slope):
    """slope 'm1': v-u constant (winds both via (1,1) steps);
    slope 'p1': v+u constant."""
    coords = [g.coordinates[i] for i in combo]
    if slope == "m1":
        vals = {(v - u) % (2 * g.L) for u, v in coords}
    else:
        vals = {(v + u) % (2 * g.L) for u, v in coords}
    return len(coords) == 2 * g.L and len(vals) == 1


def cells_at(g, k):
    return [type_pair(g, c) for c in itertools.combinations(range(g.n), k)]


def run_diamond(L):
    g = diamond_geometry(L)
    N = 2 * L * L
    print(f"== diamond L={L} (N={N})")
    ok = all(t == ("neither", "both") for k in range(2 * L)
             for t in cells_at(g, k))
    check(f"I1: n*b=C(N,k) for k<{2*L}", ok)
    ok = all(t == ("both", "neither") for k in range(N - L + 1, N + 1)
             for t in cells_at(g, k))
    ok = ok and any(t != ("both", "neither") for t in cells_at(g, N - L))
    check(f"I2: b*n=C(N,k) for k>={N-L+1}, strict at {N-L}", ok)
    cnt = 0
    bb = 0
    all_bb_are_diag = True
    for combo in itertools.combinations(range(g.n), 2 * L):
        bt, wt = type_pair(g, combo)
        if bt == "dir0" and wt == "dir0":
            cnt += 1
        if bt == "both" and wt == "both":
            bb += 1
            all_bb_are_diag &= (is_diag_line_diamond(g, combo, "m1")
                                or is_diag_line_diamond(g, combo, "p1"))
    check(f"diamond: #black dir0-only sets at k={2*L} == L*C(2L,L) == {L*comb(2*L,L)}",
          cnt == L * comb(2 * L, L))
    check(f"diamond: b*b({2*L}) == 2L == {2*L}", bb == 2 * L)
    check("diamond: every b*b set at k=2L is a straight diagonal line",
          all_bb_are_diag)
    c = cell_counts(g, 2 * L)
    check(f"diamond: d0 cell at k={2*L} == L*C(2L,L)",
          c.get(("dir0", "dir0"), 0) == L * comb(2 * L, L))
    check(f"diamond: d1 cell at k={2*L} == L*C(2L,L)",
          c.get(("dir1", "dir1"), 0) == L * comb(2 * L, L))
    ok = True
    for k in range(2 * L, min(3 * L - 1, N + 1)):
        c = cell_counts(g, k)
        ok &= (c.get(("both", "neither"), 0) == 0)
        ok &= (comb(N, k) - c.get(("neither", "both"), 0)
               == c.get(("dir0", "dir0"), 0) + c.get(("dir1", "dir1"), 0)
               + c.get(("both", "both"), 0))
    check(f"diamond I4: C-n*b = d0+d1+b*b for {2*L} <= k < {3*L-1}", ok)


if __name__ == "__main__":
    run_axis(3)
    run_axis(4)
    run_diamond(2)
    run_diamond(3)
    print("ALL CHECKS PASS" if not failures else f"FAILURES: {failures}")
    sys.exit(1 if failures else 0)

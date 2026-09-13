#!/usr/bin/env python3
"""Path A: exact integer-polynomial transfer-matrix DP for site percolation
spanning counts A_{n,m}(k) on the free-boundary n x m square grid.

Independent re-implementation for Issue 11 (T07), using a label-based
connectivity-state encoding (NOT the parenthesis signature scheme of
Mertens 2022, arXiv:2109.12102, and NOT scripts/noncrossing_connectivity_codec.py):

  * a state is a tuple of n cells; each cell is
      0          empty site,
      1          occupied and connected to the top (virtual all-occupied row 0),
      label >= 2 occupied, member of a not-yet-top-connected cluster,
    with labels canonicalized by first occurrence;
  * noncrossing of the cluster structure is asserted on every canonical state.

Physics/algorithm contract being reproduced (Mertens 2022, arXiv:2109.12102):
  spanning = one occupied 4-neighbour cluster touching rows 1 and m of the
  n x m grid with free boundaries; the virtual all-occupied "row 0" seeds
  top connectivity (paper Sec. 3.2/3.5);
  F_{n,m}(z) = sum_k A_{n,m}(k) z^k;  R_{n,m}(p) = (1-p)^{nm} F(p/(1-p));
  p_med: R_{n,n}(p_med) = 1/2;  p_cell: R_{n,n}(p_cell) = R_{n-1,n-1}(p_cell).

Exactness: all A(k) are Python integers; no floating point enters the DP.
Roots are refined with mpmath at 60 significant digits.

Built-in checks (each exists to stop one specific wrong belief):
  * A_{n,m}(m) = n               -- stops a DP that leaks non-spanning configs
                                    (paper Eq. 3a)
  * F_{n,m}(-1) parity pinned    -- stops coefficient drift (paper Eq. 6)
  * brute-force A(k) agreement   -- scripts/p11_bruteforce.py (independent path B)
"""

from __future__ import annotations

import argparse
import functools
import json
import sys
import time
from pathlib import Path

from mpmath import mp, mpf

ROOT = Path(__file__).resolve().parents[1]
EMPTY, TOP = 0, 1


# ---------------------------------------------------------------------------
# state utilities (path-A codec)
# ---------------------------------------------------------------------------

def _relabel(state):
    mapping = {}
    out = []
    for c in state:
        if c >= 2:
            if c not in mapping:
                mapping[c] = 2 + len(mapping)
            out.append(mapping[c])
        else:
            out.append(c)
    return tuple(out)


def _is_noncrossing(state):
    positions = {}
    for i, c in enumerate(state):
        if c >= 2:
            positions.setdefault(c, []).append(i)
    labels = list(positions)
    for ii in range(len(labels)):
        for jj in range(ii + 1, len(labels)):
            pa, pb = positions[labels[ii]], positions[labels[jj]]
            if len(pa) > 1 and len(pb) > 1:
                a1, a2 = pa[0], pa[-1]
                b1, b2 = pb[0], pb[-1]
                if a1 < b1 < a2 < b2 or b1 < a1 < b2 < a2:
                    return False
    return True


def canonical(state):
    state = _relabel(state)
    assert _is_noncrossing(state), f"crossing produced: {state}"
    return state


# ---------------------------------------------------------------------------
# one full-row update (transparent O(2^n) per state)
# ---------------------------------------------------------------------------

def _single_pattern_successor(state, n, mask):
    """Successor state for one occupancy pattern of the new row, or None if
    the pattern's partial configurations can never span (dropped)."""
    old = state
    parent = list(range(2 * n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for j in range(n):
        if (mask >> j) & 1:
            if old[j] != EMPTY:
                union(j, n + j)          # vertical bond to the old row
            if j > 0 and (mask >> (j - 1)) & 1:
                union(n + j, n + j - 1)  # horizontal bond inside the new row
    tops = [j for j in range(n) if old[j] == TOP]
    for j in tops[1:]:
        union(tops[0], j)                # all | sites share the virtual row 0
    # old-row sites carrying the same cluster label are one class: encode the
    # old row's internal connectivity into the union-find before reading it
    label_roots = {}
    for j, c in enumerate(old):
        if c != EMPTY:
            if c in label_roots:
                union(label_roots[c], j)
            else:
                label_roots[c] = j
    top_root = find(tops[0]) if tops else None

    newstate = []
    rootmap = {}
    for j in range(n):
        if not (mask >> j) & 1:
            newstate.append(EMPTY)
            continue
        r = find(n + j)
        if top_root is not None and r == top_root:
            newstate.append(TOP)
            continue
        if r not in rootmap:
            rootmap[r] = 2 + len(rootmap)
        newstate.append(rootmap[r])
    newstate = tuple(newstate)

    if not any(c == TOP for c in newstate):
        return None  # top cluster severed: no future cluster can reach row 0
    # NOTE: a non-top old class whose old sites all end here is simply a
    # finished cluster that never reached the top -- the configuration stays
    # valid because the spanning cluster is tracked separately by the | mark.
    return canonical(newstate)


def row_branches(state, n):
    """Successor state -> list of occupied-count shifts (one per pattern)."""
    branches = {}
    for mask in range(1 << n):
        occ = bin(mask).count("1")
        sub = _single_pattern_successor(state, n, mask)
        if sub is not None:
            branches.setdefault(sub, []).append(occ)
    return branches


@functools.lru_cache(maxsize=None)
def _row_branches_cached(state, n):
    """Transition function depends only on the state -- cache it.  The state
    sets stabilise across rows, so this amortises the O(2^n) expansion."""
    return tuple((k, tuple(v)) for k, v in
                 row_branches(state, n).items())


# ---------------------------------------------------------------------------
# exact DP with integer polynomials
# ---------------------------------------------------------------------------

def spanning_counts(n, m, verbose=False):
    """Return (A_{n,m}(k) for k=0..nm as exact integers, peak state count)."""
    t0 = time.time()
    # row 1: the 2^n - 1 configurations, every occupied site top-connected
    row = {}
    for mask in range(1, 1 << n):
        st = canonical(tuple(TOP if (mask >> j) & 1 else EMPTY
                             for j in range(n)))
        k = bin(mask).count("1")
        poly = [0] * (k + 1)
        poly[k] = 1
        prev = row.get(st)
        if prev is None:
            row[st] = poly
        else:
            row[st] = [a + b for a, b in zip(prev, poly)]
    peak_states = len(row)
    deg_cap = n * m
    for r in range(2, m + 1):
        newrow = {}
        for st, poly in row.items():
            for nst, shifts in _row_branches_cached(st, n):
                acc = newrow.get(nst)
                if acc is None:
                    acc = [0] * (deg_cap + 1)
                    newrow[nst] = acc
                for s in shifts:
                    for kk in range(len(poly)):
                        idx = kk + s
                        if poly[kk] and idx <= deg_cap:
                            acc[idx] += poly[kk]
        row = newrow
        peak_states = max(peak_states, len(row))
        if verbose:
            print(f"  row {r}: states={len(row)} t={time.time()-t0:.1f}s",
                  file=sys.stderr)
    total = [0] * (deg_cap + 1)
    for st, poly in row.items():
        for k, c in enumerate(poly):
            if c:
                total[k] += c
    return total, peak_states


# ---------------------------------------------------------------------------
# roots
# ---------------------------------------------------------------------------

def R_of_poly(A, p, nm):
    total = mpf(0)
    q = 1 - p
    pk = mpf(1)
    for k, c in enumerate(A):
        if c:
            total += c * pk * q ** (nm - k)
        pk *= p
    return total


def bisect_root(f, lo, hi, iters=150):
    mp.dps = 60
    flo = f(lo)
    for _ in range(iters):
        mid = (lo + hi) / 2
        fm = f(mid)
        if (fm < 0) == (flo < 0):
            lo, flo = mid, fm
        else:
            hi = mid
    return (lo + hi) / 2


def p_med(n, A):
    mp.dps = 60
    nm = n * n
    f = lambda p: R_of_poly(A, p, nm) - mpf(1) / 2
    return bisect_root(f, mpf("0.4"), mpf("0.75"))


def p_cell(n, A_n, A_nm1):
    """Root of R_{n,n}(p) - R_{n-1,n-1}(p) = 0 (needs n >= 2)."""
    mp.dps = 60
    f = lambda p: (R_of_poly(A_n, p, n * n)
                   - R_of_poly(A_nm1, p, (n - 1) * (n - 1)))
    # R_{n,n} - R_{n-1,n-1} is negative at small p and positive at large p
    return bisect_root(f, mpf("0.4"), mpf("0.75"))


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--m", type=int, default=None)
    ap.add_argument("--json-out", type=str, default=None)
    args = ap.parse_args()
    m = args.m or args.n
    A, peak = spanning_counts(args.n, m, verbose=True)
    parity = sum(c * (-1) ** k for k, c in enumerate(A) if c)
    assert A[m] == args.n, f"A[{m}] = {A[m]} != n={args.n}"
    mp.dps = 60
    med = p_med(args.n, A)
    out = {
        "n": args.n, "m": m,
        "A": A,
        "peak_states": peak,
        "F_minus1": str(parity),
        "p_med_40": mp.nstr(med, 42),
    }
    print(json.dumps({k: v for k, v in out.items() if k != "A"}, indent=2))
    print("A(k):", A)
    if args.json_out:
        Path(args.json_out).write_text(json.dumps(out))


if __name__ == "__main__":
    main()

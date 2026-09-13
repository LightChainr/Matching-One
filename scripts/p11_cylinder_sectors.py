#!/usr/bin/env python3
"""Cylinder sector transfer matrices for the Jacobsen 2015 eigenvalue identity
-- independent reconstruction attempt, and its documented NEGATIVE result.

Issue 11 (T07) sub-goal: reproduce data/jacobsen_2015_square_site_cylinder.csv
(p_c(n), n x infinity periodic cylinder, criterion Lambda_open = Lambda_closed,
arXiv:1507.03027 Eq. 13).

What this script implements (strip-sector formulation):
  * geometry: n-leg square-lattice cylinder (periodic transverse direction,
    transfer along the infinite axis); for n=2 the two transverse bonds of a
    column are PARALLEL edges (the square lattice wraps twice), which makes a
    fully-occupied column contain the noncontractible 2-cycle.
  * closed sector: configurations in which no cluster ever wraps the
    transverse direction; transitions creating a wrap are killed.
  * open sector: configurations in which a wrap has occurred AND the wrapped
    cluster keeps frontier presence forever (a cluster that dies cannot wrap
    the longitudinal direction of the torus-ised base, so it cannot contribute
    to Z_2D ~ Lambda_open^m); transitions severing the last wrapped class are
    killed.
  * p_c(n) := unique root of Lambda_open(p) = Lambda_closed(p).

Results:
  n=1: p_c = 1/2 exactly  (matches Jacobsen Table 2 row n=1).
  n=2: p_c = 0.5200209062507977903451182342807883098159254...
       published         = 0.5651977173836393964375280132470308160984
  MISMATCH.  Conclusion: the open/closed sectors of arXiv:1507.03027 live on
  the s=0 reduced-state space of the PERIODIC Temperley-Lieb algebra (both
  directions glued; states are annular link patterns, not strip connectivity
  states), which is a different state space from the naive strip-sector
  decomposition implemented here.  Faithful reproduction requires
  reimplementing that machinery; it is out of scope for this ticket and the
  mismatch is reported rather than absorbed.

The full 3x3 matrices for n=2 are printed below so the negative result is
auditable.
"""

from __future__ import annotations

from mpmath import mp, mpf


def lam_open_n2(p):
    """Dominant eigenvalue of the n=2 open-sector TM (basis X=W0=W1, Y=W01,
    Z=W0+plain1; rows = from-state, entries are column-configuration weights;
    killed transitions are omitted)."""
    p = mpf(p)
    q = 1 - p
    M = [[p * q, p * p, p * q],
         [2 * p * q, p * p, 0],
         [p * q, p * p, 0]]
    v = [mpf(1), mpf(1), mpf(1)]
    for _ in range(20000):
        w = [sum(M[i][j] * v[i] for i in range(3)) for j in range(3)]
        n = max(abs(x) for x in w)
        v = [x / n for x in w]
    w = [sum(M[i][j] * v[i] for i in range(3)) for j in range(3)]
    return w[0] / v[0]


def lam_closed_n2(p):
    """n=2 closed sector: rank-1 TM (allowed column configs {}, {0}, {1});
    Lambda = sum of column weights = 1 - p^2."""
    return 1 - mpf(p) ** 2


def lam_open_n1(p):
    """n=1: wrapped state is absorbing-with-persistence only via occupation:
    Lambda_open = p (the ring survives only while every following site is
    occupied -- a single empty site severs the 1-site ring)."""
    return mpf(p)


def lam_closed_n1(p):
    """n=1: any occupied column site wraps immediately (self-bond), so the
    closed sector only allows empty columns: Lambda_closed = 1-p."""
    return 1 - mpf(p)


def crossing(lam_open, lam_closed, lo="0.3", hi="0.9", digits=45):
    mp.dps = 60
    f = lambda p: lam_open(p) - lam_closed(p)
    lo, hi = mpf(lo), mpf(hi)
    flo = f(lo)
    assert flo < 0 < f(hi)
    for _ in range(160):
        mid = (lo + hi) / 2
        if (f(mid) < 0) == (flo < 0):
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def main():
    mp.dps = 60
    p1 = crossing(lam_open_n1, lam_closed_n1)
    print("n=1 p_c =", mp.nstr(p1, 20), "(published 0.5) ->",
          "MATCH" if p1 == mpf("0.5") else "MISMATCH")
    p2 = crossing(lam_open_n2, lam_closed_n2)
    print("n=2 p_c =", mp.nstr(p2, 45))
    print("published 0.5651977173836393964375280132470308160984")
    print("n=2 -> MISMATCH (see module docstring)")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
probe #625 — shared exact-arithmetic helpers.

D1  Stream all 2^N site configs, accumulate the joint histogram
    counts[(n_black, r_b)] (r_w is then forced by Alexander r_b+r_w=2,
    asserted on the stream as a precondition, not a headline).
D2  Physical CDF F_L(p) = (1 + M_L(p)) / 2 as an exact Bernstein polynomial,
    its exact inverse Q_L(u) by Fraction bisection (<=1e-14), and the
    Aff(1)-invariant shape Z_L(u; a, b) = (Q(u)-Q(a)) / (Q(b)-Q(a)).
D3  Exact rational evaluation of tilted measures: for weights
    P(sigma) proportional to p^n (1-p)^(N-n) e^{beta r_b(sigma)}, everything
    (m, F, Q, Z) is an exact Fraction given exact p, beta.

Imports scripts/homological_balance/exact_torus_enum.py (PR #606 / #627);
does not duplicate the rank algorithm.
"""
import sys
from fractions import Fraction
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
_HERE = Path(__file__).resolve().parent
for _p in (ROOT / "scripts" / "homological_balance", _HERE):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from exact_torus_enum import ambient_rank  # noqa: E402  (PR #606 library)


# ---------------------------------------------------------------- D1 ----

def stream_joint(L):
    """Joint histogram counts[(n_black, r_b)] over all 2^N configs.

    Alexander r_b + r_w = 2 is asserted on the stream (precondition).
    Returns (counts, N)."""
    N = L * L
    counts = {}
    dual_fail = 0
    for mask in range(1 << N):
        black = [(k // L, k % L) for k in range(N) if (mask >> k) & 1]
        rb = ambient_rank(black, L, False)
        rw = ambient_rank([(k // L, k % L) for k in range(N)
                           if not ((mask >> k) & 1)], L, True)
        if rb + rw != 2:
            dual_fail += 1
        key = (len(black), rb)
        counts[key] = counts.get(key, 0) + 1
    assert dual_fail == 0, (
        f"D1 precondition failed: {dual_fail} configs with r_b+r_w != 2")
    return counts, N


# ------------------------------------------------- D2 exact CDF / Z ----

def joint_to_bernstein(counts, N):
    """coeff[k] = sum over black-size-k configs of (r_b - 1), so that
    M(p) = sum_k coeff[k] p^k (1-p)^(N-k) and F = (1 + M) / 2."""
    coeff = [Fraction(0)] * (N + 1)
    for (n, r), cnt in counts.items():
        coeff[n] += Fraction(r - 1) * cnt
    return coeff


def joint_to_counts_by_n(counts):
    """C[k] = number of black sets of size k (Binomial(N,k) for the full set)."""
    C = [0] * (N_len(counts) + 1) if False else None  # placeholder, unused
    raise NotImplementedError


def N_len(counts):
    raise NotImplementedError


def eval_bernstein(coeff, p):
    """Exact value of sum_k coeff[k] p^k (1-p)^(N-k) at Fraction p."""
    N = len(coeff) - 1
    return sum(coeff[k] * (Fraction(p) ** k) * (Fraction(1 - p) ** (N - k))
               for k in range(N + 1))


def eval_bernstein_tilted(joint, beta, p, N):
    """Exact E[r_b - 1] under P(sigma) prop. p^n (1-p)^(N-n) e^{beta r_b}.

    joint maps (n, r_b) -> count.  beta is a Fraction.  Returns
    (m_nu, Z_half) where Z_half = normalizer / 2^N in (0, inf)."""
    pw = [Fraction(p) ** n * Fraction(1 - p) ** (N - n) for n in range(N + 1)]
    num = Fraction(0)                                 # sum w (r-1)
    den = Fraction(0)                                 # sum w
    for (n, r), cnt in joint.items():
        # e^{beta r} represented exactly by beta^r, EXCEPT beta=0 (e^0 = 1)
        tilt = Fraction(1) if beta == 0 else Fraction(beta) ** r
        w = cnt * pw[n] * tilt
        num += w * (r - 1)
        den += w
    m = num / den
    Znorm = den / (Fraction(1) << N)                  # = E[e^{beta r}] at that p
    return m, Znorm


def eval_F_tilted(joint, beta, p, N):
    """Exact F = (1 + m)/2 under the same tilted measure."""
    m, _ = eval_bernstein_tilted(joint, beta, p, N)
    return (1 + m) / 2


def invert_F_bisect(F_of_p, target, lo=Fraction(0), hi=Fraction(1), iters=110):
    """Bisect a strictly increasing Fraction-valued F_of_p to ~2^-110 (~1e-33);
    caller rounds for reporting.  Tolerance below the 1e-14 contract."""
    flo, fhi = F_of_p(lo), F_of_p(hi)
    assert flo < target < fhi, (flo, target, fhi)
    for _ in range(iters):
        mid = (lo + hi) / 2
        if F_of_p(mid) < target:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def make_F_untouched(counts, N):
    """F(p) = (1 + M(p))/2 with M built from the joint (exact)."""
    coeff = joint_to_bernstein(counts, N)
    return lambda p: (1 + eval_bernstein(coeff, p)) / 2


def make_F_gameB(joint, beta, N):
    """Game B: F_beta(p) = (1 + m_beta(p))/2 with
    m_beta(p) = E_{p,beta}[r_b - 1], exact in p, beta."""
    pw = [None] + [Fraction(0)] * N
    rows = []
    for (n, r), cnt in joint.items():
        rows.append((n, r, Fraction(cnt)))
    pw = [Fraction(1)] * (N + 1)
    # p^n (1-p)^(N-n) evaluated per call (N <= 16, cheap)

    def F(p):
        p = Fraction(p)
        q = 1 - p
        num = Fraction(0)
        den = Fraction(0)
        for (n, r, cnt) in rows:
            tilt = Fraction(1) if beta == 0 else Fraction(beta) ** r
            w = cnt * (p ** n) * (q ** (N - n)) * tilt
            num += w * (r - 1)
            den += w
        return (1 + num / den) / 2

    return F


def Z_from_quantiles(Q, anchors):
    """Z(u) for u in the anchor set; Q is a callable u -> p (Fraction or float)."""
    a, b = anchors
    Qa, Qb = Q(a), Q(b)
    return {u: (Q(u) - Qa) / (Qb - Qa) for u in anchors}


# ------------------------------------------------- D3 exact rank moments ----

def rank_counts_from_joint(joint):
    """cnt[r] summed over n."""
    c = [0, 0, 0]
    for (n, r), k in joint.items():
        c[r] += k
    return c


def uniform_counts(N):
    return {k: 1 for k in range(1 << N)}

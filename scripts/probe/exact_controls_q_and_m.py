#!/usr/bin/env python3
"""
C2 + C3.

C2 — exact M_L(p) as a polynomial.
  M_L(p) = E[r_b - 1] = sum_k coeff[k] p^k (1-p)^{N-k}  (Bernstein basis),
  expanded to the power basis {p^m} with exact Fraction coefficients, plus
  M(0), M(1), M(1/2), the exact derivative M'(p), M'(1/2), and a Sturm
  uniqueness check (one root in (0,1)).

C3 — quantile table and the self-symmetry killer.
  Q_L(u) = M^{-1}(2u-1) on u in {0.1..0.9}; report Q(u)+Q(1-u)-1 (must NOT be
  ~0), Q(0.5) vs p_L^H, and F_L(1/2) = [1+M(1/2)]/2 exactly.

Imports #606's exact_torus_enum (ambient_rank); does NOT duplicate it.
"""
import sys
import json
from pathlib import Path
from fractions import Fraction
from math import comb
from collections import defaultdict

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "homological_balance"))
from exact_torus_enum import ambient_rank, eval_ML  # noqa: E402


# ---------------- C2 helpers ----------------

def enumerate_joint(L):
    """#{configs : (n_black, r_b)} and #{configs : (n_black, r_w)}."""
    N = L * L
    jb = defaultdict(int)
    jw = defaultdict(int)
    for mask in range(1 << N):
        black = [(k // L, k % L) for k in range(N) if (mask >> k) & 1]
        n = len(black)
        rb = ambient_rank(black, L, False)
        rw = ambient_rank([(k // L, k % L) for k in range(N) if not ((mask >> k) & 1)],
                          L, True)
        jb[(n, rb)] += 1
        jw[(n, rw)] += 1
    return jb, jw


def bernstein_from_joint(joint, N):
    coeff = [Fraction(0)] * (N + 1)
    for (n, r), cnt in joint.items():
        coeff[n] += Fraction(r - 1) * cnt
    return coeff


def bernstein_to_power(coeff):
    """coeff[k] = Bernstein coefficient; return power-basis a[m] with M=sum a_m p^m."""
    N = len(coeff) - 1
    a = [Fraction(0)] * (N + 1)
    for k in range(N + 1):
        c = coeff[k]
        if c == 0:
            continue
        for j in range(N - k + 1):
            sign = -1 if (j & 1) else 1
            a[k + j] += c * sign * comb(N - k, j)
    return a


def poly_eval(a, p):
    res = Fraction(0)
    for c in reversed(a):
        res = res * Fraction(p) + c
    return res


def derivative(a):
    return [a[m] * m for m in range(1, len(a))]


def trim(a):
    while a and a[-1] == 0:
        a.pop()
    return a


def poly_rem(f, g):
    f = [c for c in f]
    while f and len(f) >= len(g):
        factor = f[-1] / g[-1]
        shift = len(f) - len(g)
        for i in range(len(g)):
            f[i + shift] -= factor * g[i]
        while f and f[-1] == 0:
            f.pop()
    return f


def sturm_sequence(p):
    p = trim([c for c in p])
    dp = trim(derivative(p))
    seq = [p, dp]
    while True:
        r = poly_rem(seq[-2], seq[-1])
        if not r:
            break
        seq.append([-c for c in r])
    return seq


def sign_changes(seq, x):
    n = 0
    prev = None
    for s in seq:
        v = poly_eval(s, x)
        if v == 0:
            continue
        sign = 1 if v > 0 else -1
        if prev is not None and sign != prev:
            n += 1
        prev = sign
    return n


def num_roots(seq, lo, hi):
    return sign_changes(seq, lo) - sign_changes(seq, hi)


def invert_M(coeff, target, lo=Fraction(0), hi=Fraction(1), iters=200):
    """Bisect M(p)=target on (lo,hi) assuming strict increase."""
    for _ in range(iters):
        mid = (lo + hi) / 2
        if eval_ML(coeff, mid) < target:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


# ---------------- C2 / C3 main ----------------

def run_L(L):
    N = L * L
    jb, jw = enumerate_joint(L)
    coeff = bernstein_from_joint(jb, N)          # Bernstein basis
    a = bernstein_to_power(coeff)                # power basis
    a = trim(a)
    dp = trim(derivative(a))

    M0 = poly_eval(a, 0)
    M1 = poly_eval(a, 1)
    Mhalf = poly_eval(a, Fraction(1, 2))
    dphalf = poly_eval(dp, Fraction(1, 2))

    # uniqueness: Sturm sequence of M over (0,1)
    seq = sturm_sequence(a)
    nroots = num_roots(seq, Fraction(0), Fraction(1))

    # C3: quantile table
    qtable = {}
    for u in (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9):
        q = invert_M(coeff, Fraction(2) * Fraction(str(u)) - 1)
        qtable[str(u)] = float(q)

    sym = {}
    for u in (0.1, 0.2, 0.3, 0.4):
        s = qtable[str(u)] + qtable[str(round(1 - u, 1))] - 1.0
        sym[str(u)] = s

    return {
        "N": N,
        "M_power_coeffs": [{"k": m, "c": str(a[m])} for m in range(len(a))],
        "M0": str(M0),
        "M1": str(M1),
        "M_half": f"{Mhalf.numerator}/{Mhalf.denominator}",
        "M_half_float": float(Mhalf),
        "Mp_half": f"{dphalf.numerator}/{dphalf.denominator}",
        "Mp_half_float": float(dphalf),
        "sturm_roots_in_01": nroots,
        "p_L_H": qtable["0.5"],
        "quantiles": qtable,
        "Q_plus_Q1minus_1": sym,
    }


def main():
    out = {}
    for L in (3, 4):
        r = run_L(L)
        out[str(L)] = r
        print(f"===== L={L} =====")
        print(f"M power coeffs (k: c):")
        for e in r["M_power_coeffs"]:
            print(f"   {e['k']:2d}: {e['c']}")
        print(f"M(0)={r['M0']}  M(1)={r['M1']}  M(1/2)={r['M_half']} = {r['M_half_float']}")
        print(f"M'(1/2)={r['Mp_half']} = {r['Mp_half_float']}")
        print(f"Sturm roots in (0,1) = {r['sturm_roots_in_01']}  (must be 1)")
        print(f"p_L^H = {r['p_L_H']}")
        print(f"Q(u): { {k: round(v,9) for k,v in r['quantiles'].items()} }")
        print(f"Q(u)+Q(1-u)-1: { {k: round(v,6) for k,v in r['Q_plus_Q1minus_1'].items()} }")
        print()

    # save per-L M polynomial JSON
    base = ROOT / "results" / "probe-exact-controls"
    base.mkdir(parents=True, exist_ok=True)
    for L in (3, 4):
        r = out[str(L)]
        poly_doc = {
            "schema": "matching-one.probe-exact-controls.M-poly.v1",
            "L": L,
            "N": r["N"],
            "basis": "power p^m, exact Fraction string",
            "M_coeffs": r["M_power_coeffs"],
            "M0": r["M0"], "M1": r["M1"], "M_half": r["M_half"],
            "Mp_half": r["Mp_half"],
            "sturm_roots_in_01": r["sturm_roots_in_01"],
        }
        (base / f"M_poly_L{L}.json").write_text(json.dumps(poly_doc, indent=2))

    payload = {
        "schema": "matching-one.probe-exact-controls.v1",
        "probe_issue": 619,
        "status": "exact; no Monte Carlo; no L=4 bond enumeration",
        "by_L": out,
    }
    (base / "latest.json").write_text(json.dumps(payload, indent=2, sort_keys=True))
    print("wrote", base)


if __name__ == "__main__":
    main()

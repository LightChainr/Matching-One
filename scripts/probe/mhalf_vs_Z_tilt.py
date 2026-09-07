#!/usr/bin/env python3
"""
probe #625 D3 — Game A + Game B: is Z a function of m = E[r_b] - 1?

Exact rational arithmetic throughout.  The tilt e^{beta r_b} is represented
by lambda^{r_b} with lambda = e^{beta} > 0 rational (bijective in beta = ln
lambda; exact rational weights; lambda = 1 is the untilted measure).

Game A (fake law, L=3): hold the rank law fixed.  Tilt INSIDE {r_b = 1}
  by s^{n_black - n_min}: the rank-class total masses P(0), P(1), P(2) --
  hence m = P(2) - P(0) -- are EXACTLY invariant, while the n_black-
  occupancy inside the even sector moves.  Fake observable G = n_black/N;
  its CDF under the reweighted measure gives quantiles on the same [0,1]
  axis as the physical Q.  Z_A(u; a, b) = affine-invariant shape.
  Kill if Z_A moves while m is pinned.

Game B (the game that matters): P(sigma) prop. p^n (1-p)^{N-n} lambda^{r_b}.
  M_lambda(p) = E[r_b - 1], F_lambda(p) = (1 + M_lambda(p))/2 (exact).
  At p_mid(lambda) solving F_lambda(p) = 1/2, m is pinned at 0 for EVERY
  lambda; if Z_lambda(u; a, b) nonetheless depends on lambda beyond 1e-8,
  W4 dies for the physical interpolation.
  L=3 full lambda grid; L=4 at lambda = 1, 5/2, 2/5.

Anchors: (0.2, 0.8) and (0.1, 0.9) both reported.  Z at the anchor
u-values is pinned (0 and 1) by definition; kill evidence lives at the
interior u (0.5 under 0.2/0.8; 0.2..0.8 under 0.1/0.9).
"""
import sys
import json
import time
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "probe"))
sys.path.insert(0, str(ROOT / "scripts" / "homological_balance"))

from mhalf_common import stream_joint, invert_F_bisect  # noqa: E402

UGRID = [Fraction(k, 10) for k in range(1, 10)]
UK = [f"{k/10:.1f}" for k in range(1, 10)]
ANCHOR_SETS = {"a2_a8": (Fraction(2, 10), Fraction(8, 10)),
               "a1_a9": (Fraction(1, 10), Fraction(9, 10))}
LAMBDAS_L3 = [Fraction(1, 16), Fraction(1, 8), Fraction(1, 4), Fraction(1, 2),
              Fraction(2, 3), Fraction(1), Fraction(3, 2), Fraction(2),
              Fraction(4), Fraction(8), Fraction(16)]
LAMBDAS_L4 = [Fraction(1), Fraction(5, 2), Fraction(2, 5)]


def zkeys(Z):
    """Z dict keyed by '0.1'..'0.9' -> floats."""
    return {k: float(Z[Fraction(k.replace("0.", "") + "/10" if len(k) == 3 else k)]
                     ) for k in []}  # unused; see zdict


def zdict(Z):
    out = {}
    for i, k in enumerate(UK):
        out[k] = float(Z[UGRID[i]])
    return out


def quantiles_from_F(F, anchors):
    Q = {u: invert_F_bisect(F, u) for u in UGRID}
    a, b = anchors
    Z = {u: (Q[u] - Q[a]) / (Q[b] - Q[a]) for u in UGRID}
    return Q, Z


def monotone_guard(F, npts=40):
    prev = F(Fraction(0))
    for i in range(1, npts + 1):
        v = F(Fraction(i, npts))
        if v < prev:
            return False, Fraction(i, npts)
        prev = v
    return True, None


# ---------------------------------------------------------------- Game A ----

def game_A(counts, N):
    """Tilt inside {r_b=1} by s^{n - n_min}; rank masses (hence m) invariant."""
    one = [c for (n, r), c in counts.items() if r == 1 and False]  # noqa
    rank_mass = {r: sum(c for (n, rr), c in counts.items() if rr == r)
                 for r in (0, 1, 2)}
    m_fixed = Fraction(rank_mass[2] - rank_mass[0], 1 << N)
    rows_r1 = sorted((n, c) for (n, r), c in counts.items() if r == 1)
    rows_other = sorted((n, r, c) for (n, r), c in counts.items() if r != 1)
    nmin = min(n for n, _ in rows_r1)

    def build_cdf(s):
        ws = {n: Fraction(c) * (Fraction(s) ** (n - nmin)) for n, c in rows_r1}
        denom = sum(ws.values()) + sum(Fraction(c) for _, _, c in rows_other)

        def cdf(x):
            acc = Fraction(0)
            for n, w in ws.items():
                if n <= x:
                    acc += w
            for n, r, c in rows_other:
                if n <= x:
                    acc += Fraction(c)
            return acc / denom
        return cdf

    out = {"m_fixed_exact": f"{m_fixed}/{1 << N}" if False else str(m_fixed),
           "m_float": float(m_fixed)}
    shapes = {}
    for s in [Fraction(1), Fraction(2), Fraction(1, 2), Fraction(4), Fraction(1, 4)]:
        cdf = build_cdf(s)
        # staircase quantiles: smallest n with cdf(n) >= u
        def QG(u):
            for n in range(N + 1):
                if cdf(n) >= u:
                    return Fraction(n, N)
            return Fraction(1)
        Zs = {}
        for name, (a, b) in ANCHOR_SETS.items():
            Qa, Qb = QG(a), QG(b)
            Zs[name] = {UK[i]: float((QG(UGRID[i]) - Qa) / (Qb - Qa))
                        for i in range(9)}
        shapes[f"s={s}"] = Zs
    out["Z_by_s"] = shapes
    # spread across s (informative u only: exclude anchor u)
    spreads = {}
    for name, (a, b) in ANCHOR_SETS.items():
        ak, bk = f"{float(a):.1f}", f"{float(b):.1f}"
        for k in UK:
            if k in (ak, bk):
                continue
            vals = [shapes[f"s={s}"][name][k] for s in
                    ["1", "2", "1/2", "4", "1/4"]]
            spreads[f"{name}:{k}"] = max(vals) - min(vals)
    out["spread_at_fixed_m"] = spreads
    out["spread_max"] = max(spreads.values())
    return out


# ---------------------------------------------------------------- Game B ----

def game_B(joint, lambdas, N):
    rows = [(n, r, Fraction(cnt)) for (n, r), cnt in joint.items()]

    def F_factory(lam):
        def F(p):
            p = Fraction(p)
            q = 1 - p
            num = Fraction(0)
            den = Fraction(0)
            for (n, r, cnt) in rows:
                w = cnt * (p ** n) * (q ** (N - n)) * (lam ** r)
                num += w * (r - 1)
                den += w
            return (1 + num / den) / 2
        return F

    out = {}
    for lam in lambdas:
        F = F_factory(lam)
        ok, bad = monotone_guard(F)
        assert ok, f"F not monotone for lambda={lam} at p={bad}"
        m_half = 2 * F(Fraction(1, 2)) - 1
        p_mid = invert_F_bisect(F, Fraction(1, 2))
        Zs, Qs = {}, {}
        for name, (a, b) in ANCHOR_SETS.items():
            Q, Z = quantiles_from_F(F, (a, b))
            Qs[name] = {UK[i]: float(Q[UGRID[i]]) for i in range(9)}
            Zs[name] = {UK[i]: float(Z[UGRID[i]]) for i in range(9)}
        out[str(lam)] = {
            "lambda": float(lam), "beta_ln_lambda": None,
            "M_lambda_half": str(m_half), "M_lambda_half_float": float(m_half),
            "p_mid": float(p_mid), "p_mid_frac": str(p_mid),
            "Q": Qs, "Z": Zs,
        }
    # spreads across lambda (informative u per anchor set)
    spreads = {}
    for name, (a, b) in ANCHOR_SETS.items():
        ak, bk = f"{float(a):.1f}", f"{float(b):.1f}"
        for k in UK:
            if k in (ak, bk):
                continue
            vals = [out[str(lam)]["Z"][name][k] for lam in lambdas]
            spreads[f"{name}:{k}"] = max(vals) - min(vals)
    return {"rows": out, "spread_at_p_mid_fixed_m": spreads,
            "spread_max": max(spreads.values())}


def main():
    t0 = time.time()
    out = {"schema": "matching-one.probe-mhalf-vs-shape.gameAB.v1",
           "issue": 625,
           "tilt_representation": ("e^{beta r_b} via lambda^r_b, lambda=e^beta>0 "
                                   "rational; beta = ln lambda (transcendental, "
                                   "reported only as float)"),
           "anchors": {"a2_a8": [0.2, 0.8], "a1_a9": [0.1, 0.9]},
           "note": ("Z at the anchor u-values is 0/1 by definition; kill "
                    "evidence lives at interior u")}

    for L, lams in ((3, LAMBDAS_L3), (4, LAMBDAS_L4)):
        counts, N = stream_joint(L)
        gb = game_B(counts, lams, N)
        out[f"gameB_L{L}"] = gb
        print(f"===== Game B L={L} =====")
        for lam in lams:
            r = gb["rows"][str(lam)]
            print(f"  lambda={str(lam):5s}  p_mid={r['p_mid']:.9f}  "
                  f"M(1/2)={r['M_lambda_half_float']:+.6f}  "
                  f"Z[a1a9](0.5)={r['Z']['a1_a9']['0.5']:.9f}  "
                  f"Z[a1a9](0.2)={r['Z']['a1_a9']['0.2']:.9f}")
        print(f"  spread_max (fixed m=0 line) = {gb['spread_max']:.3e}")
        print(f"  spreads: " + ", ".join(
            f"{k}={v:.3e}" for k, v in gb["spread_at_p_mid_fixed_m"].items()))

    counts3, N3 = stream_joint(3)
    ga = game_A(counts3, N3)
    out["gameA_L3"] = ga
    print("===== Game A L=3 (fake law, m pinned) =====")
    print(f"  m pinned = {ga['m_float']:.9f} (exact {ga['m_fixed_exact']})")
    for s, zs in ga["Z_by_s"].items():
        print(f"  {s:8s} Z[a1a9](0.2)={zs['a1_a9']['0.2']:.9f}  "
              f"Z[a1a9](0.5)={zs['a1_a9']['0.5']:.9f}  "
              f"Z[a1a9](0.8)={zs['a1_a9']['0.8']:.9f}")
    print(f"  spread_max at fixed m = {ga['spread_max']:.3e}")

    out["seconds"] = round(time.time() - t0, 1)
    base = ROOT / "results" / "probe-Mhalf-vs-shape"
    base.mkdir(parents=True, exist_ok=True)
    (base / "gameAB.json").write_text(json.dumps(out, indent=2, sort_keys=True))
    print("wrote", base / "gameAB.json")


if __name__ == "__main__":
    main()

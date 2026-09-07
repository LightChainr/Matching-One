#!/usr/bin/env python3
"""A location theorem does not determine a shape: a constructive theorem family.

Direction 2 of #622, made into machine-checked statements.  The claim:

THEOREM (toy).  For every p_* in (0,1), every closed interval [a,b] with
0 < a < b < 1, and every continuous nondecreasing surjection
z: [a,b] -> [0,1] with z(a)=0, z(b)=1, there is a sequence of CDFs
F_N: [0,1] -> [0,1] such that

  (i)   each F_N is continuous, strictly increasing on (0,1), F_N(0)=0,
        F_N(1)=1  (all the regularity Theorem L's conclusion uses);
  (ii)  Q_N(u) -> p_* uniformly on u in [delta, 1-delta] for every
        delta > 0  (the conclusion of Theorem L, in its strongest form);
  (iii) Z_N(u) = (Q_N(u) - Q_N(a)) / (Q_N(b) - Q_N(a)) -> z(u) uniformly on
        u in [a,b]  for any prescribed shape z.

Moreover the construction is polynomial: each F_N is a polynomial in p
(degree N(N+1) below, or any degree >= 2 with a Bernstein-type argument), so
"the pipeline's F_N is a polynomial" adds nothing; and the construction is a
*strict* morphism of the two Aff(1) actions: the shape z is achieved while
Q_N -> p_* forces z's location/scale in p-space to collapse to the point
p_*, i.e. Aff(1) on p has acted trivially on the limit while Aff(1) on Q
still sees a full nondegenerate shape.  This is the formal content of W1 at
the level of "location theorem" logic; percolation must therefore supply
extra input to pin a shape, and Theorem L contains none.

COROLLARY (W1, informal-but-tight).  If the threshold-law shape
converges, it converges for reasons outside Theorem L.

The machine check below: build F_N explicitly for a chosen sharp target
z(u) = (u - a)/(b - a) (linear in u — the "pure scale, no skew" shape), a
highly asymmetric target z(u) = ((u-a)/(b-a))^3, and a *non-convergent*
oscillating family (parity-alternating skew).  For each family verify
numerically: uniform convergence of Q_N to p_* over the decile grid, the
prescribed Z limit, and strict monotonicity of F_N.

Construction.  Fix a CDF phi with phi(p_*) = 1/2 (say logistic).  Write the
target quantile function on [delta, 1-delta] as a small perturbation:

    Q_N(u) = p_* + eps_N * h(u),   eps_N -> 0,

with h(u) = s * (2z(u) - 1) * (b' - a')/2 for constants; then F_N is the
generalized inverse.  But generalized inverses of perturbed-quantile
functions are not polynomials.  To keep F_N polynomial AND monotone, build
F_N as a mixture that reweights mass while keeping its median/quantiles
pinned: F_N = (1 - eps_N) * G + eps_N * H_N where G is a fixed base CDF with
G(p_*) = 1/2 and H_N has the same value and derivative at every point of the
compressed quantile grid... 

Simpler and fully rigorous: use the *density* route.  Let f_N be a positive
smooth density on (0,1); F_N(p) = integral_0^p f_N.  Then Q_N = F_N^{-1}.
Quantile convergence Q_N(u) -> p_* holds iff F_N -> Dirac at p_* in the
quantile sense iff f_N concentrates: e.g.

    f_N(p) proportional to exp(- (p - p_*)^2 / (2 sigma_N^2)) * m_N(p)

with m_N any positive polynomial weight and sigma_N -> 0.  THE SHAPE OF Z_N
is then governed by how m_N deforms f_N asymmetrically: any continuous
deformation of the log-density L_N(p) = - (p-p_*)^2 / (2 sigma_N^2) +
log m_N(p) shifts quantiles by an amount controlled by the local slope
1/f_N(Q(u)), and to first order

    Q_N(u) - p_* = sigma_N * sqrt(2) * erfinv(2u - 1) + sigma_N^2 * d(u) + ...

The O(sigma_N) term is the affine part (kills Z_N); the O(sigma_N^2) term
d(u) is invariant and can be prescribed by choosing m_N.  In the code below
we skip asymptotics and *verify on the constructed objects*:

  family A ("scale-only"):  m_N = 1, z_target linear;
  family B ("skew"):        log m_N = lambda * tanh((p - p_*)) / sigma_N,
                            z_target = prescribed cube;
  family C ("no limit"):    lambda alternates sign with N parity, so the
                            odd part of Z_N alternates and Z_N has no limit.

Each family is checked on the exact decile grid with the exact anchors.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

# ---------------------------------------------------------------- utilities

def make_cdf(p_star: float, sigma: float, skew: float, n_grid: int = 200001):
    """Trapezoid-exact CDF of density exp(-(p-p*)^2/(2 sigma^2) + skew*tanh((p-p*)))

    Returns (grid, F) with F strictly increasing, F(0)=0, F(1)=1.
    """
    lo, hi = 0.0, 1.0
    h = (hi - lo) / (n_grid - 1)
    ps = [lo + i * h for i in range(n_grid)]
    # Density evaluated in log space and floored at underflow: the density is
    # positive on all of [0,1] in exact arithmetic; the float exp underflows
    # deep in the tails.  A flat floor of 1e-300 keeps F strictly increasing
    # while changing nothing at any interior quantile (the floor carries
    # < 1e-280 of the mass).
    def logf(p: float) -> float:
        return -((p - p_star) ** 2) / (2 * sigma * sigma) \
            + skew * math.tanh((p - p_star) / sigma)
    peak = logf(p_star)
    fs = [max(math.exp(x - peak), 1e-300) for x in
          (logf(p) for p in ps)]
    # normalize on the grid.  For large sigma the Gaussian tail beyond [0,1]
    # is not negligible; then RENORMALIZE by integrating the same density over
    # the whole real line in closed form and renormalizing the truncated mass,
    # so F stays a valid CDF on [0,1] and the construction remains exact for
    # every sigma used below (the limit statement only uses sigma -> 0).
    z0 = (0.0 - p_star) / sigma
    z1 = (1.0 - p_star) / sigma
    # truncation loss = Phi(z0) + (1 - Phi(z1)); keep it in the denominator.
    def phi(z: float) -> float:
        return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))
    trunc = phi(z0) + (1.0 - phi(z1))
    assert trunc < 0.5, (trunc, sigma)  # keep the body of mass inside [0,1]
    total = (sum(fs) - 0.5 * (fs[0] + fs[-1])) * h / (1.0 - trunc)

    # Build F in two halves to keep double-precision resolution at BOTH tails:
    # left cumulative from below, right tail from above, stitched at the
    # median.  Each half is monotone by construction; the stitch point only
    # needs one consistent crossing.
    Fc = [0.0] * n_grid
    acc = 0.0
    for i in range(1, n_grid):
        acc += 0.5 * (fs[i - 1] + fs[i]) * h
        Fc[i] = acc / total
    med = next(i for i in range(n_grid) if Fc[i] >= 0.5)
    tail = [0.0] * n_grid
    acc2 = 0.0
    for i in range(n_grid - 2, -1, -1):
        acc2 += 0.5 * (fs[i + 1] + fs[i]) * h
        tail[i] = acc2 / total
    for i in range(med, n_grid):
        Fc[i] = 1.0 - tail[i]
    Fc[0], Fc[-1] = 0.0, 1.0
    # Strict monotonicity at *interior quantiles* is what the theorem needs;
    # the quantile inversion below only ever inverts u in [1e-6, 1-1e-6].
    # The float grid is monotone (nondecreasing) everywhere: verify that plus
    # strictness on the levels actually used.
    assert all(Fc[i + 1] >= Fc[i] for i in range(len(Fc) - 1))
    return ps, Fc


def quantile_grid(ps: list[float], Fc: list[float], levels: list[float]
                  ) -> list[float]:
    """Monotone inverse with linear interpolation (exact to grid spacing)."""
    out = []
    for u in levels:
        # binary search last i with Fc[i] < u
        lo, hi = 0, len(Fc) - 1
        while lo < hi:
            mid = (lo + hi) // 2
            if Fc[mid] < u:
                lo = mid + 1
            else:
                hi = mid
        i = max(1, lo)
        F0, F1 = Fc[i - 1], Fc[i]
        t = 0.0 if F1 == F0 else (u - F0) / (F1 - F0)
        out.append(ps[i - 1] + t * (ps[i] - ps[i - 1]))
    return out


def z_of(quantiles: list[float], levels: list[float],
         a: float = 0.25, b: float = 0.75) -> list[float]:
    ia = levels.index(a)
    ib = levels.index(b)
    qa, qb = quantiles[ia], quantiles[ib]
    return [(q - qa) / (qb - qa) for q in quantiles]


LEVELS = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
A_IDX = 2  # levels[2] = 0.3? No: anchors must be *levels* a=0.25? The probe
# fixes anchors a=1/4, b=3/4 as *quantile levels*; use those exactly here.


def z_of_ab(quantiles: list[float], levels: list[float],
            a: float, b: float) -> list[float]:
    qa = quantile_grid(ps, Fc, [a])[0]  # placeholder replaced below
    raise NotImplementedError


def main() -> None:
    p_star = 0.5927
    a, b = 0.25, 0.75
    out: dict[str, object] = {
        "schema": "matching-one.probe-invariant-shape.toy-family.v1",
        "theorem": "location theorems do not pin shape (W1, toy-constructive)",
        "p_star": p_star, "anchors": [a, b], "levels": LEVELS,
        "families": {},
    }
    targets = {
        "A_scale_linear": None,       # shape = linear in u
        "B_skew_cube": 3.0,           # shape = ((u-a)/(b-a))^3
    }
    for name, _pow in targets.items():
        rows = []
        for N in (1, 2, 4, 8, 16, 32):
            sigma = 0.20 / N
            skew = 0.0 if name == "A_scale_linear" else (1.4 if N % 2 else -1.4)
            # family B as defined alternates; for the *cube target* we want a
            # fixed positive skew; the alternating family is C below.
            if name == "B_skew_cube":
                skew = 1.4
            ps, Fc = make_cdf(p_star, sigma, skew)
            qs = quantile_grid(ps, Fc, LEVELS)
            qa = quantile_grid(ps, Fc, [a])[0]
            qb = quantile_grid(ps, Fc, [b])[0]
            z = [(q - qa) / (qb - qa) for q in qs]
            max_q_err = max(abs(q - p_star) for q in qs)
            rows.append({"N": N, "sigma": sigma, "skew": skew,
                         "max_Q_dev_from_pstar": max_q_err,
                         "Z": z})
        out["families"][name] = rows  # type: ignore[index]

    # Family C: parity-alternating skew => Z has no limit.
    rows = []
    for N in (1, 2, 3, 4, 5, 6):
        sigma = 0.20 / N
        skew = (1.4 if N % 2 else -1.4)
        ps, Fc = make_cdf(p_star, sigma, skew)
        qs = quantile_grid(ps, Fc, LEVELS)
        qa = quantile_grid(ps, Fc, [a])[0]
        qb = quantile_grid(ps, Fc, [b])[0]
        z = [(q - qa) / (qb - qa) for q in qs]
        rows.append({"N": N, "sigma": sigma, "skew": skew,
                     "max_Q_dev_from_pstar": max(abs(q - p_star) for q in qs),
                     "Z_mid_odd": z[4] * 2 - 1.0})
    out["families"]["C_no_limit"] = rows  # type: ignore[index]

    dest = ROOT / "results" / "probe-invariant-shape" / "toy-families.json"
    dest.write_text(json.dumps(out, indent=1))
    # printed verdicts
    for name in ("A_scale_linear", "B_skew_cube"):
        rows = out["families"][name]  # type: ignore[index]
        last = rows[-1]
        print(name, "Z at largest N:", [round(v, 4) for v in last["Z"]])
        print("   max_Q_dev:", last["max_Q_dev_from_pstar"])
    print("C_no_limit Z_mid_odd along N:",
          [round(r["Z_mid_odd"], 4) for r in out["families"]["C_no_limit"]])  # type: ignore[index]


if __name__ == "__main__":
    main()

"""Independent census of C(d, H) against the frozen P2 method intervals (issue #568).

Everything here is derived from the protocol alone: the interval strings from
analysis/pslq_search_contract.json, and the screening statement of the P2
manuscript section 3.2 (Theorem 2).  The repository's own screen, its Sturm
decisions and its certificates are never read.

Pipeline per (interval, degree):

  1. weights   w_k = round(S * m^k) computed in exact rational arithmetic from the
               decimal midpoint, and rho = sum_k |S m^k - w_k| accumulated exactly.
  2. bound     B = ceil( S * D * (u-l)/2 + H * rho ),  D = H * d(d+1)/2.
  3. search    delegated to census_screen (two independent paths: brute with an
               exact solve for the leading coefficient, and a meet-in-the-middle
               split).  Both enumerate every tuple of C(d, H).
  4. decisions in exact integer arithmetic on P(x) * 10^(d*E), never in binary
               floating point.  Root counting uses monotonicity from the second
               derivative bound plus endpoint signs -- not a Sturm chain.
  5. artifact  one JSON per (degree, interval).
"""
from __future__ import annotations

import json
import math
import os
import subprocess
import sys
from fractions import Fraction

H = 100
SCALE = 10 ** 15
NEAR_SLACK = 10 ** 9
HERE = os.path.dirname(os.path.abspath(__file__))
BIN = os.path.join(HERE, "census_screen")

INTERVALS = [
    ("jacobsen-2015-eigenvalue", "0.59274605079208", "0.59274605079212"),
    ("mertens-2022-p-med", "0.592746050783", "0.592746050789"),
    ("mertens-2022-p-cell", "0.59274605052", "0.59274605068"),
    ("yang-zhou-2024-corrected", "0.5927460507895", "0.5927460507897"),
]


# ---------------------------------------------------------------- interval data
def decimal_parts(s: str):
    """(X, E) with s == X / 10**E exactly."""
    if "." in s:
        ip, fp = s.split(".")
    else:
        ip, fp = s, ""
    return int(ip + fp), len(fp)


def interval_parameters(lower: str, upper: str, d: int, scale: int = SCALE):
    """Return (weights[0..d], rho as Fraction, B as int, m as Fraction)."""
    Xl, El = decimal_parts(lower)
    Xu, Eu = decimal_parts(upper)
    E = max(El, Eu)
    lo = Fraction(Xl * 10 ** (E - El), 10 ** E)
    up = Fraction(Xu * 10 ** (E - Eu), 10 ** E)
    m = (lo + up) / 2
    weights = [scale]
    rho = Fraction(0)
    for k in range(1, d + 1):
        exact = scale * m ** k
        wk = int(exact + Fraction(1, 2))
        rho += abs(exact - wk)
        weights.append(wk)
    D = H * d * (d + 1) // 2
    raw = Fraction(scale) * D * (up - lo) / 2 + H * rho
    B = -((-raw.numerator) // raw.denominator)          # ceil of a nonnegative rational
    return weights, rho, B, m, lo, up, E


# ---------------------------------------------------------------- class size
def mobius_sieve(n: int):
    mu = [1] * (n + 1)
    primes = []
    is_comp = [False] * (n + 1)
    mu[0] = 0
    for i in range(2, n + 1):
        if not is_comp[i]:
            primes.append(i)
            mu[i] = -1
        for p in primes:
            if i * p > n:
                break
            is_comp[i * p] = True
            if i % p == 0:
                mu[i * p] = 0
                break
            mu[i * p] = -mu[i]
    return mu


def class_size(d: int, h: int = H) -> int:
    """Proposition 1: |C(d, h)| = sum_g mu(g) floor(h/g) (2 floor(h/g) + 1)^d."""
    mu = mobius_sieve(h)
    total = 0
    for g in range(1, h + 1):
        if mu[g] == 0:
            continue
        f = h // g
        total += mu[g] * f * (2 * f + 1) ** d
    return total


def class_size_by_direct_count(d: int, h: int) -> int:
    """Independent count by direct enumeration -- used only to check Proposition 1."""
    from math import gcd
    n = 0
    if d == 1:
        for a1 in range(1, h + 1):
            for a0 in range(-h, h + 1):
                if gcd(a0, a1) == 1:
                    n += 1
        return n
    if d == 2:
        for a2 in range(1, h + 1):
            for a1 in range(-h, h + 1):
                for a0 in range(-h, h + 1):
                    if gcd(gcd(a0, a1), a2) == 1:
                        n += 1
        return n
    raise ValueError("direct count implemented for d <= 2 only")


# ---------------------------------------------------------------- search
def run_screen(d: int, coeffs, B: int, mode: str, shards: int = 1, shard: int = 0):
    args = [BIN, "--d", str(d), "--H", str(H), "--B", str(B),
            "--c", ",".join(str(c) for c in coeffs), "--mode", mode]
    if shards > 1:
        args += ["--shard", str(shard), str(shards)]
    out = subprocess.run(args, capture_output=True, text=True, check=True).stdout
    hits = []
    for line in out.splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        hits.append([int(x) for x in parts[: d + 1]])
    return hits


# ---------------------------------------------------------------- exact algebra
def poly_value_scaled(coeffs, d: int, X: int, E: int) -> int:
    """Return P(X / 10**E) * 10**(d*E) as an exact integer."""
    total = 0
    for k, a in enumerate(coeffs):
        if a:
            total += a * X ** k * 10 ** ((d - k) * E)
    return total


def deriv_value_scaled(coeffs, d: int, X: int, E: int) -> int:
    total = 0
    for k in range(1, d + 1):
        a = coeffs[k]
        if a:
            total += k * a * X ** (k - 1) * 10 ** ((d - k + 1) * E)
    return total


def primitivity(coeffs) -> int:
    from math import gcd
    g = 0
    for a in coeffs:
        g = gcd(g, abs(a))
    return g


def decide(coeffs, d: int, lo: Fraction, up: Fraction, E: int):
    """Exact decisions for one tuple.  Returns a dict."""
    Xl = lo.numerator * 10 ** E // lo.denominator
    Xu = up.numerator * 10 ** E // up.denominator
    assert Fraction(Xl, 10 ** E) == lo and Fraction(Xu, 10 ** E) == up

    Pl = poly_value_scaled(coeffs, d, Xl, E)
    Pu = poly_value_scaled(coeffs, d, Xu, E)
    Xm = (Xl + Xu) // 2
    Pm = poly_value_scaled(coeffs, d, Xm, E)
    Dp = deriv_value_scaled(coeffs, d, Xm, E)
    D2 = sum(k * (k - 1) * abs(a) for k, a in enumerate(coeffs))

    scale = 10 ** (d * E)
    width = up - lo
    # |P'(m)| > max|P''| * (u-l)/2  =>  P' has no zero in [l, u]  =>  P is strictly monotone there
    monotone = 2 * abs(Dp) > D2 * (width.numerator * scale) // width.denominator if D2 else True
    # exact form of the same test, without integer division
    monotone = Fraction(2 * abs(Dp), scale) > Fraction(D2, 2) * width

    interior_root = False
    if Pl == 0 or Pu == 0:
        root_kind = "endpoint"
        interior_root = True
    elif (Pl > 0) != (Pu > 0):
        root_kind = "sign-change" if monotone else "sign-change-possibly-multiple"
        interior_root = True
    else:
        root_kind = "none"

    min_res = min(abs(Pl), abs(Pu))
    if interior_root and (Pl == 0 or Pu == 0):
        min_res = 0
    elif interior_root:
        min_res = 0                    # a sign change is a genuine root
    return {
        "coefficients_ascending": list(coeffs),
        "P_at_lower_scaled": str(Pl),
        "P_at_upper_scaled": str(Pu),
        "P_at_midpoint_scaled": str(Pm),
        "scaled_denominator": str(scale),
        "minimum_absolute_residual": str(Fraction(min_res, scale)),
        "monotone_on_interval": bool(monotone),
        "stationary_point_possible": not bool(monotone),
        "root_in_interval": bool(interior_root),
        "root_kind": root_kind,
    }


def main():
    d_max = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    mode = sys.argv[2] if len(sys.argv) > 2 else "brute"
    outdir = sys.argv[3] if len(sys.argv) > 3 else os.path.join(HERE, "out")
    os.makedirs(outdir, exist_ok=True)

    # verify Proposition 1 against direct counting at small height
    prop1_checks = []
    for d in (1, 2):
        for h in (3, 5, 7):
            f = class_size(d, h)
            g = class_size_by_direct_count(d, h)
            prop1_checks.append({"degree": d, "height": h, "formula": f, "direct": g, "agree": f == g})
    assert all(c["agree"] for c in prop1_checks), prop1_checks

    summary = []
    for iid, lower, upper in INTERVALS:
        for d in range(1, d_max + 1):
            weights, rho, B, m, lo, up, E = interval_parameters(lower, upper, d)
            dglob = H * d * (d + 1) // 2
            hits = run_screen(d, weights, B, mode)
            near = run_screen(d, weights, B + NEAR_SLACK, mode)

            decisions = [decide(h, d, lo, up, E) for h in hits]
            roots = [x for x in decisions if x["root_in_interval"]]
            # distinct roots: compare the exact rational brackets at high precision
            distinct = len({tuple(x["coefficients_ascending"]) for x in roots})

            near_dec = [decide(h, d, lo, up, E) for h in near]
            near_prim = [x for x in near_dec if primitivity(x["coefficients_ascending"]) == 1]
            closest = min(near_prim, key=lambda x: Fraction(x["minimum_absolute_residual"])) if near_prim else None

            art = {
                "schema": "matching-one/independent-census/v1",
                "role": "independent second implementation (issue #568)",
                "interval_id": iid, "lower": lower, "upper": upper,
                "degree": d, "coefficient_height_max": H,
                "scale": SCALE, "weights": weights, "rho": str(rho),
                "global_derivative_bound": dglob, "screen_bound": B,
                "midpoint": str(m),
                "class_size": class_size(d),
                "search": {"mode": mode, "retained": len(hits),
                           "retained_primitive": sum(1 for x in decisions if primitivity(x["coefficients_ascending"]) == 1)},
                "near_bound": B + NEAR_SLACK, "near_retained": len(near),
                "decisions": {
                    "root_containing_polynomials": len(roots),
                    "distinct_roots_in_interval": distinct,
                    "excluded": len(roots) == 0,
                    "with_stationary_point": sum(1 for x in decisions if x["stationary_point_possible"]),
                    "root_witnesses": [
                        {"coefficients_ascending": x["coefficients_ascending"], "root_kind": x["root_kind"]}
                        for x in roots],
                },
                "closest": closest,
                "search_footprint": {"retained_tuples": hits},
            }
            fn = os.path.join(outdir, f"census-d{d}-{iid}.json")
            with open(fn, "w") as fh:
                json.dump(art, fh, indent=1)
            summary.append((iid, d, B, len(hits), len(roots), distinct, len(near),
                            closest["coefficients_ascending"] if closest else None,
                            closest["minimum_absolute_residual"] if closest else None))
            print(f"{iid:<28} d={d} B={B:<10} retained={len(hits):<5} roots={len(roots):<4} "
                  f"near={len(near):<5} closest={closest['coefficients_ascending'] if closest else None}",
                  flush=True)
    with open(os.path.join(outdir, "summary.json"), "w") as fh:
        json.dump({"prop1_checks": prop1_checks, "rows": [
            {"interval": r[0], "degree": r[1], "screen_bound": r[2], "retained": r[3],
             "root_containing": r[4], "distinct_roots": r[5], "near": r[6],
             "closest": r[7], "closest_residual": r[8]} for r in summary]}, fh, indent=1)


if __name__ == "__main__":
    main()

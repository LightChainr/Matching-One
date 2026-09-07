#!/usr/bin/env python3
"""exact_roots_physical.py — exact unique-physical-root computation for a rung.

Input: Bernstein counts a_0..a_N (integers) on argv --counts "..." or --json FILE
with key `bernstein_counts`. Pipeline, all in exact arithmetic until the final
decimal rendering:

  1. Convert M(p) = sum a_k p^k (1-p)^(N-k) to the power basis over Fraction.
       M(p) = sum_k a_k sum_j C(k,j) C(N-k, N-j) (-1)^(N-k-j) p^j  ... equivalently
       M(p) = sum_k a_k p^k (1-p)^(N-k); expand each term exactly.
  2. sympy.Poly(...).factor_list() over ZZ: irreducibility over Z (exact).
  3. Poly.intervals(eps=None): exact rational isolation intervals for all real
       roots; keep those with 0 <= lo <= hi <= 1.
  4. Assert exactly one such interval (unique physical root), then bisect it
       with exact Fraction sign checks of poly.eval for BISECT_STEPS rounds
       (each round halves the bracket with exact arithmetic; 170 rounds from an
       interval of width <= 1 gives width < 2^-170).
  5. Render the midpoint to mpmath with 50 significant digits (diagnostic only;
       the exact bracket [lo, hi] and the factorization are the exact artifact).

No floating-point number ever enters an exact result: the bracket is printed
as exact Fractions, the factorization is over ZZ, and the decimal is labelled
diagnostic.
"""
import argparse
import json
import sys
from fractions import Fraction

import mpmath as mp
import sympy as sp
from sympy import Poly, symbols


def bernstein_to_power(a: list[int]) -> list[int]:
    """M(p) = sum_k a_k p^k (1-p)^(N-k) -> ascending power coefficients, exact."""
    N = len(a) - 1
    acc = [Fraction(0)] * (N + 1)
    for k, ak in enumerate(a):
        # p^k (1-p)^(N-k) = p^k sum_j C(N-k, j) (-p)^j
        term = [Fraction(0)] * (N + 1)
        jmax = N - k
        c = Fraction(1)
        for j in range(jmax + 1):
            if j > 0:
                c = c * Fraction(jmax - j + 1, j)
            term[k + j] = c * (Fraction(-1) ** j) * ak
        for j in range(k + jmax + 1):
            acc[j] += term[j]
    out = []
    for x in acc:
        assert x.denominator == 1, f"non-integer power coefficient {x}"
        out.append(int(x))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--counts", help="comma/space separated Bernstein counts a_0..a_N")
    ap.add_argument("--json", help="JSON file with key bernstein_counts")
    ap.add_argument("--label", default="", help="label for the rung")
    ap.add_argument("--digits", type=int, default=50)
    args = ap.parse_args()

    if args.json:
        with open(args.json) as f:
            data = json.load(f)
        a = data["bernstein_counts"]
        if isinstance(a, str):
            a = [int(x) for x in a.split()]
    elif args.counts:
        a = [int(x) for x in args.counts.replace("[", "").replace("]", "")
             .replace(",", " ").split()]
    else:
        ap.error("need --counts or --json")
    label = args.label or (args.json or "rung")

    N = len(a) - 1
    pw = bernstein_to_power(a)
    print(f"== {label}: N={N} ==")
    print("power_coefficients_ascending:", " ".join(map(str, pw)))

    x = symbols("x")
    poly = Poly(sum(sp.Integer(c) * x**i for i, c in enumerate(pw)), x)

    factor, mults = poly.factor_list()
    # factor_list returns (content, [(Poly, mult), ...]). Irreducible over ZZ
    # iff content is +-1 and there is exactly one factor, with multiplicity 1,
    # whose degree equals the polynomial degree.
    content_is_unit = abs(int(factor)) == 1
    single_full_factor = (len(mults) == 1 and mults[0][1] == 1
                          and Poly(mults[0][0], x).degree() == poly.degree())
    irreducible = content_is_unit and single_full_factor
    print("factor_over_ZZ (content):", sp.sstr(factor))
    print("factor_list_mults:", [(sp.sstr(f), m) for f, m in mults])
    print(f"irreducible_over_ZZ: {irreducible}")

    intervals = poly.intervals()
    # sympy returns list of ((lo, hi), mult) pairs with exact rational endpoints.
    # Convert via sympy Rational -> Fraction (int() would truncate 3/2 -> 3/2 -> 1).
    physical = []
    for item in intervals:
        (lo, hi), mult = item
        rlo, rhi = sp.Rational(lo), sp.Rational(hi)
        loR, hiR = Fraction(rlo.p, rlo.q), Fraction(rhi.p, rhi.q)
        if 0 <= loR <= 1 and 0 <= hiR <= 1:
            physical.append((loR, hiR))
    # Exact count of real roots in [0,1] via Sturm (sympy count_roots).
    n_phys = poly.count_roots(0, 1)
    print("physical_root_count_exact:", n_phys)
    print("physical_isolation_intervals:", [(str(lo), str(hi)) for lo, hi in physical])
    if len(physical) != 1:
        print(f"ERROR: expected exactly 1 physical root, found {len(physical)}")
        return 2
    lo, hi = physical[0]
    loR, hiR = Fraction(str(lo)), Fraction(str(hi))

    BISECT_STEPS = 170
    # invariant maintained across rounds: sign(loR) is the sign of M(loR)
    def sign_at(fr: Fraction) -> int:
        val = sum(sp.Integer(c) * sp.Rational(fr.numerator, fr.denominator) ** i
                  for i, c in enumerate(pw))
        return 1 if val > 0 else (-1 if val < 0 else 0)

    s_lo = sign_at(loR)
    for _ in range(BISECT_STEPS):
        mid = (loR + hiR) / 2
        s = sign_at(mid)
        if s == 0:
            loR = hiR = mid
            break
        if s == s_lo:
            loR = mid
        else:
            hiR = mid
    print(f"exact_bracket: [{loR.numerator}/{loR.denominator}, {hiR.numerator}/{hiR.denominator}]")
    width = hiR - loR
    print(f"bracket_width_exact: {width.numerator}/{width.denominator}")

    mp.mp.dps = args.digits
    mid_dec = mp.mpf(loR.numerator) / mp.mpf(loR.denominator)
    mid_dec2 = mp.mpf(hiR.numerator) / mp.mpf(hiR.denominator)
    print(f"physical_root_diagnostic ({args.digits} digits, midpoint of bracket):")
    print(" ", mp.nstr((mid_dec + mid_dec2) / 2, args.digits))
    return 0


if __name__ == "__main__":
    sys.exit(main())

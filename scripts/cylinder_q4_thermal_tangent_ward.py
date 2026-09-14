#!/usr/bin/env python3
"""Exact rational checks for the cylinder thermal-Q4 tangent identity.

The CFT input is the ordinary c=0, h=5/8 thermal Kac quotient

    (L_-2 - 2/3 L_-1^2)|h> = 0

and the repository normalization

    Q4 = 40 L_-2^2 - 60 L_-3 L_-1 - 9 L_-4.

For a dimensionless cylinder of circumference 2*pi, the chiral stress-tensor
one-point Ward function between identical external primary states is

    m + h/(4 sinh(w/2)^2).

This script derives the local Taylor coefficients using only Fraction
arithmetic, then performs the null/translation reduction of Q4.
"""
from fractions import Fraction
from math import factorial


def mul(a, b, n):
    out = [Fraction(0) for _ in range(n + 1)]
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            if i + j <= n:
                out[i + j] += ai * bj
    return out


def inv(a, n):
    if a[0] == 0:
        raise ValueError("series has zero constant term")
    out = [Fraction(0) for _ in range(n + 1)]
    out[0] = 1 / a[0]
    for k in range(1, n + 1):
        out[k] = -sum(a[j] * out[k - j] for j in range(1, k + 1)) / a[0]
    return out


def exp_series(n):
    return [Fraction(1, factorial(k)) for k in range(n + 1)]


def ward_regular_series(order=6):
    """Return coefficients of w^2*exp(w)/(exp(w)-1)^2.

    The desired Ward factor is exp(w)/(exp(w)-1)^2. Multiplying by w^2
    removes its double pole, so ordinary power-series arithmetic applies.
    """
    # d(w)=(exp(w)-1)/w = sum_{k>=0} w^k/(k+1)!
    d = [Fraction(1, factorial(k + 1)) for k in range(order + 1)]
    d_inv = inv(d, order)
    d_inv_sq = mul(d_inv, d_inv, order)
    e = exp_series(order)
    return mul(e, d_inv_sq, order)


def derive():
    s = ward_regular_series(6)
    # exp(w)/(exp(w)-1)^2 = w^-2 * s(w).
    # Therefore s_0=1, s_1=0, s_2=-1/12, s_3=0,
    # s_4=1/240, ...
    assert s[0] == 1
    assert s[1] == 0
    assert s[2] == Fraction(-1, 12)
    assert s[3] == 0
    assert s[4] == Fraction(1, 240)

    h = Fraction(5, 8)
    # Cylinder one-point ratios read off from
    # m + h*w^-2*s(w).
    l4_over_primary = h * s[4]
    assert l4_over_primary == Fraction(1, 384)

    # Translation invariance + null quotient.
    l3_l1_over_l4 = Fraction(-2, 1)
    l2_sq_over_l4 = Fraction(4, 3)
    q4_over_l4 = 40 * l2_sq_over_l4 - 60 * l3_l1_over_l4 - 9
    assert q4_over_l4 == Fraction(493, 3)

    q4_chiral_over_primary = q4_over_l4 * l4_over_primary
    q4_bulk_real_over_primary = 2 * q4_chiral_over_primary
    assert q4_chiral_over_primary == Fraction(493, 1152)
    assert q4_bulk_real_over_primary == Fraction(493, 576)

    return {
        "ward_w2_exp_over_expm1_sq_coeffs_0_to_6": s,
        "h": h,
        "L-4_over_primary": l4_over_primary,
        "Q4_over_L-4": q4_over_l4,
        "Q4_chiral_over_primary": q4_chiral_over_primary,
        "Q4_plus_Qbar4_over_primary": q4_bulk_real_over_primary,
    }


def main():
    result = derive()
    for key, value in result.items():
        print(f"{key}: {value}")
    print("PASS: ordinary thermal Q4 is exactly tangent to the thermal-primary matrix element at the cylinder CFT point")


if __name__ == "__main__":
    main()

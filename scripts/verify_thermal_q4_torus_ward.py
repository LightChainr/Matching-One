#!/usr/bin/env python3
"""Exact rational checker for the thermal Q4 torus Ward reduction.

This script does not evaluate any percolation data.  It checks only the rational
Virasoro/null-state algebra used in notes/thermal-q4-torus-ward-identity.md.

Input identities:

    [L_-1, L_-3] = 2 L_-4
    (L_-2 - 2/3 L_-1^2) |h=5/8> = 0
    <L_-1 Psi>_torus = 0
    <L_-4 phi>_torus = (g2*h/20) <phi>_torus

Repository quasiprimary normalization:

    Q4 = 40 L_-2^2 - 60 L_-3 L_-1 - 9 L_-4.
"""

from fractions import Fraction


def derive_coefficients():
    """Return exact coefficients relative to C=<L_-4 phi>."""
    h = Fraction(5, 8)

    # 0 = <L_-1 L_-3 phi>
    #   = 2 <L_-4 phi> + <L_-3 L_-1 phi>.
    l3_l1_over_c = Fraction(-2, 1)

    # Null relation and translation invariance:
    # <L_-2^2 phi> = (2/3)<L_-2 L_-1^2 phi>
    #                 = -(2/3)<L_-3 L_-1 phi>.
    l2_sq_over_c = -Fraction(2, 3) * l3_l1_over_c

    q4_over_c = (
        40 * l2_sq_over_c
        - 60 * l3_l1_over_c
        - 9
    )

    # Brehm--Runkel: C/<phi> = g2*h/20.
    c_over_g2_primary = h / 20
    q4_over_g2_primary = q4_over_c * c_over_g2_primary

    return {
        "h": h,
        "L-3L-1_over_C": l3_l1_over_c,
        "L-2^2_over_C": l2_sq_over_c,
        "Q4_over_C": q4_over_c,
        "C_over_g2_primary": c_over_g2_primary,
        "Q4_over_g2_primary": q4_over_g2_primary,
    }


def self_test():
    c = derive_coefficients()
    assert c["L-3L-1_over_C"] == Fraction(-2, 1)
    assert c["L-2^2_over_C"] == Fraction(4, 3)
    assert c["Q4_over_C"] == Fraction(493, 3)
    assert c["C_over_g2_primary"] == Fraction(1, 32)
    assert c["Q4_over_g2_primary"] == Fraction(493, 96)
    return c


def main():
    c = self_test()
    for key, value in c.items():
        print(f"{key}: {value}")
    print("PASS: <Q4 phi>/<phi> = (493/96) g2(tau)")


if __name__ == "__main__":
    main()

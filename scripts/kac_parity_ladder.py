#!/usr/bin/env python3
"""Exact arithmetic audit for the conditional c=0 Potts Kac-parity ladder.

Exact statements checked here:
- c=0 critical Potts Kac weights h_(r,s)=((2r-3s)^2-1)/24;
- dimensions x=2h for diagonal fields V_<1,s>;
- finite-size correction exponent omega=x-2 for a dimensionless torus observable;
- relative exponents to the leading candidate matching-odd H4 field x=21/4;
- root exponent after division by M'~L^(3/4).

Conditional statement (NOT proved here): if the continuum matching involution is an
OPE/interchiral automorphism with V_<1,2> odd and the usual adjacent Kac fusion,
then eta_s=(-1)^(s-1).  The script labels this separately as a hypothesis.
"""

from fractions import Fraction


def h_kac(r: int, s: int) -> Fraction:
    return Fraction((2 * r - 3 * s) ** 2 - 1, 24)


def x_diag(r: int, s: int) -> Fraction:
    return 2 * h_kac(r, s)


def parity_hypothesis(s: int) -> int:
    return 1 if (s - 1) % 2 == 0 else -1


def fmt(q: Fraction) -> str:
    return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"


def main() -> int:
    leading_x = Fraction(21, 4)
    leading_M_power = leading_x - 2
    thermal_y = Fraction(3, 4)

    rows = []
    for s in range(1, 6):
        x = x_diag(1, s)
        omega = x - 2
        relative_q = omega - leading_M_power
        root_power = omega + thermal_y
        rows.append((s, x, parity_hypothesis(s), omega, relative_q, root_power))

    print("conditional V_<1,s> matching-parity ladder")
    print("s  x       eta  M exponent x-2   q vs T4   root exponent")
    for s, x, eta, omega, q, root in rows:
        print(f"{s:<2} {fmt(x):<7} {eta:+d}   {fmt(omega):<14} {fmt(q):<9} {fmt(root)}")

    assert x_diag(1, 1) == 0
    assert x_diag(1, 2) == Fraction(5, 4)
    assert x_diag(1, 3) == 4
    assert x_diag(1, 4) == Fraction(33, 4)
    assert leading_M_power == Fraction(13, 4)
    assert leading_M_power + thermal_y == 4

    v13_omega = x_diag(1, 3) - 2
    assert v13_omega == 2
    assert 4 + v13_omega == 6

    v14_omega = x_diag(1, 4) - 2
    assert v14_omega == Fraction(25, 4)
    assert v14_omega - leading_M_power == 3
    assert v14_omega + thermal_y == 7

    print("\nexact focal checks: PASS")
    print("H4*H0 -> H4; H4*H4 -> H0+H8 (trigonometric selection)")
    print("matching parities are conditional on the unproved continuum involution hypothesis")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

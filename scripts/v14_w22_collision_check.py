#!/usr/bin/env python3
"""Exact algebra for the V_<1,4> / W(2,2)-top collision at percolation.

Use t=beta^2 and critical-Potts weights

    h_{r,s}=(c-1)/24 + 1/4 (r beta - s/beta)^2.

The common central-charge term cancels from the difference between h_{1,4}
and h_{2,-2}.  Percolation has t=2/3 and Q=4 cos^2(pi t).
"""
from fractions import Fraction
import math


def delta_h(t: Fraction) -> Fraction:
    """h_{1,4}-h_{2,-2}, with t=beta^2."""
    return Fraction(1, 4) * (-3 * t - 16 + Fraction(12, 1) / t)


def d_delta_x_dt(t: Fraction) -> Fraction:
    """Derivative of x_{1,4}-x_{2,-2}=2 delta_h with respect to t."""
    return Fraction(1, 2) * (-3 - Fraction(12, 1) / (t * t))


def main() -> None:
    t = Fraction(2, 3)
    assert delta_h(t) == 0
    assert d_delta_x_dt(t) == -15

    # Zero condition: 3 t^2 + 16 t - 12 = 0.
    # Physical roots are t=2/3 and t=-6; only the first is positive.
    assert 3 * t * t + 16 * t - 12 == 0

    dQ_dt = 2 * math.pi * math.sqrt(3.0)  # at t=2/3
    dx_dQ = -15.0 / dQ_dt

    print("t_percolation = 2/3")
    print("Delta h = h_(1,4)-h_(2,-2) = 0")
    print("d_t Delta x = -15")
    print("d_t Q = 2*pi*sqrt(3)")
    print("d_Q Delta x = -15/(2*pi*sqrt(3)) =", dx_dQ)
    print("PASS: V_<1,4> collides with the W(2,2) diagonal top/bottom weight exactly at Q=1")


if __name__ == "__main__":
    main()

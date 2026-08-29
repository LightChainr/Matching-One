#!/usr/bin/env python3
"""Exact Q=1 tangent oracle for Cai's CLE boundary four-point ODE.

The script differentiates the general-kappa third-order ODE in
arXiv:2603.28161v2 and constructs the high Frobenius branch
V_{3h+1}.  All series coefficients and their kappa derivatives are
computed with exact rational dual numbers; only the displayed decimal
approximations use floating point arithmetic.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Dict, Iterable, List


@dataclass(frozen=True)
class Dual:
    """Value and first derivative with respect to kappa."""

    value: Fraction
    derivative: Fraction = Fraction(0)

    @staticmethod
    def coerce(value: object) -> "Dual":
        if isinstance(value, Dual):
            return value
        return Dual(Fraction(value))

    def __add__(self, other: object) -> "Dual":
        rhs = self.coerce(other)
        return Dual(self.value + rhs.value, self.derivative + rhs.derivative)

    __radd__ = __add__

    def __neg__(self) -> "Dual":
        return Dual(-self.value, -self.derivative)

    def __sub__(self, other: object) -> "Dual":
        return self + (-self.coerce(other))

    def __rsub__(self, other: object) -> "Dual":
        return self.coerce(other) - self

    def __mul__(self, other: object) -> "Dual":
        rhs = self.coerce(other)
        return Dual(
            self.value * rhs.value,
            self.derivative * rhs.value + self.value * rhs.derivative,
        )

    __rmul__ = __mul__

    def reciprocal(self) -> "Dual":
        if self.value == 0:
            raise ZeroDivisionError("dual reciprocal at zero")
        return Dual(
            1 / self.value,
            -self.derivative / (self.value * self.value),
        )

    def __truediv__(self, other: object) -> "Dual":
        return self * self.coerce(other).reciprocal()

    def __rtruediv__(self, other: object) -> "Dual":
        return self.coerce(other) / self

    def __pow__(self, exponent: int) -> "Dual":
        if exponent < 0:
            return (self.reciprocal()) ** (-exponent)
        result = Dual(Fraction(1))
        for _ in range(exponent):
            result *= self
        return result


def falling(value: Dual, order: int) -> Dual:
    result = Dual(Fraction(1))
    for offset in range(order):
        result *= value - offset
    return result


def cle_operator(kappa: Dual) -> Dict[int, List[Dual]]:
    """Return polynomial coefficients a_j(lambda) of sum a_j D^j.

    Lists are in ascending powers of lambda.  This is the unscaled ODE
    in Theorem 1.3 / Eq. (1.10) of arXiv:2603.28161v2.
    """

    zero = Dual(Fraction(0))
    a3 = kappa**3 / 2
    a2 = kappa**2 * (3 * kappa - 16)
    b1 = 18 * kappa**2 - 212 * kappa + 608
    a1_constant = 3 * kappa * (kappa - 4) * (kappa - 8)
    a0_scale = 6 * (kappa - 4) * (kappa - 8) ** 2
    return {
        3: [zero, zero, a3, -2 * a3, a3],
        2: [zero, a2, -3 * a2, 2 * a2],
        1: [a1_constant, -kappa * b1, kappa * b1],
        0: [-a0_scale, 2 * a0_scale],
    }


def frobenius_high_branch(order: int) -> tuple[Dual, List[Dual]]:
    """Return r=24/kappa-2 and c_0..c_order for lambda^r sum c_n lambda^n."""

    kappa = Dual(Fraction(6), Fraction(1))
    exponent = 24 / kappa - 2
    operator = cle_operator(kappa)
    coefficients: List[Dual] = [Dual(Fraction(1))]

    for m in range(1, order + 1):
        known = Dual(Fraction(0))
        pivot = Dual(Fraction(0))
        for derivative_order, polynomial in operator.items():
            for power, polynomial_coefficient in enumerate(polynomial):
                n = m - 1 + derivative_order - power
                if n < 0 or n > m:
                    continue
                factor = polynomial_coefficient * falling(
                    exponent + n, derivative_order
                )
                if n == m:
                    pivot += factor
                else:
                    known += factor * coefficients[n]
        if pivot.value == 0:
            raise ZeroDivisionError(f"resonant Frobenius pivot at order {m}")
        coefficients.append(-known / pivot)

    return exponent, coefficients


def _fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def _poly_payload(values: Iterable[Fraction]) -> List[str]:
    return [_fraction_text(value) for value in values]


def _operator_at_six_payload() -> dict:
    operator = cle_operator(Dual(Fraction(6), Fraction(1)))
    # Cai divides the kappa=6 equation by 12.  Apply the same fixed scale to
    # both L_6 and partial_kappa L|_6.
    return {
        f"D{order}": {
            "L6": _poly_payload(item.value / 12 for item in polynomial),
            "d_kappa_L_at_6": _poly_payload(
                item.derivative / 12 for item in polynomial
            ),
        }
        for order, polynomial in sorted(operator.items(), reverse=True)
    }


def build_payload(order: int) -> dict:
    exponent, coefficients = frobenius_high_branch(order)
    kappa_q_scale = -Fraction(3, 2)  # d kappa/dQ = scale * sqrt(3)/pi
    exponent_q_scale = kappa_q_scale * exponent.derivative
    coefficient_rows = []
    for n, coefficient in enumerate(coefficients):
        regular_q_scale = kappa_q_scale * coefficient.derivative
        coefficient_rows.append(
            {
                "n": n,
                "c_n_at_Q1": _fraction_text(coefficient.value),
                "d_kappa_c_n_at_6": _fraction_text(coefficient.derivative),
                "d_Q_c_n_sqrt3_over_pi_scale": _fraction_text(regular_q_scale),
                "d_Q_c_n_decimal": float(regular_q_scale)
                * math.sqrt(3)
                / math.pi,
            }
        )

    return {
        "schema": "matching-one.p263-boundary-tangent-ode.v1",
        "source": {
            "title": "Boundary four-point connectivities of conformal loop ensembles",
            "author": "Gefei Cai",
            "arxiv": "2603.28161v2",
            "url": "https://arxiv.org/abs/2603.28161",
            "equation": "general-kappa third-order boundary four-point ODE",
        },
        "fk_cle_map": {
            "equation": "sqrt(Q)=-2 cos(4 pi/kappa)",
            "kappa_at_Q1": "6",
            "d_kappa_d_Q_at_Q1": "-3 sqrt(3)/(2 pi)",
            "h": "8/kappa-1",
            "d_h_d_Q_at_Q1": "sqrt(3)/(3 pi)",
        },
        "inhomogeneous_tangent_equation": {
            "definition": "T(lambda)=partial_Q U_kappa(lambda)|_{Q=1}",
            "equation": "L6[T]=(3 sqrt(3)/(2 pi)) (partial_kappa L)_6[U0]",
            "operator_polynomials_ascending_lambda_power": _operator_at_six_payload(),
        },
        "full_green_function_tangent": {
            "definition": "G=K(x_i)^(2h) U(lambda)",
            "equation": "partial_Q G|_1=K^(2/3)[T+(2 sqrt(3)/(3 pi)) log(K) U0]",
            "warning": "The conformal-prefactor derivative is distinct from the cross-ratio tangent.",
        },
        "high_frobenius_branch": {
            "branch": "V_{3h+1}=lambda^r sum_n c_n lambda^n, c_0=1",
            "r_at_Q1": _fraction_text(exponent.value),
            "d_kappa_r_at_6": _fraction_text(exponent.derivative),
            "d_Q_r_sqrt3_over_pi_scale": _fraction_text(exponent_q_scale),
            "d_Q_r": "sqrt(3)/pi",
            "universal_log_statement": "partial_Q V|_1 contains (sqrt(3)/pi) V_2(lambda) log(lambda)",
            "coefficients": coefficient_rows,
        },
        "amplitude_gauge": {
            "statement": "T is defined modulo a homogeneous solution when the Q-dependent normalization is not fixed.",
            "gauge_free_shape": "partial_Q log(U(lambda)/U(lambda_anchor))|_1",
            "invariant": "The leading high-branch log coefficient sqrt(3)/pi is unchanged by analytic amplitude rescaling.",
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build the exact #263 Q=1 boundary tangent ODE oracle."
    )
    parser.add_argument("--order", type=int, default=6)
    parser.add_argument("--json", type=Path)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.order < 0:
        raise SystemExit("--order must be non-negative")
    payload = build_payload(args.order)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""First-order T Tbar shape operator and the P321 ratio obstruction.

This is an algebraic oracle for the small-coupling expansion of Cardy's
real-analytic torus one-point deformation (arXiv:2201.00478, Theorem 2).
It does not assume that the P321 homology projector is transparent.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Sequence


def first_order_mode_multiplier(
    scaling_dimension: float,
    radial_frequency: float,
    angular_frequency: float,
    modulus_imaginary_part: float,
) -> float:
    """Multiplier of one real-analytic Fourier mode at first order.

    The undeformed mode is
    exp(-2*pi*x*y + 2*pi*i*p*theta).  After normalizing Cardy's deformed
    form to have unit zeroth-order limit, its alpha coefficient is

      4*pi^2*(x^2-p^2)*y^2 - 2*pi*(1+k)*x*y.
    """

    k = float(scaling_dimension)
    x = float(radial_frequency)
    p = float(angular_frequency)
    y = float(modulus_imaginary_part)
    if y <= 0:
        raise ValueError("modulus imaginary part must be positive")
    return 4.0 * math.pi**2 * (x * x - p * p) * y * y - 2.0 * math.pi * (1.0 + k) * x * y


def normalized_deformed_mode(
    scaling_dimension: float,
    radial_frequency: float,
    angular_frequency: float,
    modulus_imaginary_part: float,
    modulus_real_part: float,
    alpha: float,
) -> complex:
    """Evaluate the normalized single-mode factor in Cardy's Theorem 2."""

    k = float(scaling_dimension)
    x = float(radial_frequency)
    p = float(angular_frequency)
    y = float(modulus_imaginary_part)
    theta = float(modulus_real_part)
    a = float(alpha)
    if y <= 0 or a <= 0:
        raise ValueError("y and alpha must be positive")
    radical = math.sqrt(1.0 + 8.0 * math.pi * x * a * y + (4.0 * math.pi * p * a * y) ** 2)
    prefactor = ((1.0 + radical) ** (1.0 - k) / radical) / (2.0 ** (1.0 - k))
    exponent = -(radical - 1.0) / (2.0 * a)
    phase = complex(math.cos(2.0 * math.pi * p * theta), math.sin(2.0 * math.pi * p * theta))
    return prefactor * math.exp(exponent) * phase


def undeformed_mode(
    radial_frequency: float,
    angular_frequency: float,
    modulus_imaginary_part: float,
    modulus_real_part: float,
) -> complex:
    x = float(radial_frequency)
    p = float(angular_frequency)
    y = float(modulus_imaginary_part)
    theta = float(modulus_real_part)
    phase = complex(math.cos(2.0 * math.pi * p * theta), math.sin(2.0 * math.pi * p * theta))
    return math.exp(-2.0 * math.pi * x * y) * phase


def product_rule_obstruction() -> dict[str, Any]:
    """Return the exact scalarized product-rule identity for the operator."""

    return {
        "operator": "L_k[f]=y^2*(d_y^2+d_theta^2)f+(1+k)*y*d_y f",
        "product_identity": (
            "L_(kf+kg)[f*g]-g*L_kf[f]-f*L_kg[g]="
            "2*y^2*grad(f).grad(g)+kg*y*g*d_y(f)+kf*y*f*d_y(g)"
        ),
        "ratio_identity": (
            "L_(kf+kg)[f*g]/(f*g)-L_kf[f]/f="
            "L_kg[g]/g+2*y^2*grad(log f).grad(log g)+"
            "kg*y*d_y(log f)+kf*y*d_y(log g)"
        ),
        "p321_substitution": {
            "f": "restricted thermal one-point F_t",
            "g": "E4 factor in F_u=kappa*E4*F_t",
            "consequence": (
                "The first subleading root ratio contains derivatives of F_t; "
                "the leading E4 ratio alone does not determine D_width/C_width."
            ),
        },
    }


def build_artifact() -> dict[str, Any]:
    return {
        "schema": "matching-one/p321-ttbar-shape-obstruction/v1",
        "data_class": "exact first-order operator identity; no new Monte Carlo",
        "source": {
            "cardy": "https://arxiv.org/abs/2201.00478",
            "he_sun": "https://arxiv.org/abs/2004.07486",
        },
        "fourier_mode_multiplier": (
            "4*pi^2*(x^2-p^2)*y^2-2*pi*(1+k)*x*y"
        ),
        **product_rule_obstruction(),
        "decision": {
            "closed": (
                "An E4-only parameter-free formula for the N^-3/N^-2 root-amplitude "
                "ratio does not follow from the already established leading E4 curve."
            ),
            "next_exact_object": (
                "Compute the homology-resolved thermal one-point F_t(tau), or an "
                "equivalent defect/tube-algebra logarithmic derivative, before naming "
                "the subleading shape as ordinary identity dressing."
            ),
            "empirical_use": (
                "The deep equal-area stream may estimate D/C and test a supplied "
                "F_t model, but must not fit an arbitrary E4-only correction after reveal."
            ),
        },
        "spin_boundary": (
            "Q4 epsilon is spin four, whereas Cardy's displayed Theorem 2 is scalar. "
            "The scalarized product identity is therefore a favorable no-cancellation "
            "test, not a claimed spinful descendant formula."
        ),
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    rendered = json.dumps(build_artifact(), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

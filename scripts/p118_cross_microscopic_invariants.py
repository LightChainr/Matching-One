#!/usr/bin/env python3
"""Exact normalization bookkeeping for Issue #118.

The oracle deliberately separates two statements:

* independent microscopic couplings forbid a nonconstant invariant made from
  one bare H4 amplitude and one bare E6 amplitude;
* using two matched moduli produces a double ratio in which both couplings,
  both field normalizations, and both metric factors cancel.

It also records the projective Gram diagnostic for a logarithmic pair.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any


def raw_charge(p: int, q: int) -> tuple[int, int]:
    """Gauge charge of A4**p A6**q under independent couplings."""

    return p, q


def double_ratio(a4_j: complex, a4_k: complex, a6_j: complex, a6_k: complex) -> complex:
    """Weight-12 cross-microscopic shape double ratio."""

    return (a4_j / a4_k) ** 3 / (a6_j / a6_k) ** 2


def polynomial_null(a4_j: complex, a4_k: complex, a6_j: complex, a6_k: complex) -> complex:
    """Division-free form of double_ratio == 1."""

    return a4_j**3 * a6_k**2 - a4_k**3 * a6_j**2


def hex_child_forms(scale4: complex = 1, scale6: complex = 1) -> tuple[list[complex], list[complex]]:
    """Degree-2 hex-child E4 and E6 character vectors."""

    zeta = complex(-0.5, math.sqrt(3) / 2)
    return (
        [scale4, scale4 * zeta, scale4 * zeta**2],
        [scale6, scale6, scale6],
    )


def jordan_k(ll: float, ld: float, dd: float) -> float:
    """Projective Gram coordinate K=LL*DD/LD^2 (J in #234 is K-1)."""

    if ld == 0:
        raise ZeroDivisionError("the mixed pairing LD must be nonzero")
    return ll * dd / (ld * ld)


def rescale_gram(
    ll: float, ld: float, dd: float, local_scale: float, top_scale: float
) -> tuple[float, float, float]:
    return (
        local_scale**2 * ll,
        local_scale * top_scale * ld,
        top_scale**2 * dd,
    )


def shear_top(ll: float, ld: float, dd: float, shear: float) -> tuple[float, float, float]:
    """Apply top -> top + shear*bottom to a symmetric two-field Gram form."""

    return ll, ld + shear * ll, dd + 2 * shear * ld + shear**2 * ll


def build_artifact() -> dict[str, Any]:
    a4, a6 = hex_child_forms()
    pair_values: dict[str, dict[str, Any]] = {}
    for j in range(3):
        for k in range(j + 1, 3):
            value = double_ratio(a4[j], a4[k], a6[j], a6[k])
            pair_values[f"{j}:{k}"] = {
                "exact": "1",
                "character_exponent_mod_3": (3 * (j - k)) % 3,
                "floating_oracle_absolute_error": abs(value - 1),
            }

    return {
        "schema": "matching-one/cross-microscopic-amplitude-invariants/v1",
        "issue": 118,
        "status": "exact_gauge_counting_and_conditional_modular_prediction",
        "raw_single_modulus_no_go": {
            "gauge_group": "(C*)_H4 x (C*)_E6",
            "charge_matrix": [[1, 0], [0, 1]],
            "kernel_dimension": 0,
            "statement": "A4^p A6^q is invariant only for p=q=0.",
        },
        "cross_microscopic_double_ratio": {
            "formula": "U46(j,k)=(A4_square(j)/A4_square(k))^3/(A6_triangular(j)/A6_triangular(k))^2",
            "division_free_null": "A4_square(j)^3*A6_triangular(k)^2-A4_square(k)^3*A6_triangular(j)^2=0",
            "powers_reason": "3*weight(E4)=2*weight(E6)=12",
            "cancelled": [
                "square_H4_irrelevant_coupling",
                "triangular_E6_irrelevant_coupling",
                "both_field_normalizations",
                "both_constant_metric_factors",
            ],
            "hex_degree2_child_pair_values": pair_values,
            "frozen_target": "U46(j,k)=1 for all three child pairs",
            "conditional_bridge": "the two lattice responses must be typed as the same-frame E4 and E6 torus one-point shapes",
        },
        "log_pair": {
            "coflow_observables": ["LL", "LD", "DD"],
            "K": "LL*DD/LD^2",
            "existing_issue_234_J": "K-1",
            "rescaling_invariant": True,
            "finite_cutoff_shear_invariant": False,
            "jordan_limit": "LL->0 with LD!=0 implies K->0 and J->-1, invariant under top->top+alpha*bottom",
            "minimal_required_observable": "the same-flow mixed pairing LD; LL and DD alone cannot cancel independent field normalizations",
        },
        "claim_boundary": {
            "exact": "gauge no-go, cancellation in U46, C3 character closure, and Gram transformation laws",
            "conditional": "identification of microscopic H4/E6 responses and the continuum Jordan limit",
            "not_claimed": "universality of either bare amplitude or the finite-cutoff Issue #234 J values",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    artifact = build_artifact()
    rendered = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()

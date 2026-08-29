#!/usr/bin/env python3
"""Exact response generating function and operators for Krawtchouk jets.

The functions here operate on the score coefficients emitted by
``threshold_score_modes.py``.  They separate exact finite-N identities from
the scaling hypotheses recorded in
``predictions/hermite_krawtchouk_jet_20260829.yaml``.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from typing import Sequence

import mpmath as mp


def response_from_modes(
    coefficients: Sequence[mp.mpf], n: int, p0: mp.mpf, p: mp.mpf
) -> mp.mpf:
    """Evaluate the exact finite-N response from orthonormal score modes.

    For the positive-score convention used by ``threshold_score_modes.py``,

      E_p[H_r^(p0)] = sqrt(C(N,r)) ((p-p0)/sqrt(p0(1-p0)))^r.
    """

    if len(coefficients) > n + 1:
        raise ValueError("a degree-N response has at most N+1 score modes")
    if not 0 < p0 < 1 or not 0 <= p <= 1:
        raise ValueError("probabilities lie outside their valid intervals")
    coordinate = (p - p0) / mp.sqrt(p0 * (1 - p0))
    return mp.fsum(
        coefficient * mp.sqrt(math.comb(n, order)) * coordinate**order
        for order, coefficient in enumerate(coefficients)
    )


def scaling_derivative_jet(
    coefficients: Sequence[mp.mpf],
    n: int,
    p0: mp.mpf,
    alpha: mp.mpf,
) -> list[mp.mpf]:
    """Convert score coefficients to derivatives of the scaling function.

    If c_r ~ N^(-alpha-r/8) and z=(p-p0)N^(3/8), then the limiting
    derivative is

      d_r = sqrt(r!) N^(alpha+r/8) c_r / (p0(1-p0))^(r/2).
    """

    if n <= 0 or not 0 < p0 < 1:
        raise ValueError("invalid size or center")
    variance = p0 * (1 - p0)
    return [
        mp.sqrt(math.factorial(order))
        * mp.power(n, alpha + mp.mpf(order) / 8)
        * coefficient
        / mp.power(variance, mp.mpf(order) / 2)
        for order, coefficient in enumerate(coefficients)
    ]


def translate_jet(jet: Sequence[mp.mpf], displacement: mp.mpf) -> list[mp.mpf]:
    """Return the truncated derivative jet of F(z+displacement)."""

    return [
        mp.fsum(
            jet[order + step] * displacement**step / math.factorial(step)
            for step in range(len(jet) - order)
        )
        for order in range(len(jet))
    ]


def dilate_jet(jet: Sequence[mp.mpf], scale: mp.mpf) -> list[mp.mpf]:
    """Return the derivative jet of F(scale*z)."""

    return [value * scale**order for order, value in enumerate(jet)]


def width_normalized_jet(
    jet: Sequence[mp.mpf], scaled_rank_gap: mp.mpf
) -> list[mp.mpf]:
    """Remove a thermal-window width w_N from the derivative jet.

    If F_N(z)=A_N F(z/w_N), then d_(r,N)=A_N w_N^(-r)d_r and multiplying
    by w_N^r removes the width drift while retaining the common amplitude.
    """

    if scaled_rank_gap <= 0:
        raise ValueError("scaled rank gap must be positive")
    return [
        value * scaled_rank_gap**order for order, value in enumerate(jet)
    ]


def width_cross_residual(
    parent: Sequence[mp.mpf],
    child: Sequence[mp.mpf],
    parent_width: mp.mpf,
    child_width: mp.mpf,
) -> list[mp.mpf]:
    """Division-free residual for common shape up to width and amplitude."""

    if len(parent) != len(child) or not parent:
        raise ValueError("parent and child jets must have equal nonzero length")
    left = width_normalized_jet(parent, parent_width)
    right = width_normalized_jet(child, child_width)
    return [
        right[order] * left[0] - left[order] * right[0]
        for order in range(1, len(parent))
    ]


def cocycle_residual(
    parent: Sequence[mp.mpf],
    norm2: Sequence[mp.mpf],
    norm5: Sequence[mp.mpf],
    multiplier: mp.mpf,
) -> list[mp.mpf]:
    """Return Z_5N-c Z_2N+(c-1)Z_N component by component."""

    if not (len(parent) == len(norm2) == len(norm5)):
        raise ValueError("cocycle jets must have equal length")
    return [
        norm5[index]
        - multiplier * norm2[index]
        + (multiplier - 1) * parent[index]
        for index in range(len(parent))
    ]


def pooled_gap_convention_shift(
    endpoint_shifts: Sequence[tuple[int, int]],
) -> Fraction:
    """Mean change in G under integer Kminus/Kplus convention shifts.

    Each tuple is ``(delta_Kminus, delta_Kplus)`` for one equally weighted
    orientation.  With two orientations, a convention-only pooled shift lies
    in one-half integers; a quarter cannot be forced by the frozen rank
    off-by-one convention.
    """

    if not endpoint_shifts:
        raise ValueError("at least one orientation is required")
    return Fraction(
        sum(delta_plus - delta_minus for delta_minus, delta_plus in endpoint_shifts),
        len(endpoint_shifts),
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-order", type=int, default=6)
    args = parser.parse_args()
    if not 1 <= args.max_order <= 12:
        raise SystemExit("--max-order must lie in 1..12")
    order = args.max_order
    payload = {
        "exact_generating_function": (
            "R_N(p)=sum_r c_r sqrt(C(N,r)) "
            "((p-p0)/sqrt(p0(1-p0)))^r"
        ),
        "scaling_derivative_jet": (
            "d_r=sqrt(r!)*N^(alpha+r/8)*c_r/(p0(1-p0))^(r/2)"
        ),
        "translation_generator": [f"d_{r + 1}" for r in range(order)],
        "width_generator": [f"{r}*d_{r}" for r in range(order + 1)],
        "rank_gap_width": "w_N=N^(-5/8) E[K_plus-K_minus]",
    }
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

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
    """Convert score coefficients to the exact finite-N scaling jet.

    Substituting ``p-p0=z*N^(-3/8)`` in the exact finite-N response and
    multiplying by ``N^alpha`` gives

      d_fin_r = r! N^alpha sqrt(C(N,r)) N^(-3r/8) c_r
                / (p0(1-p0))^(r/2).
    """

    if n <= 0 or not 0 < p0 < 1:
        raise ValueError("invalid size or center")
    if len(coefficients) > n + 1:
        raise ValueError("a degree-N response has at most N+1 score modes")
    variance = p0 * (1 - p0)
    return [
        math.factorial(order)
        * mp.sqrt(math.comb(n, order))
        * mp.power(n, alpha - 3 * mp.mpf(order) / 8)
        * coefficient
        / mp.power(variance, mp.mpf(order) / 2)
        for order, coefficient in enumerate(coefficients)
    ]


def asymptotic_scaling_derivative_jet(
    coefficients: Sequence[mp.mpf],
    n: int,
    p0: mp.mpf,
    alpha: mp.mpf,
) -> list[mp.mpf]:
    """Return the large-N view obtained by replacing ``(N)_r`` with ``N^r``."""

    if n <= 0 or not 0 < p0 < 1:
        raise ValueError("invalid size or center")
    if len(coefficients) > n + 1:
        raise ValueError("a degree-N response has at most N+1 score modes")
    variance = p0 * (1 - p0)
    return [
        mp.sqrt(math.factorial(order))
        * mp.power(n, alpha + mp.mpf(order) / 8)
        * coefficient
        / mp.power(variance, mp.mpf(order) / 2)
        for order, coefficient in enumerate(coefficients)
    ]


def finite_jet_factor(n: int, order: int) -> mp.mpf:
    """Return ``d_fin/d_asym=sqrt((N)_r/N^r)`` exactly before the square root."""

    if n <= 0 or not 0 <= order <= n:
        raise ValueError("order must lie in 0..n for positive n")
    falling = math.factorial(n) // math.factorial(n - order)
    return mp.sqrt(mp.mpf(falling) / mp.power(n, order))


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
    jet: Sequence[mp.mpf], canonical_width: mp.mpf
) -> list[mp.mpf]:
    """Remove the canonical dimensionless width from the derivative jet.

    If F_N(z)=A_N F(z/w_can), then d_(r,N)=A_N w_can^(-r)d_r and
    multiplying by w_can^r removes width drift while retaining the common
    amplitude.
    """

    if canonical_width <= 0:
        raise ValueError("canonical width must be positive")
    return [
        value * canonical_width**order for order, value in enumerate(jet)
    ]


def canonical_dimensionless_width(n: int, mean_rank_gap: mp.mpf) -> mp.mpf:
    """Return ``N^(3/8) E[G]/(N+1)`` from the exact canonical area bridge."""

    if n <= 0 or mean_rank_gap < 0:
        raise ValueError("size must be positive and mean rank gap nonnegative")
    return mp.power(n, mp.mpf(3) / 8) * mean_rank_gap / (n + 1)


def rank_normalized_width(n: int, mean_rank_gap: mp.mpf) -> mp.mpf:
    """Return the rank surrogate ``N^(-5/8) E[G]``."""

    if n <= 0 or mean_rank_gap < 0:
        raise ValueError("size must be positive and mean rank gap nonnegative")
    return mp.power(n, -mp.mpf(5) / 8) * mean_rank_gap


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


def microcanonical_matching_sign(
    n: int, occupation: int, k_minus: int, k_plus: int
) -> int:
    """Return the three-valued matching sign for one threshold pair.

    Thresholds may equal ``n+1``, which represents a transition beyond the
    last Bernstein layer.  The neutral plateau is ``k_minus <= k < k_plus``.
    """

    if n < 0 or not 0 <= occupation <= n:
        raise ValueError("occupation must lie in 0..n")
    if not 0 <= k_minus <= k_plus <= n + 1:
        raise ValueError("thresholds must satisfy 0 <= Kminus <= Kplus <= n+1")
    if occupation < k_minus:
        return -1
    if occupation < k_plus:
        return 0
    return 1


def bernstein_basis_integral(n: int, occupation: int) -> Fraction:
    """Integrate ``C(n,k) p^k (1-p)^(n-k)`` exactly over ``p in [0,1]``."""

    if n < 0 or not 0 <= occupation <= n:
        raise ValueError("occupation must lie in 0..n")
    return Fraction(
        math.comb(n, occupation)
        * math.factorial(occupation)
        * math.factorial(n - occupation),
        math.factorial(n + 1),
    )


def canonical_neutral_window_area(
    n: int, k_minus: int, k_plus: int
) -> Fraction:
    """Return the exact canonical area of ``U=1-m^2`` for one permutation."""

    return sum(
        (
            bernstein_basis_integral(n, occupation)
            for occupation in range(n + 1)
            if 1
            - microcanonical_matching_sign(
                n, occupation, k_minus, k_plus
            )
            ** 2
        ),
        Fraction(0),
    )


def mean_canonical_neutral_window_area(
    n: int, threshold_pairs: Sequence[tuple[int, int]]
) -> Fraction:
    """Average the exact canonical neutral-window area over permutations."""

    if not threshold_pairs:
        raise ValueError("at least one threshold pair is required")
    return sum(
        (
            canonical_neutral_window_area(n, k_minus, k_plus)
            for k_minus, k_plus in threshold_pairs
        ),
        Fraction(0),
    ) / len(threshold_pairs)


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
        "primary_finite_scaling_jet": (
            "d_fin_r=r!*N^alpha*sqrt(C(N,r))*N^(-3r/8)*c_r/"
            "(p0(1-p0))^(r/2)"
        ),
        "asymptotic_scaling_jet": (
            "d_asym_r=sqrt(r!)*N^(alpha+r/8)*c_r/"
            "(p0(1-p0))^(r/2)"
        ),
        "finite_to_asymptotic_factor": "sqrt((N)_r/N^r)",
        "translation_generator": [f"d_{r + 1}" for r in range(order)],
        "width_generator": [f"{r}*d_{r}" for r in range(order + 1)],
        "canonical_dimensionless_width": (
            "w_can(N)=N^(3/8) E[K_plus-K_minus]/(N+1)"
        ),
        "rank_width_surrogate": (
            "w_rank(N)=N^(-5/8) E[K_plus-K_minus]="
            "(1+1/N) w_can(N)"
        ),
        "exact_neutral_area_bridge": (
            "integral_0^1 E[U_{K~Bin(N,p)}] dp="
            "E[K_plus-K_minus]/(N+1)"
        ),
    }
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

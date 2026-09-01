#!/usr/bin/env python3
"""Exact Virasoro checks for the c=0 thermal Q4 Jordan descendant.

The calculation is finite rational algebra.  It verifies the level-2 null state,
the level-4 Gram matrix of the repository Q4 basis, the nonzero Q4 norm, and
records the representation-theoretic Jordan inheritance under D=L0+Lbar0.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

C = Fraction(0)
H = Fraction(5, 8)
MU = Fraction(-5, 4)


def category(mode: int) -> int:
    return 0 if mode < 0 else 1 if mode == 0 else 2


@lru_cache(maxsize=None)
def expectation(word: tuple[int, ...]) -> Fraction:
    """Return <h| product L_mode |h> by exact Virasoro normal ordering."""
    modes = list(word)
    for index in range(len(modes) - 1):
        left, right = modes[index], modes[index + 1]
        if category(left) <= category(right):
            continue

        swapped = tuple(modes[:index] + [right, left] + modes[index + 2 :])
        value = expectation(swapped)

        # [L_m,L_n]=(m-n)L_(m+n)+c*m(m^2-1)/12 delta_(m+n,0)
        commutator_mode = left + right
        reduced = tuple(
            modes[:index] + [commutator_mode] + modes[index + 2 :]
        )
        value += Fraction(left - right) * expectation(reduced)
        if left + right == 0 and left != 0:
            central = C * Fraction(left * (left * left - 1), 12)
            value += central * expectation(tuple(modes[:index] + modes[index + 2 :]))
        return value

    # Normal ordered: negative modes kill the bra; positive modes kill the ket.
    if any(mode != 0 for mode in modes):
        return Fraction(0)
    return H ** len(modes)


def gram(left: tuple[int, ...], right: tuple[int, ...]) -> Fraction:
    bra = tuple(reversed(left))
    ket = tuple(-mode for mode in right)
    return expectation(bra + ket)


def quadratic_form(coefficients: tuple[Fraction, ...], matrix: list[list[Fraction]]) -> Fraction:
    return sum(
        coefficients[i] * matrix[i][j] * coefficients[j]
        for i in range(len(coefficients))
        for j in range(len(coefficients))
    )


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def render() -> dict:
    # Level-2 basis L_-2, L_-1^2 and the percolation thermal null vector.
    level2_basis = ((2,), (1, 1))
    level2_gram = [[gram(a, b) for b in level2_basis] for a in level2_basis]
    null_vector = (Fraction(1), Fraction(-2, 3))
    null_norm = quadratic_form(null_vector, level2_gram)

    # Q4 basis follows the repository convention exactly.
    level4_basis = ((2, 2), (3, 1), (4,))
    level4_gram = [[gram(a, b) for b in level4_basis] for a in level4_basis]
    q4 = (Fraction(40), Fraction(-60), Fraction(-9))
    q4_norm = quadratic_form(q4, level4_gram)

    parent_h = H
    parent_hbar = H
    parent_x = parent_h + parent_hbar
    descendant_h = parent_h + 4
    descendant_hbar = parent_hbar
    descendant_x = descendant_h + descendant_hbar
    descendant_spin = descendant_h - descendant_hbar

    return {
        "schema": "matching-one/q4-jordan-inheritance/v1",
        "status": "exact Virasoro algebra; lattice overlap remains conditional",
        "module": {"c": fraction_text(C), "h": fraction_text(H)},
        "level2_null": {
            "basis": ["L_-2", "L_-1^2"],
            "gram": [[fraction_text(x) for x in row] for row in level2_gram],
            "coefficients": [fraction_text(x) for x in null_vector],
            "norm": fraction_text(null_norm),
            "state": "(L_-2 - 2/3 L_-1^2)|h>"
        },
        "q4": {
            "basis": ["L_-2^2", "L_-3 L_-1", "L_-4"],
            "gram": [[fraction_text(x) for x in row] for row in level4_gram],
            "coefficients": [fraction_text(x) for x in q4],
            "norm": fraction_text(q4_norm),
            "state": "(40 L_-2^2 - 60 L_-3 L_-1 - 9 L_-4)|h>"
        },
        "jordan_inheritance": {
            "parent_weights": [fraction_text(parent_h), fraction_text(parent_hbar)],
            "parent_dimension": fraction_text(parent_x),
            "descendant_weights": [fraction_text(descendant_h), fraction_text(descendant_hbar)],
            "descendant_dimension": fraction_text(descendant_x),
            "descendant_spin": fraction_text(descendant_spin),
            "level_commutator": "[L0+Lbar0,Q4]=4 Q4",
            "bottom": "[D-21/4] Q4|eps> = 0",
            "top": "[D-21/4] Q4|eps_tilde> = Q4|eps>",
            "rank": 2
        },
        "he_normalization": {
            "parent_log_coupling_mu": fraction_text(MU),
            "unnormalized_descendant_cross_overlap_factor": fraction_text(MU * q4_norm),
            "warning": "the 4930 factor is Q4-normalization dependent; it is not a lattice amplitude prediction"
        },
        "finite_size": {
            "dimensionless_power_in_N": "N^-13/8",
            "rank2_top_form": "N^-13/8 [A + B log L] = N^-13/8 [A + (B/2) log N]",
            "gaussian_log_increment": "Delta log L = (1/2) log Norm(m)"
        }
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = render()
    text = json.dumps(payload, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    else:
        print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

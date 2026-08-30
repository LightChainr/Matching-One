#!/usr/bin/env python3
"""Thermal-Q4 aspect-ratio oracle for the common P321 projector.

The oracle is conditional on the ordinary (non-logarithmic) thermal
``Q4 epsilon`` bridge.  It does not apply the vacuum KdV operator to the
critical Pinson--Arguin probabilities.  Instead it uses the torus Ward ratio

    <Q4 epsilon>/<epsilon> = (493*pi^4/72) E4(tau)

and removes the overall descendant and lattice normalizations by dividing by
the rectangular-cylinder cusp ``E4(i infinity)=1``.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import mpmath as mp


DEFAULT_RATIOS = ("1", "16/9", "9/4", "4", "9")


@dataclass(frozen=True)
class E4Evaluation:
    value: mp.mpf
    terms: int
    absolute_tail_bound: mp.mpf


def parse_positive_rational(text: str) -> mp.mpf:
    """Parse a positive decimal or ``numerator/denominator`` at mp precision."""

    pieces = text.split("/")
    if len(pieces) == 1:
        result = mp.mpf(pieces[0])
    elif len(pieces) == 2:
        result = mp.mpf(pieces[0]) / mp.mpf(pieces[1])
    else:
        raise ValueError(f"invalid aspect ratio: {text}")
    if result <= 0:
        raise ValueError("aspect ratio must be positive")
    return result


def sigma3(n: int) -> int:
    total = 0
    divisor = 1
    while divisor * divisor <= n:
        if n % divisor == 0:
            quotient = n // divisor
            total += divisor**3
            if quotient != divisor:
                total += quotient**3
        divisor += 1
    return total


def _tail_bound(q: mp.mpf, cutoff: int) -> mp.mpf:
    """Bound ``240 sum_{n>cutoff} sigma_3(n) q^n``."""

    first_n = cutoff + 1
    ratio = q * (mp.mpf(first_n + 1) / first_n) ** 3
    if ratio >= 1:
        return mp.inf
    # sigma_3(n) <= zeta(3) n^3, followed by a geometric ratio bound.
    return 240 * mp.zeta(3) * first_n**3 * q**first_n / (1 - ratio)


def e4_imaginary(
    rho: mp.mpf, *, dps: int = 90, max_terms: int = 1_000_000
) -> E4Evaluation:
    """Evaluate normalized ``E4(i*rho)`` with an explicit tail bound."""

    if rho <= 0:
        raise ValueError("rho must be positive")
    if dps < 30:
        raise ValueError("dps must be at least 30")
    with mp.workdps(dps + 25):
        q = mp.exp(-2 * mp.pi * rho)
        tolerance = mp.power(10, -(dps + 8))
        total = mp.mpf(1)
        for n in range(1, max_terms + 1):
            total += 240 * sigma3(n) * q**n
            bound = _tail_bound(q, n)
            if bound < tolerance:
                return E4Evaluation(+total, n, +bound)
    raise ArithmeticError("E4 series did not reach the requested tail bound")


def _text(value: mp.mpf, digits: int = 70) -> str:
    return mp.nstr(value, digits)


def oracle_payload(
    ratios: Sequence[str] = DEFAULT_RATIOS, *, dps: int = 90
) -> dict[str, object]:
    with mp.workdps(dps + 25):
        parsed = [(label, parse_positive_rational(label)) for label in ratios]
        square = e4_imaginary(mp.mpf(1), dps=dps)
        records = []
        for label, rho in parsed:
            evaluation = e4_imaginary(rho, dps=dps)
            e4 = evaluation.value
            records.append(
                {
                    "rho": label,
                    "rho_decimal": _text(rho),
                    "tau": f"i*({label})",
                    "E4_i_rho": _text(e4),
                    "width_C_over_cylinder_C": _text(e4),
                    "width_C_over_square_C": _text(e4 / square.value),
                    "equal_area_C_over_square_C": _text(
                        rho**2 * e4 / square.value
                    ),
                    "equal_area_descaled_C_over_cylinder_C": _text(e4),
                    "distance_from_cylinder_cusp": _text(e4 - 1),
                    "q_series_terms": evaluation.terms,
                    "absolute_tail_bound": _text(evaluation.absolute_tail_bound),
                }
            )

        holomorphic_prefactor = mp.mpf(493) * mp.pi**4 / 72
        return {
            "schema": "matching-one/p321-thermal-q4-aspect-ratio-oracle/v1",
            "issue": 321,
            "status": "conditional_parameter_free_shape_prediction",
            "projector": "A_top=P2-P0",
            "field_hypothesis": "ordinary thermal Q4 epsilon, h=hbar=5/8",
            "ward_identity": {
                "Q4": "40 L_-2^2 - 60 L_-3 L_-1 - 9 L_-4",
                "holomorphic_ratio": "<Q4 epsilon>/<epsilon>=(493*pi^4/72) E4(tau)",
                "holomorphic_prefactor_decimal": _text(holomorphic_prefactor),
                "real_aligned_ratio": "add the antiholomorphic conjugate; its factor two is absorbed by the real lattice coupling",
            },
            "normalization_free_prediction": {
                "width_law": "p_n(rho)-pc=C_width(rho)n_width^-4+...",
                "ratio": "C_width(rho)/C_width(infinity)=E4(i*rho)",
                "cylinder_cusp": "E4(i*infinity)=1",
                "equal_area_conversion": "if N=rho*n_width^2 and p-pc=C_N(rho)N^-2, then C_width(rho)=C_N(rho)/rho^2",
                "equal_area_descaling": "C_N(rho)/rho^2 has the finite cylinder limit and is the quantity to compare with E4(i*rho)",
            },
            "square_exact_value": {
                "identity": "E4(i)=3*Gamma(1/4)^8/(64*pi^6)",
                "decimal": _text(square.value),
            },
            "records": records,
            "limit_strategy": [
                "At each fixed rho, first take n_width to infinity and extract C_width(rho).",
                "Only then take rho to infinity; E4(i*rho)=1+240 exp(-2*pi*rho)+O(exp(-4*pi*rho)).",
                "For equal-area fits, divide the N^-2 coefficient C_N(rho) by rho^2 before taking the cylinder limit.",
                "rho=9 is a numerical cusp check for this continuum formula, not proof that finite lattice data have reached the TL endpoint.",
            ],
            "operator_boundary": {
                "vacuum_KdV": "Applying D2D0 to critical Pinson-Arguin sector probabilities is Alexander-even and gives exactly zero in P2-P0.",
                "thermal_Q4": "The nonzero candidate is the descendant of the thermal insertion that defines F_t, not a descendant of the vacuum partition function.",
                "unconditional_missing_input": "a proof that the homology projector is Virasoro-transparent for the restricted thermal one-point, excluding defect/contact terms in the descendant Ward recursion",
                "absolute_missing_input": "the lattice coupling ratio g_u/g_t; it cancels from aspect-ratio ratios",
                "logarithmic_loophole": "a top Jordan-partner readout can add a second modulus function and is not covered by the plain E4 curve",
            },
            "sources": [
                "https://arxiv.org/abs/hep-th/0111193",
                "https://arxiv.org/abs/2112.01563",
                "https://arxiv.org/abs/2604.24491",
            ],
        }


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rho", nargs="+", default=list(DEFAULT_RATIOS))
    parser.add_argument("--dps", type=int, default=90)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    rendered = json.dumps(
        oracle_payload(args.rho, dps=args.dps), indent=2, sort_keys=True
    ) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

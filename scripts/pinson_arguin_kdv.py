#!/usr/bin/env python3
"""Vacuum-KdV K4 response of Q=1 Pinson--Arguin primitive sectors.

At Q=1 the random-cluster partition sum is exactly one, so the normalized
primitive-sector probability is also the restricted numerator in this
normalization.  This module applies the chiral modular derivative

    K4 = D2 D0 = (delta-E2/6) delta,  delta=q d/dq,

and derives the reflection-even real response as K4+K4bar = 2 Re(K4).
The analytic Gaussian q-series is checked against a numerical Wirtinger
derivative with tau-bar held fixed.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Iterable, Sequence

import mpmath as mp

from pinson_arguin_primitive import (
    canonical_primitive,
    dedekind_eta_abs_squared,
    primitive_probability_direct,
)


SECTORS = ((1, 0), (0, 1), (1, 1))
CONTRAST_ORDER = ("C_nontrivial_real", "Q_reflection_null", "S_scalar")
CONTRAST = (
    (1.0, -0.5, -0.5),
    (0.0, -math.sqrt(3.0) / 2.0, math.sqrt(3.0) / 2.0),
    (1.0, 1.0, 1.0),
)


def _sigma_power(n: int, power: int) -> int:
    total = 0
    divisor = 1
    while divisor * divisor <= n:
        if n % divisor == 0:
            quotient = n // divisor
            total += divisor**power
            if quotient != divisor:
                total += quotient**power
        divisor += 1
    return total


def eisenstein_e2(tau: mp.mpc, terms: int = 120) -> mp.mpc:
    q = mp.exp(2 * mp.pi * 1j * tau)
    return 1 - 24 * mp.fsum(
        _sigma_power(n, 1) * q**n for n in range(1, terms + 1)
    )


def eisenstein_e4(tau: mp.mpc, terms: int = 120) -> mp.mpc:
    q = mp.exp(2 * mp.pi * 1j * tau)
    return 1 + 240 * mp.fsum(
        _sigma_power(n, 3) * q**n for n in range(1, terms + 1)
    )


def dedekind_eta_holomorphic(
    tau: mp.mpc, *, tolerance: mp.mpf, max_terms: int = 10000
) -> mp.mpc:
    """Holomorphic eta product without a fractional-q branch ambiguity."""

    q = mp.exp(2 * mp.pi * 1j * tau)
    radius = abs(q)
    if radius >= 1:
        raise ValueError("holomorphic eta requires Im(tau)>0")
    product = mp.mpc(1)
    for n in range(1, max_terms + 1):
        product *= 1 - q**n
        if radius ** (n + 1) / (1 - radius) < tolerance:
            return mp.exp(mp.pi * 1j * tau / 12) * product
    raise ArithmeticError("holomorphic eta product did not converge")


def _coefficient(k: int) -> mp.mpf:
    return mp.cos(2 * mp.pi * k / 3) - mp.cos(mp.pi * k)


def primitive_numerator_extended(
    a: int,
    b: int,
    tau: mp.mpc,
    tau_bar: mp.mpc,
    *,
    dps: int = 80,
    gaussian_terms: int = 80,
) -> mp.mpc:
    """Complexified numerator Z_ab(tau,tau_bar) for Wirtinger derivatives."""

    a, b = canonical_primitive((a, b))
    with mp.workdps(dps + 30):
        return +_primitive_numerator_extended_current_precision(
            a, b, mp.mpc(tau), mp.mpc(tau_bar), gaussian_terms
        )


def _primitive_numerator_extended_current_precision(
    a: int, b: int, z: mp.mpc, zbar: mp.mpc, gaussian_terms: int
) -> mp.mpc:
    """Internal form that lets ``mp.diff`` raise working precision itself."""

    height = (z - zbar) / (2j)
    u = a - b * z
    ubar = a - b * zbar
    tolerance = mp.power(10, -(mp.mp.dps - 12))
    eta = dedekind_eta_holomorphic(z, tolerance=tolerance)
    eta_bar = mp.conj(
        dedekind_eta_holomorphic(mp.conj(zbar), tolerance=tolerance)
    )
    gaussian = mp.fsum(
        2
        * _coefficient(k)
        * mp.exp(-2 * mp.pi * k * k * u * ubar / (3 * height))
        for k in range(1, gaussian_terms + 1)
    )
    return mp.sqrt(2 / (3 * height)) * gaussian / (eta * eta_bar)


def primitive_k4_holomorphic_series(
    a: int,
    b: int,
    tau: mp.mpc,
    *,
    dps: int = 80,
    gaussian_terms: int = 80,
    eisenstein_terms: int = 120,
) -> mp.mpc:
    """Analytic termwise q-series for D2 D0 Z_ab."""

    a, b = canonical_primitive((a, b))
    with mp.workdps(dps + 30):
        z = mp.mpc(tau)
        height = mp.im(z)
        if height <= 0:
            raise ValueError("tau must lie in the upper half-plane")
        u = a - b * z
        ubar = mp.conj(u)
        e2 = eisenstein_e2(z, eisenstein_terms)
        e4 = eisenstein_e4(z, eisenstein_terms)
        prefactor = mp.sqrt(2 / (3 * height)) / dedekind_eta_abs_squared(
            z, dps=dps + 5
        )
        response = mp.mpc(0)
        for k in range(1, gaussian_terms + 1):
            term = (
                2
                * prefactor
                * _coefficient(k)
                * mp.exp(-2 * mp.pi * k * k * abs(u) ** 2 / (3 * height))
            )
            # L_k = delta log(term_k), with tau_bar held fixed.
            log_first = (
                1 / (8 * mp.pi * height)
                - e2 / 24
                - k * k * ubar * ubar / (6 * height * height)
            )
            log_second = (
                1 / (32 * mp.pi**2 * height**2)
                + (e4 - e2 * e2) / 288
                - k * k * ubar * ubar / (12 * mp.pi * height**3)
            )
            response += term * (
                log_first * log_first
                + log_second
                - e2 * log_first / 6
            )
        return +response


def primitive_k4_holomorphic_numeric(
    a: int,
    b: int,
    tau: mp.mpc,
    *,
    dps: int = 80,
    gaussian_terms: int = 80,
    eisenstein_terms: int = 120,
) -> mp.mpc:
    """Independent numerical Wirtinger derivative of the complexified sum."""

    with mp.workdps(dps + 30):
        z = mp.mpc(tau)
        zbar = mp.conj(z)

        def numerator(variable: mp.mpc) -> mp.mpc:
            return _primitive_numerator_extended_current_precision(
                a, b, variable, zbar, gaussian_terms
            )

        first = mp.diff(numerator, z)
        second = mp.diff(numerator, z, 2)
        delta_first = first / (2 * mp.pi * 1j)
        delta_second = second / (2 * mp.pi * 1j) ** 2
        return +delta_second - eisenstein_e2(z, eisenstein_terms) * delta_first / 6


def primitive_k4_reflection_even(
    a: int, b: int, tau: mp.mpc, *, dps: int = 80
) -> mp.mpf:
    """Response to a real aligned J4+J4bar coupling: 2 Re[D2D0 Z_ab]."""

    return +mp.mpf(2 * mp.re(primitive_k4_holomorphic_series(a, b, tau, dps=dps)))


def project_contrasts(values: Sequence[mp.mpf]) -> list[mp.mpf]:
    if len(values) != 3:
        raise ValueError("the registered C/Q/S projection requires three sectors")
    return [
        mp.fsum(CONTRAST[row][column] * values[column] for column in range(3))
        for row in range(3)
    ]


def _text(value: mp.mpf | mp.mpc, digits: int = 50) -> str:
    return mp.nstr(value, digits)


def oracle_record(
    identifier: str,
    tau: mp.mpc,
    *,
    n: int | None,
    omega1_length_squared: int | None,
    dps: int,
) -> dict[str, object]:
    series = []
    numeric = []
    probabilities = []
    for a, b in SECTORS:
        probabilities.append(primitive_probability_direct(a, b, tau, dps=dps))
        series.append(primitive_k4_holomorphic_series(a, b, tau, dps=dps))
        numeric.append(primitive_k4_holomorphic_numeric(a, b, tau, dps=dps))
    real_response = [2 * mp.re(value) for value in series]
    contrasts = project_contrasts(real_response)
    return {
        "id": identifier,
        "N": n,
        "omega1_length_squared": omega1_length_squared,
        "tau_real": _text(mp.re(tau)),
        "tau_imag": _text(mp.im(tau)),
        "paper_sector_order": [list(sector) for sector in SECTORS],
        "probability": [_text(value) for value in probabilities],
        "K4_holomorphic_series": [_text(value) for value in series],
        "K4_holomorphic_numeric": [_text(value) for value in numeric],
        "dual_path_absolute_difference": [
            _text(abs(first - second)) for first, second in zip(series, numeric)
        ],
        "reflection_even_K4_plus_K4bar": [_text(value) for value in real_response],
        "contrast_order": list(CONTRAST_ORDER),
        "reflection_even_contrasts": [_text(value) for value in contrasts],
        "finite_size_design_vector": (
            [_text(value / omega1_length_squared) for value in contrasts]
            if omega1_length_squared is not None
            else None
        ),
    }


def artifact_payload(dps: int = 80) -> dict[str, object]:
    with mp.workdps(dps + 30):
        designs = (
            ("square", mp.mpc(0, 1), None, None),
            ("hexagonal", mp.mpc(mp.mpf("0.5"), mp.sqrt(3) / 2), None, None),
            (
                "pell_Dminus2_N30",
                mp.mpc(mp.mpf("0.5"), mp.mpf(5) / 6),
                30,
                36,
            ),
            (
                "pell_Dplus1_N56",
                mp.mpc(mp.mpf("0.5"), mp.mpf(7) / 8),
                56,
                64,
            ),
        )
        return {
            "schema": "matching-one/p231-pinson-arguin-vacuum-kdv-oracle/v1",
            "issue": 231,
            "status": "deterministic_formula_oracle_not_preregistration",
            "operator": "K4=D2D0=(delta-E2/6)delta at c=0",
            "normalization": (
                "Q=1 random-cluster total partition sum is exactly one; "
                "the primitive probability is the restricted numerator"
            ),
            "real_response_derivation": (
                "for a real fixed-label numerator and a real aligned reflection-even "
                "coupling, K4bar is the complex conjugate of K4, hence K4+K4bar=2Re(K4)"
            ),
            "finite_size_law": (
                "delta P_a = g4 * |omega1|^{-2} * (K4+K4bar) Z_a; "
                "N^{-1} differs by the known Im(tau) shape factor"
            ),
            "records": [
                oracle_record(
                    identifier,
                    tau,
                    n=n,
                    omega1_length_squared=scale_squared,
                    dps=dps,
                )
                for identifier, tau, n, scale_squared in designs
            ],
            "interpretation_boundary": [
                "The overall lattice coupling g4 is not fixed by the continuum calculation.",
                "The 2Re projection assumes a real reflection-even coupling aligned with the registered omega1 frame.",
                "A modular or period-basis move must also transport the spin-4 coupling phase.",
                "The N30/N56 use is retrospective and is not new independent evidence.",
            ],
        }


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args(list(argv) if argv is not None else None)
    payload = artifact_payload(dps=args.dps)
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(text, encoding="utf-8")
    else:
        print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

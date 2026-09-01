#!/usr/bin/env python3
"""Critical Q=1 primitive homology-sector probabilities on a complex torus.

The formulas are the percolation specialization of Arguin's restricted
Fortuin--Kasteleyn torus partition functions.  The paper convention labels a
primitive type ``{a,b}`` by the physical cycle ``a*omega_1-b*omega_2``.
Repository winding vectors instead use ``u*omega_1+v*omega_2``; consequently
the exact interface map is ``(u,v) -> {u,-v}`` before canonicalizing the
unoriented primitive line.

This module evaluates continuum critical probabilities.  It does not assert
that a finite lattice probability is equal to its continuum baseline.
"""

from __future__ import annotations

import argparse
import json
from math import gcd
from pathlib import Path
from typing import Iterable, Optional, Sequence, Tuple, Union

import mpmath as mp


Pair = Tuple[int, int]


def _canonical_sign(pair: Pair) -> Pair:
    first, second = pair
    if first < 0 or (first == 0 and second < 0):
        return -first, -second
    return pair


def canonical_primitive(pair: Pair, *, saturate: bool = False) -> Pair:
    """Return the canonical unoriented primitive representative.

    Formula inputs must already be primitive, while engine winding generators
    may carry a nonprimitive multiple and are therefore explicitly saturated.
    """

    first, second = pair
    divisor = gcd(abs(first), abs(second))
    if divisor == 0:
        raise ValueError("the zero winding has no primitive homology type")
    if divisor != 1 and not saturate:
        raise ValueError("Pinson/Arguin sector labels must be coprime")
    return _canonical_sign((first // divisor, second // divisor))


def engine_to_paper(winding: Pair) -> Pair:
    """Map repository ``+omega_1,+omega_2`` winding to paper ``{a,b}``."""

    first, second = canonical_primitive(winding, saturate=True)
    return canonical_primitive((first, -second))


def paper_to_engine(sector: Pair) -> Pair:
    """Map a Pinson/Arguin ``{a,b}`` label to repository winding."""

    first, second = canonical_primitive(sector)
    return canonical_primitive((first, -second), saturate=True)


def _reduce_for_eta(tau: mp.mpc) -> tuple[mp.mpc, mp.mpf]:
    """Reduce ``tau`` and return the factor for ``|eta(original)|^2``.

    Only the absolute square is needed.  The transformations used are
    ``|eta(tau+n)|=|eta(tau)|`` and
    ``|eta(-1/tau)|^2=|tau| |eta(tau)|^2``.
    """

    reduced = mp.mpc(tau)
    factor = mp.mpf(1)
    for _ in range(100):
        nearest = int(mp.floor(mp.re(reduced) + mp.mpf("0.5")))
        reduced -= nearest
        if abs(reduced) >= 1:
            return reduced, factor
        modulus = abs(reduced)
        factor /= modulus
        reduced = -1 / reduced
    raise ArithmeticError("modular reduction of tau did not converge")


def dedekind_eta_abs_squared(
    tau: Union[complex, mp.mpc],
    *,
    dps: int = 80,
) -> mp.mpf:
    """Evaluate ``|eta(tau)|^2`` with an explicit product-tail bound."""

    if dps < 30:
        raise ValueError("dps must be at least 30")
    with mp.workdps(dps + 20):
        value = mp.mpc(tau)
        if mp.im(value) <= 0:
            raise ValueError("tau must lie in the upper half-plane")
        reduced, factor = _reduce_for_eta(value)
        q = mp.exp(2j * mp.pi * reduced)
        radius = abs(q)
        tolerance = mp.power(10, -(dps + 8))
        product = mp.mpc(1)
        for index in range(1, 1_000_001):
            product *= 1 - q**index
            remainder = radius ** (index + 1) / (
                (1 - radius) * (1 - radius ** (index + 1))
            )
            if remainder < tolerance:
                eta = mp.exp(mp.pi * 1j * reduced / 12) * product
                return +mp.mpf(factor * abs(eta) ** 2)
        raise ArithmeticError("Dedekind eta product did not reach tolerance")


def _direct_tail_bound(prefactor: mp.mpf, exponent: mp.mpf, cutoff: int) -> mp.mpf:
    """Bound the omitted two-sided Gaussian sum after ``|k|<=cutoff``."""

    first = mp.exp(-exponent * (cutoff + 1) ** 2)
    ratio = mp.exp(-exponent * (2 * cutoff + 3))
    return 4 * prefactor * first / (1 - ratio)


def primitive_probability_direct(
    a: int,
    b: int,
    tau: Union[complex, mp.mpc],
    *,
    dps: int = 80,
) -> mp.mpf:
    """Evaluate ``pi_tau({a,b})`` by the direct Gaussian winding sum."""

    a, b = canonical_primitive((a, b))
    with mp.workdps(dps + 30):
        value = mp.mpc(tau)
        height = mp.im(value)
        if height <= 0:
            raise ValueError("tau must lie in the upper half-plane")
        eta_squared = dedekind_eta_abs_squared(value, dps=dps + 10)
        cycle = a - b * value
        exponent = 2 * mp.pi * abs(cycle) ** 2 / (3 * height)
        prefactor = mp.sqrt(mp.mpf(2) / (3 * height)) / eta_squared
        tolerance = mp.power(10, -(dps + 8))
        total = mp.mpf(0)
        for cutoff in range(1, 1_000_001):
            coefficient = mp.cos(2 * mp.pi * cutoff / 3) - mp.cos(
                mp.pi * cutoff
            )
            total += 2 * mp.exp(-exponent * cutoff**2) * coefficient
            if _direct_tail_bound(prefactor, exponent, cutoff) < tolerance:
                return +mp.mpf(prefactor * total)
        raise ArithmeticError("direct Gaussian sum did not reach tolerance")


def _theta3_imaginary(t: mp.mpf, tolerance: mp.mpf) -> mp.mpf:
    if t <= 0:
        raise ValueError("theta modulus must be positive")
    if t < 1:
        return _theta3_imaginary(1 / t, tolerance * mp.sqrt(t)) / mp.sqrt(t)
    total = mp.mpf(1)
    for n in range(1, 1_000_001):
        term = 2 * mp.exp(-mp.pi * t * n**2)
        total += term
        if term < tolerance:
            return total
    raise ArithmeticError("theta3 sum did not reach tolerance")


def _theta4_imaginary(t: mp.mpf, tolerance: mp.mpf) -> mp.mpf:
    if t <= 0:
        raise ValueError("theta modulus must be positive")
    if t < 1:
        return _theta2_imaginary(1 / t, tolerance * mp.sqrt(t)) / mp.sqrt(t)
    total = mp.mpf(1)
    for n in range(1, 1_000_001):
        term = 2 * mp.exp(-mp.pi * t * n**2)
        total += -term if n % 2 else term
        if term < tolerance:
            return total
    raise ArithmeticError("theta4 sum did not reach tolerance")


def _theta2_imaginary(t: mp.mpf, tolerance: mp.mpf) -> mp.mpf:
    if t <= 0:
        raise ValueError("theta modulus must be positive")
    if t < 1:
        return _theta4_imaginary(1 / t, tolerance * mp.sqrt(t)) / mp.sqrt(t)
    total = mp.mpf(0)
    for n in range(0, 1_000_001):
        term = 2 * mp.exp(-mp.pi * t * (n + mp.mpf("0.5")) ** 2)
        total += term
        if term < tolerance:
            return total
    raise ArithmeticError("theta2 sum did not reach tolerance")


def primitive_probability_theta(
    a: int,
    b: int,
    tau: Union[complex, mp.mpc],
    *,
    dps: int = 80,
) -> mp.mpf:
    """Independently evaluate the compact Pinson/Arguin theta formula."""

    a, b = canonical_primitive((a, b))
    with mp.workdps(dps + 40):
        value = mp.mpc(tau)
        height = mp.im(value)
        if height <= 0:
            raise ValueError("tau must lie in the upper half-plane")
        cycle_abs = abs(a - b * value)
        scaled_height = height / cycle_abs**2
        tolerance = mp.power(10, -(dps + 15))
        first = _theta3_imaginary(scaled_height / 6, tolerance)
        second_argument = 3 * scaled_height / 2
        second = _theta3_imaginary(second_argument, tolerance)
        third = _theta2_imaginary(second_argument, tolerance)
        eta_squared = dedekind_eta_abs_squared(value, dps=dps + 10)
        return +mp.mpf(
            (first - second - 2 * third) / (2 * cycle_abs * eta_squared)
        )


def _probability_text(value: mp.mpf, digits: int = 50) -> str:
    return mp.nstr(value, digits)


def baseline_records(dps: int = 90) -> list[dict[str, object]]:
    """Return the #156 square, sheared, hexagonal and small-Pell oracles."""

    with mp.workdps(dps + 20):
        rho_height = mp.sqrt(3) / 2
        designs: Sequence[tuple[str, mp.mpf, mp.mpf, str]] = (
            ("square", mp.mpf(0), mp.mpf(1), "tau=i"),
            ("half_sheared", mp.mpf("0.5"), mp.mpf(1), "tau=1/2+i"),
            (
                "hexagonal_limit",
                mp.mpf("0.5"),
                rho_height,
                "tau=1/2+i*sqrt(3)/2",
            ),
            ("pell_Dminus2_N30", mp.mpf("0.5"), mp.mpf(5) / 6, "x=5,m=3"),
            ("pell_Dplus1_N56", mp.mpf("0.5"), mp.mpf(7) / 8, "x=7,m=4"),
        )
        engine_lines: Sequence[tuple[str, Pair]] = (
            ("l0", (1, 0)),
            ("l1", (0, 1)),
            ("l2", (1, -1)),
        )
        records: list[dict[str, object]] = []
        for identifier, real, imag, provenance in designs:
            tau = mp.mpc(real, imag)
            sectors = []
            for line_name, engine_winding in engine_lines:
                paper_type = engine_to_paper(engine_winding)
                probability = primitive_probability_direct(
                    *paper_type, tau, dps=dps
                )
                sectors.append(
                    {
                        "line": line_name,
                        "engine_winding": list(engine_winding),
                        "paper_type": list(paper_type),
                        "probability": _probability_text(probability),
                    }
                )
            records.append(
                {
                    "id": identifier,
                    "tau_real": _probability_text(real),
                    "tau_imag": _probability_text(imag),
                    "provenance": provenance,
                    "sectors": sectors,
                }
            )
        return records


def artifact_payload(dps: int = 90) -> dict[str, object]:
    return {
        "schema": "pinson-arguin-primitive-baseline-v1",
        "issue": 156,
        "claim_level": "continuum_exact_formula_numerical_evaluation",
        "model": "critical_Q1_FK_percolation",
        "sources": [
            {
                "author": "Louis-Pierre Arguin",
                "arxiv": "hep-th/0111193v2",
                "doi": "10.1023/A:1019979326380",
            },
            {
                "author": "Alexi Morin-Duchesne and Yvan Saint-Aubin",
                "arxiv": "0812.2925v2",
                "doi": "10.1103/PhysRevE.80.021130",
            },
        ],
        "convention": {
            "torus": "C/(Z+tau Z), Im(tau)>0",
            "paper_type": "{a,b} is the unoriented cycle a*omega1-b*omega2",
            "engine_winding": "(u,v) is u*omega1+v*omega2",
            "engine_to_paper": "(u,v)->{u,-v}, followed by primitive saturation and sign canonicalization",
            "event": "the full configuration homology image is exactly the primitive rank-1 subgroup; rank-0 and rank-2 cross are excluded",
        },
        "numerics": {
            "dps": dps,
            "direct_sum": "adaptive two-sided Gaussian sum with an explicit tail bound",
            "eta": "modularly reduced adaptive Dedekind product with a log-product tail bound",
            "independent_check": "compact theta formula with independent imaginary-modulus theta sums",
        },
        "records": baseline_records(dps=dps),
        "interpretation_boundary": [
            "These are critical continuum probabilities, not exact finite-lattice probabilities.",
            "They apply to site and bond percolation through Q=1 universality in the isotropic physical metric.",
            "They do not supply a primal/matching joint-channel identity or an off-critical baseline.",
        ],
    }


def main(argv: Optional[Iterable[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dps", type=int, default=90)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args(list(argv) if argv is not None else None)
    payload = artifact_payload(dps=args.dps)
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.json is None:
        print(text, end="")
    else:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

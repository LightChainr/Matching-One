#!/usr/bin/env python3
"""Freeze norm-five chiral Hecke phases and a tiny exact Z5 channel oracle."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from math import atan2, degrees, gcd
from pathlib import Path
import sys
from typing import Sequence

sys.path.insert(0, str(Path(__file__).resolve().parent))

from integer_period_torus import classify_configuration, integer_torus_geometry  # noqa: E402


Q = Fraction
Cyclo = tuple[Q, Q, Q, Q]  # coefficients in 1,zeta,zeta^2,zeta^3
ZERO: Cyclo = (Q(0), Q(0), Q(0), Q(0))
ONE: Cyclo = (Q(1), Q(0), Q(0), Q(0))
P = Q(2, 5)
SPINS = (4, 8, 12)


def zeta_power(power: int) -> Cyclo:
    basis = (
        (Q(1), Q(0), Q(0), Q(0)),
        (Q(0), Q(1), Q(0), Q(0)),
        (Q(0), Q(0), Q(1), Q(0)),
        (Q(0), Q(0), Q(0), Q(1)),
        (Q(-1), Q(-1), Q(-1), Q(-1)),
    )
    return basis[power % 5]


def cadd(left: Cyclo, right: Cyclo) -> Cyclo:
    return tuple(a + b for a, b in zip(left, right))  # type: ignore[return-value]


def cscale(value: Cyclo, scalar: Q | int) -> Cyclo:
    return tuple(Q(scalar) * item for item in value)  # type: ignore[return-value]


def cmul(left: Cyclo, right: Cyclo) -> Cyclo:
    result = ZERO
    for first, left_value in enumerate(left):
        for second, right_value in enumerate(right):
            result = cadd(
                result,
                cscale(zeta_power(first + second), left_value * right_value),
            )
    return result


def cconjugate(value: Cyclo) -> Cyclo:
    result = ZERO
    for power, coefficient in enumerate(value):
        result = cadd(result, cscale(zeta_power(-power), coefficient))
    return result


def ccomplex(value: Cyclo) -> complex:
    import cmath

    root = cmath.exp(2j * cmath.pi / 5)
    return sum(float(coefficient) * root**power for power, coefficient in enumerate(value))


def cyclo_record(value: Cyclo) -> dict:
    approximate = ccomplex(value)
    return {
        "basis": "1,zeta5,zeta5^2,zeta5^3",
        "coefficients": [str(item) for item in value],
        "decimal": [approximate.real, approximate.imag],
    }


def gaussian_ratio_power(spin: int) -> tuple[int, int, int]:
    """Return exact real numerator, imag numerator, denominator of q^spin."""

    real, imag = 1, 0
    for _ in range(spin):
        real, imag = 3 * real - 4 * imag, 4 * real + 3 * imag
    denominator = 5**spin
    common = gcd(gcd(abs(real), abs(imag)), denominator)
    return real // common, imag // common, denominator // common


def phase_record(spin: int) -> dict:
    real, imag, denominator = gaussian_ratio_power(spin)
    value = complex(real / denominator, imag / denominator)
    sign = "+" if imag >= 0 else "-"
    return {
        "spin": spin,
        "exact": {
            "real": f"{real}/{denominator}",
            "imag": f"{imag}/{denominator}",
            "gaussian_rational": f"({real}{sign}{abs(imag)}i)/{denominator}",
        },
        "decimal": [value.real, value.imag],
        "phase_degrees_principal": degrees(atan2(value.imag, value.real)),
        "modulus_squared_exact": "1",
    }


def geometry(a: int, b: int):
    return integer_torus_geometry(((a, -b), (b, a)), name=f"gaussian-{a}-{b}")


def labels_to_active(geom, a: int, b: int, labels: Sequence[bool]) -> list[bool]:
    n = a * a + b * b
    return [labels[(a * x + b * y) % n] for x, y in geom.coordinates]


def matching_odd_cross(geom, a: int, b: int, labels: Sequence[bool]) -> int:
    active = labels_to_active(geom, a, b, labels)
    primal, _ = classify_configuration(geom, active)
    matching, _ = classify_configuration(
        geom, [not value for value in active], matching=True
    )
    return int(primal.cross) - int(matching.cross)


def pivotal_delta(
    geom, a: int, b: int, labels: Sequence[bool], site: int
) -> int:
    absent = list(labels)
    absent[site] = False
    low = matching_odd_cross(geom, a, b, absent)
    present = list(labels)
    present[site] = True
    high = matching_odd_cross(geom, a, b, present)
    return high - low


def marked_row(geom, a: int, b: int, labels: Sequence[bool]) -> Cyclo:
    result = ZERO
    for site in range(5):
        result = cadd(
            result,
            cscale(zeta_power(-site), pivotal_delta(geom, a, b, labels, site)),
        )
    return result


def score_numerator(labels: Sequence[bool]) -> Cyclo:
    result = ZERO
    for site, state in enumerate(labels):
        result = cadd(result, cscale(zeta_power(site), int(state) - P))
    return result


def exact_marked_response(a: int, b: int) -> tuple[Cyclo, Cyclo]:
    """Return score response and an independent probability-polynomial derivative."""

    geom = geometry(a, b)
    response = ZERO
    direct_derivative = ZERO
    q = 1 - P
    for mask in range(1 << 5):
        labels = [bool(mask & (1 << site)) for site in range(5)]
        occupied = sum(labels)
        weight = P**occupied * q ** (5 - occupied)
        row = marked_row(geom, a, b, labels)
        response = cadd(
            response,
            cscale(cmul(row, score_numerator(labels)), weight / (P * q)),
        )

        # Coefficient of epsilon in prod_k [p+eps*zeta^k]^X
        # [q-eps*zeta^k]^(1-X), derived independently of the score formula.
        polynomial: list[Cyclo] = [ONE]
        for site, state in enumerate(labels):
            constant = P if state else q
            slope = zeta_power(site)
            if not state:
                slope = cscale(slope, -1)
            updated = [ZERO] * (len(polynomial) + 1)
            for degree, coefficient in enumerate(polynomial):
                updated[degree] = cadd(updated[degree], cscale(coefficient, constant))
                updated[degree + 1] = cadd(
                    updated[degree + 1], cmul(coefficient, slope)
                )
            polynomial = updated
        direct_derivative = cadd(direct_derivative, cmul(row, polynomial[1]))
    return response, direct_derivative


def build_artifact() -> dict:
    phases = [phase_record(spin) for spin in SPINS]
    exact_pairs = {
        (entry["exact"]["real"], entry["exact"]["imag"]) for entry in phases
    }
    if len(exact_pairs) != len(SPINS):
        raise AssertionError("the handed pair does not separate the frozen spins")
    pairwise_phase_separation = {}
    for first_index, first in enumerate(phases):
        first_value = complex(*first["decimal"])
        for second in phases[first_index + 1 :]:
            second_value = complex(*second["decimal"])
            separation = abs(
                degrees(atan2((first_value / second_value).imag, (first_value / second_value).real))
            )
            pairwise_phase_separation[f"H{first['spin']}_vs_H{second['spin']}"] = separation

    tiny = {}
    responses = []
    for label, (a, b) in {"2+i": (2, 1), "2-i": (2, -1)}.items():
        response, derivative = exact_marked_response(a, b)
        if response != derivative or response == ZERO:
            raise AssertionError(f"tiny marked response failed for {label}")
        tiny[label] = {
            "response": cyclo_record(response),
            "direct_symbolic_derivative": cyclo_record(derivative),
        }
        responses.append(response)
    if responses[1] != cconjugate(responses[0]):
        raise AssertionError("reflection transport failed in the tiny oracle")

    return {
        "schema": "matching-one.norm5-chiral-hecke-phase.v1",
        "issues": [226, 244],
        "status": "exact_prediction_plus_tiny_channel_oracle",
        "handed_pair": {
            "multipliers": ["2+i", "2-i"],
            "ratio": "q=(2+i)/(2-i)=(3+4i)/5",
            "reflection": "(x,y)->(x,-y)",
            "deck_transport": "k_plus -> k_minus",
            "character_transport": "chi_plus(k)=zeta5^k -> conjugate chi_minus(k)=zeta5^-k",
        },
        "response_design": {
            "score": "S_chi=sum_(j,k) zeta5^k (X_(j,k)-p)",
            "opposite_character_row": "O_chibar=sum_(j,k) zeta5^-k Delta_(j,k)[matching-odd H4/local-pivotal readout]",
            "matrix_element": "R_m=E[O_chibar S_chi]/[p(1-p)]",
            "handed_ratio": "R_(2+i)/R_(2-i) after exact reflection transport",
            "covariance": "estimate Re/Im of both handed responses in one shared-randomness batch vector",
        },
        "hecke_eigenfield_predictions": {
            "formula": "A_s(2+i)/A_s(2-i)=((2+i)/(2-i))^s",
            "values": phases,
            "pairwise_principal_phase_separation_degrees": pairwise_phase_separation,
            "minimum_separation_degrees": min(pairwise_phase_separation.values()),
            "all_three_exactly_distinct": True,
            "single_handed_pair_discriminates": "yes, conditional on one-field dominance and the frozen reflection/character transport",
            "unknown_common_complex_normalization": "cancels in the handed ratio",
        },
        "tiny_exact_oracle": {
            "parent_order": 1,
            "child_order": 5,
            "p": str(P),
            "enumeration": "all 2^5 configurations for each hand",
            "observable": "Z5-opposite-character pivotal projection of matching-odd cross wrap",
            "hands": tiny,
            "reflection_conjugacy_exact": True,
            "scientific_scope": "validates charged channel and transport only; N=5 does not test the Hecke spin phase",
        },
        "interpretation_boundary": {
            "result": "H4, H8 and H12 give three distinct unit phases for one handed norm-5 pair",
            "mechanism_hypothesis": "a stable measured phase near one target identifies the dominant spin-s Hecke eigenfield",
            "not_implied": "the tiny real oracle response is not evidence for any of H4/H8/H12",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    artifact = build_artifact()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    print(json.dumps(artifact, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Exact p-biased Boolean noise-semigroup gate for Issue #227.

The oracle keeps only level-summed Fourier spectra.  It independently checks
them by enumerating noisy configuration pairs, and relates the rho=1
autocorrelation derivative to exact pivotal mass on N=5 and N=10 tori.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction
import json
from pathlib import Path
from typing import Iterable, Sequence

from c4_self_matching_exact import c4_self_matching_torus
from integer_period_torus import (
    IntegerTorusGeometry,
    classify_configuration,
    gaussian_integer_torus,
)


P = Fraction(2, 5)
RHO_GRID = (Fraction(0), Fraction(1, 4), Fraction(1, 2), Fraction(3, 4), Fraction(1))
INDICATORS = (
    "primal_cross",
    "primal_direction0",
    "primal_direction1",
    "matching_complement_cross",
    "matching_complement_direction0",
    "matching_complement_direction1",
)
DERIVED = ("orientation_difference", "matching_odd_cross")
OBSERVABLES = INDICATORS + DERIVED
CROSS_PAIRS = (
    ("primal_direction0", "primal_direction1", "orientation_cross_spectrum"),
    ("primal_cross", "matching_complement_cross", "matching_cross_spectrum"),
    ("orientation_difference", "matching_odd_cross", "H4_matching_odd_cross_spectrum"),
)


def popcount(value: int) -> int:
    """Return the Hamming weight on every supported Python (3.9+)."""

    return bin(value).count("1")


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def active_from_mask(mask: int, n: int) -> tuple[bool, ...]:
    return tuple(bool(mask & (1 << vertex)) for vertex in range(n))


def truth_tables(geometry: IntegerTorusGeometry) -> dict[str, list[int]]:
    tables = {name: [] for name in OBSERVABLES}
    full = (1 << geometry.n) - 1
    for mask in range(1 << geometry.n):
        active = active_from_mask(mask, geometry.n)
        complement = active_from_mask(full ^ mask, geometry.n)
        primal, _ = classify_configuration(geometry, active, matching=False)
        matching, _ = classify_configuration(geometry, complement, matching=True)
        row = {
            "primal_cross": int(primal.cross),
            "primal_direction0": int(primal.direction_0),
            "primal_direction1": int(primal.direction_1),
            "matching_complement_cross": int(matching.cross),
            "matching_complement_direction0": int(matching.direction_0),
            "matching_complement_direction1": int(matching.direction_1),
        }
        row["orientation_difference"] = (
            row["primal_direction0"] - row["primal_direction1"]
        )
        row["matching_odd_cross"] = (
            row["primal_cross"] - row["matching_complement_cross"]
        )
        for name in OBSERVABLES:
            tables[name].append(row[name])
    return tables


def mask_probability(mask: int, n: int, p: Fraction) -> Fraction:
    occupied = popcount(mask)
    return p**occupied * (1 - p) ** (n - occupied)


def biased_fourier_moments(
    tables: dict[str, list[int]], n: int, p: Fraction
) -> tuple[dict[str, list[Fraction]], dict[str, Fraction]]:
    """Return unnormalized p-biased moments E[f prod_i(X_i-p)]."""

    moments = {name: [Fraction(0) for _ in range(1 << n)] for name in OBSERVABLES}
    means = {name: Fraction(0) for name in OBSERVABLES}
    q = 1 - p
    for mask in range(1 << n):
        probability = mask_probability(mask, n, p)
        subset_products = [Fraction(1)] * (1 << n)
        for subset in range(1, 1 << n):
            low = subset & -subset
            vertex = low.bit_length() - 1
            centered_bit = q if mask & low else -p
            subset_products[subset] = subset_products[subset ^ low] * centered_bit
        for name in OBSERVABLES:
            value = tables[name][mask]
            if not value:
                continue
            means[name] += probability * value
            for subset, product in enumerate(subset_products):
                moments[name][subset] += probability * value * product
    # Centered observables have zero empty-set coefficient.  Nonempty moments
    # are unchanged because every nontrivial basis function has zero mean.
    for name in OBSERVABLES:
        moments[name][0] = Fraction(0)
    return moments, means


def level_cross_spectrum(
    first: Sequence[Fraction],
    second: Sequence[Fraction],
    n: int,
    p: Fraction,
) -> list[Fraction]:
    pq = p * (1 - p)
    levels = [Fraction(0) for _ in range(n + 1)]
    for subset in range(1, 1 << n):
        level = popcount(subset)
        levels[level] += first[subset] * second[subset] / pq**level
    return levels


def evaluate_generating(coefficients: Sequence[Fraction], rho: Fraction) -> Fraction:
    return sum(value * rho**level for level, value in enumerate(coefficients))


def derivative_at_one(coefficients: Sequence[Fraction]) -> Fraction:
    return sum(level * value for level, value in enumerate(coefficients))


def direct_noisy_pair_histograms(
    tables: dict[str, list[int]], n: int
) -> dict[tuple[str, str], dict[tuple[int, int, int], int]]:
    """Enumerate x,y pairs and retain integer sufficient statistics."""

    pairs = [(name, name) for name in OBSERVABLES]
    pairs.extend((first, second) for first, second, _ in CROSS_PAIRS)
    histograms = {pair: defaultdict(int) for pair in pairs}
    for first_mask in range(1 << n):
        for second_mask in range(1 << n):
            both_one = popcount(first_mask & second_mask)
            different = popcount(first_mask ^ second_mask)
            both_zero = n - both_one - different
            key = (both_one, both_zero, different)
            for first, second in pairs:
                product = tables[first][first_mask] * tables[second][second_mask]
                if product:
                    histograms[(first, second)][key] += product
    return {pair: dict(histogram) for pair, histogram in histograms.items()}


def direct_noisy_covariance(
    histogram: dict[tuple[int, int, int], int],
    mean_first: Fraction,
    mean_second: Fraction,
    p: Fraction,
    rho: Fraction,
) -> Fraction:
    q = 1 - p
    pq = p * q
    joint_11 = p * p + pq * rho
    joint_00 = q * q + pq * rho
    joint_different = pq * (1 - rho)
    expectation = sum(
        count
        * joint_11**both_one
        * joint_00**both_zero
        * joint_different**different
        for (both_one, both_zero, different), count in histogram.items()
    )
    return expectation - mean_first * mean_second


def pivotal_mass_exact(table: Sequence[int], n: int, p: Fraction) -> Fraction:
    """Return the unsigned total pivotal mass of a Boolean indicator."""

    total = Fraction(0)
    for vertex in range(n):
        low = 1 << vertex
        for mask in range(1 << n):
            if mask & low:
                continue
            difference = table[mask | low] - table[mask]
            if difference not in (-1, 0, 1):
                raise ArithmeticError("indicator edge difference is not Boolean")
            occupied = popcount(mask)
            # mask has vertex fixed to zero, hence n-1 free coordinates.
            total += difference * difference * p**occupied * (1 - p) ** (
                n - 1 - occupied
            )
    return total


def _spectrum_payload(coefficients: Sequence[Fraction]) -> dict:
    return {
        "level_order": list(range(len(coefficients))),
        "coefficients": [fraction_text(value) for value in coefficients],
        "generating_function": "sum_k coefficient[k]*rho^k",
        "value_at_rho_1": fraction_text(sum(coefficients)),
        "derivative_at_rho_1": fraction_text(derivative_at_one(coefficients)),
    }


def analyze_geometry(geometry: IntegerTorusGeometry, p: Fraction) -> dict:
    n = geometry.n
    if n > 10:
        raise ValueError("exact noise oracle is limited to N<=10")
    tables = truth_tables(geometry)
    moments, means = biased_fourier_moments(tables, n, p)
    auto = {
        name: level_cross_spectrum(moments[name], moments[name], n, p)
        for name in OBSERVABLES
    }
    cross = {
        label: level_cross_spectrum(moments[first], moments[second], n, p)
        for first, second, label in CROSS_PAIRS
    }
    histograms = direct_noisy_pair_histograms(tables, n)

    direct_checks = {}
    pair_specs = [(name, name, f"auto:{name}") for name in OBSERVABLES]
    pair_specs.extend(
        (first, second, f"cross:{label}")
        for first, second, label in CROSS_PAIRS
    )
    for first, second, label in pair_specs:
        coefficients = auto[first] if first == second else cross[label.split(":", 1)[1]]
        values = []
        for rho in RHO_GRID:
            direct = direct_noisy_covariance(
                histograms[(first, second)], means[first], means[second], p, rho
            )
            fourier = evaluate_generating(coefficients, rho)
            if direct != fourier:
                raise ArithmeticError(f"direct/Fourier noise mismatch for {label}")
            values.append(
                {
                    "rho": fraction_text(rho),
                    "direct_noisy_pair_covariance": fraction_text(direct),
                    "fourier_generating_value": fraction_text(fourier),
                    "difference": "0",
                }
            )
        direct_checks[label] = values

    pivotal = {}
    pq = p * (1 - p)
    for name in INDICATORS:
        mass = pivotal_mass_exact(tables[name], n, p)
        derivative = derivative_at_one(auto[name])
        if derivative != pq * mass:
            raise ArithmeticError(f"rho derivative/pivotal mismatch for {name}")
        pivotal[name] = {
            "total_unsigned_pivotal_mass": fraction_text(mass),
            "raw_autocorrelation_derivative": fraction_text(derivative),
            "derivative_divided_by_p_times_q": fraction_text(derivative / pq),
            "exact_identity": "C_f'(1)=p*(1-p)*sum_i P(i pivotal)",
        }

    return {
        "geometry": geometry.name,
        "N": n,
        "configurations": 1 << n,
        "p": fraction_text(p),
        "observable_definitions": {
            "orientation_difference": "primal_direction0-primal_direction1",
            "matching_odd_cross": "primal_cross-matching_complement_cross",
            "matching_complement": "matching-lattice wrapping of the white complement, represented as a decreasing function of the black field",
        },
        "means": {name: fraction_text(value) for name, value in means.items()},
        "autocorrelation_generating_functions": {
            name: _spectrum_payload(values) for name, values in auto.items()
        },
        "cross_spectra": {
            label: {
                "observables": [first, second],
                **_spectrum_payload(cross[label]),
            }
            for first, second, label in CROSS_PAIRS
        },
        "direct_noisy_pair_cross_checks": {
            "rho_grid": [fraction_text(value) for value in RHO_GRID],
            "all_exact": True,
            "pairs": direct_checks,
        },
        "rho_one_pivotal_checks": pivotal,
    }


def render() -> dict:
    geometries = (
        gaussian_integer_torus(2, 1),
        c4_self_matching_torus(3, 1),
    )
    return {
        "schema": "matching-one.exact-boolean-noise-semigroup.v1",
        "issue": 227,
        "claim_level": "exact_finite_Boolean_identity",
        "p_biased_basis": "chi_S=prod_{i in S}(X_i-p)/sqrt(p*(1-p))",
        "noise_operator": "keep each bit with probability rho; otherwise resample Bernoulli(p)",
        "geometries": [analyze_geometry(geometry, P) for geometry in geometries],
        "interpretation_boundary": [
            "The rho curve is one generating function; grid points are not independent evidence.",
            "Krawtchouk occupation-count modes are not the spatial Boolean Fourier spectrum.",
            "The rho=1 derivative equals p*(1-p) times pivotal mass before normalization.",
            "Tiny exact spectral organization does not imply a large-N exponent or continuum spectral sample.",
            "No Jordan or logarithmic mechanism is inferred from these finite identities.",
        ],
    }


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(list(argv) if argv is not None else None)
    text = json.dumps(render(), indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(text, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

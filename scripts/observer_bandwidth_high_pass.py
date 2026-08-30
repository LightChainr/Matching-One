#!/usr/bin/env python3
"""Exact finite-noise high-pass filter for the p-biased Walsh clock."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Mapping

from observer_bandwidth_k_centered_euler import euler_observer_values
from observer_bandwidth_product_walsh import (
    apply_noise,
    centered_basis,
    fraction,
    walsh_degree_coefficients,
)


def popcount(mask: int) -> int:
    return bin(mask).count("1")


def filter_coefficients(rho: Fraction, degree_cutoff: int) -> list[Fraction]:
    coefficients = [Fraction(1)]
    for degree in range(degree_cutoff + 1):
        factor = rho ** (-degree)
        updated = [Fraction(0)] * (len(coefficients) + 1)
        for index, coefficient in enumerate(coefficients):
            updated[index] += coefficient
            updated[index + 1] -= factor * coefficient
        coefficients = updated
    return coefficients


def spectral_multiplier(rho: Fraction, degree_cutoff: int, walsh_degree: int) -> Fraction:
    return product(
        Fraction(1) - rho ** (walsh_degree - removed_degree)
        for removed_degree in range(degree_cutoff + 1)
    )


def product(values) -> Fraction:
    result = Fraction(1)
    for value in values:
        result *= value
    return result


def apply_high_pass(
    values: Mapping[int, Fraction],
    n: int,
    p: Fraction,
    rho: Fraction,
    degree_cutoff: int,
) -> dict[int, Fraction]:
    coefficients = filter_coefficients(rho, degree_cutoff)
    levels = [apply_noise(values, n, p, rho**power) for power in range(len(coefficients))]
    return {
        mask: sum(
            coefficient * level[mask]
            for coefficient, level in zip(coefficients, levels)
        )
        for mask in range(1 << n)
    }


def raw_monomial_values(n: int, subset: int) -> dict[int, Fraction]:
    return {
        mask: Fraction(1 if mask & subset == subset else 0) for mask in range(1 << n)
    }


def build_report(manifest: Mapping[str, object]) -> dict[str, object]:
    rho = fraction(manifest["rho"])
    cutoff = int(manifest["degree_cutoff"])
    n = int(manifest["basis_cube_n"])
    p = fraction(manifest["basis_cube_p"])
    coefficients = filter_coefficients(rho, cutoff)
    if len(coefficients) != cutoff + 2:
        raise AssertionError("high-pass representation does not use d+2 noise levels")
    multipliers = {
        degree: spectral_multiplier(rho, cutoff, degree)
        for degree in range(cutoff + 9)
    }
    if any(multipliers[degree] for degree in range(cutoff + 1)):
        raise AssertionError("high-pass filter did not annihilate a low Walsh degree")
    if any(not 0 < multipliers[degree] < 1 for degree in range(cutoff + 1, cutoff + 9)):
        raise AssertionError("high-pass multiplier is not a positive contraction above cutoff")

    monomial_checks = 0
    for subset in range(1 << n):
        if popcount(subset) > cutoff:
            continue
        filtered = apply_high_pass(raw_monomial_values(n, subset), n, p, rho, cutoff)
        if any(filtered.values()):
            raise AssertionError("degree-at-most-d raw monomial survived high-pass")
        monomial_checks += len(filtered)

    degree_five_subset = (1 << n) - 1
    degree_five = {
        mask: centered_basis(mask, degree_five_subset, p) for mask in range(1 << n)
    }
    filtered_five = apply_high_pass(degree_five, n, p, rho, cutoff)
    expected_multiplier = multipliers[cutoff + 1]
    if any(
        filtered_five[mask] != expected_multiplier * degree_five[mask]
        for mask in range(1 << n)
    ):
        raise AssertionError("degree-five control has the wrong high-pass multiplier")
    if not any(filtered_five.values()):
        raise AssertionError("degree-five control vanished")

    euler_values, incidence = euler_observer_values(int(manifest["euler_torus_L"]))
    euler_n = incidence["sites"]
    euler_energy = walsh_degree_coefficients(euler_values, euler_values, euler_n, p)
    if any(euler_energy[degree] for degree in range(cutoff + 1, euler_n + 1)):
        raise AssertionError("Euler observer contains Walsh degree above four")
    filtered_euler_energy = sum(
        multipliers[degree] ** 2 * energy
        for degree, energy in euler_energy.items()
    )
    if filtered_euler_energy != 0:
        raise AssertionError("Euler positive control survived the exact high-pass")

    return {
        "schema": manifest["schema"],
        "status": "exact_high_pass_controls_verified",
        "rho": str(rho),
        "degree_cutoff": cutoff,
        "noise_levels": [str(rho**power) for power in range(len(coefficients))],
        "linear_combination_coefficients": [str(value) for value in coefficients],
        "coefficient_l1_norm": str(sum(abs(value) for value in coefficients)),
        "multipliers": {str(degree): str(value) for degree, value in multipliers.items()},
        "annihilated_raw_monomial_point_checks": monomial_checks,
        "degree_five_control": {
            "multiplier": str(expected_multiplier),
            "nonzero_output_points": sum(bool(value) for value in filtered_five.values()),
        },
        "euler_control": {
            "L": int(manifest["euler_torus_L"]),
            "incidence": incidence,
            "degree_energies": {str(degree): str(value) for degree, value in euler_energy.items()},
            "filtered_l2_energy": str(filtered_euler_energy),
        },
        "boundary": manifest["boundary"],
    }


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=root / "analysis/observer_bandwidth_high_pass_manifest.json",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = build_report(json.loads(args.manifest.read_text(encoding="utf-8")))
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

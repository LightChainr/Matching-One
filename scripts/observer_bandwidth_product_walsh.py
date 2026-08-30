#!/usr/bin/env python3
"""Exact p-biased Walsh bandwidth oracle for a finite Bernoulli cube."""

from __future__ import annotations

import argparse
from fractions import Fraction
import itertools
import json
from pathlib import Path
from typing import Mapping


def fraction(value: object) -> Fraction:
    return Fraction(str(value))


def product_weight(mask: int, n: int, p: Fraction) -> Fraction:
    ones = mask.bit_count()
    return p**ones * (1 - p) ** (n - ones)


def centered_basis(mask: int, subset: int, p: Fraction) -> Fraction:
    value = Fraction(1)
    bit = 0
    remaining = subset
    while remaining:
        if remaining & 1:
            value *= Fraction((mask >> bit) & 1) - p
        remaining >>= 1
        bit += 1
    return value


def expectation(values: Mapping[int, Fraction], n: int, p: Fraction) -> Fraction:
    return sum(
        product_weight(mask, n, p) * values[mask] for mask in range(1 << n)
    )


def center(values: Mapping[int, Fraction], n: int, p: Fraction) -> dict[int, Fraction]:
    mean = expectation(values, n, p)
    return {mask: values[mask] - mean for mask in range(1 << n)}


def noise_transition(
    source: int, target: int, n: int, p: Fraction, rho: Fraction
) -> Fraction:
    probability = Fraction(1)
    for site in range(n):
        x = (source >> site) & 1
        y = (target >> site) & 1
        resample = p if y else 1 - p
        probability *= (rho if x == y else 0) + (1 - rho) * resample
    return probability


def apply_noise(
    values: Mapping[int, Fraction], n: int, p: Fraction, rho: Fraction
) -> dict[int, Fraction]:
    return {
        source: sum(
            noise_transition(source, target, n, p, rho) * values[target]
            for target in range(1 << n)
        )
        for source in range(1 << n)
    }


def direct_covariance(
    observer: Mapping[int, Fraction],
    source: Mapping[int, Fraction],
    n: int,
    p: Fraction,
    rho: Fraction,
) -> Fraction:
    observer_centered = center(observer, n, p)
    source_centered = center(source, n, p)
    noised_source = apply_noise(source_centered, n, p, rho)
    return expectation(
        {
            mask: observer_centered[mask] * noised_source[mask]
            for mask in range(1 << n)
        },
        n,
        p,
    )


def walsh_degree_coefficients(
    observer: Mapping[int, Fraction],
    source: Mapping[int, Fraction],
    n: int,
    p: Fraction,
) -> dict[int, Fraction]:
    observer_centered = center(observer, n, p)
    source_centered = center(source, n, p)
    coefficients = {degree: Fraction(0) for degree in range(1, n + 1)}
    for subset in range(1, 1 << n):
        degree = subset.bit_count()
        observer_projection = expectation(
            {
                mask: observer_centered[mask] * centered_basis(mask, subset, p)
                for mask in range(1 << n)
            },
            n,
            p,
        )
        source_projection = expectation(
            {
                mask: source_centered[mask] * centered_basis(mask, subset, p)
                for mask in range(1 << n)
            },
            n,
            p,
        )
        norm = (p * (1 - p)) ** degree
        coefficients[degree] += observer_projection * source_projection / norm
    return coefficients


def evaluate_multilinear(
    terms: Mapping[int, Fraction], n: int
) -> dict[int, Fraction]:
    return {
        mask: sum(
            coefficient for subset, coefficient in terms.items() if mask & subset == subset
        )
        for mask in range(1 << n)
    }


def source_fixture(n: int) -> dict[int, Fraction]:
    return {
        mask: Fraction((mask * mask + 3 * mask.bit_count() + 7) % 17 - 8)
        for mask in range(1 << n)
    }


def subset_mask(indices: list[int]) -> int:
    return sum(1 << index for index in indices)


def verify_basis_eigenvalues(n: int, p: Fraction, rhos: list[Fraction]) -> int:
    checks = 0
    for rho, subset in itertools.product(rhos, range(1 << n)):
        values = {
            mask: centered_basis(mask, subset, p) for mask in range(1 << n)
        }
        transformed = apply_noise(values, n, p, rho)
        eigenvalue = rho ** subset.bit_count()
        for mask in range(1 << n):
            if transformed[mask] != eigenvalue * values[mask]:
                raise AssertionError("p-biased centered basis is not a noise eigenfunction")
            checks += 1
    return checks


def build_report(manifest: Mapping[str, object]) -> dict[str, object]:
    n = int(manifest["n"])
    p = fraction(manifest["p"])
    cutoff = int(manifest["degree_cutoff"])
    rhos = [fraction(value) for value in manifest["rho_values"]]
    terms = {
        subset_mask(list(row["subset"])): fraction(row["coefficient"])
        for row in manifest["observer_terms"]
    }
    if max(mask.bit_count() for mask in terms) > cutoff:
        raise ValueError("observer fixture exceeds declared degree cutoff")
    observer = evaluate_multilinear(terms, n)
    source = source_fixture(n)
    coefficients = walsh_degree_coefficients(observer, source, n, p)
    if any(coefficients[degree] for degree in range(cutoff + 1, n + 1)):
        raise AssertionError("observer has covariance above its degree cutoff")
    comparisons = []
    for rho in rhos:
        direct = direct_covariance(observer, source, n, p, rho)
        spectral = sum(rho**degree * value for degree, value in coefficients.items())
        if direct != spectral:
            raise AssertionError("direct and Walsh covariance disagree")
        comparisons.append(
            {
                "rho": str(rho),
                "direct_covariance": str(direct),
                "spectral_covariance": str(spectral),
                "residual": str(direct - spectral),
            }
        )
    return {
        "schema": manifest["schema"],
        "status": "exact_identity_verified",
        "n": n,
        "p": str(p),
        "degree_cutoff": cutoff,
        "basis_point_checks": verify_basis_eigenvalues(n, p, rhos),
        "degree_coefficients": {
            str(degree): str(value) for degree, value in coefficients.items()
        },
        "active_nonconstant_degrees": [
            degree for degree, value in coefficients.items() if value
        ],
        "comparisons": comparisons,
        "boundary": manifest["boundary"],
    }


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=root / "analysis/observer_bandwidth_product_walsh_manifest.json",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    report = build_report(manifest)
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

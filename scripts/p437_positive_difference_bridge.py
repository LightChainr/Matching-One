#!/usr/bin/env python3
"""Exact Hermite--Genocchi bridge and minimal fifth-difference kernel (#437)."""
from __future__ import annotations
import argparse
from fractions import Fraction
import json
from math import comb, factorial
from pathlib import Path

NODES = tuple(Fraction(1, 2 ** level) for level in range(6))
H5 = Fraction(9765, 32768)


def product(values):
    result = Fraction(1)
    for value in values:
        result *= value
    return result


def divided_difference_coefficients():
    return [H5 / product(node - other for j, other in enumerate(NODES) if j != i)
            for i, node in enumerate(NODES)]


def complete_homogeneous(degree):
    values = [Fraction(1)] + [Fraction(0)] * degree
    for node in NODES:
        for order in range(1, degree + 1):
            values[order] += node * values[order - 1]
    return values[degree]


def dirichlet_mixture_moment(degree):
    # T=sum W_l t_l, W~Dirichlet(1,...,1), six coordinates.
    return Fraction(factorial(degree) * factorial(5), factorial(degree + 5)) * complete_homogeneous(degree)


def spectral_multiplier(degree):
    return product(1 - Fraction(2) ** (r - degree) for r in range(5))


def mixed_derivative(mask, sites, observable):
    """Rademacher D_S=2^-|S| alternating vertex sum; never resamples S."""
    sites = tuple(sites)
    if len(sites) != 5 or len(set(sites)) != 5 or min(sites) < 0:
        raise ValueError("exactly five distinct nonnegative bond indices required")
    base = mask & ~sum(1 << site for site in sites)
    value = 0
    for subset in range(32):
        vertex = base | sum(1 << site for j, site in enumerate(sites) if (subset >> j) & 1)
        sign = (-1) ** (5 - bin(subset).count("1"))
        value += sign * observable(vertex)
    return value / 32


def paired_difference(mask_a, mask_b, sites, observable, proposal_probability):
    """Unbiased *signed* integrand if T/noise/proposal follow the note.

    Caller supplies q(S|X,Y,T)>0. This routine does not silently invent a
    pivotal proposal or assume that zero-probability omitted sets vanish.
    """
    if not 0 < proposal_probability <= 1:
        raise ValueError("proposal probability must be known and positive")
    left = mixed_derivative(mask_a, sites, observable)
    right = mixed_derivative(mask_b, sites, observable)
    return float(H5) * left.conjugate() * right / proposal_probability


def certificate():
    degrees = []
    for j in range(21):
        hg = (H5 * Fraction(factorial(j), factorial(j - 5) * factorial(5)) *
              dirichlet_mixture_moment(j - 5)) if j >= 5 else Fraction(0)
        exact = spectral_multiplier(j)
        if hg != exact:
            raise AssertionError("Hermite--Genocchi multiplier identity failed")
        degrees.append({"degree": j, "highpass_multiplier": str(exact), "HG_multiplier": str(hg)})
    count = comb(224, 5)
    return {"schema": "matching-one/p437-positive-difference-bridge/v1",
            "noise_nodes": [str(node) for node in NODES], "h5": str(H5),
            "divided_difference_coefficients": [str(value) for value in divided_difference_coefficients()],
            "multiplier_identity": degrees,
            "sparse_degree5_uniform_subset_obstruction": {
                "number_of_five_subsets": count,
                "nonzero_probability": str(Fraction(1, count)),
                "nonzero_response": str(H5 * count),
                "mean": str(H5), "variance": str(H5 * H5 * (count - 1)),
                "variance_float": float(H5 * H5 * (count - 1)),
                "relative_variance": count - 1,
                "conclusion": "pointwise annihilation does not make uniform five-set sampling efficient"},
            "pair_product_boundary": "degree6 mode: D_SF=x6, so product=-1 when outside bits differ; common-ancestor expectation is a square, individual pair is not",
            "next_minimal_gate": "known nonuniform inclusion probabilities, full support or proved zero omissions, then same estimand variance pilot; no new MC here"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(certificate(), indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()

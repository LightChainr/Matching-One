#!/usr/bin/env python3
"""Exact Galois-group certificates for the finite matching polynomials reached so far.

For selected exact square-site matching polynomials we combine:

1. an irreducibility-over-Q certificate (hence transitive Galois group);
2. a squarefree finite-field factorization with cycle type containing one
   2-cycle and one large prime cycle;
3. an elementary permutation-group argument.

The modular factorization is verified here without CAS factorization: the
prespecified monic factors must multiply to the monic reduction of the source
polynomial and each factor is checked irreducible by Rabin's criterion.

Group argument used by the generated certificate:

- Dedekind/Frobenius gives an element with the recorded squarefree factor-degree
  cycle type.
- A suitable odd power kills the other odd cycles and leaves a transposition.
- A suitable power kills the other cycles and leaves the large prime cycle r.
- Since r exceeds both the size and number of blocks for every nontrivial block
  system of degree n, that r-cycle rules out imprimitivity.
- In a primitive permutation group, the conjugates of one transposition form a
  connected graph; their edge transpositions generate S_n. Hence G=S_n.

This proves full symmetric Galois groups for the listed finite polynomials. It
makes no statement about the infinite-volume threshold.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Iterable

from certify_axis_matching_irreducibility import (
    ROOT,
    CERTIFICATE_PRIMES as AXIS_IRREDUCIBILITY_PRIMES,
    load_axis_coefficients,
    monic,
    poly_mul,
    rabin_irreducible,
)
from certify_diamond_matching_irreducibility import (
    CERTIFICATE_PRIMES as DIAMOND_IRREDUCIBILITY_PRIMES,
    load_diamond_coefficients,
)


CASES = [
    {
        "geometry": "axis",
        "L": 3,
        "factor_prime": 13,
        "factors": [
            [10, 8, 1],
            [1, 7, 6, 8, 6, 10, 7, 1],
        ],
        "large_prime_cycle": 7,
    },
    {
        "geometry": "axis",
        "L": 4,
        "factor_prime": 331,
        "factors": [
            [86, 1],
            [41, 211, 1],
            [201, 265, 256, 33, 121, 57, 63, 292, 194, 71, 180, 223, 139, 1],
        ],
        "large_prime_cycle": 13,
    },
    {
        "geometry": "axis",
        "L": 5,
        "factor_prime": 863,
        "factors": [
            [765, 412, 1],
            [568, 204, 247, 336, 420, 774, 317, 72, 668, 606, 158, 413,
             59, 665, 507, 502, 162, 695, 627, 431, 28, 732, 401, 1],
        ],
        "large_prime_cycle": 23,
    },
    {
        "geometry": "diamond",
        "L": 2,
        "factor_prime": 127,
        "factors": [
            [121, 1],
            [72, 64, 1],
            [61, 90, 99, 108, 69, 1],
        ],
        "large_prime_cycle": 5,
    },
    {
        "geometry": "diamond",
        "L": 3,
        "factor_prime": 241,
        "factors": [
            [111, 131, 1],
            [211, 237, 102, 1],
            [152, 212, 67, 59, 22, 66, 160, 232, 19, 168, 114, 144, 240, 1],
        ],
        "large_prime_cycle": 13,
    },
    {
        "geometry": "gaussian-3-1",
        "L": 0,
        "factor_prime": 13,
        "factors": [
            [3, 1],
            [2, 4, 1],
            [1, 2, 5, 1, 4, 0, 9, 1],
        ],
        "large_prime_cycle": 7,
    },
]


def product_mod(factors: Iterable[list[int]], prime: int) -> list[int]:
    product = [1]
    for factor in factors:
        product = poly_mul(product, factor, prime)
    return monic(product, prime)


def proper_block_max(n: int) -> int:
    maxima = []
    for block_size in range(2, n):
        if n % block_size == 0:
            block_count = n // block_size
            if block_count > 1:
                maxima.append(max(block_size, block_count))
    return max(maxima, default=1)


def lcm(values: Iterable[int]) -> int:
    result = 1
    for value in values:
        result = math.lcm(result, int(value))
    return result


def _gaussian_3_1_coefficients() -> list[int]:
    payload = json.loads(
        (
            ROOT
            / "results"
            / "exact-axis-l5-frontier"
            / "gaussian_3_1_target.json"
        ).read_text(encoding="utf-8")
    )
    geometry = payload["geometry"]
    if geometry["a"] != 3 or geometry["b"] != 1 or geometry["N"] != 10:
        raise RuntimeError("unexpected Gaussian target metadata")
    return [int(value) for value in payload["power_coefficients_ascending"]]


def source_coefficients(geometry: str, L: int) -> list[int]:
    if geometry == "axis":
        return load_axis_coefficients()[L]
    if geometry == "diamond":
        return load_diamond_coefficients()[L]
    if geometry == "gaussian-3-1":
        return _gaussian_3_1_coefficients()
    raise ValueError(geometry)


def irreducibility_prime(geometry: str, L: int) -> int:
    if geometry == "axis":
        return AXIS_IRREDUCIBILITY_PRIMES[L]
    if geometry == "diamond":
        return DIAMOND_IRREDUCIBILITY_PRIMES[L]
    if geometry == "gaussian-3-1":
        return 31
    raise ValueError(geometry)


def certify_case(case: dict[str, object]) -> dict[str, object]:
    geometry = str(case["geometry"])
    L = int(case["L"])
    q = int(case["factor_prime"])
    factors = [[int(v) for v in factor] for factor in case["factors"]]
    coefficients = source_coefficients(geometry, L)
    n = len(coefficients) - 1

    irreducibility = rabin_irreducible(coefficients, irreducibility_prime(geometry, L))
    if not irreducibility["irreducible"]:
        raise AssertionError(f"missing Q-irreducibility prerequisite: {geometry} L={L}")

    source_reduction = monic(coefficients, q)
    factor_product = product_mod(factors, q)
    if source_reduction != factor_product:
        raise AssertionError(f"factor product mismatch: {geometry} L={L} mod {q}")

    factor_rows = []
    for factor in factors:
        result = rabin_irreducible(factor, q)
        if not result["irreducible"]:
            raise AssertionError(f"certificate factor reducible: {geometry} L={L} mod {q}")
        factor_rows.append(
            {
                "degree": len(factor) - 1,
                "monic_factor_ascending": factor,
                "irreducible_mod_prime": True,
            }
        )

    degrees = sorted(row["degree"] for row in factor_rows)
    if degrees.count(2) != 1:
        raise AssertionError("cycle type must contain exactly one 2-cycle")
    if any(degree % 2 == 0 for degree in degrees if degree != 2):
        raise AssertionError("non-2 cycle lengths must be odd so an odd power isolates transposition")

    long_cycle = int(case["large_prime_cycle"])
    if long_cycle not in degrees:
        raise AssertionError("declared large prime cycle absent")

    transposition_power = lcm(degree for degree in degrees if degree != 2)
    if transposition_power % 2 != 1:
        raise AssertionError("transposition-isolating power must be odd")

    other_for_long = [degree for degree in degrees if degree != long_cycle]
    long_cycle_power = lcm(other_for_long)
    if math.gcd(long_cycle_power, long_cycle) != 1:
        raise AssertionError("long-cycle isolating power destroys the long cycle")

    block_bound = proper_block_max(n)
    if long_cycle <= block_bound:
        raise AssertionError(
            f"large cycle r={long_cycle} does not rule out all block systems; bound={block_bound}"
        )

    return {
        "geometry": geometry,
        "L": L,
        "N": n,
        "integer_degree": n,
        "irreducible_over_Q": True,
        "irreducibility_prime": irreducibility_prime(geometry, L),
        "transitive_galois_action": True,
        "factorization_prime": q,
        "factorization_squarefree": True,
        "factor_degrees": degrees,
        "factors": factor_rows,
        "dedekind_cycle_type": degrees,
        "transposition_power": transposition_power,
        "large_prime_cycle": long_cycle,
        "large_cycle_isolation_power": long_cycle_power,
        "maximum_nontrivial_block_size_or_count": block_bound,
        "primitive": True,
        "contains_transposition": True,
        "galois_group": f"S_{n}",
        "not_solvable_by_radicals": n >= 5,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()

    certificates = [certify_case(case) for case in CASES]
    payload = {
        "schema": "matching-one/finite-matching-galois-certificate/v1",
        "argument": (
            "Q-irreducible => transitive; squarefree modular cycle type supplies a "
            "transposition and a large prime cycle; large cycle excludes every "
            "proper block system; primitive + transposition => full symmetric group."
        ),
        "certificates": certificates,
    }
    for row in certificates:
        print(
            f"{row['geometry']} L={row['L']} degree={row['integer_degree']} "
            f"mod={row['factorization_prime']} cycle_type={row['factor_degrees']} "
            f"G={row['galois_group']}"
        )
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

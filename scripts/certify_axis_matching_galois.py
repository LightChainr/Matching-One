#!/usr/bin/env python3
"""Frobenius certificates for Galois groups of committed axis matching polynomials.

The exact power-basis polynomials are those recorded in
``results/exact_small_matching_polynomials.md``.  This module does **not**
read the unmerged axis L=5 frontier (PR #84).  Certificates use only
Dedekind–Frobenius factorizations at primes that do not divide the
leading coefficient and at which the reduction is square-free.

Locked conclusions:

- L=2 (degree 4, even in ``p^2``): Gal = C4, not S4.  The cyclic group is
  the ``p -> -p`` lattice artifact of a biquadratic, not a closed-form
  mechanism for ``p_c``.
- L=3 (degree 9): primitive (contains an 8-cycle) and contains a
  transposition, hence Gal = S9.
- L=4 (degree 16): primitive (contains a 15-cycle) and contains a
  transposition, hence Gal = S16.

A full ``S_n`` finite-cell Galois group excludes a radicals/solvable-group
expression for that finite physical root.  It does **not** speak to
algebraicity or transcendence of the infinite-volume threshold.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

ROOT = Path(__file__).resolve().parents[1]

# Ascending integer coefficients from results/exact_small_matching_polynomials.md.
AXIS_POWER_POLYNOMIALS: Dict[int, List[int]] = {
    2: [-1, 0, 4, 0, -2],
    3: [-1, 0, 0, 6, 0, 0, 0, -18, 18, -4],
    4: [-1, 0, 0, 0, 8, 0, 32, -64, 172, -704, 1104, -608, -56, 128, 16, -32, 6],
}

# Dedekind primes used in the group-theory deductions.  Each reduction is
# square-free and does not drop degree.
FROBENIUS_CERTIFICATES: Dict[int, Dict[str, object]] = {
    2: {
        "irreducible_prime": 3,
        "irreducible_partition": (4,),
        "group": "C4",
        "reason": "biquadratic_splitting_field_degree_4_with_4_cycles",
    },
    3: {
        "irreducible_prime": 5,
        "irreducible_partition": (9,),
        "primitive_prime": 23,
        "primitive_partition": (8, 1),
        "transposition_prime": 11,
        "transposition_partition": (5, 2, 1, 1),
        "transposition_power": 5,
        "group": "S9",
        "reason": "primitive_plus_transposition",
    },
    4: {
        "irreducible_prime": 5,
        "irreducible_partition": (16,),
        "primitive_prime": 19,
        "primitive_partition": (15, 1),
        "transposition_prime": 31,
        "transposition_partition": (11, 3, 2),
        "transposition_power": 33,
        "group": "S16",
        "reason": "primitive_plus_transposition",
    },
}


def degree(coefficients: Sequence[int]) -> int:
    poly = list(coefficients)
    while poly and poly[-1] == 0:
        poly.pop()
    if not poly:
        raise ValueError("zero polynomial has no degree")
    return len(poly) - 1


def mod_p(coefficients: Sequence[int], prime: int) -> List[int]:
    reduced = [int(value) % prime for value in coefficients]
    while reduced and reduced[-1] == 0:
        reduced.pop()
    return reduced or [0]


def make_monic(coefficients: Sequence[int], prime: int) -> List[int]:
    poly = mod_p(coefficients, prime)
    if poly == [0]:
        return [0]
    scale = pow(poly[-1], prime - 2, prime)
    return [(value * scale) % prime for value in poly]


def add_mod(left: Sequence[int], right: Sequence[int], prime: int) -> List[int]:
    width = max(len(left), len(right))
    total = [0] * width
    for index, value in enumerate(left):
        total[index] = (total[index] + value) % prime
    for index, value in enumerate(right):
        total[index] = (total[index] + value) % prime
    while total and total[-1] == 0:
        total.pop()
    return total or [0]


def sub_mod(left: Sequence[int], right: Sequence[int], prime: int) -> List[int]:
    return add_mod(left, [(prime - value) % prime for value in right], prime)


def mul_mod(left: Sequence[int], right: Sequence[int], prime: int) -> List[int]:
    if left == [0] or right == [0]:
        return [0]
    total = [0] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            total[i + j] = (total[i + j] + x * y) % prime
    while total and total[-1] == 0:
        total.pop()
    return total or [0]


def divmod_mod(
    numerator: Sequence[int], denominator: Sequence[int], prime: int
) -> Tuple[List[int], List[int]]:
    denom = list(denominator)
    while denom and denom[-1] == 0:
        denom.pop()
    if not denom or denom == [0]:
        raise ZeroDivisionError("polynomial division by zero")
    inv = pow(denom[-1], prime - 2, prime)
    remainder = list(numerator)
    while remainder and remainder[-1] == 0:
        remainder.pop()
    quotient = [0] * max(1, len(remainder) - len(denom) + 1)
    while remainder and len(remainder) >= len(denom):
        coef = (remainder[-1] * inv) % prime
        shift = len(remainder) - len(denom)
        if shift >= len(quotient):
            quotient.extend([0] * (shift - len(quotient) + 1))
        quotient[shift] = coef
        for index, value in enumerate(denom):
            remainder[index + shift] = (remainder[index + shift] - coef * value) % prime
        while remainder and remainder[-1] == 0:
            remainder.pop()
    while quotient and quotient[-1] == 0:
        quotient.pop()
    return quotient or [0], remainder or [0]


def rem_mod(numerator: Sequence[int], denominator: Sequence[int], prime: int) -> List[int]:
    return divmod_mod(numerator, denominator, prime)[1]


def gcd_mod(left: Sequence[int], right: Sequence[int], prime: int) -> List[int]:
    a, b = list(left), list(right)
    while b != [0]:
        a, b = b, rem_mod(a, b, prime)
    if a == [0]:
        return [0]
    return make_monic(a, prime)


def pow_mod(base: Sequence[int], exponent: int, modulus: Sequence[int], prime: int) -> List[int]:
    result = [1]
    power = rem_mod(base, modulus, prime)
    while exponent:
        if exponent & 1:
            result = rem_mod(mul_mod(result, power, prime), modulus, prime)
        power = rem_mod(mul_mod(power, power, prime), modulus, prime)
        exponent >>= 1
    return result


def derivative_mod(coefficients: Sequence[int], prime: int) -> List[int]:
    if len(coefficients) <= 1:
        return [0]
    deriv = [((index + 1) * coefficients[index + 1]) % prime for index in range(len(coefficients) - 1)]
    while deriv and deriv[-1] == 0:
        deriv.pop()
    return deriv or [0]


def is_squarefree_mod(coefficients: Sequence[int], prime: int) -> bool:
    poly = make_monic(coefficients, prime)
    return gcd_mod(poly, derivative_mod(poly, prime), prime) == [1]


def factor_degree_partition(coefficients: Sequence[int], prime: int) -> Tuple[int, ...]:
    """Return the Frobenius cycle type via distinct-degree factorization."""

    reduced = mod_p(coefficients, prime)
    if reduced == [0] or len(reduced) - 1 != degree(coefficients):
        raise ValueError(f"p={prime} divides the leading coefficient")
    if not is_squarefree_mod(reduced, prime):
        raise ValueError(f"p={prime} ramifies (reduction is not square-free)")
    poly = make_monic(reduced, prime)
    parts: List[int] = []
    remaining = poly
    value = [0, 1]  # x
    index = 1
    while remaining != [1] and 2 * index <= degree(remaining):
        value = pow_mod(value, prime, remaining, prime)
        factor = gcd_mod(sub_mod(value, [0, 1], prime), remaining, prime)
        if factor != [1]:
            block = degree(factor)
            if block % index != 0:
                raise AssertionError("distinct-degree factor degree is not a multiple of i")
            parts.extend([index] * (block // index))
            remaining, _ = divmod_mod(remaining, factor, prime)
            remaining = make_monic(remaining, prime)
            value = rem_mod(value, remaining, prime)
        index += 1
    if remaining != [1]:
        parts.append(degree(remaining))
    parts.sort(reverse=True)
    if sum(parts) != degree(coefficients):
        raise AssertionError("factor-degree partition does not sum to the polynomial degree")
    return tuple(parts)


def q_primitive(coefficients: Sequence[int]) -> List[Fraction]:
    poly = [Fraction(value) for value in coefficients]
    while poly and poly[-1] == 0:
        poly.pop()
    if not poly:
        return [Fraction(0)]
    leading = poly[-1]
    return [value / leading for value in poly]


def q_mod(numerator: Sequence[Fraction], denominator: Sequence[Fraction]) -> List[Fraction]:
    remainder = [Fraction(value) for value in numerator]
    denom = q_primitive(denominator)
    while remainder and len(remainder) >= len(denom):
        coef = remainder[-1]
        shift = len(remainder) - len(denom)
        for index, value in enumerate(denom):
            remainder[index + shift] -= coef * value
        while remainder and remainder[-1] == 0:
            remainder.pop()
    return remainder or [Fraction(0)]


def q_gcd(left: Sequence[int], right: Sequence[int]) -> List[Fraction]:
    a, b = q_primitive(left), q_primitive(right)
    while b != [Fraction(0)]:
        a, b = b, q_primitive(q_mod(a, b))
    return a


def pairwise_gcds() -> Dict[str, str]:
    sizes = sorted(AXIS_POWER_POLYNOMIALS)
    rows = {}
    for i, left in enumerate(sizes):
        for right in sizes[i + 1 :]:
            gcd = q_gcd(AXIS_POWER_POLYNOMIALS[left], AXIS_POWER_POLYNOMIALS[right])
            rows[f"L{left}_L{right}"] = "1" if gcd == [Fraction(1)] else str(gcd)
    return rows


def l2_splitting_field_degree() -> Dict[str, object]:
    """C4 certificate from the biquadratic tower, independent of Frobenius."""

    poly = AXIS_POWER_POLYNOMIALS[2]
    # f(p) = -2 p^4 + 4 p^2 - 1.  Even: g(u)= -2 u^2 + 4 u - 1, u=p^2.
    # disc(g) = 16 - 8 = 8, not a square, so Q(u)=Q(sqrt(2)).
    disc_g = 4 * 4 - 4 * (-2) * (-1)
    return {
        "even_in_p2": poly[1] == 0 and poly[3] == 0,
        "quadratic_in_u_discriminant": disc_g,
        "quadratic_in_u_disc_is_square": disc_g > 0 and int(disc_g ** 0.5) ** 2 == disc_g,
        "splitting_field_degree": 4,
        "transitive": True,
        "has_4_cycles": True,
        "group": "C4",
    }


def certify_size(L: int) -> Dict[str, object]:
    poly = AXIS_POWER_POLYNOMIALS[L]
    spec = FROBENIUS_CERTIFICATES[L]
    partitions = {}
    irreducible_p = int(spec["irreducible_prime"])
    partitions[irreducible_p] = factor_degree_partition(poly, irreducible_p)
    if partitions[irreducible_p] != tuple(spec["irreducible_partition"]):
        raise AssertionError(f"L={L}: irreducible Frobenius partition drifted")
    payload: Dict[str, object] = {
        "L": L,
        "N": L * L,
        "degree": degree(poly),
        "coefficients_ascending": list(poly),
        "leading_coefficient": poly[-1],
        "irreducible_over_Q": True,
        "irreducible_witness": {
            "p": irreducible_p,
            "partition": list(partitions[irreducible_p]),
        },
        "group": spec["group"],
        "reason": spec["reason"],
        "not_a_statement_about_infinite_pc": True,
    }
    if L == 2:
        payload["splitting_field"] = l2_splitting_field_degree()
        payload["not_Sn"] = True
        return payload
    primitive_p = int(spec["primitive_prime"])
    transposition_p = int(spec["transposition_prime"])
    partitions[primitive_p] = factor_degree_partition(poly, primitive_p)
    partitions[transposition_p] = factor_degree_partition(poly, transposition_p)
    if partitions[primitive_p] != tuple(spec["primitive_partition"]):
        raise AssertionError(f"L={L}: (n-1)-cycle partition drifted")
    if partitions[transposition_p] != tuple(spec["transposition_partition"]):
        raise AssertionError(f"L={L}: transposition-bearing partition drifted")
    payload["primitive_witness"] = {
        "p": primitive_p,
        "partition": list(partitions[primitive_p]),
        "criterion": "transitive_plus_(n-1)_cycle_is_primitive",
    }
    payload["transposition_witness"] = {
        "p": transposition_p,
        "partition": list(partitions[transposition_p]),
        "power": spec["transposition_power"],
        "criterion": "odd_power_of_the_2_cycle_factor_is_a_transposition",
    }
    payload["jordan_criterion"] = "primitive_permutation_group_containing_a_transposition_is_Sn"
    return payload


def run_suite() -> Dict[str, object]:
    certificates = [certify_size(L) for L in (2, 3, 4)]
    return {
        "schema": "axis matching Galois certificates v1",
        "excludes_L5": True,
        "excludes_PR84": True,
        "pairwise_gcds": pairwise_gcds(),
        "certificates": certificates,
        "interpretation": (
            "Finite-cell physical roots at L=3,4 are not solvable by radicals. "
            "This is silent about the infinite-volume threshold."
        ),
    }


def render_report(payload: Dict[str, object]) -> str:
    lines = [
        "# Axis matching-polynomial Galois certificates",
        "",
        "Source: `scripts/certify_axis_matching_galois.py`.",
        "Claim level: C5 for L=2 `C4`; C5 for L=3 `S9` and L=4 `S16` via Dedekind–Frobenius",
        "plus the Jordan transposition criterion. Silent about infinite-volume `p_c`.",
        "Axis L=5 / PR #84 is deliberately excluded.",
        "",
        "## Pairwise gcds over Q",
        "",
        "```text",
    ]
    for key, value in payload["pairwise_gcds"].items():
        lines.append(f"{key}: {value}")
    lines.extend(
        [
            "```",
            "",
            "## Certificates",
            "",
            "| L | degree | group | irreducible p | primitive p | transposition p |",
            "|---:|---:|---|---:|---:|---:|",
        ]
    )
    for row in payload["certificates"]:
        primitive = row.get("primitive_witness", {}).get("p", "—")
        trans = row.get("transposition_witness", {}).get("p", "—")
        irr = row["irreducible_witness"]["p"]
        lines.append(
            f"| {row['L']} | {row['degree']} | {row['group']} | {irr} | {primitive} | {trans} |"
        )
    lines.extend(
        [
            "",
            "## What this does not establish",
            "",
            "- algebraicity or transcendence of square-site `p_c`;",
            "- Galois groups of diamond polynomials or of axis L=5;",
            "- a finite-cell closed form for the infinite threshold.",
            "",
        ]
    )
    return "\n".join(lines)


def json_ready(value: object) -> object:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    return value


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    payload = run_suite()
    for row in payload["certificates"]:
        print(f"L={row['L']} deg={row['degree']} Gal={row['group']}")
    print("pairwise gcds", payload["pairwise_gcds"])
    if args.report is not None:
        args.report.write_text(render_report(payload), encoding="utf-8")
        print("wrote " + str(args.report))
    if args.json is not None:
        args.json.write_text(json.dumps(json_ready(payload), indent=2) + "\n", encoding="utf-8")
        print("wrote " + str(args.json))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

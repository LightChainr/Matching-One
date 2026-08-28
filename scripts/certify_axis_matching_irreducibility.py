#!/usr/bin/env python3
"""Dependency-free irreducibility certificates for exact axis matching polynomials.

For a primitive integer polynomial f, an irreducible same-degree reduction
modulo one prime proves f irreducible over Q (Gauss lemma).  This script loads
the exact axis L=2..5 coefficient arrays already committed in the exact-zero
and L=5 frontier results, then verifies prespecified finite-field certificates
with Rabin's irreducibility criterion.

No CAS factorization is used by the verifier.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE_PRIMES = {2: 3, 3: 5, 4: 5, 5: 19}


def trim(poly: list[int]) -> list[int]:
    out = poly[:]
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def inv_mod(value: int, prime: int) -> int:
    return pow(value % prime, prime - 2, prime)


def monic(poly: Iterable[int], prime: int) -> list[int]:
    out = trim([int(value) % prime for value in poly])
    if out == [0]:
        raise ValueError("zero polynomial has no monic normalization")
    scale = inv_mod(out[-1], prime)
    return trim([(value * scale) % prime for value in out])


def poly_add(a: list[int], b: list[int], prime: int) -> list[int]:
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(n):
        out[i] = ((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)) % prime
    return trim(out)


def poly_sub(a: list[int], b: list[int], prime: int) -> list[int]:
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(n):
        out[i] = ((a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)) % prime
    return trim(out)


def poly_mul(a: list[int], b: list[int], prime: int) -> list[int]:
    out = [0] * (len(a) + len(b) - 1)
    for i, av in enumerate(a):
        for j, bv in enumerate(b):
            out[i + j] = (out[i + j] + av * bv) % prime
    return trim(out)


def poly_divmod(a: list[int], b: list[int], prime: int) -> tuple[list[int], list[int]]:
    numerator = trim([value % prime for value in a])
    denominator = trim([value % prime for value in b])
    if denominator == [0]:
        raise ZeroDivisionError("polynomial division by zero")
    if len(numerator) < len(denominator):
        return [0], numerator

    quotient = [0] * (len(numerator) - len(denominator) + 1)
    inv_lead = inv_mod(denominator[-1], prime)
    remainder = numerator[:]
    while remainder != [0] and len(remainder) >= len(denominator):
        shift = len(remainder) - len(denominator)
        coefficient = remainder[-1] * inv_lead % prime
        quotient[shift] = coefficient
        for j, value in enumerate(denominator):
            remainder[shift + j] = (remainder[shift + j] - coefficient * value) % prime
        remainder = trim(remainder)
    return trim(quotient), remainder


def poly_mod(a: list[int], modulus: list[int], prime: int) -> list[int]:
    return poly_divmod(a, modulus, prime)[1]


def poly_mul_mod(a: list[int], b: list[int], modulus: list[int], prime: int) -> list[int]:
    return poly_mod(poly_mul(a, b, prime), modulus, prime)


def poly_pow_mod(base: list[int], exponent: int, modulus: list[int], prime: int) -> list[int]:
    result = [1]
    factor = poly_mod(base, modulus, prime)
    power = int(exponent)
    while power:
        if power & 1:
            result = poly_mul_mod(result, factor, modulus, prime)
        factor = poly_mul_mod(factor, factor, modulus, prime)
        power >>= 1
    return result


def poly_gcd(a: list[int], b: list[int], prime: int) -> list[int]:
    left = trim([value % prime for value in a])
    right = trim([value % prime for value in b])
    while right != [0]:
        left, right = right, poly_divmod(left, right, prime)[1]
    return monic(left, prime)


def prime_divisors(n: int) -> list[int]:
    remaining = n
    divisors: list[int] = []
    candidate = 2
    while candidate * candidate <= remaining:
        if remaining % candidate == 0:
            divisors.append(candidate)
            while remaining % candidate == 0:
                remaining //= candidate
        candidate += 1
    if remaining > 1:
        divisors.append(remaining)
    return divisors


def frobenius_x(power: int, modulus: list[int], prime: int) -> list[int]:
    """Return x^(prime^power) mod modulus without constructing the huge exponent."""

    value = [0, 1]
    for _ in range(power):
        value = poly_pow_mod(value, prime, modulus, prime)
    return value


def rabin_irreducible(integer_coefficients: list[int], prime: int) -> dict[str, object]:
    reduced = monic(integer_coefficients, prime)
    degree = len(reduced) - 1
    if degree <= 0:
        return {"irreducible": False, "reason": "degree <= 0"}

    x = [0, 1]
    full = frobenius_x(degree, reduced, prime)
    full_condition = poly_sub(full, x, prime) == [0]

    gcd_checks = []
    for divisor in prime_divisors(degree):
        exponent_level = degree // divisor
        probe = poly_sub(frobenius_x(exponent_level, reduced, prime), x, prime)
        gcd_value = poly_gcd(reduced, probe, prime)
        gcd_checks.append(
            {
                "degree_divisor": divisor,
                "frobenius_level": exponent_level,
                "gcd": gcd_value,
                "gcd_is_one": gcd_value == [1],
            }
        )

    return {
        "irreducible": full_condition and all(row["gcd_is_one"] for row in gcd_checks),
        "degree": degree,
        "prime": prime,
        "monic_reduction_ascending": reduced,
        "x_qn_equals_x": full_condition,
        "gcd_checks": gcd_checks,
    }


def load_axis_coefficients() -> dict[int, list[int]]:
    pilot_path = ROOT / "results" / "exact-zero-map-pilot" / "zero_map.json"
    frontier_path = ROOT / "results" / "exact-axis-l5-frontier" / "raw" / "axis_L5_cpp.json"
    pilot = json.loads(pilot_path.read_text(encoding="utf-8"))
    frontier = json.loads(frontier_path.read_text(encoding="utf-8"))

    coefficients: dict[int, list[int]] = {}
    for row in pilot["polynomials"]:
        if row["geometry"] == "axis" and row["L"] in (2, 3, 4):
            coefficients[int(row["L"])] = [int(value) for value in row["power_coefficients_ascending"]]
    coefficients[5] = [int(value) for value in frontier["power_coefficients_ascending"]]
    if sorted(coefficients) != [2, 3, 4, 5]:
        raise RuntimeError(f"missing exact axis coefficients: {sorted(coefficients)}")
    return coefficients


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path, help="optional output certificate JSON")
    args = parser.parse_args()

    coefficients = load_axis_coefficients()
    payload = {
        "schema": "matching-one/axis-irreducibility-certificate/v1",
        "method": "Gauss lemma + Rabin irreducibility over a prespecified finite field",
        "certificates": [],
    }

    for L in sorted(coefficients):
        coeffs = coefficients[L]
        primitive_content = math.gcd(*[abs(value) for value in coeffs])
        prime = CERTIFICATE_PRIMES[L]
        result = rabin_irreducible(coeffs, prime)
        same_degree = coeffs[-1] % prime != 0
        passed = primitive_content == 1 and same_degree and bool(result["irreducible"])
        row = {
            "L": L,
            "N": L * L,
            "integer_degree": len(coeffs) - 1,
            "integer_content": primitive_content,
            "certificate_prime": prime,
            "same_degree_mod_prime": same_degree,
            "finite_field": result,
            "irreducible_over_Q": passed,
            "physical_root_algebraic_degree_if_root_is_selected": (L * L if passed else None),
        }
        payload["certificates"].append(row)
        print(
            f"L={L} degree={len(coeffs)-1} prime={prime} "
            f"irreducible_mod_p={result['irreducible']} irreducible_Q={passed}"
        )
        if not passed:
            raise SystemExit(f"certificate failed at L={L}")

    payload["all_irreducible_over_Q"] = True
    payload["pairwise_gcd_consequence"] = (
        "Degrees 4,9,16,25 are distinct and every polynomial is irreducible over Q; "
        "therefore every pair has gcd 1 in Q[p]."
    )

    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

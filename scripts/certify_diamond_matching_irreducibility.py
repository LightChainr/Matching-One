#!/usr/bin/env python3
"""Finite-field irreducibility certificates for exact diamond matching polynomials.

This is the diamond companion to ``certify_axis_matching_irreducibility.py``.
It reuses the dependency-free Rabin verifier and the exact coefficients already
committed by the tiny zero-map pilot.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from certify_axis_matching_irreducibility import ROOT, rabin_irreducible


CERTIFICATE_PRIMES = {1: 3, 2: 3, 3: 79}


def load_diamond_coefficients() -> dict[int, list[int]]:
    payload = json.loads(
        (ROOT / "results" / "exact-zero-map-pilot" / "zero_map.json").read_text(
            encoding="utf-8"
        )
    )
    coefficients: dict[int, list[int]] = {}
    for row in payload["polynomials"]:
        if row["geometry"] == "diamond" and row["L"] in (1, 2, 3):
            coefficients[int(row["L"])] = [
                int(value) for value in row["power_coefficients_ascending"]
            ]
    if sorted(coefficients) != [1, 2, 3]:
        raise RuntimeError(f"missing exact diamond coefficients: {sorted(coefficients)}")
    return coefficients


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()

    coefficients = load_diamond_coefficients()
    payload: dict[str, object] = {
        "schema": "matching-one/diamond-irreducibility-certificate/v1",
        "method": "Gauss lemma + Rabin irreducibility over a prespecified finite field",
        "certificates": [],
    }

    for L in sorted(coefficients):
        values = coefficients[L]
        prime = CERTIFICATE_PRIMES[L]
        content = math.gcd(*[abs(value) for value in values])
        result = rabin_irreducible(values, prime)
        same_degree = values[-1] % prime != 0
        passed = content == 1 and same_degree and bool(result["irreducible"])
        row = {
            "geometry": "diamond",
            "L": L,
            "N": 2 * L * L,
            "integer_degree": len(values) - 1,
            "integer_content": content,
            "certificate_prime": prime,
            "same_degree_mod_prime": same_degree,
            "finite_field": result,
            "irreducible_over_Q": passed,
            "physical_root_algebraic_degree_if_root_is_selected": (
                2 * L * L if passed else None
            ),
        }
        payload["certificates"].append(row)
        print(
            f"diamond L={L} degree={len(values)-1} prime={prime} "
            f"irreducible_Q={passed}"
        )
        if not passed:
            raise SystemExit(f"diamond certificate failed at L={L}")

    payload["all_irreducible_over_Q"] = True
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

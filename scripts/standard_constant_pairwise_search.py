#!/usr/bin/env python3
"""Certified bounded pairwise searches for the six frozen standard constants."""

from __future__ import annotations

import argparse
from bisect import bisect_left
from fractions import Fraction
from functools import lru_cache
from math import gcd
import hashlib
import json
import mpmath as mp
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = ROOT / "analysis" / "pslq_search_contract.json"
DEFAULT_OUTPUT = ROOT / "results" / "pslq-standard-constant-pairwise" / "latest.json"
SCHEMA = "matching-one/standard-constant-pairwise-search/v1"

CONSTANT_BOUNDS = {
    "pi": ("3.14159265358979323846264338327950288419716939937510", "3.14159265358979323846264338327950288419716939937511"),
    "e": ("2.71828182845904523536028747135266249775724709369995", "2.71828182845904523536028747135266249775724709369996"),
    "log2": ("0.69314718055994530941723212145817656807550013436025", "0.69314718055994530941723212145817656807550013436026"),
    "sqrt2": ("1.41421356237309504880168872420969807856967187537694", "1.41421356237309504880168872420969807856967187537695"),
    "sqrt3": ("1.73205080756887729352744634150587236694280525381038", "1.73205080756887729352744634150587236694280525381039"),
    "sqrt5": ("2.23606797749978969640917366873127623544061835961152", "2.23606797749978969640917366873127623544061835961153"),
}
IV_CONSTANTS = {
    "pi": lambda: mp.iv.pi,
    "e": lambda: mp.iv.e,
    "log2": lambda: mp.iv.log(2),
    "sqrt2": lambda: mp.iv.sqrt(2),
    "sqrt3": lambda: mp.iv.sqrt(3),
    "sqrt5": lambda: mp.iv.sqrt(5),
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def _ceil(value: Fraction) -> int:
    return -((-value.numerator) // value.denominator)


def _floor(value: Fraction) -> int:
    return value.numerator // value.denominator


def _scaled_interval(coefficient: int, lower: Fraction, upper: Fraction) -> tuple[Fraction, Fraction]:
    return (coefficient * lower, coefficient * upper) if coefficient >= 0 else (coefficient * upper, coefficient * lower)


def _distance(lower: Fraction, upper: Fraction) -> Fraction:
    if lower <= 0 <= upper:
        return Fraction(0)
    return min(abs(lower), abs(upper))


def primitive_relation_count(height: int) -> int:
    return sum(
        1
        for b in range(1, height + 1)
        for c in range(-height, height + 1)
        if c
        for a in range(-height, height + 1)
        if gcd(gcd(abs(a), b), abs(c)) == 1
    )


def verify_constant_enclosures() -> None:
    mp.iv.dps = 80
    for constant_id, (lower, upper) in CONSTANT_BOUNDS.items():
        enclosure = mp.iv.mpf([lower, upper])
        _require(IV_CONSTANTS[constant_id]() in enclosure, f"{constant_id} enclosure is not outward")


def _candidate_as(allowed: Sequence[int], target_lower: Fraction, target_upper: Fraction) -> tuple[int, ...]:
    start = bisect_left(allowed, _ceil(target_lower))
    candidates = set()
    for index in (start - 1, start, start + 1):
        if 0 <= index < len(allowed):
            candidates.add(allowed[index])
    if start < len(allowed) and allowed[start] <= _floor(target_upper):
        candidates.add(allowed[start])
    return tuple(sorted(candidates))


def run_search(contract: Mapping[str, Any]) -> tuple[int, list[dict[str, Any]]]:
    stage = contract["search_stages"]["standard_constant_pairwise"]
    height = stage["coefficient_height_max"]
    _require(height == 100 and stage["basis_template"] == ["1", "p_c", "constant"], "pairwise contract drift")
    library = [row["id"] for row in stage["library"]]
    _require(library == list(CONSTANT_BOUNDS), "constant library/order drift")
    intervals = [(row, Fraction(row["lower"]), Fraction(row["upper"])) for row in contract["intervals"]]
    relation_count = primitive_relation_count(height)
    results = []
    for constant_id in library:
        constant_lower, constant_upper = map(Fraction, CONSTANT_BOUNDS[constant_id])
        best = [None] * len(intervals)
        containing = [0] * len(intervals)
        for b in range(1, height + 1):
            for c in range(-height, height + 1):
                if c == 0:
                    continue
                common = gcd(b, abs(c))
                allowed = [a for a in range(-height, height + 1) if gcd(abs(a), common) == 1]
                c_lower, c_upper = _scaled_interval(c, constant_lower, constant_upper)
                for index, (_, p_lower, p_upper) in enumerate(intervals):
                    base_lower = b * p_lower + c_lower
                    base_upper = b * p_upper + c_upper
                    for a in _candidate_as(allowed, -base_upper, -base_lower):
                        residual_lower, residual_upper = a + base_lower, a + base_upper
                        distance = _distance(residual_lower, residual_upper)
                        coefficients = (a, b, c)
                        candidate = (distance, coefficients, residual_lower, residual_upper)
                        if distance == 0:
                            containing[index] += 1
                        if best[index] is None or candidate[:2] < best[index][:2]:
                            best[index] = candidate
        for index, ((row, _, _), winner, zero_count) in enumerate(zip(intervals, best, containing)):
            assert winner is not None
            distance, coefficients, residual_lower, residual_upper = winner
            results.append({
                "constant_id": constant_id,
                "interval_id": row["id"],
                "source_id": row["source_id"],
                "primitive_relations_checked": relation_count,
                "zero_containing_residuals": zero_count,
                "excluded": zero_count == 0,
                "closest_relation": {
                    "coefficients_for_1_p_constant": list(coefficients),
                    "height": max(abs(value) for value in coefficients),
                    "minimum_absolute_residual": _text(distance),
                    "residual_interval": [_text(residual_lower), _text(residual_upper)],
                    "p_slope_abs": coefficients[1],
                    "root_condition_to_additive_residual": _text(Fraction(1, coefficients[1])),
                },
            })
    return relation_count, results


@lru_cache(maxsize=2)
def build_result(contract_path: Path = DEFAULT_CONTRACT) -> dict[str, Any]:
    contract_bytes = contract_path.read_bytes()
    contract = json.loads(contract_bytes)
    provenance = contract["provenance"]
    digest = hashlib.sha256((ROOT / provenance["path"]).read_bytes()).hexdigest()
    _require(digest == provenance["sha256"], "provenance digest drift")
    verify_constant_enclosures()
    checked, results = run_search(contract)
    return {
        "schema": SCHEMA,
        "issue": 1,
        "status": "standard_constant_pairwise_search_complete",
        "contract_sha256": hashlib.sha256(contract_bytes).hexdigest(),
        "provenance_sha256": digest,
        "constant_enclosures": {key: list(value) for key, value in CONSTANT_BOUNDS.items()},
        "constant_enclosures_verified_with_mpmath_iv": True,
        "search": {
            "basis": ["1", "p_c", "constant"],
            "coefficient_height_max": 100,
            "primitive_only": True,
            "sign_normalization": "p_coefficient_positive",
            "constant_coefficient_nonzero": True,
            "relations_per_constant_interval": checked,
            "method_interval_count": len(contract["intervals"]),
            "constant_count": len(CONSTANT_BOUNDS),
            "arithmetic": "exact Fraction arithmetic over fixed outward decimal enclosures",
        },
        "results": results,
        "conclusion": {
            "all_frozen_pairwise_bases_excluded": all(row["excluded"] for row in results),
            "zero_containing_residual_total": sum(row["zero_containing_residuals"] for row in results),
        },
        "claim_boundary": {
            "included": "frozen pairwise standard-constant bases at height 100 over each method interval",
            "excluded": "expanded libraries, multi-constant bases, degree-3/4 polynomial exclusions, near-hit promotion, closed forms, or transcendence",
            "parent_issue": "remain open",
        },
    }


def validate_result(result: Mapping[str, Any], contract_path: Path = DEFAULT_CONTRACT) -> Mapping[str, Any]:
    expected = build_result(contract_path)
    _require(result == expected, "standard-constant result does not exactly reproduce")
    return {"schema": SCHEMA, "status": "valid", "rows": len(expected["results"]), "all_excluded": expected["conclusion"]["all_frozen_pairwise_bases_excluded"]}


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--validate", type=Path)
    args = parser.parse_args(argv)
    if args.validate:
        print(json.dumps(validate_result(json.loads(args.validate.read_text()), args.contract), indent=2, sort_keys=True))
        return 0
    rendered = json.dumps(build_result(args.contract), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

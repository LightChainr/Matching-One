#!/usr/bin/env python3
"""Certified exclusion of primitive polynomials C(1..6, 3) on one frozen interval.

This is a **second, independent implementation** of the census in
``scripts/degree6_low_height_exclusion.py``, kept alongside it rather than
replacing it.  It was written separately from the frozen issue #559 protocol,
and the two differ where it matters: this one screens at the interval
*midpoint* (``|P(m)| > D*(u-l)/2`` rules out a root), the other screens at both
*endpoints* (``|P(l)| > D*(u-l)``).  Their enumerations of ``C(d,3)`` are also
written independently.  ``scripts/degree6_implementation_agreement.py`` compares
the two artifact sets cell by cell and is what makes the agreement a checked
claim rather than a remark.

Both implementations import the repository's ``exact_polynomial_root_certificate``
unchanged, so the Sturm path is shared, not replicated.  That shared code
contributes nothing to the null reported here: the screen retains zero
candidates at every degree on every interval, so root isolation never runs
during the census.  It does run in the planted-root controls, where the
replication is therefore partial.

Committed as received apart from this note; no line of its logic was altered.

Implements the frozen issue #559 protocol:

  hypothesis: no primitive integer polynomial of degree <= 6 and coefficient
  height <= 3 has a real root inside any of the four frozen method intervals.

Protocol (per interval):
  1. Enumerate the class C(d,3) for d = 1..6: primitive, sign-normalised
     a_d >= 1, |a_i| <= 3, exact degree d.  Total 409,584 polynomials.
  2. Screen by the mean-value bound.  For the interval midpoint m and
     derivative bound D = sum_{k>=1} k|a_k| <= 63, any polynomial with a root
     in the interval satisfies |P(m)| <= D * width / 2.  The midpoint value is
     evaluated in exact rational arithmetic (no rounding, no binary float), so
     a polynomial failing the screen is certified root-free in the interval.
  3. Exact decision.  Every screen survivor is sent to
     exact_polynomial_root_certificate.isolate_roots over [lower, upper] at
     120-bit dyadic isolation with exact rational endpoint evaluation.  A root
     is reported only when the certified isolation is non-empty; exclusion is
     asserted only when the isolation is empty.

No binary floating point appears in any exclusion claim.

This driver mirrors the house pattern of scripts/degree4_interval_exclusion.py:
it exposes build_result / validate_result / run_search and a --output/--validate
CLI, and reproduces deterministically.  It differs from the degree-4 census in
scope (degree 1..6, height 3 -- the frozen #559 class) and in that its screen
is a pure-Python exact rational screen rather than the C++ fixed-point
meet-in-the-middle screen that the height-100 quartic census requires.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from functools import lru_cache
from math import gcd
import hashlib
import json
from pathlib import Path
from typing import Any, Iterator, Mapping, Optional, Sequence

try:
    from scripts.exact_polynomial_root_certificate import evaluate, isolate_roots
except ModuleNotFoundError:
    from exact_polynomial_root_certificate import evaluate, isolate_roots


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = ROOT / "analysis" / "pslq_search_contract.json"
SCHEMA = "matching-one/degree56-interval-exclusion/v1"
ISSUE = 559
DEGREE_MIN = 1
DEGREE_MAX = 6
HEIGHT = 3
ISOLATION_BITS = 120
EXPECTED_CLASS_SIZES = {1: 15, 2: 129, 3: 975, 4: 7041, 5: 49935, 6: 351489}
EXPECTED_CLASS_TOTAL = 409584


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def class_size(degree: int, height: int = HEIGHT) -> int:
    """Number of primitive sign-normalised degree-`degree` polynomials with
    |a_i| <= height.  Validated against the frozen counts for C(d,3)."""
    count = 0
    for coefficient in _enumerate_degree(degree, height):
        count += 1
    return count


def _enumerate_degree(degree: int, height: int = HEIGHT) -> Iterator[tuple[int, ...]]:
    """Yield primitive (gcd 1), sign-normalised (a_d >= 1), |a_i| <= height,
    exact-degree-`degree` polynomials as ascending-coefficient tuples."""
    if degree < 1:
        return
    lead_min, lead_max = 1, height
    rest_min, rest_max = -height, height
    # a_0..a_{d-1} in [-h..h], a_d in [1..h]
    from itertools import product

    for leading in range(lead_min, lead_max + 1):
        for rest in product(range(rest_min, rest_max + 1), repeat=degree):
            coefficients = (*rest, leading)
            common = 0
            for value in coefficients:
                common = gcd(common, abs(value))
                if common == 1:
                    break
            if common == 1:
                yield coefficients


@lru_cache(maxsize=8)
def _degree_polynomials(degree: int, height: int = HEIGHT) -> tuple[tuple[int, ...], ...]:
    return tuple(_enumerate_degree(degree, height))


def _derivative_bound(coefficients: Sequence[int]) -> int:
    """D = sum_{k>=1} k |a_k|; a certified |P'| bound on x in (0,1)."""
    return sum(coefficient_index * abs(coefficient)
               for coefficient_index, coefficient in enumerate(coefficients))


def _screen(coefficients: Sequence[int], midpoint: Fraction,
            half_width: Fraction) -> tuple[Fraction, Fraction, bool]:
    """Exact midpoint value and derivative-bound screen.

    Returns (midpoint_value, derivative_bound, keep).  keep=True means the
    polynomial is not certifiably root-free by the mean-value screen.
    """
    bound = Fraction(_derivative_bound(coefficients))
    polynomial = [Fraction(value) for value in coefficients]
    value = evaluate(polynomial, midpoint)
    keep = abs(value) <= bound * half_width
    return value, bound, keep


def _interval_parameters(lower: Fraction, upper: Fraction) -> dict[str, Any]:
    width = upper - lower
    midpoint = (lower + upper) / 2
    return {
        "lower": lower,
        "upper": upper,
        "width": width,
        "half_width": width / 2,
        "midpoint": midpoint,
    }


def _is_rootless_isolate(coefficients: Sequence[int], lower: Fraction,
                         upper: Fraction) -> list[tuple[Fraction, Fraction]]:
    polynomial = [Fraction(value) for value in coefficients]
    return isolate_roots(polynomial, lower, upper, bits=ISOLATION_BITS)


def run_search(interval: Mapping[str, Any],
               degree_min: int = DEGREE_MIN,
               degree_max: int = DEGREE_MAX,
               height: int = HEIGHT) -> dict[str, Any]:
    """Full class census of C(degree_min..degree_max, height) on one interval.

    The returned record is fully deterministic (no wall-clock fields) so that
    validate_result can require exact reproduction.
    """
    lower, upper = Fraction(interval["lower"]), Fraction(interval["upper"])
    _require(lower < upper, "interval must be nonempty")
    parameters = _interval_parameters(lower, upper)
    midpoint = parameters["midpoint"]
    half_width = parameters["half_width"]

    class_sizes: dict[int, int] = {}
    per_degree: list[dict[str, Any]] = []
    screen_survivors_total = 0
    root_containing_total = 0
    distinct_roots_total = 0
    root_witnesses: list[dict[str, Any]] = []
    closest: tuple[Fraction, tuple[int, ...]] | None = None

    for degree in range(degree_min, degree_max + 1):
        polynomials = _degree_polynomials(degree, height)
        _require(len(polynomials) == EXPECTED_CLASS_SIZES[degree],
                 f"class size drift at degree {degree}")
        class_sizes[degree] = len(polynomials)
        survivors = 0
        root_containing = 0
        distinct_roots = 0
        for coefficients in polynomials:
            _require(coefficients[-1] > 0, "sign normalisation drift")
            common = 0
            for value in coefficients:
                common = gcd(common, abs(value))
            _require(common == 1, "nonprimitive polynomial escaped enumeration")
            midpoint_value, derivative_bound, keep = _screen(
                coefficients, midpoint, half_width)
            if closest is None or abs(midpoint_value) < abs(closest[0]):
                closest = (midpoint_value, coefficients)
            if not keep:
                continue
            survivors += 1
            roots = _is_rootless_isolate(coefficients, lower, upper)
            if roots:
                root_containing += 1
                distinct_roots += len(roots)
                root_witnesses.append({
                    "degree": degree,
                    "height": max(abs(value) for value in coefficients),
                    "coefficients_ascending": list(coefficients),
                    "root_brackets": [[_text(lo), _text(hi)] for lo, hi in roots],
                    "isolation_bits": ISOLATION_BITS,
                })
        screen_survivors_total += survivors
        root_containing_total += root_containing
        distinct_roots_total += distinct_roots
        per_degree.append({
            "degree": degree,
            "class_size": class_sizes[degree],
            "screen_survivors": survivors,
            "root_containing_polynomials": root_containing,
            "distinct_roots_in_interval": distinct_roots,
        })

    _require(closest is not None, "empty polynomial class")
    return {
        "interval_id": interval["id"],
        "source_id": interval.get("source_id"),
        "lower": _text(lower),
        "upper": _text(upper),
        "width": _text(parameters["width"]),
        "midpoint": _text(midpoint),
        "class_size_total": sum(class_sizes.values()),
        "by_degree": per_degree,
        "screen_survivors_total": screen_survivors_total,
        "root_containing_polynomials": root_containing_total,
        "distinct_roots_in_interval": distinct_roots_total,
        "root_witnesses": root_witnesses,
        "excluded": root_containing_total == 0,
        "closest_polynomial_at_midpoint": {
            "coefficients_ascending": list(closest[1]),
            "degree": len(closest[1]) - 1,
            "height": max(abs(value) for value in closest[1]),
            "minimum_absolute_residual_at_midpoint": _text(abs(closest[0])),
        },
    }


@lru_cache(maxsize=8)
def build_result(interval_id: str,
                 contract_path: Path = DEFAULT_CONTRACT) -> dict[str, Any]:
    raw = contract_path.read_bytes()
    contract = json.loads(raw)
    rows = [row for row in contract["intervals"] if row["id"] == interval_id]
    _require(len(rows) == 1, "interval id is not uniquely frozen")
    provenance = contract["provenance"]
    digest = hashlib.sha256((ROOT / provenance["path"]).read_bytes()).hexdigest()
    _require(digest == provenance["sha256"], "provenance digest drift")
    interval_result = run_search(rows[0])
    return {
        "schema": SCHEMA,
        "issue": ISSUE,
        "status": "degree56_interval_exact_exclusion_complete",
        "contract_sha256": hashlib.sha256(raw).hexdigest(),
        "provenance_sha256": digest,
        "search": {
            "degree_min": DEGREE_MIN,
            "degree_max": DEGREE_MAX,
            "coefficient_height_max": HEIGHT,
            "primitive_coefficients_only": True,
            "sign_normalization": "a_d_positive",
            "screen": "exact rational midpoint evaluation with mean-value derivative bound",
            "exact_decision": "120-bit Sturm isolation with exact rational endpoints",
            "arithmetic": "exact rational throughout; no binary float in any exclusion claim",
        },
        "interval_result": interval_result,
        "claim_boundary": {
            "included": (f"certified exclusion on the frozen {interval_id} method interval "
                         f"over the declared finite class C(1..6,{HEIGHT})"),
            "excluded": ("other method intervals, higher degree/height, library expansion, "
                         "near-hit promotion, p-values, closed forms, cross-interval "
                         "resolution, or transcendence"),
            "parent_issue": "remain open",
        },
    }


def validate_result(result: Mapping[str, Any], interval_id: str) -> Mapping[str, Any]:
    expected = build_result(interval_id)
    _require(result == expected, "degree-56 interval result does not exactly reproduce")
    row = expected["interval_result"]
    return {
        "schema": SCHEMA,
        "status": "valid",
        "interval_id": interval_id,
        "class_size_total": row["class_size_total"],
        "excluded": row["excluded"],
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("interval_id")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--validate", type=Path)
    args = parser.parse_args(argv)
    if args.validate:
        print(json.dumps(validate_result(json.loads(args.validate.read_text()), args.interval_id),
                         indent=2, sort_keys=True))
        return 0
    rendered = json.dumps(build_result(args.interval_id), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Certified meet-in-the-middle exclusion of frozen quartics on one interval."""

from __future__ import annotations

import argparse
from fractions import Fraction
from functools import lru_cache
from math import gcd
import hashlib
import json
from pathlib import Path
import re
import subprocess
import tempfile
from typing import Any, Mapping, Optional, Sequence

try:
    from scripts.exact_polynomial_root_certificate import derivative, evaluate, isolate_roots
except ModuleNotFoundError:
    from exact_polynomial_root_certificate import derivative, evaluate, isolate_roots


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = ROOT / "analysis" / "pslq_search_contract.json"
CPP_SOURCE = ROOT / "scripts" / "degree4_fixed_point_screen.cpp"
SCHEMA = "matching-one/degree4-interval-exclusion/v1"
HEIGHT = 100
SCALE = 10**15
NEAR_MARGIN_SCALED = 10**9  # 1e-6 at SCALE.


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def _ceil(value: Fraction) -> int:
    return -((-value.numerator) // value.denominator)


def _nearest_integer(value: Fraction) -> int:
    _require(value >= 0, "fixed-point powers must be nonnegative")
    return (2 * value.numerator + value.denominator) // (2 * value.denominator)


def _mobius(limit: int) -> list[int]:
    values = [1] * (limit + 1)
    prime = [True] * (limit + 1)
    for p in range(2, limit + 1):
        if not prime[p]:
            continue
        for multiple in range(p, limit + 1, p):
            prime[multiple] = False if multiple != p else prime[multiple]
            values[multiple] *= -1
        square = p * p
        for multiple in range(square, limit + 1, square):
            values[multiple] = 0
    return values


def primitive_quartic_count(height: int = HEIGHT) -> int:
    mu = _mobius(height)
    return sum(
        mu[d] * (height // d) * (2 * (height // d) + 1) ** 4
        for d in range(1, height + 1)
    )


def output_path(interval_id: str) -> Path:
    return ROOT / "results" / f"pslq-degree4-{interval_id}" / "latest.json"


def _screen_parameters(lower: Fraction, upper: Fraction) -> dict[str, Any]:
    midpoint = (lower + upper) / 2
    exact_scaled = [Fraction(SCALE) * midpoint**power for power in range(1, 5)]
    weights = [_nearest_integer(value) for value in exact_scaled]
    rounding_error = sum((abs(value - weight) for value, weight in zip(exact_scaled, weights)), Fraction())
    global_derivative_bound = HEIGHT * sum(range(1, 5))
    root_error = Fraction(SCALE * global_derivative_bound) * (upper - lower) / 2 + HEIGHT * rounding_error
    root_bound = _ceil(root_error)
    return {
        "midpoint": midpoint,
        "weights": weights,
        "weight_rounding_error_bound_scaled": rounding_error,
        "global_derivative_bound": global_derivative_bound,
        "root_filter_bound_scaled": root_bound,
        "near_filter_bound_scaled": root_bound + NEAR_MARGIN_SCALED,
    }


def _run_screen(parameters: Mapping[str, Any]) -> tuple[list[tuple[tuple[int, ...], int, bool]], str]:
    with tempfile.TemporaryDirectory(prefix="degree4-screen-", dir=ROOT) as directory:
        binary = Path(directory) / "degree4_fixed_point_screen"
        subprocess.run(
            ["g++", "-std=c++17", "-O3", "-DNDEBUG", str(CPP_SOURCE), "-o", str(binary)],
            check=True,
            capture_output=True,
            text=True,
        )
        completed = subprocess.run(
            [
                str(binary),
                *[str(value) for value in parameters["weights"]],
                str(parameters["root_filter_bound_scaled"]),
                str(parameters["near_filter_bound_scaled"]),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
    candidates = []
    for line in completed.stdout.splitlines():
        values = [int(value) for value in line.split()]
        _require(len(values) == 7, "malformed C++ candidate row")
        coefficients = tuple(values[:5])
        candidates.append((coefficients, values[5], bool(values[6])))
    match = re.fullmatch(r"near_candidates=(\d+) root_filter_candidates=(\d+)\n?", completed.stderr)
    _require(match is not None, "malformed C++ screen summary")
    _require(int(match.group(1)) == len(candidates), "C++ candidate count drift")
    _require(int(match.group(2)) == sum(row[2] for row in candidates), "C++ root-filter count drift")
    return candidates, completed.stderr.strip()


def run_search(interval: Mapping[str, Any]) -> dict[str, Any]:
    lower, upper = Fraction(interval["lower"]), Fraction(interval["upper"])
    _require(lower < upper, "interval must be nonempty")
    parameters = _screen_parameters(lower, upper)
    candidates, screen_summary = _run_screen(parameters)
    root_filter_candidates = sum(row[2] for row in candidates)
    root_containing = 0
    distinct_roots = 0
    root_witnesses: list[dict[str, Any]] = []
    stationary_candidates = 0
    best: tuple[Fraction, tuple[int, ...], tuple[Fraction, Fraction]] | None = None

    for coefficients, fixed_residual, root_filter in candidates:
        _require(coefficients[-1] > 0, "quartic sign normalization drift")
        common = 0
        for value in coefficients:
            common = gcd(common, abs(value))
        _require(common == 1, "nonprimitive quartic escaped the screen")
        _require(abs(fixed_residual) <= parameters["near_filter_bound_scaled"], "candidate escaped near bound")
        if root_filter:
            _require(abs(fixed_residual) <= parameters["root_filter_bound_scaled"], "root-filter flag drift")

        polynomial = [Fraction(value) for value in coefficients]
        stationary = isolate_roots(derivative(polynomial), lower, upper, bits=120)
        stationary_candidates += bool(stationary)
        _require(not stationary, "near candidate has an internal stationary point; exact algebraic range path required")
        roots = isolate_roots(polynomial, lower, upper, bits=120)
        root_containing += bool(roots)
        distinct_roots += len(roots)
        if roots:
            root_witnesses.append(
                {
                    "coefficients_ascending": list(coefficients),
                    "root_brackets": [[_text(lo), _text(hi)] for lo, hi in roots],
                    "isolation_bits": 120,
                }
            )
        endpoint_values = (evaluate(polynomial, lower), evaluate(polynomial, upper))
        minimum = Fraction(0) if roots else min(abs(value) for value in endpoint_values)
        candidate = (minimum, coefficients, endpoint_values)
        if best is None or candidate[:2] < best[:2]:
            best = candidate

    assert best is not None
    _require(best[0] < Fraction(1, 10**6), "near filter is too narrow to certify the global closest quartic")
    # If a quartic had an interval root, the MVT plus fixed-point rounding
    # bounds force it into root_filter_candidates, all of which receive exact
    # root decisions above.  Likewise, any quartic
    # closer than 1e-6 must enter the near set; the winning residual is below
    # that margin, so the global closest witness is certified.
    return {
        "interval_id": interval["id"],
        "source_id": interval["source_id"],
        "lower": interval["lower"],
        "upper": interval["upper"],
        "primitive_quartics_covered": primitive_quartic_count(),
        "fixed_point_scale": SCALE,
        "fixed_point_weights_for_powers_1_to_4": parameters["weights"],
        "global_derivative_bound": parameters["global_derivative_bound"],
        "weight_rounding_error_bound_scaled": _text(parameters["weight_rounding_error_bound_scaled"]),
        "root_filter_bound_scaled": parameters["root_filter_bound_scaled"],
        "near_filter_bound_scaled": parameters["near_filter_bound_scaled"],
        "near_candidates_exactly_checked": len(candidates),
        "root_filter_candidates": root_filter_candidates,
        "near_candidates_with_stationary_point": stationary_candidates,
        "root_containing_polynomials": root_containing,
        "distinct_roots_in_interval": distinct_roots,
        "root_witnesses": root_witnesses,
        "excluded": root_containing == 0,
        "screen_summary": screen_summary,
        "closest_polynomial": {
            "coefficients_ascending": list(best[1]),
            "height": max(abs(value) for value in best[1]),
            "minimum_absolute_residual": _text(best[0]),
            "polynomial_endpoint_values": [_text(value) for value in best[2]],
            "global_closest_certified_by_margin": "1e-6",
        },
    }


@lru_cache(maxsize=8)
def build_result(interval_id: str, contract_path: Path = DEFAULT_CONTRACT) -> dict[str, Any]:
    raw = contract_path.read_bytes()
    contract = json.loads(raw)
    stage = contract["search_stages"]["algebraic_polynomial"]
    _require(stage["degree_max"] >= 4 and stage["coefficient_height_max"] == HEIGHT, "quartic contract drift")
    rows = [row for row in contract["intervals"] if row["id"] == interval_id]
    _require(len(rows) == 1, "interval id is not uniquely frozen")
    provenance = contract["provenance"]
    digest = hashlib.sha256((ROOT / provenance["path"]).read_bytes()).hexdigest()
    _require(digest == provenance["sha256"], "provenance digest drift")
    return {
        "schema": SCHEMA,
        "issue": 1,
        "status": "degree4_interval_exact_census_complete",
        "contract_sha256": hashlib.sha256(raw).hexdigest(),
        "provenance_sha256": digest,
        "search": {
            "degree": 4,
            "coefficient_height_max": HEIGHT,
            "primitive_coefficients_only": True,
            "sign_normalization": "a4_positive",
            "screen": "certified fixed-point meet-in-the-middle followed by exact rational/Sturm decisions",
        },
        "interval_result": run_search(rows[0]),
        "claim_boundary": {
            "included": f"degree-4 height-100 exclusion on {interval_id} only",
            "excluded": "other method intervals, higher bounds, library expansion, near-hit promotion, p-values, closed forms, or transcendence",
            "parent_issue": "remain open",
        },
    }


def validate_result(result: Mapping[str, Any], interval_id: str) -> Mapping[str, Any]:
    expected = build_result(interval_id)
    _require(result == expected, "degree-4 interval result does not exactly reproduce")
    row = expected["interval_result"]
    return {
        "schema": SCHEMA,
        "status": "valid",
        "interval_id": interval_id,
        "primitive_quartics_covered": row["primitive_quartics_covered"],
        "excluded": row["excluded"],
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("interval_id")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--validate", type=Path)
    args = parser.parse_args(argv)
    if args.validate:
        print(json.dumps(validate_result(json.loads(args.validate.read_text()), args.interval_id), indent=2, sort_keys=True))
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

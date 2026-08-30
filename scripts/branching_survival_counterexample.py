#!/usr/bin/env python3
"""Exact N=16 trace-equivalent but branching-distinct survival witness."""

from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from itertools import combinations
import json
from math import comb
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, Sequence

try:
    from scripts.integer_period_torus import IntegerTorusGeometry, axis_integer_torus
    from scripts.rank_one_survival_certificate import (
        RankCache,
        fraction_text,
        trigger_layers,
        vacant_vertices,
    )
except ModuleNotFoundError:
    from integer_period_torus import IntegerTorusGeometry, axis_integer_torus  # type: ignore
    from rank_one_survival_certificate import (  # type: ignore
        RankCache,
        fraction_text,
        trigger_layers,
        vacant_vertices,
    )


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "analysis" / "branching_survival_counterexample.json"
SCHEMA = "matching-one/branching-survival-counterexample/v1"

A_COORDINATES = (
    (0, 0), (1, 0), (2, 0), (3, 0),
    (1, 1), (3, 1), (0, 3), (1, 3),
)
B_COORDINATES = (
    (0, 0), (1, 0), (2, 0),
    (0, 1), (1, 1), (2, 1), (3, 1), (0, 3),
)


def mask_from_coordinates(
    geometry: IntegerTorusGeometry, coordinates: Iterable[tuple[int, int]]
) -> int:
    mask = 0
    for coordinate in coordinates:
        if (
            not isinstance(coordinate, tuple)
            or len(coordinate) != 2
            or not all(isinstance(value, int) for value in coordinate)
        ):
            raise ValueError("coordinates must be integer pairs")
        vertex = geometry.vertex(coordinate)
        bit = 1 << vertex
        if mask & bit:
            raise ValueError("duplicate quotient coordinate")
        mask |= bit
    return mask


def row_major_mask(coordinates: Iterable[tuple[int, int]], length: int) -> int:
    if not isinstance(length, int) or length <= 0:
        raise ValueError("length must be positive")
    mask = 0
    for x, y in coordinates:
        if not 0 <= x < length or not 0 <= y < length:
            raise ValueError("coordinate outside the row-major fundamental domain")
        bit = 1 << (x + length * y)
        if mask & bit:
            raise ValueError("duplicate row-major coordinate")
        mask |= bit
    return mask


def complete_survival_counts(cache: RankCache, mask: int) -> tuple[int, ...]:
    if cache.rank(mask) != 1:
        raise ValueError("survival counts require a rank-one state")
    vacant = vacant_vertices(cache.geometry.n, mask)
    return tuple(
        sum(
            cache.rank(mask | sum(1 << vertex for vertex in added)) == 1
            for added in combinations(vacant, horizon)
        )
        for horizon in range(len(vacant) + 1)
    )


def complete_survival_probabilities(cache: RankCache, mask: int) -> tuple[Fraction, ...]:
    counts = complete_survival_counts(cache, mask)
    q = len(counts) - 1
    return tuple(Fraction(count, comb(q, horizon)) for horizon, count in enumerate(counts))


def successor_h2_distribution(cache: RankCache, mask: int) -> tuple[int, dict[int, int]]:
    if cache.rank(mask) != 1:
        raise ValueError("successor distribution requires a rank-one state")
    absorbed = 0
    safe_h2: Counter[int] = Counter()
    for vertex in vacant_vertices(cache.geometry.n, mask):
        child = mask | (1 << vertex)
        if cache.rank(child) == 2:
            absorbed += 1
        else:
            singleton, _ = trigger_layers(cache, child)
            safe_h2[len(singleton)] += 1
    return absorbed, dict(sorted(safe_h2.items()))


def branch_success_direct(cache: RankCache, mask: int) -> tuple[int, int]:
    """Enumerate a common insertion followed by one insertion in each clone."""

    if cache.rank(mask) != 1:
        raise ValueError("branch experiment requires a rank-one state")
    common_choices = vacant_vertices(cache.geometry.n, mask)
    numerator = 0
    denominator = 0
    for common in common_choices:
        child = mask | (1 << common)
        branch_choices = vacant_vertices(cache.geometry.n, child)
        denominator += len(branch_choices) ** 2
        if cache.rank(child) == 2:
            continue
        for first in branch_choices:
            for second in branch_choices:
                if (
                    cache.rank(child | (1 << first)) == 1
                    and cache.rank(child | (1 << second)) == 1
                ):
                    numerator += 1
    return numerator, denominator


def branch_success_from_h2(cache: RankCache, mask: int) -> Fraction:
    absorbed, safe_h2 = successor_h2_distribution(cache, mask)
    q = len(vacant_vertices(cache.geometry.n, mask))
    if q < 2:
        raise ValueError("branch experiment requires at least two vacancies")
    if absorbed + sum(safe_h2.values()) != q:
        raise ArithmeticError("successor distribution does not partition insertions")
    numerator = sum(count * (q - 1 - h2) ** 2 for h2, count in safe_h2.items())
    return Fraction(numerator, q * (q - 1) ** 2)


def witness_record(
    cache: RankCache, coordinates: Sequence[tuple[int, int]], expected_row_mask: int
) -> dict[str, Any]:
    mask = mask_from_coordinates(cache.geometry, coordinates)
    rank, line = cache.rank_and_line(mask)
    counts = complete_survival_counts(cache, mask)
    probabilities = complete_survival_probabilities(cache, mask)
    absorbed, safe_h2 = successor_h2_distribution(cache, mask)
    direct_numerator, direct_denominator = branch_success_direct(cache, mask)
    formula = branch_success_from_h2(cache, mask)
    if Fraction(direct_numerator, direct_denominator) != formula:
        raise ArithmeticError("direct and H2 branching calculations disagree")
    declared_mask = row_major_mask(coordinates, 4)
    if declared_mask != expected_row_mask:
        raise ArithmeticError("coordinate witness does not match its declared mask")
    return {
        "coordinates": [list(coordinate) for coordinate in coordinates],
        "row_major_mask": declared_mask,
        "k": len(coordinates),
        "rank": rank,
        "primitive_line": list(line) if line is not None else None,
        "complete_survival_counts": list(counts),
        "complete_survival_probabilities": [fraction_text(value) for value in probabilities],
        "successor_distribution": {
            "already_absorbed": absorbed,
            "safe_by_h2": {str(h2): count for h2, count in safe_h2.items()},
        },
        "branch_success_count": direct_numerator,
        "branch_total_count": direct_denominator,
        "branch_success_probability": fraction_text(formula),
    }


def build_artifact() -> dict[str, Any]:
    geometry = axis_integer_torus(4)
    cache = RankCache(geometry)
    first = witness_record(cache, A_COORDINATES, 12463)
    second = witness_record(cache, B_COORDINATES, 4343)
    if first["complete_survival_counts"] != second["complete_survival_counts"]:
        raise ArithmeticError("witnesses are not trace-equivalent")
    branch_gap = Fraction(first["branch_success_count"], first["branch_total_count"]) - Fraction(
        second["branch_success_count"], second["branch_total_count"]
    )
    return {
        "schema": SCHEMA,
        "issue": 429,
        "status": "exact_trace_equivalent_branching_counterexample",
        "geometry": {
            "name": "axis",
            "length": 4,
            "site_count": geometry.n,
            "label_rule": "row-major label x+4*y on coordinates modulo four",
        },
        "witnesses": {"A": first, "B": second},
        "comparison": {
            "complete_survival_counts_equal": True,
            "unbranched_two_step_survival": first["complete_survival_probabilities"][2],
            "branch_success_gap": fraction_text(branch_gap),
            "distinguishing_test": "one common insertion, clone, then one independent insertion in each clone",
        },
        "claim_boundary": {
            "included": "coordinate-defined N=16 full-survival and branching counterexample",
            "excluded": "prefix-history table, six-quotient refinement census, predictive algebra, approximate or scaled process claims",
            "parent_issue": "remain open",
        },
    }


def validate_artifact(artifact: Mapping[str, Any]) -> dict[str, Any]:
    regenerated = build_artifact()
    if artifact != regenerated:
        raise ValueError("branching survival artifact does not exactly reproduce")
    return {
        "schema": SCHEMA,
        "status": "valid_exact_branching_survival_counterexample",
        "survival_counts": regenerated["witnesses"]["A"]["complete_survival_counts"],
        "branch_success_gap": regenerated["comparison"]["branch_success_gap"],
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--validate", type=Path)
    args = parser.parse_args(argv)
    if args.validate:
        artifact = json.loads(args.validate.read_text(encoding="utf-8"))
        print(json.dumps(validate_artifact(artifact), indent=2, sort_keys=True))
        return 0
    rendered = json.dumps(build_artifact(), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

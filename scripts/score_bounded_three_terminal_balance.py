#!/usr/bin/env python3
"""Score exact three-terminal balance roots across the bounded gadget census."""

from __future__ import annotations

import argparse
from decimal import Decimal, getcontext
from fractions import Fraction
from functools import reduce
import hashlib
import json
from math import comb, gcd
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, Sequence

try:
    from scripts.bounded_gadget_census import is_connected, vertex_degrees
    from scripts.gadget_graph_canonical import (
        Graph,
        canonical_graph,
        decode_graph,
        enumerate_graphs,
    )
    from scripts.terminal_partition_canonical import full_symmetric_group
    from scripts.terminal_reliability_polynomial import enumerate_reliability
except ModuleNotFoundError:  # Direct ``python scripts/...`` execution.
    from bounded_gadget_census import is_connected, vertex_degrees
    from gadget_graph_canonical import Graph, canonical_graph, decode_graph, enumerate_graphs
    from terminal_partition_canonical import full_symmetric_group
    from terminal_reliability_polynomial import enumerate_reliability


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results" / "issue13-three-terminal-balance-screen" / "latest.json"
SCHEMA = "matching-one/bounded-three-terminal-balance-screen/v1"
SOURCE_COMMIT = "e7d63b4c08cfb63a9c027959087c8beaa6228fdf"
TARGET = Decimal("0.5927460507896")
ROOT_BITS = 192
ALL_CONNECTED = (0, 0, 0)
ALL_SEPARATE = (0, 1, 2)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _digest(payload: Mapping[str, Any]) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def connected_orbit_catalog() -> Mapping[str, tuple[int, Graph]]:
    """Return canonical encoding -> (labeled orbit multiplicity, representative)."""

    vertex_count = 4
    terminal_count = 3
    group = full_symmetric_group(terminal_count)
    buckets: dict[str, list[Graph]] = {}
    for graph in enumerate_graphs(vertex_count, terminal_count):
        key = canonical_graph(vertex_count, terminal_count, graph, group)
        buckets.setdefault(key, []).append(graph)
    return {
        key: (len(members), decode_graph(key)[2])
        for key, members in sorted(buckets.items())
        if is_connected(vertex_count, decode_graph(key)[2])
    }


def bernstein_to_power(counts: Sequence[int]) -> tuple[int, ...]:
    """Expand sum_k counts[k] p^k (1-p)^(m-k), low degree first."""

    edge_count = len(counts) - 1
    coefficients = [0] * (edge_count + 1)
    for opened, count in enumerate(counts):
        for extra in range(edge_count - opened + 1):
            coefficients[opened + extra] += count * comb(edge_count - opened, extra) * ((-1) ** extra)
    return tuple(coefficients)


def primitive_polynomial(coefficients: Iterable[int]) -> tuple[int, ...]:
    values = list(coefficients)
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    divisor = reduce(gcd, (abs(value) for value in values if value), 0)
    _require(divisor > 0, "zero polynomial has no primitive normalization")
    values = [value // divisor for value in values]
    if values[0] > 0:
        values = [-value for value in values]
    return tuple(values)


def evaluate_power(coefficients: Sequence[int], p: Fraction) -> Fraction:
    value = Fraction(0)
    for coefficient in reversed(coefficients):
        value = value * p + coefficient
    return value


def isolate_unit_root(coefficients: Sequence[int], bits: int = ROOT_BITS) -> tuple[Fraction, Fraction]:
    """Isolate the monotone all-connected/all-separated crossing by dyadic bisection."""

    low, high = Fraction(0), Fraction(1)
    _require(evaluate_power(coefficients, low) < 0, "balance must be negative at p=0")
    _require(evaluate_power(coefficients, high) > 0, "balance must be positive at p=1")
    for _ in range(bits):
        middle = (low + high) / 2
        if evaluate_power(coefficients, middle) < 0:
            low = middle
        else:
            high = middle
    return low, high


def _fraction_text(value: Fraction) -> str:
    return "%d/%d" % (value.numerator, value.denominator)


def _decimal(value: Fraction) -> Decimal:
    return Decimal(value.numerator) / Decimal(value.denominator)


def evaluate_bernstein_decimal(counts: Sequence[int], p: Decimal) -> Decimal:
    edge_count = len(counts) - 1
    return sum(
        Decimal(count) * p**opened * (Decimal(1) - p) ** (edge_count - opened)
        for opened, count in enumerate(counts)
    )


def score_candidate(encoding: str, multiplicity: int, graph: Graph) -> dict[str, Any]:
    gadget = {
        "id": encoding,
        "vertex_count": 4,
        "terminal_count": 3,
        "edges": [list(edge) for edge in graph],
    }
    reliability = enumerate_reliability(gadget)
    edge_count = len(graph)
    zero = (0,) * (edge_count + 1)
    connected_counts = reliability.get(ALL_CONNECTED, zero)
    separate_counts = reliability.get(ALL_SEPARATE, zero)
    balance_bernstein = tuple(a - b for a, b in zip(connected_counts, separate_counts))
    balance_power = primitive_polynomial(bernstein_to_power(balance_bernstein))
    low, high = isolate_unit_root(balance_power)
    midpoint = (_decimal(low) + _decimal(high)) / 2
    all_connected_probability = evaluate_bernstein_decimal(connected_counts, midpoint)
    all_separate_probability = evaluate_bernstein_decimal(separate_counts, midpoint)
    proper_mass = Decimal(1) - all_connected_probability - all_separate_probability
    internal_degree = vertex_degrees(4, graph)[3]
    return {
        "canonical_graph_encoding": encoding,
        "gadget_sha256": _digest(gadget),
        "edges": gadget["edges"],
        "edge_count": edge_count,
        "labeled_orbit_multiplicity": multiplicity,
        "internal_degree": internal_degree,
        "primary_subset": internal_degree >= 3,
        "nonzero_terminal_partitions": len(reliability),
        "all_connected_bernstein_counts": list(connected_counts),
        "all_separate_bernstein_counts": list(separate_counts),
        "balance_bernstein_counts": list(balance_bernstein),
        "primitive_balance_power_coefficients_low_to_high": list(balance_power),
        "root_isolation": {
            "bits": ROOT_BITS,
            "lower": _fraction_text(low),
            "upper": _fraction_text(high),
        },
        "root_decimal": format(midpoint, ".50f"),
        "absolute_distance_to_square_site_reference": format(abs(midpoint - TARGET), ".50f"),
        "proper_partial_partition_mass_at_root": format(proper_mass, ".50f"),
        "balance_probability_gap_at_reported_midpoint": format(
            all_connected_probability - all_separate_probability, ".6E"
        ),
    }


def build_artifact() -> dict[str, Any]:
    getcontext().prec = 100
    catalog = connected_orbit_catalog()
    _require(len(catalog) == 11, "connected census drift")
    rows = [score_candidate(key, multiplicity, graph) for key, (multiplicity, graph) in catalog.items()]
    rows.sort(
        key=lambda row: (
            Decimal(row["absolute_distance_to_square_site_reference"]),
            row["edge_count"],
            row["canonical_graph_encoding"],
        )
    )
    for rank, row in enumerate(rows, 1):
        row["all_connected_rank"] = rank
    primary = [row for row in rows if row["primary_subset"]]
    _require(len(primary) == 4, "primary degree-filtered census drift")
    for rank, row in enumerate(primary, 1):
        row["primary_rank"] = rank
    retained = primary[0]
    closest = rows[0]
    return {
        "schema": SCHEMA,
        "issue": 13,
        "status": "exact_bounded_three_terminal_balance_screen",
        "source_commit": SOURCE_COMMIT,
        "protocol": "predictions/p13_three_terminal_balance_screen_20260830.yaml",
        "probability_model": "one shared independent bond probability p",
        "formal_balance": "P(000)-P(012)=0 in restricted-growth-string notation",
        "square_site_reference": format(TARGET, "f"),
        "candidate_counts": {
            "connected_orbits": len(rows),
            "internal_degree_at_least_3_orbits": len(primary),
        },
        "ranking": rows,
        "decision": {
            "primary_retained_for_structural_followup": retained["canonical_graph_encoding"],
            "primary_retained_root": retained["root_decimal"],
            "primary_retained_distance": retained["absolute_distance_to_square_site_reference"],
            "closest_connected_candidate": closest["canonical_graph_encoding"],
            "closest_connected_root": closest["root_decimal"],
            "closest_connected_is_primary": closest["primary_subset"],
            "exact_square_site_claims_certified": 0,
        },
        "interpretation": {
            "exact": (
                "all eleven connected canonical orbits have exact reliability polynomials and "
                "a unique formal all-connected/all-separate balance root in the open unit interval"
            ),
            "screen": (
                "the frozen target-distance ranking chooses one degree-filtered candidate for "
                "planarity, tiling, and dual-cell follow-up"
            ),
            "elimination": (
                "none of the bounded candidates is certified as an exact square-site construction "
                "by a balance-root match alone"
            ),
        },
        "claim_boundary": {
            "included": (
                "exact terminal reliability and formal three-terminal balance roots for the complete "
                "connected three-terminal plus one-internal-vertex census"
            ),
            "excluded": (
                "a square-site threshold formula, percolation bound, periodic tiling, self-dual "
                "embedding, critical manifold, or convergent gadget family"
            ),
            "parent_issue": "remain open",
        },
    }


def validate_artifact(artifact: Mapping[str, Any]) -> dict[str, Any]:
    expected = build_artifact()
    _require(artifact == expected, "balance-screen artifact does not exactly reproduce")
    _require(artifact.get("schema") == SCHEMA, "unknown schema")
    _require(artifact.get("decision", {}).get("exact_square_site_claims_certified") == 0, "claim drift")
    _require(artifact.get("claim_boundary", {}).get("parent_issue") == "remain open", "parent drift")
    return {
        "schema": SCHEMA,
        "status": "valid_exact_bounded_three_terminal_balance_screen",
        "connected_orbits": artifact["candidate_counts"]["connected_orbits"],
        "primary_orbits": artifact["candidate_counts"]["internal_degree_at_least_3_orbits"],
        "retained": artifact["decision"]["primary_retained_for_structural_followup"],
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--validate", type=Path)
    args = parser.parse_args(argv)
    if args.validate is not None:
        artifact = json.loads(args.validate.read_text(encoding="utf-8"))
        print(json.dumps(validate_artifact(artifact), indent=2, sort_keys=True))
        return 0
    artifact = build_artifact()
    rendered = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

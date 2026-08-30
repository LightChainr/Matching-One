#!/usr/bin/env python3
"""Exactly enumerate terminal-connectivity reliability polynomials."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from math import comb
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, Sequence, Tuple

try:
    from scripts.terminal_partition_canonical import blocks_to_rgs, enumerate_rgs
except ModuleNotFoundError:  # Direct ``python scripts/...`` execution.
    from terminal_partition_canonical import blocks_to_rgs, enumerate_rgs


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results" / "terminal-reliability" / "star4" / "latest.json"
SCHEMA = "matching-one/terminal-reliability-polynomial/v1"
Edge = Tuple[int, int]


STAR4 = {
    "id": "four-terminal-independent-bond-star",
    "vertex_count": 5,
    "terminal_count": 4,
    "edges": [[0, 4], [1, 4], [2, 4], [3, 4]],
    "edge_probability": "one shared symbolic p; edges mutually independent",
}


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def validate_gadget(gadget: Mapping[str, Any]) -> tuple[int, int, Tuple[Edge, ...]]:
    vertex_count = gadget.get("vertex_count")
    terminal_count = gadget.get("terminal_count")
    _require(type(vertex_count) is int and vertex_count >= 2, "invalid vertex count")
    _require(type(terminal_count) is int and 2 <= terminal_count <= vertex_count, "invalid terminal count")
    edges = []
    for raw_edge in gadget.get("edges", []):
        _require(isinstance(raw_edge, list) and len(raw_edge) == 2, "edges must be endpoint pairs")
        u, v = raw_edge
        _require(type(u) is int and type(v) is int, "edge endpoints must be integers")
        _require(0 <= u < vertex_count and 0 <= v < vertex_count, "edge endpoint out of range")
        _require(u != v, "self-loops are forbidden")
        edges.append((min(u, v), max(u, v)))
    _require(bool(edges), "gadget must contain at least one edge")
    _require(len(edges) == len(set(edges)), "parallel/duplicate edges are forbidden")
    _require(len(edges) <= 24, "exhaustive evaluator is limited to 24 edges")
    return vertex_count, terminal_count, tuple(sorted(edges))


def _terminal_partition(vertex_count: int, terminal_count: int, open_edges: Iterable[Edge]) -> tuple[int, ...]:
    parent = list(range(vertex_count))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        root_a, root_b = find(a), find(b)
        if root_a != root_b:
            parent[root_b] = root_a

    for u, v in open_edges:
        union(u, v)
    blocks: dict[int, list[int]] = {}
    for terminal in range(terminal_count):
        blocks.setdefault(find(terminal), []).append(terminal)
    return blocks_to_rgs(blocks.values(), terminal_count)


def enumerate_reliability(gadget: Mapping[str, Any]) -> Mapping[tuple[int, ...], tuple[int, ...]]:
    """Return exact Bernstein counts indexed by terminal RGS and open-edge count."""

    vertex_count, terminal_count, edges = validate_gadget(gadget)
    counts = {partition: [0] * (len(edges) + 1) for partition in enumerate_rgs(terminal_count)}
    for mask in range(1 << len(edges)):
        open_edges = tuple(edge for index, edge in enumerate(edges) if mask & (1 << index))
        partition = _terminal_partition(vertex_count, terminal_count, open_edges)
        counts[partition][len(open_edges)] += 1
    return {partition: tuple(values) for partition, values in counts.items() if any(values)}


def evaluate_at(counts: Sequence[int], p: Fraction) -> Fraction:
    edge_count = len(counts) - 1
    return sum(count * p**opened * (1 - p) ** (edge_count - opened) for opened, count in enumerate(counts))


def _gadget_digest(gadget: Mapping[str, Any]) -> str:
    payload = json.dumps(gadget, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def build_star4_result() -> dict[str, Any]:
    counts = enumerate_reliability(STAR4)
    edge_count = len(STAR4["edges"])
    rows = []
    for partition, coefficients in sorted(counts.items()):
        probability_half = evaluate_at(coefficients, Fraction(1, 2))
        rows.append(
            {
                "terminal_partition_rgs": list(partition),
                "bernstein_counts_by_open_edge_count": list(coefficients),
                "probability_at_one_half": "%d/%d" % (probability_half.numerator, probability_half.denominator),
            }
        )
    normalization = [sum(row["bernstein_counts_by_open_edge_count"][k] for row in rows) for k in range(edge_count + 1)]
    expected_normalization = [comb(edge_count, k) for k in range(edge_count + 1)]
    _require(normalization == expected_normalization, "partition probabilities do not normalize")
    return {
        "schema": SCHEMA,
        "issue": 14,
        "status": "exact_small_gadget_baseline",
        "gadget": STAR4,
        "gadget_sha256": _gadget_digest(STAR4),
        "probability_basis": "sum_k count[k] * p^k * (1-p)^(m-k)",
        "enumeration": {
            "edge_count": edge_count,
            "configurations": 1 << edge_count,
            "nonzero_terminal_partitions": len(rows),
            "terminal_partition_encoding": "restricted_growth_string",
            "arithmetic": "integer Bernstein counts and fractions.Fraction evaluations",
        },
        "terminal_partition_polynomials": rows,
        "normalization_counts": normalization,
        "endpoint_checks": {
            "p_zero_partition": [0, 1, 2, 3],
            "p_one_partition": [0, 0, 0, 0],
            "both_have_probability_one": True,
        },
        "conclusion": {
            "exact_result": "all 16 bond configurations of the four-terminal star are partitioned into exact terminal-connectivity Bernstein polynomials",
            "new_percolation_bound": False,
            "theorem_claim": False,
        },
        "claim_boundary": {
            "included": "exact terminal-connectivity polynomials for one declared four-terminal independent-bond star",
            "excluded": "baseline bound reproduction, broad enumeration, optimization, domination certification, theorem assumptions, or a new percolation bound",
            "parent_issue": "remain open",
        },
    }


def validate_result(result: Mapping[str, Any]) -> dict[str, Any]:
    expected = build_star4_result()
    _require(result == expected, "result does not exactly reproduce the declared gadget")
    _require(result.get("conclusion", {}).get("new_percolation_bound") is False, "bound overclaim")
    _require(result.get("conclusion", {}).get("theorem_claim") is False, "theorem overclaim")
    _require(result.get("claim_boundary", {}).get("parent_issue") == "remain open", "parent boundary drift")
    return {
        "schema": result["schema"],
        "status": "valid_exact_small_gadget_baseline",
        "configurations": result["enumeration"]["configurations"],
        "nonzero_terminal_partitions": result["enumeration"]["nonzero_terminal_partitions"],
        "normalization_counts": result["normalization_counts"],
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--validate", type=Path)
    args = parser.parse_args(argv)
    if args.validate is not None:
        result = json.loads(args.validate.read_text(encoding="utf-8"))
        print(json.dumps(validate_result(result), indent=2, sort_keys=True))
        return 0
    result = build_star4_result()
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

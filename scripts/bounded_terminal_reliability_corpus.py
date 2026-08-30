#!/usr/bin/env python3
"""Build exact reliability polynomials for the bounded four-terminal census."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from math import comb
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

try:
    from scripts.bounded_gadget_census import is_connected, validate_artifact, vertex_degrees
    from scripts.gadget_graph_canonical import canonical_graph, decode_graph, enumerate_graphs
    from scripts.terminal_partition_canonical import full_symmetric_group
    from scripts.terminal_reliability_polynomial import STAR4, build_star4_result, enumerate_reliability
except ModuleNotFoundError:  # Direct ``python scripts/...`` execution.
    from bounded_gadget_census import is_connected, validate_artifact, vertex_degrees
    from gadget_graph_canonical import canonical_graph, decode_graph, enumerate_graphs
    from terminal_partition_canonical import full_symmetric_group
    from terminal_reliability_polynomial import STAR4, build_star4_result, enumerate_reliability


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CENSUS = ROOT / "analysis" / "bounded_gadget_census.json"
DEFAULT_OUTPUT = ROOT / "results" / "terminal-reliability" / "bounded-four-terminal-corpus.json"
STAR_FIXTURE = ROOT / "results" / "terminal-reliability" / "star4" / "latest.json"
SCHEMA = "matching-one/bounded-terminal-reliability-corpus/v1"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _sha256_file(path: Path) -> str:
    return _sha256_bytes(path.read_bytes())


def selected_canonical_graphs() -> tuple[str, ...]:
    terminal_count = 4
    vertex_count = 5
    internal_vertex = 4
    group = full_symmetric_group(terminal_count)
    keys = {
        canonical_graph(vertex_count, terminal_count, graph, group)
        for graph in enumerate_graphs(vertex_count, terminal_count)
    }
    selected = []
    for key in sorted(keys):
        _, _, graph = decode_graph(key)
        if is_connected(vertex_count, graph) and vertex_degrees(vertex_count, graph)[internal_vertex] >= 3:
            selected.append(key)
    return tuple(selected)


def _polynomial_rows(counts: Mapping[tuple[int, ...], tuple[int, ...]]) -> list[dict[str, Any]]:
    return [
        {
            "terminal_partition_rgs": list(partition),
            "bernstein_counts_by_open_edge_count": list(coefficients),
        }
        for partition, coefficients in sorted(counts.items())
    ]


def build_candidate(canonical_encoding: str) -> dict[str, Any]:
    terminal_count, vertex_count, graph = decode_graph(canonical_encoding)
    _require((terminal_count, vertex_count) == (4, 5), "candidate left the four-terminal census")
    edge_count = len(graph)
    gadget = {
        "vertex_count": vertex_count,
        "terminal_count": terminal_count,
        "edges": [list(edge) for edge in graph],
    }
    counts = enumerate_reliability(gadget)
    polynomials = _polynomial_rows(counts)
    normalization = [
        sum(row["bernstein_counts_by_open_edge_count"][opened] for row in polynomials)
        for opened in range(edge_count + 1)
    ]
    _require(normalization == [comb(edge_count, opened) for opened in range(edge_count + 1)], "normalization drift")
    return {
        "canonical_graph_encoding": canonical_encoding,
        "canonical_graph_sha256": _sha256_bytes(canonical_encoding.encode("utf-8")),
        "edge_count": edge_count,
        "internal_degree": vertex_degrees(vertex_count, graph)[4],
        "configurations": 1 << edge_count,
        "nonzero_terminal_partitions": len(polynomials),
        "normalization_counts": normalization,
        "terminal_partition_polynomials": polynomials,
        "canonical_polynomial_signature_sha256": _sha256_bytes(_canonical_bytes(polynomials)),
    }


def _star_corpus_polynomials() -> list[dict[str, Any]]:
    return _polynomial_rows(enumerate_reliability(STAR4))


def build_result(census_path: Path = DEFAULT_CENSUS) -> dict[str, Any]:
    census_bytes = census_path.read_bytes()
    census = json.loads(census_bytes)
    validate_artifact(census)
    four_terminal = next(row for row in census["rows"] if row["terminal_count"] == 4)
    keys = selected_canonical_graphs()
    _require(len(keys) == four_terminal["connected_internal_degree_at_least_3_orbits"], "census filter count drift")
    candidates = [build_candidate(key) for key in keys]
    edge_histogram = Counter(row["edge_count"] for row in candidates)

    star_polynomials = _star_corpus_polynomials()
    star_matches = [row for row in candidates if row["terminal_partition_polynomials"] == star_polynomials]
    _require(len(star_matches) == 1, "four-star fixture must match exactly one canonical candidate")
    checked_star = json.loads(STAR_FIXTURE.read_text(encoding="utf-8"))
    _require(checked_star == build_star4_result(), "checked-in four-star fixture drift")

    signatures = Counter(row["canonical_polynomial_signature_sha256"] for row in candidates)
    corpus_digest = _sha256_bytes(_canonical_bytes(candidates))
    return {
        "schema": SCHEMA,
        "issue": 14,
        "status": "exact_bounded_candidate_reliability_corpus",
        "sources": {
            "census_path": str(census_path.relative_to(ROOT)),
            "census_sha256": _sha256_bytes(census_bytes),
            "star_fixture_path": str(STAR_FIXTURE.relative_to(ROOT)),
            "star_fixture_sha256": _sha256_file(STAR_FIXTURE),
        },
        "selection": {
            "terminal_count": 4,
            "internal_count": 1,
            "filter": "connected carrier and sole internal vertex degree at least three",
            "candidate_orbits": len(candidates),
            "edge_count_histogram": {str(k): v for k, v in sorted(edge_histogram.items())},
        },
        "enumeration": {
            "total_independent_edge_configurations": sum(row["configurations"] for row in candidates),
            "total_nonzero_terminal_partition_polynomials": sum(row["nonzero_terminal_partitions"] for row in candidates),
            "distinct_byte_identical_canonical_polynomial_signatures": len(signatures),
            "signature_multiplicity_histogram": {
                str(size): count for size, count in sorted(Counter(signatures.values()).items())
            },
            "star_fixture_matches": len(star_matches),
            "corpus_sha256": corpus_digest,
        },
        "probability_basis": "sum_k count[k] * p^k * (1-p)^(m-k)",
        "candidates": candidates,
        "conclusion": {
            "exact_result": "every selected bounded census representative has a normalized exact terminal-partition reliability polynomial",
            "planarity_certified": False,
            "self_duality_tested": False,
            "new_percolation_bound": False,
            "theorem_claim": False,
        },
        "claim_boundary": {
            "included": "exact reliability-polynomial corpus for the 27 frozen connected/internal-degree>=3 four-terminal representatives",
            "excluded": "planarity, periodic tiling, reductions, domination, optimization, self-duality, critical manifolds, baseline bounds, theorem assumptions, or a new bound",
            "parent_issue": "remain open",
        },
    }


def validate_result(result: Mapping[str, Any], census_path: Path = DEFAULT_CENSUS) -> Mapping[str, Any]:
    expected = build_result(census_path)
    _require(result == expected, "bounded reliability corpus does not exactly reproduce")
    _require(result.get("claim_boundary", {}).get("parent_issue") == "remain open", "parent boundary drift")
    return {
        "schema": result["schema"],
        "status": "valid_exact_bounded_candidate_reliability_corpus",
        "candidate_orbits": result["selection"]["candidate_orbits"],
        "total_configurations": result["enumeration"]["total_independent_edge_configurations"],
        "polynomial_rows": result["enumeration"]["total_nonzero_terminal_partition_polynomials"],
        "distinct_signatures": result["enumeration"]["distinct_byte_identical_canonical_polynomial_signatures"],
        "corpus_sha256": result["enumeration"]["corpus_sha256"],
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--census", type=Path, default=DEFAULT_CENSUS)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--validate", type=Path)
    args = parser.parse_args(argv)
    if args.validate:
        result = json.loads(args.validate.read_text(encoding="utf-8"))
        print(json.dumps(validate_result(result, args.census), indent=2, sort_keys=True))
        return 0
    result = build_result(args.census)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

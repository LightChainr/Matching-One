#!/usr/bin/env python3
"""Generate and hash deterministic small-width planar state transitions."""

from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

try:
    from scripts.noncrossing_connectivity_codec import noncrossing_states
    from scripts.planar_state_operations import (
        detach_blocks,
        detach_rgs,
        join_cyclic_adjacent_blocks,
        join_cyclic_adjacent_rgs,
    )
except ModuleNotFoundError:  # Direct ``python scripts/...`` execution.
    from noncrossing_connectivity_codec import noncrossing_states
    from planar_state_operations import (
        detach_blocks,
        detach_rgs,
        join_cyclic_adjacent_blocks,
        join_cyclic_adjacent_rgs,
    )


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "analysis" / "planar_transition_table_contract.json"
SCHEMA = "matching-one/planar-transition-table/v1"
OPERATIONS = (
    ("detach", detach_rgs, detach_blocks),
    ("join_cyclic_adjacent", join_cyclic_adjacent_rgs, join_cyclic_adjacent_blocks),
)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _record_bytes(record: Mapping[str, Any]) -> bytes:
    return (json.dumps(record, sort_keys=True, separators=(",", ":")) + "\n").encode("utf-8")


def build_width_table(width: int) -> tuple[dict[str, Any], bytes]:
    _require(type(width) is int and 1 <= width <= 8, "width must be in [1,8]")
    states = noncrossing_states(width)
    ranks = {state: rank for rank, state in enumerate(states)}
    serialized = bytearray()
    changed = Counter()
    operation_edges = {name: set() for name, _, _ in OPERATIONS}
    row_targets = [set() for _ in states]

    for source_rank, state in enumerate(states):
        for point in range(width):
            for operation, implementation, independent in OPERATIONS:
                target = implementation(state, point)
                _require(target == independent(state, point), "independent transition implementations disagree")
                _require(target in ranks, "transition target left the canonical state space")
                target_rank = ranks[target]
                record = {
                    "operation": operation,
                    "point": point,
                    "source_rank": source_rank,
                    "target_rank": target_rank,
                    "width": width,
                }
                serialized.extend(_record_bytes(record))
                changed[operation] += target_rank != source_rank
                operation_edges[operation].add((source_rank, target_rank))
                row_targets[source_rank].add(target_rank)

    row_degree_histogram = Counter(len(targets) for targets in row_targets)
    cases_per_operation = len(states) * width
    summary = {
        "width": width,
        "states": len(states),
        "cases_per_operation": cases_per_operation,
        "serialized_records": len(OPERATIONS) * cases_per_operation,
        "changed_cases": {name: changed[name] for name, _, _ in OPERATIONS},
        "unique_source_target_edges": {
            name: len(operation_edges[name]) for name, _, _ in OPERATIONS
        },
        "row_unique_target_degree_histogram": {
            str(degree): count for degree, count in sorted(row_degree_histogram.items())
        },
        "canonical_jsonl_sha256": hashlib.sha256(serialized).hexdigest(),
    }
    return summary, bytes(serialized)


def build_result(maximum_width: int = 8) -> dict[str, Any]:
    _require(type(maximum_width) is int and 1 <= maximum_width <= 8, "maximum width must be in [1,8]")
    width_summaries = []
    complete = bytearray()
    for width in range(1, maximum_width + 1):
        summary, serialized = build_width_table(width)
        width_summaries.append(summary)
        complete.extend(serialized)
    return {
        "schema": SCHEMA,
        "issue": 11,
        "status": "exact_small_width_transition_table_only",
        "operation_order": [name for name, _, _ in OPERATIONS],
        "serialization": {
            "format": "UTF-8 JSON Lines",
            "json": "sorted keys and compact separators",
            "record_order": "width, source_rank, point, operation_order",
            "state_rank_source": "lexicographically enumerated canonical noncrossing RGS",
        },
        "maximum_width": maximum_width,
        "widths": width_summaries,
        "totals": {
            "states": sum(row["states"] for row in width_summaries),
            "cases_per_operation": sum(row["cases_per_operation"] for row in width_summaries),
            "serialized_records": sum(row["serialized_records"] for row in width_summaries),
            "canonical_jsonl_sha256": hashlib.sha256(complete).hexdigest(),
        },
        "validation": {
            "independent_block_set_targets_agree": True,
            "all_targets_have_canonical_ranks": True,
            "content_hash_algorithm": "SHA-256",
        },
        "claim_boundary": {
            "included": "exact width-1..8 detach/join transition records, rank targets, counts, row degrees, and content hashes",
            "excluded": "transfer weights, row propagation, topological sectors, eigensolvers, published widths, resource extrapolation, checkpointing, or frontier extension",
            "parent_issue": "remain open",
        },
    }


def validate_result(result: Mapping[str, Any]) -> Mapping[str, Any]:
    _require(result.get("schema") == SCHEMA, "unknown schema")
    expected = build_result(result.get("maximum_width"))
    _require(result == expected, "transition-table contract does not exactly reproduce")
    _require(result.get("claim_boundary", {}).get("parent_issue") == "remain open", "parent boundary drift")
    return {
        "schema": result["schema"],
        "status": "valid_exact_small_width_transition_table",
        "maximum_width": result["maximum_width"],
        "states": result["totals"]["states"],
        "cases_per_operation": result["totals"]["cases_per_operation"],
        "serialized_records": result["totals"]["serialized_records"],
        "canonical_jsonl_sha256": result["totals"]["canonical_jsonl_sha256"],
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--maximum-width", type=int, default=8)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--validate", type=Path)
    args = parser.parse_args(argv)
    if args.validate:
        result = json.loads(args.validate.read_text(encoding="utf-8"))
        print(json.dumps(validate_result(result), indent=2, sort_keys=True))
        return 0
    result = build_result(args.maximum_width)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

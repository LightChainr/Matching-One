#!/usr/bin/env python3
"""Strictly merge deterministic corrected-reservoir scan shards.

This tool changes no scan semantics.  It validates a complete set of outputs
from ``p334_tm_corrected_reservoir_scan.py``, orders rows by their frozen
``row_index``, recomputes the aggregate summary, and emits a bounded complete-
order certificate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence


SOURCE_SCHEMA = "p334-tm-corrected-reservoir-scan-v1"
MERGED_SCHEMA = "p334-tm-corrected-reservoir-merge-v1"
SOURCE_BOUNDARY = (
    "A complete order certificate requires every deterministic shard; "
    "saturation is a bounded exact result, not the general HNF theorem."
)
MERGED_BOUNDARY = (
    "All declared deterministic shards and row indices are present.  This is "
    "a complete bounded order certificate, not the general HNF theorem."
)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _row_summary(rows: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    deficiencies = []
    saturated = 0
    for row in rows:
        audit = row.get("corrected_reservoir_orbit_graph")
        _require(isinstance(audit, dict), "row is missing corrected reservoir audit")
        deficiency = audit.get("Hall_deficiency")
        _require(isinstance(deficiency, int) and deficiency >= 0, "invalid Hall deficiency")
        saturates = audit.get("saturates")
        _require(isinstance(saturates, bool), "row saturation flag must be boolean")
        _require(saturates == (deficiency == 0), "row saturation and deficiency disagree")
        deficiencies.append(deficiency)
        saturated += int(saturates)
    failed = len(rows) - saturated
    return {
        "saturated_rows": saturated,
        "failed_rows": failed,
        "minimum_Hall_deficiency": min(deficiencies, default=None),
        "maximum_Hall_deficiency": max(deficiencies, default=None),
        "status": "counterexample_found" if failed else "all_rows_saturate",
    }


def _validate_source_summary(shard: Mapping[str, Any], rows: Sequence[Mapping[str, Any]]) -> None:
    expected = _row_summary(rows)
    declared = shard.get("summary")
    _require(isinstance(declared, dict), "shard summary must be an object")
    source_expected = {
        **expected,
        "status": (
            "counterexample_found"
            if expected["failed_rows"]
            else "all_selected_rows_saturate"
        ),
    }
    _require(declared == source_expected, "shard summary drifted from its rows")


def merge_shards(
    shards: Sequence[Mapping[str, Any]],
    *,
    provenance: Sequence[Mapping[str, str]] | None = None,
) -> dict[str, Any]:
    _require(bool(shards), "at least one shard is required")
    first = shards[0]
    _require(first.get("schema") == SOURCE_SCHEMA, "source schema drifted")
    order = first.get("order")
    _require(isinstance(order, int) and order >= 1, "order must be a positive integer")
    first_shard = first.get("shard")
    _require(isinstance(first_shard, dict), "shard descriptor must be an object")
    shard_count = first_shard.get("count")
    _require(isinstance(shard_count, int) and shard_count >= 1, "shard count must be positive")
    complete_rows = first.get("complete_order_row_count")
    _require(
        isinstance(complete_rows, int) and complete_rows >= 0,
        "complete row count must be nonnegative",
    )

    by_index: dict[int, Mapping[str, Any]] = {}
    all_rows: list[Mapping[str, Any]] = []
    for shard in shards:
        _require(shard.get("schema") == SOURCE_SCHEMA, "source schema drifted")
        _require(shard.get("order") == order, "shard orders disagree")
        descriptor = shard.get("shard")
        _require(isinstance(descriptor, dict), "shard descriptor must be an object")
        _require(descriptor.get("count") == shard_count, "shard counts disagree")
        index = descriptor.get("index")
        _require(
            isinstance(index, int) and 0 <= index < shard_count,
            "shard index is out of range",
        )
        _require(index not in by_index, "duplicate shard index")
        by_index[index] = shard
        _require(
            shard.get("complete_order_row_count") == complete_rows,
            "complete row counts disagree",
        )
        _require(
            shard.get("scientific_boundary") == SOURCE_BOUNDARY,
            "source scientific boundary drifted",
        )
        rows = shard.get("rows")
        _require(isinstance(rows, list), "shard rows must be a list")
        _require(shard.get("selected_row_count") == len(rows), "selected row count drifted")
        row_indices = []
        for row in rows:
            _require(isinstance(row, dict), "each row must be an object")
            row_index = row.get("row_index")
            _require(isinstance(row_index, int), "row index must be an integer")
            _require(
                row_index % shard_count == index,
                "row index does not belong to its deterministic shard",
            )
            row_indices.append(row_index)
        _require(
            row_indices == sorted(row_indices),
            "rows within a shard must be in deterministic order",
        )
        expected_indices = list(range(index, complete_rows, shard_count))
        _require(
            row_indices == expected_indices,
            "shard row indices do not exactly cover their deterministic residue class",
        )
        _validate_source_summary(shard, rows)
        all_rows.extend(rows)

    _require(len(by_index) == shard_count, "declared shard set is incomplete")
    _require(
        sorted(by_index) == list(range(shard_count)),
        "shard indices do not give complete coverage",
    )
    sorted_rows = sorted(all_rows, key=lambda row: row["row_index"])
    merged_indices = [row["row_index"] for row in sorted_rows]
    _require(
        merged_indices == list(range(complete_rows)),
        "merged row indices are duplicated, missing, or out of range",
    )

    if provenance is not None:
        _require(len(provenance) == len(shards), "provenance length disagrees with shards")
        provenance_by_index = []
        for shard, source in zip(shards, provenance):
            index = shard["shard"]["index"]
            provenance_by_index.append({"shard_index": index, **source})
        provenance_by_index.sort(key=lambda row: row["shard_index"])
    else:
        provenance_by_index = []

    return {
        "schema": MERGED_SCHEMA,
        "source_schema": SOURCE_SCHEMA,
        "order": order,
        "shard_count": shard_count,
        "merged_shard_indices": sorted(by_index),
        "complete_order_row_count": complete_rows,
        "rows": sorted_rows,
        "summary": _row_summary(sorted_rows),
        "input_provenance": provenance_by_index,
        "source_scientific_boundary": SOURCE_BOUNDARY,
        "scientific_boundary": MERGED_BOUNDARY,
    }


def load_and_merge(paths: Sequence[Path]) -> dict[str, Any]:
    shards = []
    provenance = []
    for path in paths:
        raw = path.read_bytes()
        shards.append(json.loads(raw))
        provenance.append(
            {
                "path": str(path),
                "sha256": hashlib.sha256(raw).hexdigest(),
            }
        )
    return merge_shards(shards, provenance=provenance)


def render_markdown(result: Mapping[str, Any]) -> str:
    summary = result["summary"]
    lines = [
        f"# P334 corrected-reservoir merged N={result['order']} certificate",
        "",
        (
            f"Merged all `{result['shard_count']}` deterministic shards and "
            f"`{result['complete_order_row_count']}` rows."
        ),
        "",
        (
            f"Status: **{summary['status']}**; saturated "
            f"`{summary['saturated_rows']}`, failed `{summary['failed_rows']}`."
        ),
        "",
        "| row | HNF | Smith | carrier | line | layer | orbit matching |",
        "|---:|---|---|---|---|---:|---:|",
    ]
    for row in result["rows"]:
        audit = row["corrected_reservoir_orbit_graph"]
        lines.append(
            "| {row_index} | `{matrix}` | `{smith}` | {carrier} | `{line}` | "
            "{layer} | {matching}/{sources} |".format(
                row_index=row["row_index"],
                matrix=row["matrix"],
                smith=row["Smith_invariants"],
                carrier=row["carrier"],
                line=row["line"],
                layer=row["lower_layer"],
                matching=audit["maximum_matching"],
                sources=audit["source_tokens"],
            )
        )
    lines.extend(["", result["scientific_boundary"], ""])
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("shards", type=Path, nargs="+")
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args()
    result = load_and_merge(args.shards)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.markdown.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    args.markdown.write_text(render_markdown(result), encoding="utf-8")


if __name__ == "__main__":
    main()

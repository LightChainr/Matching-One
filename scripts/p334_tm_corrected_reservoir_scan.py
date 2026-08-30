#!/usr/bin/env python3
"""Scan the corrected P334 transport reservoir on complete HNF shells.

The expensive compatibility graph is built only after exact translation-orbit
normalization.  Independent shards are selected by the deterministic order of
``hard_rows`` so the same frozen observable can be distributed across hosts.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from p334_tm_translation_orbit_hall import (
    descriptor,
    hard_rows,
    orbit_graph_audit,
    transverse_reservoir_targets,
)
from p334_tm_configuration_cross_switch import translation_permutations


SCHEMA = "p334-tm-corrected-reservoir-scan-v1"


def row_key(row) -> tuple:
    n, matrix, _geometry, carrier, _marks, line, lower_layer, _faces = row
    return (n, matrix, carrier, line, lower_layer)


def rows_for_order(order: int):
    return sorted(
        (row for row in hard_rows(order) if row[0] == order),
        key=row_key,
    )


def scan_order(order: int, *, shard_index: int = 0, shard_count: int = 1):
    if order < 1:
        raise ValueError("order must be positive")
    if shard_count < 1:
        raise ValueError("shard_count must be positive")
    if not 0 <= shard_index < shard_count:
        raise ValueError("shard_index must lie in [0, shard_count)")

    all_rows = rows_for_order(order)
    selected = [
        (index, row)
        for index, row in enumerate(all_rows)
        if index % shard_count == shard_index
    ]
    reports = []
    for index, row in selected:
        n, matrix, geometry, carrier, marks, line, lower_layer, faces = row
        permutations = translation_permutations(geometry)
        sources = [
            (replica, coexit, flat)
            for replica in range(4)
            for coexit in faces["D"]
            for flat in faces["F"]
        ]
        audit = orbit_graph_audit(
            marks,
            line,
            sources,
            permutations,
            n,
            lambda source: transverse_reservoir_targets(
                marks,
                line,
                source,
                permutations,
                n,
                transport=True,
            ),
        )
        reports.append(
            {
                "row_index": index,
                **descriptor(n, matrix, carrier, line, lower_layer, faces),
                "corrected_reservoir_orbit_graph": audit,
            }
        )

    failed = [
        report
        for report in reports
        if not report["corrected_reservoir_orbit_graph"]["saturates"]
    ]
    return {
        "schema": SCHEMA,
        "order": order,
        "shard": {"index": shard_index, "count": shard_count},
        "complete_order_row_count": len(all_rows),
        "selected_row_count": len(reports),
        "rows": reports,
        "summary": {
            "saturated_rows": len(reports) - len(failed),
            "failed_rows": len(failed),
            "minimum_Hall_deficiency": min(
                (
                    report["corrected_reservoir_orbit_graph"][
                        "Hall_deficiency"
                    ]
                    for report in reports
                ),
                default=None,
            ),
            "maximum_Hall_deficiency": max(
                (
                    report["corrected_reservoir_orbit_graph"][
                        "Hall_deficiency"
                    ]
                    for report in reports
                ),
                default=None,
            ),
            "status": (
                "counterexample_found"
                if failed
                else "all_selected_rows_saturate"
            ),
        },
        "scientific_boundary": (
            "A complete order certificate requires every deterministic shard; "
            "saturation is a bounded exact result, not the general HNF theorem."
        ),
    }


def render_markdown(result) -> str:
    summary = result["summary"]
    lines = [
        f"# P334 corrected-reservoir HNF scan at N={result['order']}",
        "",
        (
            f"Shard `{result['shard']['index']}/{result['shard']['count']}` "
            f"contains `{result['selected_row_count']}` of "
            f"`{result['complete_order_row_count']}` deterministic hard rows."
        ),
        "",
        f"Status: **{summary['status']}**; saturated "
        f"`{summary['saturated_rows']}`, failed `{summary['failed_rows']}`.",
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
    parser = argparse.ArgumentParser()
    parser.add_argument("--order", type=int, required=True)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args()
    result = scan_order(
        args.order,
        shard_index=args.shard_index,
        shard_count=args.shard_count,
    )
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.markdown.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    args.markdown.write_text(render_markdown(result), encoding="utf-8")


if __name__ == "__main__":
    main()

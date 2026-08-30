#!/usr/bin/env python3
"""Type-safe entrypoint for the frozen threshold-rank root-doubling score."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Callable, Mapping, Sequence

from wrapping_channels import ObservableDescriptor, map_observable


SEMANTIC_GATE = "predictions/threshold_rank_root_doubling_semantic_gate_20260830.json"
SIZE_ORDER = [65, 85, 130, 170]
LINEAGES = [[65, 130], [85, 170]]
COVARIANCE_MODES = [
    "full_cross_size_covariance",
    "diagonal_cross_size_covariance",
]


def load_semantic_gate(root: Path) -> tuple[dict, ObservableDescriptor, ObservableDescriptor, object]:
    gate = json.loads((root / SEMANTIC_GATE).read_text(encoding="utf-8"))
    if gate.get("status") != "semantic_gate_added_after_frozen_threshold_root_doubling_score":
        raise ValueError("threshold root-doubling semantic gate status changed")
    if gate.get("frozen_kernel_git_blob") != "463acbb63df93eb58792db5a7db012d8e6096daa":
        raise ValueError("threshold root-doubling frozen kernel identity changed")
    if gate.get("size_order") != SIZE_ORDER or gate.get("lineages_in_order") != LINEAGES:
        raise ValueError("threshold root-doubling lineage order changed")
    if gate.get("genealogy") != "multiplication_by_1_plus_i":
        raise ValueError("threshold root-doubling genealogy changed")
    if gate.get("paired_quantity") != "threshold_rank_root_gap":
        raise ValueError("threshold root-doubling quantity changed")
    if gate.get("target_ratio") != {
        "numerator": -1, "denominator": 4, "fitted": False
    }:
        raise ValueError("threshold root-doubling target ratio changed")
    if gate.get("covariance_modes_in_order") != COVARIANCE_MODES:
        raise ValueError("threshold root-doubling covariance modes changed")
    source = ObservableDescriptor.from_dict(gate["source_descriptor"])
    stored_child = ObservableDescriptor.from_dict(gate["stored_child_descriptor"])
    transform = map_observable(source, stored_child)
    expected = gate["exact_registered_stored_child_map"]
    if (transform.scale, transform.offset) != (
        float(expected["scale"]), float(expected["offset"])
    ) or (transform.scale, transform.offset) != (-1.0, 0.0):
        raise ValueError("threshold root-doubling child sign map changed")
    return gate, source, stored_child, transform


def _score_frozen(summary: Mapping[str, object]) -> dict:
    import score_threshold_rank_root_doubling as frozen_kernel

    return frozen_kernel.score(summary)


def score_typed(
    root: Path,
    summary: Mapping[str, object],
    *,
    runner: Callable[[Mapping[str, object]], dict] = _score_frozen,
) -> dict:
    gate, source, stored_child, transform = load_semantic_gate(root)
    result = runner(summary)
    if list(result)[:2] != gate["covariance_modes_in_order"]:
        raise ValueError("threshold root-doubling covariance mode order differs from gate")
    for mode in gate["covariance_modes_in_order"]:
        block = result.get(mode, {})
        if block.get("target_ratio") != -0.25:
            raise ValueError("threshold root-doubling result target ratio differs from gate")
        actual = [
            [row.get("parent_N"), row.get("child_N")]
            for row in block.get("lineages", [])
        ]
        if actual != gate["lineages_in_order"]:
            raise ValueError("threshold root-doubling result lineage order differs from gate")
    result["observable_semantics"] = {
        "semantic_gate": SEMANTIC_GATE,
        "semantic_gate_status": gate["status"],
        "source_descriptor": source.to_dict(),
        "stored_child_descriptor": stored_child.to_dict(),
        "applied_stored_child_transform": transform.to_dict(),
        "paired_quantity": gate["paired_quantity"],
        "genealogy": gate["genealogy"],
        "validation_order": "semantic_sign_map_before_frozen_covariance_score",
        "evidence_boundary": gate["evidence_boundary"],
    }
    return result


def main(argv: Sequence[str] | None = None) -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args(argv)
    summary = json.loads(args.summary.read_text(encoding="utf-8"))
    result = score_typed(root, summary)
    import score_threshold_rank_root_doubling as frozen_kernel

    frozen_kernel.write_outputs(result, args.json, args.csv, args.report)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Type-safe entrypoint for the frozen Issue #50 N=290 scorer."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Mapping, Sequence, Tuple

import score_issue50_n290 as frozen_kernel
from wrapping_channels import ObservableDescriptor, map_observable


SEMANTIC_MANIFEST = "predictions/issue50_n290_semantic_gate_20260830.json"


def load_semantic_gate(
    root: Path,
) -> tuple[dict, ObservableDescriptor, ObservableDescriptor, object]:
    payload = json.loads((root / SEMANTIC_MANIFEST).read_text(encoding="utf-8"))
    if payload.get("status") != "semantic_gate_added_after_frozen_issue50_score":
        raise ValueError("Issue #50 semantic gate status changed")
    if payload.get("N") != frozen_kernel.EXPECTED_N:
        raise ValueError("Issue #50 semantic gate N changed")
    if tuple(payload.get("lineage_first", ())) != frozen_kernel.EXPECTED_FIRST:
        raise ValueError("Issue #50 first lineage changed")
    if tuple(payload.get("lineage_second", ())) != frozen_kernel.EXPECTED_SECOND:
        raise ValueError("Issue #50 second lineage changed")
    source = ObservableDescriptor.from_dict(payload["source_descriptor"])
    target = ObservableDescriptor.from_dict(payload["target_descriptor"])
    transform = map_observable(source, target)
    expected = payload["exact_registered_map"]
    if (transform.scale, transform.offset) != (
        float(expected["scale"]),
        float(expected["offset"]),
    ):
        raise ValueError("registered Issue #50 map differs from the semantic gate")
    if (transform.scale, transform.offset) != (1.0, 0.0):
        raise ValueError("Issue #50 source and target contrast must use exact identity")
    return payload, source, target, transform


def score_typed(
    root: Path,
    rows: Mapping[Tuple[str, int], Mapping[str, int]],
    run: Mapping[str, object],
) -> dict[str, object]:
    gate, source, target, transform = load_semantic_gate(root)
    result = frozen_kernel.score(rows, run)
    if result.get("channel") != source.channel.value:
        raise ValueError("Issue #50 kernel channel differs from the semantic gate")
    if result.get("sector") != "matching_function":
        raise ValueError("Issue #50 kernel sector differs from matching-odd")
    if tuple(result.get("lineage_first", ())) != tuple(gate["lineage_first"]):
        raise ValueError("Issue #50 kernel first lineage differs from the semantic gate")
    if tuple(result.get("lineage_second", ())) != tuple(gate["lineage_second"]):
        raise ValueError("Issue #50 kernel second lineage differs from the semantic gate")
    result["observable_semantics"] = {
        "semantic_manifest": SEMANTIC_MANIFEST,
        "semantic_manifest_status": gate["status"],
        "source_descriptor": source.to_dict(),
        "target_descriptor": target.to_dict(),
        "applied_transform": transform.to_dict(),
        "probability_coordinate_detail": gate["probability_coordinate_detail"],
        "validation_order": "semantic_map_before_frozen_kernel_score",
        "evidence_boundary": gate["evidence_boundary"],
    }
    return result


def main(argv: Sequence[str] | None = None) -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batches", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--json", required=True, type=Path)
    parser.add_argument("--csv", required=True, type=Path)
    args = parser.parse_args(argv)
    run = frozen_kernel.validate_metadata(frozen_kernel.load_metadata(args.metadata))
    rows = frozen_kernel.read_batches(args.batches, run)
    result = score_typed(root, rows, run)
    frozen_kernel.write_outputs(result, args.json, args.csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

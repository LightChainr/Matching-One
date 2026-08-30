#!/usr/bin/env python3
"""Type-safe entrypoint for the frozen Issue #212 matching-odd synthesis."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Sequence

import score_matching_odd_synthesis as frozen_kernel
from wrapping_channels import ObservableDescriptor, map_observable


SEMANTIC_MANIFEST = "predictions/matching_odd_synthesis_semantic_gate_20260830.json"


def load_semantic_gate(root: Path) -> tuple[dict, list[dict[str, object]]]:
    payload = json.loads((root / SEMANTIC_MANIFEST).read_text(encoding="utf-8"))
    if payload.get("status") != "semantic_gate_added_after_frozen_issue212_synthesis":
        raise ValueError("Issue #212 synthesis semantic gate status changed")
    expected_specs = {spec["id"]: spec for spec in frozen_kernel.BLOCK_SPECS}
    blocks = payload.get("blocks")
    if not isinstance(blocks, list) or {block.get("id") for block in blocks} != set(expected_specs):
        raise ValueError("Issue #212 synthesis semantic block IDs changed")

    checked = []
    for block in blocks:
        expected_spec = expected_specs[block["id"]]
        if block.get("raw_data_group") != expected_spec["raw_data_group"]:
            raise ValueError(f"semantic raw-data group changed for {block['id']}")
        if block.get("h4_score") != expected_spec["h4_score"]:
            raise ValueError(f"semantic H4 score changed for {block['id']}")
        source = ObservableDescriptor.from_dict(block["source_descriptor"])
        target = ObservableDescriptor.from_dict(block["target_descriptor"])
        transform = map_observable(source, target)
        registered = block["exact_registered_map"]
        if (transform.scale, transform.offset) != (
            float(registered["scale"]),
            float(registered["offset"]),
        ):
            raise ValueError(f"registered matching-odd map differs for {block['id']}")
        if (transform.scale, transform.offset) != (1.0, 0.0):
            raise ValueError(f"matching-odd synthesis requires exact identity for {block['id']}")
        checked.append(
            {
                "id": block["id"],
                "source_descriptor": source.to_dict(),
                "target_descriptor": target.to_dict(),
                "applied_transform": transform.to_dict(),
            }
        )
    return payload, checked


def synthesize_typed(
    root: Path, ledger: dict, *, source_sha256: str | None = None
) -> dict:
    gate, checked = load_semantic_gate(root)
    result = frozen_kernel.synthesize(ledger, source_sha256=source_sha256)
    result["observable_semantics"] = {
        "semantic_manifest": SEMANTIC_MANIFEST,
        "semantic_manifest_status": gate["status"],
        "blocks": checked,
        "validation_order": "all_registered_maps_before_frozen_synthesis",
        "evidence_boundary": gate["evidence_boundary"],
    }
    return result


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("ledger", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    root = Path(__file__).resolve().parents[1]
    raw = args.ledger.read_bytes()
    result = synthesize_typed(
        root, json.loads(raw), source_sha256=hashlib.sha256(raw).hexdigest()
    )
    result["source_ledger"]["path"] = args.ledger.as_posix()
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

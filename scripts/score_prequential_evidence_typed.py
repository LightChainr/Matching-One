#!/usr/bin/env python3
"""Typed replay entrypoint for the frozen prequential evidence ledger."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
from typing import Callable, Sequence

import score_prequential_evidence as frozen
from wrapping_channels import ObservableDescriptor, map_observable


SEMANTIC_GATE = "predictions/prequential_evidence_semantic_gate_20260830.json"


def git_blob_sha(data: bytes) -> str:
    return hashlib.sha1(f"blob {len(data)}\0".encode("ascii") + data).hexdigest()


def load_semantic_gate(root: Path):
    gate = json.loads((root / SEMANTIC_GATE).read_text(encoding="utf-8"))
    if gate.get("status") != "semantic_gate_added_after_frozen_evidence_ledger":
        raise ValueError("prequential semantic gate status changed")
    if gate.get("frozen_kernel_git_blob") != "364878808083c88dae6a5d9ba4ecd7eb987e2c16":
        raise ValueError("prequential frozen kernel identity changed")
    manifest_path = root / gate["canonical_manifest"]
    if git_blob_sha(manifest_path.read_bytes()) != gate.get("canonical_manifest_git_blob"):
        raise ValueError("prequential canonical manifest identity changed")

    labels = {
        name: ObservableDescriptor.from_dict(definition["descriptor"])
        for name, definition in gate["labels"].items()
    }
    validated = []
    for block in gate["blocks"]:
        source = labels[block["source"]]
        target = labels[block["target"]]
        transform = map_observable(source, target)
        if (transform.scale, transform.offset) != (
            float(block["scale"]), float(block["offset"])
        ):
            raise ValueError(f"prequential registered map changed for {block['id']}")
        validated.append({
            "contract": block,
            "source_descriptor": source,
            "target_descriptor": target,
            "transform": transform,
            "source_response_coordinate": gate["labels"][block["source"]]["response_coordinate"],
            "target_response_coordinate": gate["labels"][block["target"]]["response_coordinate"],
        })
    return gate, validated


def validate_manifest(root: Path, manifest: dict, gate: dict) -> None:
    canonical = json.loads((root / gate["canonical_manifest"]).read_text(encoding="utf-8"))
    if manifest != canonical:
        raise ValueError("typed prequential replay requires the frozen canonical manifest")
    if len(manifest["blocks"]) != len(gate["blocks"]):
        raise ValueError("prequential block count changed")
    for actual, expected in zip(manifest["blocks"], gate["blocks"]):
        for field in ("id", "raw_data_group", "role", "status"):
            if actual.get(field) != expected[field]:
                raise ValueError(f"prequential {field} changed for {expected['id']}")
        legacy = actual.get("channel", {})
        if legacy != {
            "source": expected["source"],
            "target": expected["target"],
            "exact_map": None,
        }:
            raise ValueError(f"prequential legacy channel contract changed for {expected['id']}")


def score_manifest_typed(
    root: Path,
    manifest: dict,
    *,
    scorer: Callable[[dict], dict] = frozen.score_manifest,
) -> dict:
    gate, validated = load_semantic_gate(root)
    validate_manifest(root, manifest, gate)
    result = scorer(copy.deepcopy(manifest))
    expected_order = [block["id"] for block in gate["blocks"]]
    if [block["id"] for block in result.get("blocks", [])] != expected_order:
        raise ValueError("prequential frozen result block order changed")
    result["observable_semantics"] = {
        "semantic_gate": SEMANTIC_GATE,
        "semantic_gate_status": gate["status"],
        "blocks": [
            {
                "id": row["contract"]["id"],
                "status": row["contract"]["status"],
                "source_descriptor": row["source_descriptor"].to_dict(),
                "target_descriptor": row["target_descriptor"].to_dict(),
                "applied_transform": row["transform"].to_dict(),
                "source_response_coordinate": row["source_response_coordinate"],
                "target_response_coordinate": row["target_response_coordinate"],
            }
            for row in validated
        ],
        "validation_order": "canonical_ledger_identity_and_registered_maps_before_frozen_score",
        "protocol_history_boundary": gate["protocol_history_boundary"],
        "chronology_boundary": gate["chronology_boundary"],
        "evidence_boundary": gate["evidence_boundary"],
    }
    return result


def main(argv: Sequence[str] | None = None) -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "manifest", nargs="?", type=Path,
        default=root / "analysis/evidence_ledger_manifest.yaml",
    )
    parser.add_argument("--json-output", type=Path)
    parser.add_argument("--markdown-output", type=Path)
    args = parser.parse_args(argv)
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    result = score_manifest_typed(root, manifest)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.json_output:
        args.json_output.parent.mkdir(parents=True, exist_ok=True)
        args.json_output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    if args.markdown_output:
        args.markdown_output.parent.mkdir(parents=True, exist_ok=True)
        args.markdown_output.write_text(frozen.render_markdown(result), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

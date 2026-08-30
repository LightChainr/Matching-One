#!/usr/bin/env python3
"""Type-safe entrypoint for the frozen C4 N=26 exact-law scorer."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Mapping, Sequence

import score_c4_self_matching_n26 as frozen_kernel
from wrapping_channels import ObservableDescriptor, map_observable


SEMANTIC_MANIFEST = "predictions/c4_self_matching_n26_semantic_gate_20260830.json"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_semantic_gate(
    root: Path,
) -> tuple[dict, ObservableDescriptor, ObservableDescriptor, object]:
    payload = json.loads((root / SEMANTIC_MANIFEST).read_text(encoding="utf-8"))
    if payload.get("status") != "semantic_gate_added_after_freeze_before_typed_replay":
        raise ValueError("C4 N=26 semantic gate status changed")
    source = ObservableDescriptor.from_dict(payload["source_descriptor"])
    target = ObservableDescriptor.from_dict(payload["target_descriptor"])
    transform = map_observable(source, target)
    expected = payload["exact_registered_map"]
    if (transform.scale, transform.offset) != (
        float(expected["scale"]),
        float(expected["offset"]),
    ):
        raise ValueError("registered C4 N=26 map differs from the semantic gate")
    if (transform.scale, transform.offset) != (1.0, 0.0):
        raise ValueError("C4 N=26 frozen and target descriptors must be identical")
    return payload, source, target, transform


def score_typed(
    root: Path,
    prediction_path: Path,
    primary_path: Path,
    reproduction_path: Path,
) -> dict[str, object]:
    gate, source, target, transform = load_semantic_gate(root)
    if _sha256(prediction_path) != gate["frozen_prediction_sha256"]:
        raise ValueError("C4 N=26 prediction hash differs from the semantic gate")
    prediction = json.loads(prediction_path.read_text(encoding="utf-8"))
    if prediction.get("geometry", {}).get("wrapping_channel") != source.channel.value:
        raise ValueError("C4 N=26 prediction channel differs from its source descriptor")
    result = frozen_kernel.score(prediction_path, primary_path, reproduction_path)
    if "observable_semantics" in result:
        raise ValueError("frozen C4 N=26 output already contains observable_semantics")
    result["observable_semantics"] = {
        "semantic_manifest": SEMANTIC_MANIFEST,
        "semantic_manifest_status": gate["status"],
        "source_descriptor": source.to_dict(),
        "target_descriptor": target.to_dict(),
        "applied_transform": transform.to_dict(),
        "validation_order": "semantic_map_before_frozen_kernel_score",
        "finite_control_boundary": gate["boundary"],
    }
    return result


def main(argv: Sequence[str] | None = None) -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prediction", type=Path, required=True)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--reproduction", type=Path, required=True)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args(argv)
    payload = score_typed(root, args.prediction, args.primary, args.reproduction)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

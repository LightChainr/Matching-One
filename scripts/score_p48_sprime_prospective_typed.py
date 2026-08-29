#!/usr/bin/env python3
"""Type-safe entrypoint for the frozen prospective P48 S-prime scorer.

The numerical scorer and both 2026-08-28 model artifacts remain unchanged.
This wrapper validates the post-freeze observable contract before delegating to
the frozen kernel and annotates the required output file with that contract.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys
from typing import Mapping, Sequence

import yaml

from wrapping_channels import ObservableDescriptor, map_observable


SEMANTIC_MANIFEST = "predictions/p48_sprime_semantic_gate_20260829.yaml"


def load_semantic_gate(
    root: Path,
) -> tuple[dict, ObservableDescriptor, ObservableDescriptor, object]:
    path = root / SEMANTIC_MANIFEST
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("P48 semantic manifest must be a mapping")
    if payload.get("status") != "semantic_gate_added_after_freeze_before_typed_score":
        raise ValueError("P48 semantic gate status changed")

    source = ObservableDescriptor.from_dict(payload["source_descriptor"])
    target = ObservableDescriptor.from_dict(payload["target_descriptor"])
    transform = map_observable(source, target)
    expected = payload["exact_registered_map"]
    if transform.scale != float(expected["scale"]) or transform.offset != float(
        expected["offset"]
    ):
        raise ValueError("registered P48 map no longer matches the semantic manifest")
    if transform.scale != 1.0 or transform.offset != 0.0:
        raise ValueError("P48 normalized cross/even map must be exact identity")
    return payload, source, target, transform


def find_output_path(arguments: Sequence[str]) -> Path:
    for index, argument in enumerate(arguments):
        if argument == "--output":
            if index + 1 >= len(arguments):
                raise ValueError("--output requires a path")
            return Path(arguments[index + 1])
        if argument.startswith("--output="):
            return Path(argument.split("=", 1)[1])
    raise ValueError("typed P48 scorer requires --output")


def annotate_output(
    output_path: Path,
    semantic_manifest: Mapping[str, object],
    source: ObservableDescriptor,
    target: ObservableDescriptor,
    transform: object,
) -> None:
    payload = json.loads(output_path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("frozen P48 scorer output must be a JSON object")
    if "observable_semantics" in payload:
        raise ValueError("P48 scorer output already contains observable_semantics")
    payload["observable_semantics"] = {
        "semantic_manifest": SEMANTIC_MANIFEST,
        "semantic_manifest_status": semantic_manifest["status"],
        "frozen_scoring_manifest": semantic_manifest["frozen_scoring_manifest"],
        "frozen_model_artifact": semantic_manifest["frozen_model_artifact"],
        "source_descriptor": source.to_dict(),
        "target_descriptor": target.to_dict(),
        "applied_transform": transform.to_dict(),
        "validation_order": "semantic_map_before_frozen_kernel_score",
    }
    output_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    manifest, source, target, transform = load_semantic_gate(root)
    output_path = find_output_path(sys.argv[1:])

    import score_p48_sprime_prospective as frozen_kernel

    return_code = frozen_kernel.main()
    if return_code != 0:
        return return_code
    annotate_output(output_path, manifest, source, target, transform)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

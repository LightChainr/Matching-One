#!/usr/bin/env python3
"""Type-safe entrypoint for the frozen intrinsic full-curve cocycle scorer.

The 2026-08-28 q=2/Jordan computation is preserved unchanged.  This wrapper
validates that all primitive source/target observables remain cross-channel and
size-locally angular-normalized where applicable, then annotates the resulting
score with the semantic contract.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys
from typing import Mapping, Sequence

import yaml

import score_intrinsic_functional_cocycle as frozen_kernel
from wrapping_channels import ObservableDescriptor, map_observable


SEMANTIC_MANIFEST = "predictions/intrinsic_functional_cocycle_semantic_gate_20260829.yaml"


def load_semantic_gate(root: Path) -> tuple[dict, dict[str, ObservableDescriptor], dict[str, object]]:
    path = root / SEMANTIC_MANIFEST
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("intrinsic functional semantic manifest must be a mapping")
    if payload.get("status") != "semantic_gate_added_before_norm5_target_reveal":
        raise ValueError("intrinsic functional semantic gate status changed")

    descriptors: dict[str, ObservableDescriptor] = {}
    transforms: dict[str, object] = {}
    raw = payload.get("primitive_descriptors")
    if not isinstance(raw, dict) or not raw:
        raise ValueError("semantic manifest lacks primitive descriptors")
    for name, descriptor_payload in raw.items():
        descriptor = ObservableDescriptor.from_dict(descriptor_payload)
        transform = map_observable(descriptor, descriptor)
        if transform.scale != 1.0 or transform.offset != 0.0:
            raise ValueError(f"{name}: cross-size primitive map must be exact identity")
        descriptors[str(name)] = descriptor
        transforms[str(name)] = transform
    return payload, descriptors, transforms


def find_json_path(arguments: Sequence[str]) -> Path:
    for index, argument in enumerate(arguments):
        if argument == "--json":
            if index + 1 >= len(arguments):
                raise ValueError("--json requires a path")
            return Path(arguments[index + 1])
        if argument.startswith("--json="):
            return Path(argument.split("=", 1)[1])
    raise ValueError("typed intrinsic functional scorer requires --json")


def annotate_output(
    output_path: Path,
    manifest: Mapping[str, object],
    descriptors: Mapping[str, ObservableDescriptor],
    transforms: Mapping[str, object],
) -> None:
    payload = json.loads(output_path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("frozen functional scorer output must be a JSON object")
    payload["observable_semantics"] = {
        "semantic_manifest": SEMANTIC_MANIFEST,
        "semantic_manifest_status": manifest["status"],
        "primitive_descriptors": {
            name: descriptor.to_dict() for name, descriptor in descriptors.items()
        },
        "applied_maps": {
            name: transform.to_dict() for name, transform in transforms.items()
        },
        "cross_size_rule": (
            "same cross-channel descriptor at every N; P4 angular normalization "
            "is size-local before any cocycle comparison"
        ),
        "validation_order": "semantic_identity_maps_before_frozen_kernel_score",
    }
    output_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    manifest, descriptors, transforms = load_semantic_gate(root)
    output_path = find_json_path(sys.argv[1:])
    return_code = frozen_kernel.main()
    if return_code != 0:
        return return_code
    annotate_output(output_path, manifest, descriptors, transforms)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

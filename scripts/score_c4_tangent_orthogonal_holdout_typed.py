#!/usr/bin/env python3
"""Type-safe entrypoint for the frozen C4 tangent orthogonal holdout score."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

import score_c4_tangent_orthogonal_holdout as frozen_kernel
from wrapping_channels import ObservableDescriptor, map_observable


SEMANTIC_MANIFEST = "predictions/c4_tangent_orthogonal_semantic_gate_20260830.json"


def load_semantic_gate(
    root: Path,
) -> tuple[dict, ObservableDescriptor, ObservableDescriptor, object]:
    payload = json.loads((root / SEMANTIC_MANIFEST).read_text(encoding="utf-8"))
    if payload.get("status") != "semantic_gate_added_after_frozen_c4_tangent_holdout":
        raise ValueError("C4 tangent semantic gate status changed")
    if payload.get("frozen_kernel_git_blob") != "b58cbaade09aafa1378f2c70399146efaa67dd6f":
        raise ValueError("C4 tangent frozen kernel identity changed")
    if (payload.get("source_N"), payload.get("target_N")) != (
        frozen_kernel.SOURCE_N,
        frozen_kernel.TARGET_N,
    ):
        raise ValueError("C4 tangent source/target sizes changed")
    if payload.get("primary_channel") != "cross":
        raise ValueError("C4 tangent primary channel changed")
    if payload.get("response_coordinates_in_order") != ["t", "lambda"]:
        raise ValueError("C4 tangent response-coordinate order changed")
    exponent = payload.get("thermal_size_map", {}).get("exponent_in_N", {})
    if exponent != {"numerator": 3, "denominator": 8, "fitted": False}:
        raise ValueError("C4 tangent thermal exponent changed")

    source = ObservableDescriptor.from_dict(payload["source_descriptor"])
    target = ObservableDescriptor.from_dict(payload["target_descriptor"])
    transform = map_observable(source, target)
    expected = payload["exact_registered_map"]
    if (transform.scale, transform.offset) != (
        float(expected["scale"]),
        float(expected["offset"]),
    ):
        raise ValueError("registered C4 tangent map differs from the semantic gate")
    if (transform.scale, transform.offset) != (1.0, 0.0):
        raise ValueError("C4 tangent source/target map must be exact identity")
    return payload, source, target, transform


def render_typed(root: Path, source_path: Path, target_path: Path) -> dict:
    gate, source, target, transform = load_semantic_gate(root)
    result = frozen_kernel.render(source_path, target_path)
    if result.get("primary_channel") != gate["primary_channel"]:
        raise ValueError("C4 tangent result channel differs from semantic gate")
    if (result.get("source_N"), result.get("target_N")) != (
        gate["source_N"],
        gate["target_N"],
    ):
        raise ValueError("C4 tangent result sizes differ from semantic gate")
    required = {
        "source_t",
        "source_lambda",
        "source_c",
        "target_t",
        "target_lambda",
        "target_c",
        "orthogonal_residual",
        "thermal_scaling_residual",
    }
    if set(result.get("point", {})) != required:
        raise ValueError("C4 tangent response-coordinate payload changed")
    result["observable_semantics"] = {
        "semantic_manifest": SEMANTIC_MANIFEST,
        "semantic_manifest_status": gate["status"],
        "source_descriptor": source.to_dict(),
        "target_descriptor": target.to_dict(),
        "applied_transform": transform.to_dict(),
        "response_coordinates_in_order": gate["response_coordinates_in_order"],
        "source_projection": gate["source_projection"],
        "target_projection": gate["target_projection"],
        "thermal_size_map": gate["thermal_size_map"],
        "batch_relation": gate["batch_relation"],
        "validation_order": "semantic_map_before_frozen_kernel_score",
        "evidence_boundary": gate["evidence_boundary"],
    }
    return result


def main(argv: Sequence[str] | None = None) -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--target", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    payload = render_typed(root, args.source, args.target)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Type-safe entrypoint for the chronological frozen P48 S-prime scorer."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Mapping, Sequence

import score_p48_sprime_frozen as frozen_kernel
from wrapping_channels import ObservableDescriptor, map_observable


SEMANTIC_MANIFEST = "predictions/p48_sprime_frozen_semantic_gate_20260830.json"


def load_semantic_gate(
    root: Path,
) -> tuple[dict, ObservableDescriptor, ObservableDescriptor, object]:
    payload = json.loads((root / SEMANTIC_MANIFEST).read_text(encoding="utf-8"))
    if payload.get("status") != "semantic_gate_added_after_frozen_chronological_p48_score":
        raise ValueError("chronological P48 semantic gate status changed")
    if payload.get("frozen_kernel_git_blob") != "18046c546b87e0feebbf43cb6060820630b281c4":
        raise ValueError("chronological P48 frozen kernel identity changed")
    if tuple(payload.get("target_sizes", ())) != (185, 265):
        raise ValueError("chronological P48 target sizes changed")
    if float(payload.get("leading_power_in_N")) != 1.25:
        raise ValueError("chronological P48 leading power changed")
    if payload.get("models_in_scoring_order") != [
        "pure_N^-5/4",
        "zero_effect",
        "q2_even_scalar",
        "rank2_jordan_log",
    ]:
        raise ValueError("chronological P48 model order changed")

    source = ObservableDescriptor.from_dict(payload["source_descriptor"])
    target = ObservableDescriptor.from_dict(payload["target_descriptor"])
    transform = map_observable(source, target)
    expected = payload["exact_registered_map"]
    if (transform.scale, transform.offset) != (
        float(expected["scale"]),
        float(expected["offset"]),
    ):
        raise ValueError("registered chronological P48 map differs from the semantic gate")
    if (transform.scale, transform.offset) != (1.0, 0.0):
        raise ValueError("chronological P48 descriptor map must be exact identity")
    return payload, source, target, transform


def score_typed(
    root: Path,
    target_payload: Mapping[str, object],
    scoring_manifest: Mapping[str, object],
) -> dict:
    gate, source, target, transform = load_semantic_gate(root)
    result = frozen_kernel.score(target_payload, scoring_manifest)
    if result.get("sizes") != gate["target_sizes"]:
        raise ValueError("chronological P48 result sizes differ from the semantic gate")
    if result.get("scoring_order") != gate["models_in_scoring_order"]:
        raise ValueError("chronological P48 result model order differs from the semantic gate")
    result["observable_semantics"] = {
        "semantic_manifest": SEMANTIC_MANIFEST,
        "semantic_manifest_status": gate["status"],
        "frozen_scoring_manifest": gate["frozen_scoring_manifest"],
        "shared_descriptor_evidence": gate["shared_descriptor_evidence"],
        "source_descriptor": source.to_dict(),
        "target_descriptor": target.to_dict(),
        "applied_transform": transform.to_dict(),
        "validation_order": "semantic_map_before_frozen_kernel_score",
        "evidence_boundary": gate["evidence_boundary"],
    }
    return result


def main(argv: Sequence[str] | None = None) -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path, required=True)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=root / "predictions/p48_sprime_scoring_manifest_20260828.yaml",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    result = score_typed(
        root,
        frozen_kernel.read_json(args.target),
        frozen_kernel.read_yaml(args.manifest),
    )
    rendered = json.dumps(result, indent=2, sort_keys=False) + "\n"
    print(rendered, end="")
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

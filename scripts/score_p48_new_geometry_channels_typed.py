#!/usr/bin/env python3
"""Typed entrypoint for the frozen P48 new-geometry four-projector score."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Mapping

import score_p48_new_geometry_channels as frozen
from wrapping_channels import ObservableDescriptor, map_observable


SEMANTIC_MANIFEST = "predictions/p48_new_geometry_channels_semantic_gate_20260830.json"
EXPECTED_PROJECTORS = {
    "P4_S": ("even", "second_minus_first", "value", 1, 1, "A_S_N1"),
    "P4_D": ("odd", "first_minus_second", "value", 13, 8, "A_D_N13_8"),
    "P4_S_prime": (
        "even", "second_minus_first", "first_p_derivative", 5, 4,
        "A_Sprime_N5_4",
    ),
    "P4_D_prime": (
        "odd", "first_minus_second", "first_p_derivative", 5, 8,
        "A_Dprime_N5_8",
    ),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_semantic_gate(root: Path) -> tuple[dict, dict[str, dict[str, object]]]:
    payload = json.loads((root / SEMANTIC_MANIFEST).read_text(encoding="utf-8"))
    if payload.get("status") != "semantic_gate_added_after_frozen_p48_new_geometry_score":
        raise ValueError("P48 new-geometry semantic gate status changed")
    if payload.get("frozen_kernel_git_blob") != "6302d49471a9e78f311d8f492e7845c2b3f5ed7e":
        raise ValueError("P48 new-geometry frozen kernel identity changed")
    if payload.get("source_training_sizes") != [65, 85, 130]:
        raise ValueError("P48 new-geometry source sizes changed")
    if payload.get("target_sizes") != [185, 265]:
        raise ValueError("P48 new-geometry target sizes changed")
    if payload.get("projector_order") != list(frozen.CHANNELS):
        raise ValueError("P48 new-geometry projector order changed")
    if payload.get("target_refit_parameters") != 0:
        raise ValueError("P48 new-geometry target refit boundary changed")

    validated: dict[str, dict[str, object]] = {}
    for projector in payload["projector_order"]:
        definition = payload.get("projectors", {}).get(projector, {})
        combination, orientation, response, numerator, denominator, scaled_key = (
            EXPECTED_PROJECTORS[projector]
        )
        expected_descriptor = {
            "channel": "cross",
            "combination": combination,
            "coordinate": "p",
            "orientation_order": orientation,
            "normalization": "angular_normalized",
            "quantity": "orientation_contrast",
        }
        if definition.get("source_descriptor") != expected_descriptor:
            raise ValueError("{} source descriptor changed".format(projector))
        if definition.get("target_descriptor") != expected_descriptor:
            raise ValueError("{} target descriptor changed".format(projector))
        if definition.get("response_coordinate") != response:
            raise ValueError("{} response coordinate changed".format(projector))
        if definition.get("normalization_power_in_N") != {
            "numerator": numerator,
            "denominator": denominator,
        }:
            raise ValueError("{} normalization power changed".format(projector))
        if definition.get("scaled_key") != scaled_key:
            raise ValueError("{} scaled key changed".format(projector))
        specification = frozen.CHANNELS[projector]
        if specification["scaled_key"] != scaled_key or specification["power"] != (
            numerator / denominator
        ):
            raise ValueError("{} frozen-kernel specification changed".format(projector))

        source = ObservableDescriptor.from_dict(definition["source_descriptor"])
        target = ObservableDescriptor.from_dict(definition["target_descriptor"])
        transform = map_observable(source, target)
        expected_map = definition.get("exact_registered_map", {})
        if (transform.scale, transform.offset) != (
            float(expected_map.get("scale")),
            float(expected_map.get("offset")),
        ) or (transform.scale, transform.offset) != (1.0, 0.0):
            raise ValueError("{} source-to-target map changed".format(projector))
        validated[projector] = {
            "source_descriptor": source,
            "target_descriptor": target,
            "transform": transform,
            "response_coordinate": response,
            "normalization_power_in_N": definition["normalization_power_in_N"],
            "scaled_key": scaled_key,
        }
    return payload, validated


def score_typed(
    root: Path,
    source_summary: dict,
    target_payloads: Mapping[int, dict],
) -> dict:
    gate, validated = load_semantic_gate(root)
    result = frozen.score(source_summary, dict(target_payloads))
    if result.get("source_training_sizes") != gate["source_training_sizes"]:
        raise ValueError("P48 score source sizes differ from semantic gate")
    if result.get("target_sizes") != gate["target_sizes"]:
        raise ValueError("P48 score target sizes differ from semantic gate")
    if result.get("target_refit_parameters") != gate["target_refit_parameters"]:
        raise ValueError("P48 score refit boundary differs from semantic gate")
    if list(result.get("channels", {})) != gate["projector_order"]:
        raise ValueError("P48 score projector order differs from semantic gate")
    result["observable_semantics"] = {
        "semantic_manifest": SEMANTIC_MANIFEST,
        "semantic_manifest_status": gate["status"],
        "projector_order": gate["projector_order"],
        "projectors": {
            projector: {
                "source_descriptor": values["source_descriptor"].to_dict(),
                "target_descriptor": values["target_descriptor"].to_dict(),
                "applied_transform": values["transform"].to_dict(),
                "response_coordinate": values["response_coordinate"],
                "normalization_power_in_N": values["normalization_power_in_N"],
                "scaled_key": values["scaled_key"],
            }
            for projector, values in validated.items()
        },
        "target_relation": gate["target_relation"],
        "source_error_relation": gate["source_error_relation"],
        "covariance_contract": gate["covariance_contract"],
        "evidence_boundary": gate["evidence_boundary"],
        "validation_order": "semantic_maps_before_frozen_four_projector_score",
    }
    return result


def load_canonical_inputs(
    gate: Mapping[str, object], source_path: Path, n185_path: Path, n265_path: Path
) -> tuple[dict, dict[int, dict]]:
    paths = {"source": source_path, "185": n185_path, "265": n265_path}
    for key, path in paths.items():
        expected = gate["canonical_inputs"][key]["sha256"]
        if sha256(path) != expected:
            raise ValueError("P48 canonical {} input identity changed".format(key))
    return (
        frozen._load_json(source_path),
        {185: frozen._load_json(n185_path), 265: frozen._load_json(n265_path)},
    )


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    gate, _ = load_semantic_gate(root)
    canonical = gate["canonical_inputs"]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=root / canonical["source"]["path"])
    parser.add_argument("--n185", type=Path, default=root / canonical["185"]["path"])
    parser.add_argument("--n265", type=Path, default=root / canonical["265"]["path"])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    source, targets = load_canonical_inputs(gate, args.source, args.n185, args.n265)
    payload = score_typed(root, source, targets)
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    else:
        print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

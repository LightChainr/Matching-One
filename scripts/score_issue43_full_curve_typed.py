#!/usr/bin/env python3
"""Typed entrypoint for the frozen Issue #43 full-curve score."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Mapping

import score_issue43_full_curve as frozen
from wrapping_channels import ObservableDescriptor, map_observable


SEMANTIC_MANIFEST = "predictions/issue43_full_curve_semantic_gate_20260830.json"
EXPECTED_DESCRIPTORS = {
    "DeltaM": {
        "channel": "cross",
        "combination": "odd",
        "coordinate": "p",
        "orientation_order": "first_minus_second",
        "normalization": "raw",
        "quantity": "orientation_contrast",
    },
    "DeltaS": {
        "channel": "cross",
        "combination": "even",
        "coordinate": "p",
        "orientation_order": "first_minus_second",
        "normalization": "raw",
        "quantity": "orientation_contrast",
    },
}


def load_semantic_gate(root: Path) -> tuple[dict, dict[str, dict[str, object]]]:
    payload = json.loads((root / SEMANTIC_MANIFEST).read_text(encoding="utf-8"))
    if payload.get("status") != "semantic_gate_added_after_frozen_issue43_full_curve_score":
        raise ValueError("Issue43 full-curve semantic gate status changed")
    kernels = payload.get("frozen_kernels", {})
    expected_blobs = {
        "base": ("scripts/score_issue43_full_curve.py", "067b0b613b9c2f412bc75aa0361d488ed68e44d3"),
        "production_lock": ("scripts/score_issue43_full_curve_locked.py", "022ac040596eb4518773982975edd101ff3d0a26"),
    }
    for name, (path, blob) in expected_blobs.items():
        if kernels.get(name) != {"path": path, "git_blob": blob}:
            raise ValueError("Issue43 frozen {} kernel identity changed".format(name))
    if payload.get("p_ref") != frozen.P_REF:
        raise ValueError("Issue43 frozen p_ref changed")
    if (payload.get("source_N"), payload.get("target_N")) != (185, 265):
        raise ValueError("Issue43 source/target sizes changed")
    if payload.get("sector_order") != ["DeltaM", "DeltaS"]:
        raise ValueError("Issue43 sector order changed")
    if payload.get("prediction_sha256") != frozen.PREDICTION_SHA256:
        raise ValueError("Issue43 prediction identity changed")

    validated: dict[str, dict[str, object]] = {}
    sectors = payload.get("sectors", {})
    for sector in payload["sector_order"]:
        definition = sectors.get(sector, {})
        expected_descriptor = EXPECTED_DESCRIPTORS[sector]
        if definition.get("source_descriptor") != expected_descriptor:
            raise ValueError("Issue43 {} source descriptor changed".format(sector))
        if definition.get("target_descriptor") != expected_descriptor:
            raise ValueError("Issue43 {} target descriptor changed".format(sector))
        source = ObservableDescriptor.from_dict(definition["source_descriptor"])
        target = ObservableDescriptor.from_dict(definition["target_descriptor"])
        transform = map_observable(source, target)
        expected_map = definition.get("exact_registered_map", {})
        if (transform.scale, transform.offset) != (
            float(expected_map.get("scale")),
            float(expected_map.get("offset")),
        ) or (transform.scale, transform.offset) != (1.0, 0.0):
            raise ValueError("Issue43 {} cross-size map is not exact identity".format(sector))
        validated[sector] = {
            "source_descriptor": source,
            "target_descriptor": target,
            "transform": transform,
        }
    return payload, validated


def annotate_result(
    result: dict,
    gate: Mapping[str, object],
    validated: Mapping[str, Mapping[str, object]],
    operational_entrypoint: str,
) -> dict:
    if result.get("p_ref") != gate["p_ref"]:
        raise ValueError("Issue43 score p_ref differs from semantic gate")
    if list(result.get("scores", {})) != gate["sector_order"]:
        raise ValueError("Issue43 score sector order differs from semantic gate")
    result["observable_semantics"] = {
        "semantic_manifest": SEMANTIC_MANIFEST,
        "semantic_manifest_status": gate["status"],
        "operational_entrypoint": operational_entrypoint,
        "source_N": gate["source_N"],
        "target_N": gate["target_N"],
        "p_ref": gate["p_ref"],
        "sector_order": gate["sector_order"],
        "sectors": {
            sector: {
                "source_descriptor": values["source_descriptor"].to_dict(),
                "target_descriptor": values["target_descriptor"].to_dict(),
                "applied_transform": values["transform"].to_dict(),
            }
            for sector, values in validated.items()
        },
        "source_error_relation": gate["source_error_relation"],
        "target_relation": gate["target_relation"],
        "evidence_boundary": gate["evidence_boundary"],
        "validation_order": "semantic_maps_before_frozen_full_curve_score",
    }
    if operational_entrypoint == "production_lock":
        result["observable_semantics"]["production_lock"] = gate["production_lock"]
    return result


def analyze_typed(
    root: Path,
    runs: Mapping[int, Mapping[str, object]],
    prediction_path: Path,
    *,
    operational_entrypoint: str = "base",
) -> dict:
    gate, validated = load_semantic_gate(root)
    result = frozen.analyze(runs, prediction_path)
    return annotate_result(result, gate, validated, operational_entrypoint)


def main(*, operational_entrypoint: str = "base") -> int:
    root = Path(__file__).resolve().parents[1]
    # Fail before parsing or reconstructing target data if semantics drifted.
    load_semantic_gate(root)
    args = frozen.parse_args()
    runs = {}
    for hist_name, moments_name, metadata_name in args.run:
        metadata = frozen.validate_metadata(Path(metadata_name))
        n = int(metadata["N"])
        if n in runs:
            raise ValueError("duplicate N run")
        records = frozen.read_histograms(Path(hist_name), metadata)
        frozen.validate_moments(Path(moments_name), metadata, records)
        runs[n] = {"metadata": metadata, "sectors": frozen.reconstruct(records)}
    result = analyze_typed(
        root,
        runs,
        args.predictions,
        operational_entrypoint=operational_entrypoint,
    )
    frozen.write_outputs(result, args.json, args.csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

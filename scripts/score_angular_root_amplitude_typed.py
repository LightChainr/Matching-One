#!/usr/bin/env python3
"""Typed entrypoint for the frozen P45 angular root-amplitude scorer."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Callable, Mapping, Sequence

import mpmath as mp

import score_angular_root_amplitude as frozen
from wrapping_channels import ObservableDescriptor, map_observable


SEMANTIC_GATE = "predictions/angular_root_amplitude_semantic_gate_20260830.json"


def load_semantic_gate(root: Path):
    gate = json.loads((root / SEMANTIC_GATE).read_text(encoding="utf-8"))
    if gate.get("status") != "semantic_gate_added_after_frozen_p45_root_score":
        raise ValueError("angular-root semantic gate status changed")
    if gate.get("frozen_kernel_git_blob") != "a92e575de3a33621d4196ca45c22abdf83bd4111":
        raise ValueError("angular-root frozen kernel identity changed")
    if gate.get("frozen_prediction_git_blob") != "e657a7f051f74171b70c30bbd2753afbfdb07476":
        raise ValueError("angular-root prediction identity changed")
    if gate.get("primary_sizes") != [65, 85] or gate.get("target_refit_parameters") != 0:
        raise ValueError("angular-root size or refit boundary changed")
    if gate.get("metric_order") != list(frozen.METRICS):
        raise ValueError("angular-root metric order changed")
    if gate.get("response_transform") != "A_p=-N^2*(first_root-second_root)/signed_delta_cos4":
        raise ValueError("angular-root response transform changed")

    source = ObservableDescriptor.from_dict(gate["source_descriptor"])
    target = ObservableDescriptor.from_dict(gate["target_descriptor"])
    transforms = {}
    for size in gate["primary_sizes"]:
        design = gate["designs"][str(size)]
        transform = map_observable(
            source,
            target,
            source_angular_factor=float(design["signed_delta_cos4"]),
        )
        if not math.isclose(
            transform.scale,
            float(design["registered_map_scale"]),
            rel_tol=0.0,
            abs_tol=1.0e-15,
        ) or transform.offset != 0.0:
            raise ValueError("angular-root registered map differs from semantic gate")
        transforms[size] = transform
    return gate, source, target, transforms


def validate_record_designs(records: Mapping[tuple, Mapping[str, object]], gate: dict) -> None:
    sizes = sorted({int(key[0]) for key in records})
    if sizes != gate["primary_sizes"]:
        raise ValueError("angular-root record sizes differ from semantic gate")
    for size in sizes:
        design = gate["designs"][str(size)]
        for orientation in ("first", "second"):
            geometries = {
                (int(row["a"]), int(row["b"]))
                for key, row in records.items()
                if int(key[0]) == size and key[1] == orientation
            }
            if geometries != {tuple(design[orientation])}:
                raise ValueError("angular-root orientation geometry changed")


def score_typed(
    root: Path,
    records: Mapping[tuple, Mapping[str, object]],
    p: mp.mpf,
    prediction: float,
    prediction_se: float,
    *,
    scorer: Callable[[dict, mp.mpf, float, float], dict] = frozen.score,
) -> dict:
    gate, source, target, transforms = load_semantic_gate(root)
    validate_record_designs(records, gate)
    expected_prediction = gate["frozen_prediction"]
    if mp.mpf(p) != mp.mpf(gate["p_ref"]):
        raise ValueError("angular-root p_ref differs from semantic gate")
    if prediction != float(expected_prediction["value"]) or prediction_se != float(
        expected_prediction["source_standard_error"]
    ):
        raise ValueError("angular-root frozen prediction changed")
    result = scorer(dict(records), p, prediction, prediction_se)
    if result.get("schema") != "frozen angular-normalized root amplitude score v1":
        raise ValueError("angular-root frozen result schema changed")
    if result.get("sizes") != gate["primary_sizes"]:
        raise ValueError("angular-root frozen result sizes changed")
    frozen_prediction = result.get("frozen_prediction", {})
    if (
        frozen_prediction.get("value") != prediction
        or frozen_prediction.get("source_standard_error") != prediction_se
    ):
        raise ValueError("angular-root frozen result prediction changed")
    result["observable_semantics"] = {
        "semantic_gate": SEMANTIC_GATE,
        "semantic_gate_status": gate["status"],
        "source_descriptor": source.to_dict(),
        "target_descriptor": target.to_dict(),
        "applied_transforms_by_size": {
            str(size): transforms[size].to_dict() for size in gate["primary_sizes"]
        },
        "response_coordinate": gate["response_coordinate"],
        "response_transform": gate["response_transform"],
        "batch_contract": gate["batch_contract"],
        "covariance_contract": gate["covariance_contract"],
        "target_refit_parameters": gate["target_refit_parameters"],
        "evidence_boundary": gate["evidence_boundary"],
        "validation_order": "registered_angular_maps_and_frozen_contract_before_score",
    }
    return result


def main(argv: Sequence[str] | None = None) -> int:
    root = Path(__file__).resolve().parents[1]
    gate, _, _, _ = load_semantic_gate(root)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--histograms", nargs="+", required=True, type=Path)
    parser.add_argument("--moments", nargs="+", required=True, type=Path)
    parser.add_argument("--metadata", nargs="+", required=True, type=Path)
    parser.add_argument("--p", default=gate["p_ref"])
    parser.add_argument("--dps", type=int, default=80)
    parser.add_argument("--prediction", type=float, default=gate["frozen_prediction"]["value"])
    parser.add_argument(
        "--prediction-se", type=float,
        default=gate["frozen_prediction"]["source_standard_error"],
    )
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--source-sha256", required=True)
    parser.add_argument("--binary-sha256", required=True)
    parser.add_argument("--json", required=True, type=Path)
    parser.add_argument("--csv", required=True, type=Path)
    args = parser.parse_args(argv)
    mp.mp.dps = args.dps
    records = frozen.merge_inputs(args.histograms, args.moments)
    result = score_typed(root, records, mp.mpf(args.p), args.prediction, args.prediction_se)
    provenance = frozen.validate_metadata(args.metadata, result["sizes"], args.source_commit)
    result["provenance"] = provenance | {
        "source_sha256": args.source_sha256,
        "binary_sha256": args.binary_sha256,
        "tracked_source_clean_before_and_after_build": True,
        "cross_size_rng_policy": gate["covariance_contract"],
    }
    args.json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    frozen.write_csv(args.csv, result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Type-safe entrypoint for the frozen Issue #43 secondary ledger."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Callable, Mapping, Sequence

from wrapping_channels import ObservableDescriptor, map_observable


SEMANTIC_GATE = "predictions/issue43_secondary_semantic_gate_20260830.json"
STAGE_ORDER = [
    "original_x21_H4_two_sector",
    "x17_over_4_H4_adversarial_radial",
    "zero_effect",
    "predeclared_shared_H4_plus_H12",
    "issue72_P48_S_prime",
]
FIXED_STATUSES = [
    "ALREADY_SCORED_BY_PRIMARY",
    "SCORED_FROZEN_NO_REFIT",
    "ALREADY_SCORED_BY_PRIMARY",
    "NOT_SCORABLE",
]
EXPECTED_DESCRIPTORS = {
    "DeltaM": {
        "channel": "cross", "combination": "odd", "coordinate": "p",
        "orientation_order": "first_minus_second", "normalization": "raw",
        "quantity": "orientation_contrast",
    },
    "DeltaS": {
        "channel": "cross", "combination": "even", "coordinate": "p",
        "orientation_order": "first_minus_second", "normalization": "raw",
        "quantity": "orientation_contrast",
    },
    "P4_S_prime": {
        "channel": "cross", "combination": "even", "coordinate": "p",
        "orientation_order": "first_minus_second",
        "normalization": "angular_normalized",
        "quantity": "orientation_contrast",
    },
}


def load_semantic_gate(root: Path) -> tuple[dict, dict[str, dict[str, object]]]:
    gate = json.loads((root / SEMANTIC_GATE).read_text(encoding="utf-8"))
    if gate.get("status") != "semantic_gate_added_after_frozen_issue43_secondary_ledger":
        raise ValueError("Issue43 secondary semantic gate status changed")
    if gate.get("frozen_kernel_git_blob") != "be08124b54bf27c0ccf88ee9f43b48d66a791b4d":
        raise ValueError("Issue43 secondary frozen kernel identity changed")
    if gate.get("target_sizes") != [185, 265]:
        raise ValueError("Issue43 secondary target sizes changed")
    if gate.get("stage_order") != STAGE_ORDER:
        raise ValueError("Issue43 secondary stage order changed")
    if gate.get("excluded_model") != "V_<1,3>_N^-4/3":
        raise ValueError("Issue43 secondary excluded model changed")
    expected_hashes = {
        "primary_prediction": "a370e79a10854341fac3ee75e8c518dbf3533e8c077cba2c2ec1018178144f44",
        "x17_competitor": "941af010cc146c76e26985ecf3edf58f0df28d987fc79c03725ebc21f64964f5",
        "p48_models": "0d44228ae117f94cb1f99d1e2727eb47390aae950c3ae70c21dd8bc5a09454ae",
    }
    if gate.get("artifact_sha256") != expected_hashes:
        raise ValueError("Issue43 secondary artifact identities changed")
    validated: dict[str, dict[str, object]] = {}
    if gate.get("observable_descriptors") != EXPECTED_DESCRIPTORS:
        raise ValueError("Issue43 secondary observable descriptor changed")
    for name, payload in gate["observable_descriptors"].items():
        source = ObservableDescriptor.from_dict(payload)
        target = ObservableDescriptor.from_dict(payload)
        transform = map_observable(source, target)
        expected = gate["exact_registered_map"]
        if (transform.scale, transform.offset) != (
            float(expected["scale"]), float(expected["offset"])
        ) or (transform.scale, transform.offset) != (1.0, 0.0):
            raise ValueError("Issue43 secondary registered map changed for " + name)
        validated[name] = {
            "source_descriptor": source,
            "target_descriptor": target,
            "transform": transform,
        }
    if list(validated) != ["DeltaM", "DeltaS", "P4_S_prime"]:
        raise ValueError("Issue43 secondary observable order changed")
    return gate, validated


def _run_frozen(
    primary: Mapping[str, object],
    x17_artifact: Path,
    p48_artifact: Path,
    p48_target: Path | None,
) -> dict:
    import score_issue43_secondary as frozen_kernel

    return frozen_kernel.score(primary, x17_artifact, p48_artifact, p48_target)


def score_typed(
    root: Path,
    primary: Mapping[str, object],
    x17_artifact: Path,
    p48_artifact: Path,
    p48_target: Path | None = None,
    *,
    runner: Callable[[Mapping[str, object], Path, Path, Path | None], dict] = _run_frozen,
) -> dict:
    gate, validated = load_semantic_gate(root)
    result = runner(primary, x17_artifact, p48_artifact, p48_target)
    if result.get("protocol") != "Issue #43 frozen secondary scoring ledger":
        raise ValueError("Issue43 secondary frozen protocol changed")
    if result.get("sizes") != gate["target_sizes"]:
        raise ValueError("Issue43 secondary frozen sizes changed")
    if result.get("stage_order") != gate["stage_order"]:
        raise ValueError("Issue43 secondary frozen stage order changed")
    stages = result.get("stages", [])
    if len(stages) != 5 or [stage.get("name") for stage in stages] != STAGE_ORDER:
        raise ValueError("Issue43 secondary frozen stages changed")
    if [stage.get("status") for stage in stages[:4]] != FIXED_STATUSES:
        raise ValueError("Issue43 secondary frozen stage status changed")
    if stages[4].get("status") not in {
        "READY_AWAITING_DERIVATIVE_TARGET", "SCORED_FROZEN_NO_REFIT"
    }:
        raise ValueError("Issue43 secondary P48 stage status changed")
    excluded = result.get("excluded_models")
    if excluded != [{
        "name": gate["excluded_model"],
        "status": "EXCLUDED_INVALIDATED_WRONG_KAC_BRANCH",
        "scored": False,
    }]:
        raise ValueError("Issue43 secondary excluded-model contract changed")
    result["observable_semantics"] = {
        "semantic_gate": SEMANTIC_GATE,
        "semantic_gate_status": gate["status"],
        "stage_order": gate["stage_order"],
        "stage_status_contract": gate["stage_status_contract"],
        "observable_maps": {
            name: {
                "source_descriptor": values["source_descriptor"].to_dict(),
                "target_descriptor": values["target_descriptor"].to_dict(),
                "applied_transform": values["transform"].to_dict(),
            }
            for name, values in validated.items()
        },
        "raw_data_boundary": gate["raw_data_boundary"],
        "evidence_boundary": gate["evidence_boundary"],
        "validation_order": "semantic_maps_before_frozen_secondary_ledger",
    }
    return result


def main(argv: Sequence[str] | None = None) -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary-score", type=Path, required=True)
    parser.add_argument("--p48-target", type=Path)
    parser.add_argument("--x17-artifact", type=Path, default=root / "predictions/x17_spin4_competitor_20260828.yaml")
    parser.add_argument("--p48-artifact", type=Path, default=root / "predictions/p48_sprime_correction_20260828.yaml")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    import score_issue43_secondary as frozen_kernel

    result = score_typed(
        root,
        frozen_kernel.read_json(args.primary_score),
        args.x17_artifact,
        args.p48_artifact,
        args.p48_target,
    )
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    print(rendered, end="")
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

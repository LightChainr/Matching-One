#!/usr/bin/env python3
"""Type-safe entrypoint for the frozen P49 full-curve doubling score."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Callable, Sequence

from wrapping_channels import ObservableDescriptor, map_observable


SEMANTIC_GATE = "predictions/p49_fullcurve_doubling_semantic_gate_20260830.json"
SIZES = [65, 85, 130, 170]
LINEAGES = [[65, 130], [85, 170]]
LEVELS = [0.0, 0.025, 0.05]
PROJECTORS = ["P4_S", "P4_D", "P4_S_prime", "P4_D_prime"]
PROJECTOR_SPECS = {
    "P4_S": ("even", "value", 1, 1, "second_minus_first", "first_minus_second"),
    "P4_D": ("odd", "value", 13, 8, "first_minus_second", "second_minus_first"),
    "P4_S_prime": (
        "even", "first_p_derivative", 5, 4,
        "second_minus_first", "first_minus_second",
    ),
    "P4_D_prime": (
        "odd", "first_p_derivative", 5, 8,
        "first_minus_second", "second_minus_first",
    ),
}


def _git_blob_sha(data: bytes) -> str:
    return hashlib.sha1(f"blob {len(data)}\0".encode("ascii") + data).hexdigest()


def _validate_file(path: Path, expected: dict, label: str) -> None:
    data = path.read_bytes()
    if _git_blob_sha(data) != expected["git_blob"]:
        raise ValueError(f"P49 canonical {label} git blob changed")
    if hashlib.sha256(data).hexdigest() != expected["sha256"]:
        raise ValueError(f"P49 canonical {label} sha256 changed")


def _descriptor(combination: str, orientation_order: str) -> ObservableDescriptor:
    return ObservableDescriptor.from_dict({
        "channel": "cross",
        "combination": combination,
        "coordinate": "p",
        "orientation_order": orientation_order,
        "normalization": "angular_normalized",
        "quantity": "orientation_contrast",
    })


def load_semantic_gate(root: Path) -> tuple[dict, dict[str, object]]:
    gate = json.loads((root / SEMANTIC_GATE).read_text(encoding="utf-8"))
    if gate.get("status") != "semantic_gate_added_after_frozen_p49_fullcurve_score":
        raise ValueError("P49 full-curve semantic gate status changed")
    if gate.get("frozen_kernel_git_blob") != "c68db6e90035561892be4eb88eac05b8623cb89d":
        raise ValueError("P49 full-curve frozen kernel identity changed")
    if gate.get("size_order") != SIZES or gate.get("lineages_in_order") != LINEAGES:
        raise ValueError("P49 full-curve lineage order changed")
    if gate.get("levels_in_order") != LEVELS:
        raise ValueError("P49 full-curve intrinsic level order changed")
    if gate.get("projector_order") != PROJECTORS:
        raise ValueError("P49 full-curve projector order changed")

    matching = gate["matching_transfer"]
    source = ObservableDescriptor.from_dict(matching["source_descriptor"])
    stored_child = ObservableDescriptor.from_dict(matching["stored_child_descriptor"])
    matching_map = map_observable(source, stored_child)
    expected = matching["exact_registered_stored_child_map"]
    if (matching_map.scale, matching_map.offset) != (
        float(expected["scale"]), float(expected["offset"])
    ) or (matching_map.scale, matching_map.offset) != (-1.0, 0.0):
        raise ValueError("P49 matching child-order map changed")

    projectors: dict[str, dict[str, object]] = {}
    for name in PROJECTORS:
        definition = gate["projectors"][name]
        combination, response, numerator, denominator, parent_order, child_order = (
            PROJECTOR_SPECS[name]
        )
        if (
            definition.get("combination"),
            definition.get("response_coordinate"),
            definition.get("normalization_power_in_N"),
            definition.get("source_orientation_order"),
            definition.get("stored_child_orientation_order"),
        ) != (
            combination,
            response,
            {"numerator": numerator, "denominator": denominator},
            parent_order,
            child_order,
        ):
            raise ValueError(f"P49 {name} projector contract changed")
        parent = _descriptor(combination, parent_order)
        child = _descriptor(combination, child_order)
        transform = map_observable(parent, child)
        expected_projector = gate["projector_registered_map"]
        if (transform.scale, transform.offset) != (
            float(expected_projector["scale"]),
            float(expected_projector["offset"]),
        ) or (transform.scale, transform.offset) != (1.0, 0.0):
            raise ValueError(f"P49 {name} normalized projector map changed")
        projectors[name] = {
            "source_descriptor": parent,
            "stored_child_descriptor": child,
            "transform": transform,
            "response_coordinate": response,
            "normalization_power_in_N": definition["normalization_power_in_N"],
        }
    return gate, {
        "matching_source": source,
        "matching_stored_child": stored_child,
        "matching_map": matching_map,
        "projectors": projectors,
    }


def _run_frozen(histograms: Sequence[Path]) -> dict:
    import score_p49_fullcurve_doubling as frozen

    result = frozen.calculate(frozen.merge_inputs(histograms))
    result["provenance"] = {str(path): frozen.sha256(path) for path in histograms}
    return result


def score_typed(
    root: Path,
    histograms: Sequence[Path],
    *,
    runner: Callable[[Sequence[Path]], dict] = _run_frozen,
) -> dict:
    gate, validated = load_semantic_gate(root)
    if len(histograms) != len(gate["canonical_inputs"]):
        raise ValueError("P49 requires the four canonical histograms")
    for index, (path, expected) in enumerate(
        zip(histograms, gate["canonical_inputs"])
    ):
        _validate_file(path, expected, f"histogram {index}")
    result = runner(histograms)
    if result.get("format_version") != 1 or result.get("batches") != 100:
        raise ValueError("frozen P49 format or batch contract differs from semantic gate")
    if [int(key) for key in result.get("sizes", {})] != SIZES:
        raise ValueError("frozen P49 size order differs from semantic gate")
    if list(result.get("lineages", {})) != ["65->130", "85->170"]:
        raise ValueError("frozen P49 lineage order differs from semantic gate")
    if list(result.get("joint_scores", {})) != gate["joint_score_order"]:
        raise ValueError("frozen P49 joint score order differs from semantic gate")
    replication = result.get("P48_Sprime_fresh_seed_replication", {})
    if list(replication.get("target_sizes", [])) != gate["sprime_replication"]["target_sizes"]:
        raise ValueError("frozen P49 S-prime target order differs from semantic gate")
    if list(replication.get("models_in_frozen_order", {})) != gate["sprime_replication"]["model_order"]:
        raise ValueError("frozen P49 S-prime model order differs from semantic gate")
    if replication.get("classification") != gate["sprime_replication"]["classification"]:
        raise ValueError("frozen P49 S-prime refit boundary differs from semantic gate")
    expected_provenance = {
        row["path"]: row["sha256"] for row in gate["canonical_inputs"]
    }
    if result.get("provenance") != expected_provenance:
        raise ValueError("frozen P49 provenance differs from semantic gate")
    result["observable_semantics"] = {
        "semantic_gate": SEMANTIC_GATE,
        "semantic_gate_status": gate["status"],
        "matching_transfer": {
            "source_descriptor": validated["matching_source"].to_dict(),
            "stored_child_descriptor": validated["matching_stored_child"].to_dict(),
            "applied_stored_child_transform": validated["matching_map"].to_dict(),
            "response_families": gate["matching_transfer"]["response_families"],
        },
        "projector_order": gate["projector_order"],
        "projectors": {
            name: {
                "source_descriptor": values["source_descriptor"].to_dict(),
                "stored_child_descriptor": values["stored_child_descriptor"].to_dict(),
                "applied_transform": values["transform"].to_dict(),
                "response_coordinate": values["response_coordinate"],
                "normalization_power_in_N": values["normalization_power_in_N"],
            }
            for name, values in validated["projectors"].items()
        },
        "lineage_targets": gate["lineage_targets"],
        "batch_contract": gate["batch_contract"],
        "validation_order": "semantic_maps_and_canonical_inputs_before_frozen_fullcurve_score",
        "semantic_boundary": gate["semantic_boundary"],
        "evidence_boundary": gate["evidence_boundary"],
    }
    return result


def main(argv: Sequence[str] | None = None) -> int:
    root = Path(__file__).resolve().parents[1]
    gate = json.loads((root / SEMANTIC_GATE).read_text(encoding="utf-8"))
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--histograms", nargs=4, type=Path,
        default=[root / row["path"] for row in gate["canonical_inputs"]],
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    payload = score_typed(root, args.histograms)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
        print(args.output)
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Type-safe entrypoint for the frozen P231 vacuum-KdV sector score."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Callable, Sequence

from wrapping_channels import ObservableDescriptor, map_observable


SEMANTIC_GATE = "predictions/p231_vacuum_kdv_sector_semantic_gate_20260830.json"
DESIGNS = ["pell_Dminus2_N30", "pell_Dplus1_N56"]
COORDINATES = ["C_nontrivial_real", "Q_reflection_null", "S_scalar"]
JOINT_ORDER = [f"{design}_{coordinate}" for design in DESIGNS for coordinate in COORDINATES]


def load_semantic_gate(root: Path) -> tuple[dict, ObservableDescriptor, ObservableDescriptor, object]:
    gate = json.loads((root / SEMANTIC_GATE).read_text(encoding="utf-8"))
    if gate.get("status") != "semantic_gate_added_after_frozen_p231_vacuum_kdv_score":
        raise ValueError("P231 KdV semantic gate status changed")
    if gate.get("frozen_kernel_git_blob") != "7f4166167d4a37a93912d4816d12c1447923ff59":
        raise ValueError("P231 KdV frozen kernel identity changed")
    if gate.get("designs_in_order") != DESIGNS:
        raise ValueError("P231 KdV design order changed")
    if gate.get("coordinates_in_order") != COORDINATES:
        raise ValueError("P231 KdV coordinate order changed")
    if gate.get("joint_order") != JOINT_ORDER:
        raise ValueError("P231 KdV joint order changed")
    if gate.get("theory_normalization") != "finite_size_design_vector_per_unit_g4":
        raise ValueError("P231 KdV theory normalization changed")
    if gate.get("covariance_contract") != "two_independent_3x3_design_blocks":
        raise ValueError("P231 KdV covariance contract changed")
    if gate.get("non_scalar_C_only_indices") != [0, 3]:
        raise ValueError("P231 KdV C-only projection changed")
    source = ObservableDescriptor.from_dict(gate["source_descriptor"])
    target = ObservableDescriptor.from_dict(gate["target_descriptor"])
    transform = map_observable(source, target)
    expected = gate["exact_registered_map"]
    if (transform.scale, transform.offset) != (
        float(expected["scale"]), float(expected["offset"])
    ) or (transform.scale, transform.offset) != (1.0, 0.0):
        raise ValueError("P231 KdV registered map changed")
    return gate, source, target, transform


def _build_frozen_score(pilot_path: Path, oracle_path: Path) -> dict:
    import score_p231_vacuum_kdv_sector as frozen_kernel

    return frozen_kernel.build_score(pilot_path, oracle_path)


def score_typed(
    root: Path,
    pilot_path: Path,
    oracle_path: Path,
    *,
    runner: Callable[[Path, Path], dict] = _build_frozen_score,
) -> dict:
    gate, source, target, transform = load_semantic_gate(root)
    result = runner(pilot_path, oracle_path)
    if result.get("joint_order") != gate["joint_order"]:
        raise ValueError("frozen P231 KdV joint order differs from semantic gate")
    if len(result.get("observed", [])) != 6 or len(result.get("theory_vector_per_unit_g4", [])) != 6:
        raise ValueError("frozen P231 KdV vector width differs from semantic gate")
    covariance = result.get("covariance", [])
    if len(covariance) != 6 or any(len(row) != 6 for row in covariance):
        raise ValueError("frozen P231 KdV covariance shape differs from semantic gate")
    governance = result.get("governance", {})
    if governance.get("new_independent_evidence") is not False:
        raise ValueError("frozen P231 KdV evidence boundary differs from semantic gate")
    result["observable_semantics"] = {
        "semantic_gate": SEMANTIC_GATE,
        "semantic_gate_status": gate["status"],
        "source_descriptor": source.to_dict(),
        "target_descriptor": target.to_dict(),
        "applied_transform": transform.to_dict(),
        "designs_in_order": gate["designs_in_order"],
        "coordinates_in_order": gate["coordinates_in_order"],
        "theory_normalization": gate["theory_normalization"],
        "covariance_contract": gate["covariance_contract"],
        "non_scalar_C_only_indices": gate["non_scalar_C_only_indices"],
        "validation_order": "semantic_map_before_frozen_CQS_score",
        "evidence_boundary": gate["evidence_boundary"],
    }
    return result


def main(argv: Sequence[str] | None = None) -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot-result", type=Path, required=True)
    parser.add_argument("--oracle", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    payload = score_typed(root, args.pilot_result, args.oracle)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

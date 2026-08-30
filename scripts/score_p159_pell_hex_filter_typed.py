#!/usr/bin/env python3
"""Type-safe entrypoint for the frozen P159 Pell/hex bridge score."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Callable, Sequence

from wrapping_channels import ObservableDescriptor, map_observable


SEMANTIC_GATE = "predictions/p159_pell_hex_filter_semantic_gate_20260830.json"
DESIGNS = ["pell_Dminus2_N30", "pell_Dplus1_N56"]
CATEGORIES = [
    "rank0", "l0", "l1", "l2", "rank1_other", "rank2",
    "invariant_failure",
]
TARGET_LINES = [[1, 0], [0, 1], [1, -1]]
CONTRAST_ORDER = ["C_nontrivial_real", "Q_reflection_null", "S_scalar"]
TRANSFORM = [
    [1.0, -0.5, -0.5],
    [0.0, -math.sqrt(3.0) / 2.0, math.sqrt(3.0) / 2.0],
    [1.0, 1.0, 1.0],
]


def _git_blob_sha(data: bytes) -> str:
    return hashlib.sha1(f"blob {len(data)}\0".encode("ascii") + data).hexdigest()


def _validate_file(path: Path, expected: dict, label: str) -> None:
    data = path.read_bytes()
    if _git_blob_sha(data) != expected["git_blob"]:
        raise ValueError(f"P159 canonical {label} git blob changed")
    if hashlib.sha256(data).hexdigest() != expected["sha256"]:
        raise ValueError(f"P159 canonical {label} sha256 changed")


def load_semantic_gate(
    root: Path,
) -> tuple[dict, ObservableDescriptor, ObservableDescriptor, object]:
    gate = json.loads((root / SEMANTIC_GATE).read_text(encoding="utf-8"))
    if gate.get("status") != "semantic_gate_added_after_frozen_p159_score":
        raise ValueError("P159 semantic gate status changed")
    if gate.get("frozen_kernel_git_blob") != "52dc0a79130ca25e5ced80d176ce28c34b2203ae":
        raise ValueError("P159 frozen kernel identity changed")
    if gate.get("designs_in_order") != DESIGNS:
        raise ValueError("P159 design order changed")
    if gate.get("category_order") != CATEGORIES:
        raise ValueError("P159 category order changed")
    line_contract = gate.get("primitive_line_contract", {})
    if line_contract.get("target_lines_in_order") != TARGET_LINES:
        raise ValueError("P159 primitive line order changed")
    if line_contract.get("transport_Dminus2_to_Dplus1") != [[1, 0], [0, 1]]:
        raise ValueError("P159 basis transport changed")
    if line_contract.get("positive_rho_C3_action") != [[0, -1], [1, 1]]:
        raise ValueError("P159 positive-rho C3 action changed")
    character = gate.get("character_contract", {})
    if character.get("contrast_order") != CONTRAST_ORDER:
        raise ValueError("P159 contrast order changed")
    if character.get("transform_rows") != TRANSFORM:
        raise ValueError("P159 C/Q/S transform changed")
    source = ObservableDescriptor.from_dict(gate["source_descriptor"])
    target = ObservableDescriptor.from_dict(gate["target_descriptor"])
    transform = map_observable(source, target)
    expected = gate["exact_registered_map"]
    if (transform.scale, transform.offset) != (
        float(expected["scale"]), float(expected["offset"])
    ) or (transform.scale, transform.offset) != (1.0, 0.0):
        raise ValueError("P159 registered topology-envelope map changed")
    return gate, source, target, transform


def _build_frozen_score(batch_path: Path, source_path: Path) -> dict:
    import score_p159_pell_hex_filter as frozen_kernel

    return frozen_kernel.build_score(batch_path, source_path)


def score_typed(
    root: Path,
    batch_path: Path,
    source_path: Path,
    *,
    runner: Callable[[Path, Path], dict] = _build_frozen_score,
) -> dict:
    gate, source, target, transform = load_semantic_gate(root)
    inputs = gate["canonical_inputs"]
    _validate_file(batch_path, inputs["batches"], "batches")
    _validate_file(source_path, inputs["source_result"], "source result")
    result = runner(batch_path, source_path)
    result_inputs = result.get("inputs", {})
    if (
        result_inputs.get("batches_sha256") != inputs["batches"]["sha256"]
        or result_inputs.get("source_result_sha256")
        != inputs["source_result"]["sha256"]
    ):
        raise ValueError("frozen P159 reported input hashes differ from semantic gate")
    result["inputs"] = {
        "batches": inputs["batches"]["path"],
        "batches_sha256": inputs["batches"]["sha256"],
        "source_result": inputs["source_result"]["path"],
        "source_result_sha256": inputs["source_result"]["sha256"],
    }
    result_contract = gate["result_contract"]
    for key in ("schema", "issue", "analysis_class"):
        if result.get(key) != result_contract[key]:
            raise ValueError(f"frozen P159 {key} differs from semantic gate")
    if list(result.get("design_scores", {})) != DESIGNS:
        raise ValueError("frozen P159 design order differs from semantic gate")
    if any(
        result["design_scores"].get(design, {}).get("contrast_order") != CONTRAST_ORDER
        for design in DESIGNS
    ):
        raise ValueError("frozen P159 contrast order differs from semantic gate")
    basis = result.get("basis_transport", {})
    line_contract = gate["primitive_line_contract"]
    for key in (
        "normalized_period_basis", "transport_Dminus2_to_Dplus1",
        "positive_rho_C3_action", "ordered_unoriented_line_cycle",
    ):
        if basis.get(key) != line_contract[key]:
            raise ValueError(f"frozen P159 {key} differs from semantic gate")
    decisions = result.get("decision", {})
    for key in (
        "primitive_character_bridge", "ordinary_H4_simple_zero_bridge",
        "square_site_H4_promotion",
    ):
        if decisions.get(key) != result_contract[key]:
            raise ValueError(f"frozen P159 decision {key} differs from semantic gate")
    if result.get("governance") != gate["governance_contract"]:
        raise ValueError("frozen P159 governance differs from semantic gate")
    result["observable_semantics"] = {
        "semantic_gate": SEMANTIC_GATE,
        "semantic_gate_status": gate["status"],
        "source_descriptor": source.to_dict(),
        "target_descriptor": target.to_dict(),
        "applied_transform": transform.to_dict(),
        "typed_scope": "full_configuration_rank_positive_topology_envelope_only",
        "designs_in_order": gate["designs_in_order"],
        "primitive_line_contract": line_contract,
        "character_contract": gate["character_contract"],
        "validation_order": "semantic_gate_and_canonical_inputs_before_frozen_score",
        "semantic_boundary": gate["semantic_boundary"],
        "evidence_boundary": gate["evidence_boundary"],
    }
    return result


def main(argv: Sequence[str] | None = None) -> int:
    root = Path(__file__).resolve().parents[1]
    inputs = json.loads((root / SEMANTIC_GATE).read_text(encoding="utf-8"))["canonical_inputs"]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--batches", type=Path, default=root / inputs["batches"]["path"]
    )
    parser.add_argument(
        "--source-result", type=Path, default=root / inputs["source_result"]["path"]
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    payload = score_typed(root, args.batches, args.source_result)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

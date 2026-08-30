#!/usr/bin/env python3
"""Type-safe entrypoint for the frozen v14 scalar-root projector."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Callable, Sequence

from wrapping_channels import ObservableDescriptor, map_observable


SEMANTIC_GATE = "predictions/v14_scalar_root_projector_semantic_gate_20260830.json"
DEFAULT_LINEAGES = ((65, 130), (85, 170))
FIXED_BETA = 3.5
PROJECTOR = {
    "name": "H4_null_orientation_independent_scalar_root",
    "formula": "(cos4_first*root_second-cos4_second*root_first)/(cos4_first-cos4_second)",
    "requires_nonzero_delta_cos4": True,
}


def load_semantic_gate(root: Path) -> tuple[dict, ObservableDescriptor, ObservableDescriptor, object]:
    gate = json.loads((root / SEMANTIC_GATE).read_text(encoding="utf-8"))
    if gate.get("status") != "semantic_gate_added_after_frozen_v14_scalar_root_projector":
        raise ValueError("v14 scalar-root semantic gate status changed")
    if gate.get("frozen_kernel_git_blob") != "2e4bb607777c6a75a24beeed412640b977622844":
        raise ValueError("v14 scalar-root frozen kernel identity changed")
    if (gate.get("channel"), gate.get("sector"), gate.get("response_coordinate")) != (
        "direction_1", "matching_function", "implicit_matching_root"
    ):
        raise ValueError("v14 scalar-root channel contract changed")
    if gate.get("sizes_in_order") != [65, 85, 130, 170]:
        raise ValueError("v14 scalar-root size order changed")
    if gate.get("lineages_in_order") != [[65, 130], [85, 170]]:
        raise ValueError("v14 scalar-root lineage order changed")
    if gate.get("orientation_order") != ["first", "second"]:
        raise ValueError("v14 scalar-root orientation order changed")
    if gate.get("projector") != PROJECTOR:
        raise ValueError("v14 scalar-root projector changed")
    if gate.get("fixed_beta_in_N") != {"numerator": 7, "denominator": 2}:
        raise ValueError("v14 scalar-root beta changed")
    if gate.get("doubling_q") != "2^(-7/2)":
        raise ValueError("v14 scalar-root doubling ratio changed")
    if gate.get("covariance_contract") != (
        "synchronized_delete_one_batches_with_full_cross_size_covariance"
    ):
        raise ValueError("v14 scalar-root covariance contract changed")
    source = ObservableDescriptor.from_dict(gate["source_descriptor"])
    target = ObservableDescriptor.from_dict(gate["target_descriptor"])
    transform = map_observable(source, target)
    expected = gate["exact_registered_map"]
    if (transform.scale, transform.offset) != (
        float(expected["scale"]), float(expected["offset"])
    ) or (transform.scale, transform.offset) != (1.0, 0.0):
        raise ValueError("v14 scalar-root registered map changed")
    return gate, source, target, transform


def _run_frozen(paths: Sequence[Path], beta: float, lineages: Sequence[tuple[int, int]]) -> dict:
    import score_v14_scalar_root_projector as frozen_kernel

    records = frozen_kernel.merge_inputs(paths)
    return frozen_kernel.calculate(records, beta=beta, lineages=lineages)


def score_typed(
    root: Path,
    histograms: Sequence[Path],
    *,
    beta: float = FIXED_BETA,
    lineages: Sequence[tuple[int, int]] = DEFAULT_LINEAGES,
    runner: Callable[[Sequence[Path], float, Sequence[tuple[int, int]]], dict] = _run_frozen,
) -> dict:
    gate, source, target, transform = load_semantic_gate(root)
    if not math.isclose(float(beta), FIXED_BETA, rel_tol=0.0, abs_tol=1e-15):
        raise ValueError("v14 scalar-root runtime beta differs from semantic gate")
    normalized_lineages = tuple(tuple(pair) for pair in lineages)
    if normalized_lineages != DEFAULT_LINEAGES:
        raise ValueError("v14 scalar-root runtime lineages differ from semantic gate")
    payload = runner(histograms, beta, normalized_lineages)
    if payload.get("format_version") != 1:
        raise ValueError("v14 scalar-root frozen payload version changed")
    hypothesis = payload.get("hypothesis", {})
    if not math.isclose(float(hypothesis.get("beta_in_N", -1.0)), FIXED_BETA):
        raise ValueError("v14 scalar-root frozen hypothesis changed")
    if list(payload.get("sizes", {})) != [65, 85, 130, 170]:
        raise ValueError("v14 scalar-root frozen size order changed")
    if list(payload.get("lineages", {})) != ["65->130", "85->170"]:
        raise ValueError("v14 scalar-root frozen lineage order changed")
    if payload.get("two_lineage_consistency") is None:
        raise ValueError("v14 scalar-root consistency contract changed")
    payload["observable_semantics"] = {
        "semantic_gate": SEMANTIC_GATE,
        "semantic_gate_status": gate["status"],
        "source_descriptor": source.to_dict(),
        "target_descriptor": target.to_dict(),
        "applied_transform": transform.to_dict(),
        "response_coordinate": gate["response_coordinate"],
        "projector": gate["projector"],
        "fixed_beta_in_N": gate["fixed_beta_in_N"],
        "lineages_in_order": gate["lineages_in_order"],
        "covariance_contract": gate["covariance_contract"],
        "validation_order": "semantic_map_before_implicit_root_projection",
        "evidence_boundary": gate["evidence_boundary"],
    }
    return payload


def main(argv: Sequence[str] | None = None) -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--histograms", nargs="+", required=True, type=Path)
    parser.add_argument("--beta", type=float, default=FIXED_BETA)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args(argv)
    payload = score_typed(root, args.histograms, beta=args.beta)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.csv.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    import score_v14_scalar_root_projector as frozen_kernel

    frozen_kernel.write_csv(args.csv, payload)
    args.report.write_text(frozen_kernel.report(payload), encoding="utf-8")
    print(args.json)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

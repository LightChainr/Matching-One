#!/usr/bin/env python3
"""Typed entrypoint for the frozen adjacent-axis annihilator scorer."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Callable, Sequence

import score_axis_pair_annihilator as frozen
from wrapping_channels import ObservableDescriptor, map_observable


SEMANTIC_GATE = "predictions/axis_pair_annihilator_semantic_gate_20260830.json"


def load_semantic_gate(root: Path, kernel_key: str = "base"):
    gate = json.loads((root / SEMANTIC_GATE).read_text(encoding="utf-8"))
    if gate.get("status") != "semantic_gate_added_after_frozen_axis_pair_score":
        raise ValueError("axis-pair semantic gate status changed")
    expected_blobs = {
        "base": "dcb0093cc1bd9b7379d24d92d7be094fea0e62cb",
        "stable": "152b4d62d1c2bbf3f7eaabc1f3214498f2d82e95",
    }
    if gate.get("frozen_kernel_git_blobs") != expected_blobs:
        raise ValueError("axis-pair frozen kernel identity changed")
    if kernel_key not in expected_blobs:
        raise ValueError("unknown axis-pair kernel key")
    if gate.get("matching_reconstruction") != "K_minus_tail + K_plus_tail - 1":
        raise ValueError("axis-pair matching reconstruction changed")
    if gate.get("response_coordinates") != {
        "fixed_p": "F_L(p_ref)=L^(13/4) M_L(p_ref)-(L-1)^(13/4) M_(L-1)(p_ref)",
        "root_crosscheck": "zero of F_L(p)",
    }:
        raise ValueError("axis-pair response coordinates changed")
    if gate.get("p_ref_default") != "0.592746050790":
        raise ValueError("axis-pair default p_ref changed")
    if gate.get("candidate_q_default_order") != [2.0, 3.0, 4.0, 6.0]:
        raise ValueError("axis-pair candidate order changed")
    if gate.get("candidate_map") != "w_ann=4+q":
        raise ValueError("axis-pair root-power map changed")
    source = ObservableDescriptor.from_dict(gate["source_descriptor"])
    target = ObservableDescriptor.from_dict(gate["target_descriptor"])
    transform = map_observable(source, target)
    expected = gate["exact_registered_map"]
    if (transform.scale, transform.offset) != (
        float(expected["scale"]), float(expected["offset"])
    ) or (transform.scale, transform.offset) != (1.0, 0.0):
        raise ValueError("axis-pair descriptor map must be exact identity")
    return gate, source, target, transform


def calculate_typed(
    root: Path,
    paths: Sequence[Path],
    p_ref: float,
    train_max_L: int,
    q_candidates: Sequence[float],
    *,
    kernel_key: str = "base",
    calculator: Callable[[Sequence[Path], float, int, Sequence[float]], dict] = frozen.calculate,
) -> dict:
    gate, source, target, transform = load_semantic_gate(root, kernel_key)
    result = calculator(paths, p_ref, train_max_L, q_candidates)
    if result.get("p_ref") != p_ref or result.get("train_max_L") != train_max_L:
        raise ValueError("axis-pair frozen result coordinate or split changed")
    if result.get("primary_candidate_order") != list(q_candidates):
        raise ValueError("axis-pair frozen result candidate order changed")
    result["observable_semantics"] = {
        "semantic_gate": SEMANTIC_GATE,
        "semantic_gate_status": gate["status"],
        "frozen_kernel": kernel_key,
        "source_descriptor": source.to_dict(),
        "target_descriptor": target.to_dict(),
        "applied_transform": transform.to_dict(),
        "pair_geometry": gate["pair_geometry"],
        "response_coordinates": gate["response_coordinates"],
        "fit_boundary": gate["fit_boundary"],
        "evidence_boundary": gate["evidence_boundary"],
        "validation_order": "semantic_identity_and_axis_contract_before_frozen_score",
    }
    return result


def main(
    argv: Sequence[str] | None = None,
    *,
    kernel_key: str = "base",
    calculator=None,
    reporter=None,
) -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("histograms", nargs="+", type=Path)
    parser.add_argument("--p-ref", type=float, default=0.592746050790)
    parser.add_argument("--train-max-L", type=int, required=True)
    parser.add_argument("--q", nargs="+", type=float, default=[2.0, 3.0, 4.0, 6.0])
    parser.add_argument("--json", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args(argv)
    payload = calculate_typed(
        root,
        args.histograms,
        args.p_ref,
        args.train_max_L,
        args.q,
        kernel_key=kernel_key,
        calculator=calculator or frozen.calculate,
    )
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    args.report.write_text((reporter or frozen.report)(payload), encoding="utf-8")
    print(args.json)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Type-safe entrypoint for the frozen v14 fixed-p scalar projector."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Callable, Sequence

from wrapping_channels import ObservableDescriptor, map_observable


SEMANTIC_GATE = "predictions/v14_fixedp_scalar_projector_semantic_gate_20260830.json"
P_REF = 0.592746050790
PROJECTOR_WEIGHTS = [
    "-cos4_second/(cos4_first-cos4_second)",
    "cos4_first/(cos4_first-cos4_second)",
]


def load_semantic_gate(root: Path) -> tuple[dict, ObservableDescriptor, ObservableDescriptor, object]:
    gate = json.loads((root / SEMANTIC_GATE).read_text(encoding="utf-8"))
    if gate.get("status") != "semantic_gate_added_after_frozen_v14_fixedp_scalar_projector":
        raise ValueError("v14 fixed-p scalar semantic gate status changed")
    if gate.get("frozen_kernel_git_blob") != "c7b745b60614f315955ec5bda5fd0ae152c33f24":
        raise ValueError("v14 fixed-p scalar frozen kernel identity changed")
    if (gate.get("channel"), gate.get("sector"), gate.get("response_coordinate")) != (
        "direction_1", "matching_function", "fixed_p"
    ):
        raise ValueError("v14 fixed-p scalar channel contract changed")
    if not math.isclose(float(gate.get("p_ref")), P_REF, rel_tol=0.0, abs_tol=1e-15):
        raise ValueError("v14 fixed-p scalar p_ref changed")
    if gate.get("projector") != {
        "name": "H4_null_orientation_independent_scalar",
        "weights": PROJECTOR_WEIGHTS,
        "requires_nonzero_delta_cos4": True,
    }:
        raise ValueError("v14 fixed-p scalar projector changed")
    if gate.get("normalization_power_in_N") != {
        "numerator": 25, "denominator": 8
    }:
        raise ValueError("v14 fixed-p scalar normalization changed")
    if gate.get("covariance_rule") != (
        "reconstruct_cov12_from_var_first_var_second_and_var_difference"
    ):
        raise ValueError("v14 fixed-p scalar covariance rule changed")
    source = ObservableDescriptor.from_dict(gate["source_descriptor"])
    target = ObservableDescriptor.from_dict(gate["target_descriptor"])
    transform = map_observable(source, target)
    expected = gate["exact_registered_map"]
    if (transform.scale, transform.offset) != (
        float(expected["scale"]), float(expected["offset"])
    ) or (transform.scale, transform.offset) != (1.0, 0.0):
        raise ValueError("v14 fixed-p scalar registered map changed")
    return gate, source, target, transform


def _load_frozen_rows(path: Path, channel: str) -> list[dict]:
    import score_v14_fixedp_scalar_projector as frozen_kernel

    return frozen_kernel.load_rows(path, channel)


def score_typed(
    root: Path,
    analysis_csv: Path,
    *,
    channel: str = "direction_1",
    p_ref: float = P_REF,
    runner: Callable[[Path, str], list[dict]] = _load_frozen_rows,
) -> dict:
    gate, source, target, transform = load_semantic_gate(root)
    if channel != gate["channel"]:
        raise ValueError("v14 fixed-p scalar runtime channel differs from semantic gate")
    if not math.isclose(float(p_ref), float(gate["p_ref"]), rel_tol=0.0, abs_tol=1e-15):
        raise ValueError("v14 fixed-p scalar runtime p_ref differs from semantic gate")
    rows = runner(analysis_csv, channel)
    if not rows or [row.get("N") for row in rows] != sorted(row.get("N") for row in rows):
        raise ValueError("v14 fixed-p scalar row order differs from semantic gate")
    required = {
        "N", "first", "second", "cos4_first", "cos4_second", "delta_cos4",
        "M_first", "M_second", "M_scalar_H4_null", "M_scalar_se",
        "M_scalar_z", "within_pair_correlation", "N25_8_scaled_scalar",
        "N25_8_scaled_se",
    }
    if any(set(row) != required for row in rows):
        raise ValueError("v14 fixed-p scalar row schema differs from semantic gate")
    payload = {
        "format_version": 1,
        "classification": "retrospective discovery/power diagnostic",
        "p_ref": p_ref,
        "channel": channel,
        "hypothesis": "V_<1,4> scalar: M0(pc) proportional to N^(-25/8)",
        "rows": rows,
        "rules": [
            "Do not infer absence from a non-significant scalar projector.",
            "Do not fit a radial exponent when all per-size scalar z scores are underpowered.",
            "Use smaller N/high statistics or the leading-H4 annihilator for production tests.",
        ],
    }
    payload["observable_semantics"] = {
        "semantic_gate": SEMANTIC_GATE,
        "semantic_gate_status": gate["status"],
        "source_descriptor": source.to_dict(),
        "target_descriptor": target.to_dict(),
        "applied_transform": transform.to_dict(),
        "response_coordinate": gate["response_coordinate"],
        "projector": gate["projector"],
        "normalization_power_in_N": gate["normalization_power_in_N"],
        "covariance_rule": gate["covariance_rule"],
        "validation_order": "semantic_map_before_frozen_row_projection",
        "evidence_boundary": gate["evidence_boundary"],
    }
    return payload


def main(argv: Sequence[str] | None = None) -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("analysis_csv", type=Path)
    parser.add_argument("--channel", default="direction_1")
    parser.add_argument("--p-ref", type=float, default=P_REF)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args(argv)
    payload = score_typed(
        root, args.analysis_csv, channel=args.channel, p_ref=args.p_ref
    )
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(text, encoding="utf-8")
    else:
        print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

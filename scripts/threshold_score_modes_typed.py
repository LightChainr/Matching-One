#!/usr/bin/env python3
"""Type-safe entrypoint for the frozen threshold Krawtchouk score modes."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Callable, Mapping, Sequence

from wrapping_channels import ObservableDescriptor, map_observable


SEMANTIC_GATE = "predictions/threshold_score_modes_semantic_gate_20260830.json"
DEFAULT_MAX_ORDER = 6
EVIDENCE_GUARD = (
    "Mode 0 and mode 1 are exact coordinate views of the existing value and "
    "first derivative. They are not independent evidence blocks."
)
EXPECTED_DESCRIPTORS = {
    "S": {
        "channel": "cross", "combination": "even", "coordinate": "p",
        "orientation_order": "first_minus_second",
        "normalization": "angular_normalized",
        "quantity": "orientation_contrast",
    },
    "D": {
        "channel": "cross", "combination": "odd", "coordinate": "p",
        "orientation_order": "first_minus_second",
        "normalization": "angular_normalized",
        "quantity": "orientation_contrast",
    },
}
EXACT_VIEW_KEYS = [
    "P4_S_mode_0_minus_direct_P4_S",
    "P4_D_mode_0_minus_direct_P4_D",
    "P4_S_prime_from_mode1_minus_direct",
    "P4_D_prime_from_mode1_minus_direct",
]


def load_semantic_gate(root: Path) -> tuple[dict, dict[str, dict[str, object]]]:
    gate = json.loads((root / SEMANTIC_GATE).read_text(encoding="utf-8"))
    if gate.get("status") != "semantic_gate_added_after_frozen_threshold_score_modes":
        raise ValueError("threshold score-mode semantic gate status changed")
    if gate.get("frozen_kernel_git_blob") != "0e8157bad50fb265efb3a0ca9cde1c795f34d995":
        raise ValueError("threshold score-mode frozen kernel identity changed")
    if gate.get("coordinate") != "eta=log(p/(1-p))":
        raise ValueError("threshold score-mode coordinate changed")
    if gate.get("basis") != (
        "orthonormal Bin(N,p0) Krawtchouk; H1=(K-Np)/sqrt(Np(1-p))"
    ):
        raise ValueError("threshold score-mode basis changed")
    if (gate.get("default_max_order"), gate.get("allowed_max_order")) != (
        DEFAULT_MAX_ORDER, [1, 12]
    ):
        raise ValueError("threshold score-mode order contract changed")
    if gate.get("orientation_order") != ["first", "second"]:
        raise ValueError("threshold score-mode orientation order changed")
    if gate.get("angular_denominator") != "cos4_first-cos4_second":
        raise ValueError("threshold score-mode angular denominator changed")
    if gate.get("sector_order") != ["S", "D"]:
        raise ValueError("threshold score-mode sector order changed")
    if gate.get("sector_descriptors") != EXPECTED_DESCRIPTORS:
        raise ValueError("threshold score-mode descriptor changed")
    if gate.get("exact_view_identity_keys") != EXACT_VIEW_KEYS:
        raise ValueError("threshold score-mode exact-view contract changed")
    validated: dict[str, dict[str, object]] = {}
    for sector in gate["sector_order"]:
        source = ObservableDescriptor.from_dict(gate["sector_descriptors"][sector])
        target = ObservableDescriptor.from_dict(gate["sector_descriptors"][sector])
        transform = map_observable(source, target)
        expected = gate["exact_registered_map"]
        if (transform.scale, transform.offset) != (
            float(expected["scale"]), float(expected["offset"])
        ) or (transform.scale, transform.offset) != (1.0, 0.0):
            raise ValueError("threshold score-mode registered map changed")
        validated[sector] = {
            "source_descriptor": source,
            "target_descriptor": target,
            "transform": transform,
        }
    return gate, validated


def _run_frozen(histograms: Sequence[Path], max_order: int, dps: int) -> dict:
    import mpmath as mp
    import threshold_score_modes as frozen_kernel

    mp.mp.dps = dps
    data = {}
    for path in histograms:
        block = frozen_kernel.read(path)
        overlap = set(data).intersection(block)
        if overlap:
            raise ValueError("duplicate histogram keys")
        data.update(block)
    return {
        "schema": "matching-one/threshold-krawtchouk-score-modes/v1",
        "coordinate": "eta=log(p/(1-p))",
        "basis": (
            "orthonormal Bin(N,p0) Krawtchouk; "
            "H1=(K-Np)/sqrt(Np(1-p))"
        ),
        "max_order": max_order,
        "inputs": [
            {"path": str(path), "sha256": frozen_kernel.sha256(path)}
            for path in histograms
        ],
        "by_N": {
            str(n): frozen_kernel.analyze_size(data, n, max_order)
            for n in sorted({key[0] for key in data})
        },
        "evidence_guard": EVIDENCE_GUARD,
    }


def score_typed(
    root: Path,
    histograms: Sequence[Path],
    *,
    max_order: int = DEFAULT_MAX_ORDER,
    dps: int = 50,
    runner: Callable[[Sequence[Path], int, int], dict] = _run_frozen,
) -> dict:
    gate, validated = load_semantic_gate(root)
    lower, upper = gate["allowed_max_order"]
    if not lower <= max_order <= upper:
        raise ValueError("threshold score-mode runtime order is outside semantic gate")
    if dps <= 0:
        raise ValueError("threshold score-mode precision must be positive")
    result = runner(histograms, max_order, dps)
    if result.get("schema") != "matching-one/threshold-krawtchouk-score-modes/v1":
        raise ValueError("threshold score-mode frozen schema changed")
    if (result.get("coordinate"), result.get("basis")) != (
        gate["coordinate"], gate["basis"]
    ):
        raise ValueError("threshold score-mode frozen coordinate changed")
    if result.get("max_order") != max_order:
        raise ValueError("threshold score-mode frozen order changed")
    if result.get("evidence_guard") != EVIDENCE_GUARD:
        raise ValueError("threshold score-mode evidence guard changed")
    expected_modes = [
        *(f"P4_S_mode_{order}" for order in range(max_order + 1)),
        *(f"P4_D_mode_{order}" for order in range(max_order + 1)),
    ]
    for key, row in result.get("by_N", {}).items():
        if str(row.get("N")) != key:
            raise ValueError("threshold score-mode size key changed")
        if row.get("mode_order") != expected_modes:
            raise ValueError("threshold score-mode frozen mode order changed")
        if list(row.get("exact_view_identities", {})) != EXACT_VIEW_KEYS:
            raise ValueError("threshold score-mode exact-view identities changed")
        tower = row.get("parity_tower_scaled", {})
        if len(tower) != 2 * (max_order + 1):
            raise ValueError("threshold score-mode parity tower changed")
    result["observable_semantics"] = {
        "semantic_gate": SEMANTIC_GATE,
        "semantic_gate_status": gate["status"],
        "sector_order": gate["sector_order"],
        "sector_maps": {
            sector: {
                "source_descriptor": values["source_descriptor"].to_dict(),
                "target_descriptor": values["target_descriptor"].to_dict(),
                "applied_transform": values["transform"].to_dict(),
            }
            for sector, values in validated.items()
        },
        "intrinsic_center": gate["intrinsic_center"],
        "orientation_order": gate["orientation_order"],
        "angular_denominator": gate["angular_denominator"],
        "mode_views": gate["mode_views"],
        "parity_tower": gate["parity_tower"],
        "covariance_contract": gate["covariance_contract"],
        "validation_order": "semantic_maps_before_frozen_mode_projection",
        "evidence_boundary": gate["evidence_boundary"],
    }
    return result


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("histograms", nargs="+", type=Path)
    parser.add_argument("--max-order", type=int, default=DEFAULT_MAX_ORDER)
    parser.add_argument("--dps", type=int, default=50)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    root = Path(__file__).resolve().parents[1]
    payload = score_typed(
        root, args.histograms, max_order=args.max_order, dps=args.dps
    )
    rendered = json.dumps(payload, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
        print(args.output)
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

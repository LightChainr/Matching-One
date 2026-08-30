#!/usr/bin/env python3
"""Validate the frozen P50 full-curve quantity-family contract."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Mapping

from wrapping_channels import ObservableDescriptor, map_observable


CONTRACT = "analysis/p50_fullcurve_quantity_contract.json"
EXPECTED_FEATURE_ORDER = [
    "X_even_0.0", "X_even_0.025", "X_even_0.05", "mean_slope",
    "root_gap_lineage", "P4_S", "P4_D", "P4_S_prime", "P4_D_prime",
]
EXPECTED_SCORING_ORDER = [
    "score the existing P50/P49 primary full-curve DeltaM lineage residual first",
    "score the raw asymptotic slope ratio 2^(3/8) as the historical baseline",
    "score this frozen scalar-plus-H4 1/N slope correction without refitting",
    "score the induced finite-slope root ratio",
    "only afterward fit extra powers/logarithms or orientation terms",
]
EXPECTED_DESCRIPTORS = {
    "thermal_even": {
        "channel": "cross", "combination": "even", "coordinate": "p",
        "orientation_order": "first_minus_second", "normalization": "raw",
        "quantity": "orientation_contrast",
    },
    "matching_function": {
        "channel": "cross", "combination": "matching", "coordinate": "p",
        "orientation_order": "none", "normalization": "raw",
        "quantity": "value",
    },
    "P4_S": {
        "channel": "cross", "combination": "even", "coordinate": "p",
        "orientation_order": "first_minus_second",
        "normalization": "angular_normalized",
        "quantity": "orientation_contrast",
    },
    "P4_D": {
        "channel": "cross", "combination": "odd", "coordinate": "p",
        "orientation_order": "first_minus_second",
        "normalization": "angular_normalized",
        "quantity": "orientation_contrast",
    },
}


def git_blob(path: Path) -> str:
    data = path.read_bytes()
    return hashlib.sha1(f"blob {len(data)}\0".encode("ascii") + data).hexdigest()


def load_contract(root: Path) -> tuple[dict, dict[str, dict[str, object]]]:
    contract = json.loads((root / CONTRACT).read_text(encoding="utf-8"))
    if contract.get("status") != "quantity_contract_frozen_before_typed_wrapper":
        raise ValueError("P50 quantity-contract status changed")
    if contract.get("frozen_kernel") != {
        "path": "scripts/score_p50_fullcurve_n290.py",
        "git_blob": "5008449033a85d2c81bf4ea3f025fa61217b6c4a",
    }:
        raise ValueError("P50 frozen kernel identity changed")
    if contract.get("frozen_prediction") != {
        "path": "predictions/p49_slope_two_sector_145_290_20260828.yaml",
        "git_blob": "df4d2f0c3ca5f3e3380906fd2b6636507574d108",
        "status": "frozen_before_fullcurve_N290_reveal",
    }:
        raise ValueError("P50 frozen prediction identity changed")
    if contract.get("sizes_in_order") != [145, 290]:
        raise ValueError("P50 size order changed")
    if contract.get("representations") != {
        "145": {"first": [12, 1], "second": [9, 8]},
        "290": {"first": [13, 11], "second": [17, 1]},
    }:
        raise ValueError("P50 representation order changed")
    if contract.get("rng_relation") != "independent_parent_and_child_streams":
        raise ValueError("P50 RNG relation changed")
    if contract.get("lineage_sign") != {"145": 1.0, "290": 1.0}:
        raise ValueError("P50 lineage sign changed")
    if contract.get("feature_order") != EXPECTED_FEATURE_ORDER:
        raise ValueError("P50 feature order changed")
    if contract.get("scoring_order") != EXPECTED_SCORING_ORDER:
        raise ValueError("P50 scoring order changed")
    if contract.get("topology_anchors") != EXPECTED_DESCRIPTORS:
        raise ValueError("P50 topology anchor changed")
    validated: dict[str, dict[str, object]] = {}
    for name, payload in contract["topology_anchors"].items():
        source = ObservableDescriptor.from_dict(payload)
        target = ObservableDescriptor.from_dict(payload)
        transform = map_observable(source, target)
        expected = contract["exact_registered_map"]
        if (transform.scale, transform.offset) != (
            float(expected["scale"]), float(expected["offset"])
        ) or (transform.scale, transform.offset) != (1.0, 0.0):
            raise ValueError("P50 registered topology map changed")
        validated[name] = {
            "source_descriptor": source,
            "target_descriptor": target,
            "transform": transform,
        }
    return contract, validated


def validate_repository_files(root: Path, contract: Mapping[str, object]) -> None:
    for key in ("frozen_kernel", "frozen_prediction"):
        record = contract[key]
        path = root / record["path"]
        if git_blob(path) != record["git_blob"]:
            raise ValueError("P50 repository file identity changed: " + str(path))


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    contract, validated = load_contract(root)
    validate_repository_files(root, contract)
    print(json.dumps({
        "status": "PASS",
        "contract": CONTRACT,
        "validated_topology_anchors": list(validated),
        "feature_count": len(contract["feature_order"]),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

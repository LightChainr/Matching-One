#!/usr/bin/env python3
"""Fail-closed applicability gate for observer-bandwidth theorems."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Mapping


def classify(contract: Mapping[str, object]) -> dict[str, object]:
    """Classify one declared observer/clock contract without inference."""

    measure = contract.get("base_measure")
    clock = contract.get("clock")
    reasons = []

    if not contract.get("complete_source_transformed", False):
        reasons.append("complete_source_transform_missing")
    if contract.get("adaptive_mark", False) and not contract.get(
        "conditional_degree_derived", False
    ):
        reasons.append("adaptive_mark_degree_not_derived")
    if contract.get("degenerate_quotient", False) and not contract.get(
        "incidence_polynomial_supplied", False
    ):
        reasons.append("degenerate_quotient_incidence_missing")
    if contract.get("conditional_mean_required", False) and not contract.get(
        "conditional_mean_available", False
    ):
        reasons.append("conditional_mean_missing")
    if contract.get("coupled_noise_required", False) and not contract.get(
        "coupled_noise_levels_available", False
    ):
        reasons.append("coupled_noise_statistics_missing")

    if measure == "bernoulli_product":
        if clock != "product_resample":
            reasons.append("product_resample_clock_missing")
        if contract.get("marked_birth", False):
            reasons.append("palm_conditioning_not_removed")
        accepted = "accept_product_walsh"
        theorem = "product_walsh"
    elif measure == "fixed_k_slice":
        if clock != "occupied_empty_swap":
            reasons.append("johnson_swap_clock_missing")
        n = contract.get("n")
        k = contract.get("k")
        if isinstance(n, int) and isinstance(k, int) and k in (0, n):
            if not contract.get("endpoint_handling", False):
                reasons.append("slice_endpoint_handling_missing")
        accepted = "accept_johnson_slice"
        theorem = "johnson_slice"
    elif measure == "palm_conditioned":
        theorem = "custom_kernel_only"
        if not contract.get("actual_reversible_kernel"):
            reasons.append("palm_reversible_kernel_missing")
        accepted = "accept_custom_kernel_only"
    else:
        theorem = None
        accepted = "reject"
        reasons.append("unsupported_base_measure")

    return {
        "accepted": not reasons,
        "status": accepted if not reasons else "reject",
        "theorem": theorem if not reasons else None,
        "reasons": sorted(set(reasons)),
    }


def build_report(manifest: Mapping[str, object]) -> dict[str, object]:
    fixtures = []
    for fixture in manifest["fixtures"]:
        result = classify(fixture["contract"])
        if result != fixture["expected"]:
            raise AssertionError("fixture classification drift")
        fixtures.append({"name": fixture["name"], **result})

    p267 = []
    for layer in manifest["p267_layers"]:
        result = classify(layer["contract"])
        if result != layer["expected"]:
            raise AssertionError("P267 classification drift")
        p267.append({"name": layer["name"], **result})

    return {
        "schema": manifest["schema"],
        "status": "fail_closed_contract_verified",
        "fixture_count": len(fixtures),
        "fixtures": fixtures,
        "p267_provenance": manifest["p267_provenance"],
        "p267_layers": p267,
        "boundary": manifest["boundary"],
    }


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=root / "analysis/observer_bandwidth_semantic_gate_manifest.json",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = build_report(json.loads(args.manifest.read_text(encoding="utf-8")))
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

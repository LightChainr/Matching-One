#!/usr/bin/env python3
"""Validate the bounded Issue 1 search contract without running a search."""

from __future__ import annotations

import argparse
from decimal import Decimal, InvalidOperation
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = ROOT / "analysis" / "pslq_search_contract.json"
FORBIDDEN_KEYS = frozenset(
    {
        "preferred_point",
        "target_value",
        "combined_interval",
        "exclusions",
        "near_relations",
        "retained_relations",
        "search_results",
    }
)
EXPECTED_STAGES = frozenset(
    {"algebraic_polynomial", "standard_constant_pairwise", "lattice_native_candidates"}
)


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _decimal(value: Any, label: str) -> Decimal:
    _require(isinstance(value, str) and value.strip() == value, "%s must be an exact string" % label)
    try:
        parsed = Decimal(value)
    except InvalidOperation as exc:
        raise ValueError("%s is not a decimal" % label) from exc
    _require(parsed.is_finite(), "%s must be finite" % label)
    return parsed


def _walk_forbidden(value: Any, path: str = "$") -> None:
    if isinstance(value, Mapping):
        bad = sorted(FORBIDDEN_KEYS.intersection(value))
        _require(not bad, "%s contains forbidden result/point keys: %s" % (path, ",".join(bad)))
        for key, child in value.items():
            _walk_forbidden(child, "%s.%s" % (path, key))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            _walk_forbidden(child, "%s[%d]" % (path, index))


def _quoted_estimate(source: Mapping[str, Any], value_text: str) -> Mapping[str, Any]:
    matches = [row for row in source.get("quoted_estimates", []) if row.get("value_text") == value_text]
    _require(len(matches) == 1, "value_text must select exactly one canonical quoted estimate")
    return matches[0]


def _validate_library(stage: Mapping[str, Any], label: str) -> int:
    _require(stage.get("library_frozen_before_search") is True, "%s library is not frozen" % label)
    library = stage.get("library")
    _require(isinstance(library, list) and library, "%s library must be nonempty" % label)
    ids = [row.get("id") for row in library]
    expressions = [row.get("expression") for row in library]
    _require(all(isinstance(item, str) and item for item in ids), "%s ids must be nonempty" % label)
    _require(all(isinstance(item, str) and item for item in expressions), "%s expressions must be nonempty" % label)
    _require(len(ids) == len(set(ids)), "%s has duplicate ids" % label)
    _require(len(expressions) == len(set(expressions)), "%s has duplicate expressions" % label)
    return len(library)


def validate_contract(
    contract: Mapping[str, Any],
    provenance: Mapping[str, Any],
    provenance_bytes: bytes,
) -> dict[str, Any]:
    """Return a compact audit summary or raise ValueError on contract drift."""

    _walk_forbidden(contract)
    _require(contract.get("schema") == "matching-one/pslq-search-contract/v1", "unknown schema")
    _require(contract.get("issue") == 1, "wrong issue")
    _require(contract.get("status") == "protocol_only_no_search_results", "contract is not protocol-only")

    provenance_contract = contract.get("provenance", {})
    digest = hashlib.sha256(provenance_bytes).hexdigest()
    _require(provenance_contract.get("path") == "data/literature_threshold_sources.json", "wrong provenance path")
    _require(provenance_contract.get("sha256") == digest, "provenance SHA-256 mismatch")

    policy = contract.get("input_policy", {})
    _require(policy.get("method_specific_intervals_only") is True, "method-specific interval policy required")
    _require(policy.get("combine_inconsistent_intervals") is False, "interval combination must be disabled")
    _require(policy.get("rounded_estimate_is_definition") is False, "rounded estimate cannot define p_c")
    _require(_decimal(policy.get("comparison_resolution_floor"), "resolution floor") > 0, "invalid resolution floor")

    sources = {source["id"]: source for source in provenance.get("sources", [])}
    intervals = contract.get("intervals")
    _require(isinstance(intervals, list) and len(intervals) >= 3, "at least three method-specific intervals required")
    interval_ids = []
    selected = []
    for index, interval in enumerate(intervals):
        label = "interval[%d]" % index
        interval_id = interval.get("id")
        _require(isinstance(interval_id, str) and interval_id, "%s id is required" % label)
        interval_ids.append(interval_id)
        source_id = interval.get("source_id")
        _require(source_id in sources, "%s has unknown source_id" % label)
        estimate = _quoted_estimate(sources[source_id], interval.get("value_text"))
        _require(interval.get("central_value") == estimate.get("central_value"), "%s central value drift" % label)
        _require(interval.get("quoted_uncertainty") == estimate.get("quoted_uncertainty"), "%s uncertainty drift" % label)
        center = _decimal(interval.get("central_value"), "%s central_value" % label)
        uncertainty = _decimal(interval.get("quoted_uncertainty"), "%s quoted_uncertainty" % label)
        lower = _decimal(interval.get("lower"), "%s lower" % label)
        upper = _decimal(interval.get("upper"), "%s upper" % label)
        _require(uncertainty > 0, "%s uncertainty must be positive" % label)
        _require(lower == center - uncertainty, "%s lower endpoint drift" % label)
        _require(upper == center + uncertainty, "%s upper endpoint drift" % label)
        _require(lower < upper, "%s interval is empty" % label)
        _require(interval.get("confidence_homogenized") is False, "%s confidence was homogenized" % label)
        selected.append({"id": interval_id, "source_id": source_id, "lower": str(lower), "upper": str(upper)})
    _require(len(interval_ids) == len(set(interval_ids)), "duplicate interval ids")

    stages = contract.get("search_stages", {})
    _require(set(stages) == EXPECTED_STAGES, "search stages changed")
    algebraic = stages["algebraic_polynomial"]
    _require(algebraic.get("degree_min") == 1, "algebraic degree must start at one")
    degree_max = algebraic.get("degree_max")
    height = algebraic.get("coefficient_height_max")
    _require(isinstance(degree_max, int) and 1 <= degree_max <= 8, "unsafe algebraic degree bound")
    _require(isinstance(height, int) and 1 <= height <= 10000, "unsafe algebraic height bound")
    _require(algebraic.get("primitive_coefficients_only") is True, "primitive coefficients required")
    _require(algebraic.get("search_each_interval_separately") is True, "intervals must be searched separately")
    standard_count = _validate_library(stages["standard_constant_pairwise"], "standard constants")
    lattice_count = _validate_library(stages["lattice_native_candidates"], "lattice-native candidates")

    arithmetic = contract.get("arithmetic", {})
    initial = arithmetic.get("initial_decimal_digits")
    confirmation = arithmetic.get("confirmation_decimal_digits")
    _require(isinstance(initial, int) and initial >= 50, "initial arithmetic precision is too low")
    _require(isinstance(confirmation, int) and confirmation >= 2 * initial, "confirmation precision must double")
    _require(arithmetic.get("binary_float_exclusion_claims_allowed") is False, "binary-float exclusions forbidden")
    _require("does not contain zero" in arithmetic.get("zero_test", ""), "certified zero test missing")

    controls = contract.get("false_positive_controls", {})
    _require(controls.get("interval_stability_points") == ["lower", "midpoint", "upper"], "interval stability grid changed")
    _require(controls.get("repeat_at_confirmation_precision") is True, "precision repeat required")
    _require(controls.get("report_height_and_conditioning") is True, "conditioning report required")
    synthetic = controls.get("synthetic_random_constants", {})
    _require(isinstance(synthetic.get("count"), int) and synthetic["count"] >= 100, "too few synthetic controls")
    _require(isinstance(synthetic.get("seed"), int), "synthetic seed must be frozen")
    _require(controls.get("matching_partner_is_independent_evidence") is False, "matching partner is not independent evidence")
    _require(bool(controls.get("exact_percolation_control", {}).get("expected_relation")), "exact control missing")

    result_policy = contract.get("result_policy", {})
    _require(result_policy.get("no_results_are_contained_in_this_contract") is True, "results leaked into protocol")
    _require(result_policy.get("bounded_exclusion_implies_transcendence") is False, "transcendence overclaim")
    _require(contract.get("claim_boundary", {}).get("parent_issue") == "remain open", "parent issue boundary changed")
    return {
        "schema": contract["schema"],
        "status": "valid_protocol_only",
        "provenance_sha256": digest,
        "method_specific_interval_count": len(intervals),
        "intervals": selected,
        "algebraic_degree_max": degree_max,
        "algebraic_coefficient_height_max": height,
        "standard_constant_count": standard_count,
        "lattice_native_candidate_count": lattice_count,
        "initial_decimal_digits": initial,
        "confirmation_decimal_digits": confirmation,
        "contains_search_results": False,
    }


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument("--format", choices=("json", "text"), default="text")
    args = parser.parse_args(argv)
    contract = load_json(args.contract)
    provenance_path = ROOT / contract["provenance"]["path"]
    provenance_bytes = provenance_path.read_bytes()
    summary = validate_contract(contract, json.loads(provenance_bytes), provenance_bytes)
    if args.format == "json":
        print(json.dumps(summary, indent=2, sort_keys=True))
    else:
        print(
            "valid protocol-only contract: %d method intervals; degree<=%d; height<=%d; %d/%d frozen constants"
            % (
                summary["method_specific_interval_count"],
                summary["algebraic_degree_max"],
                summary["algebraic_coefficient_height_max"],
                summary["standard_constant_count"],
                summary["lattice_native_candidate_count"],
            )
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

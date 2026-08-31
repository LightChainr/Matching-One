#!/usr/bin/env python3
"""Rank declared next-row intervals by exact certificate-boundary separation per cost."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results/model-certificates/framework/next-row-ranking/latest.json"
SCHEMA = "matching-one/exact-certificate-boundary-ranking/v1"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _interval(values: Sequence[Any], label: str) -> tuple[Fraction, Fraction]:
    _require(isinstance(values, list) and len(values) == 2, f"{label} must have two endpoints")
    lower, upper = (Fraction(value) for value in values)
    _require(lower <= upper, f"{label} endpoints are reversed")
    return lower, upper


def interval_separation(left: tuple[Fraction, Fraction], right: tuple[Fraction, Fraction]) -> Fraction:
    if left[1] < right[0]:
        return right[0] - left[1]
    if right[1] < left[0]:
        return left[0] - right[1]
    return Fraction()


def rank_candidates(candidates: Sequence[Mapping[str, Any]]) -> list[Mapping[str, Any]]:
    _require(isinstance(candidates, list) and candidates, "candidate list cannot be empty")
    seen_ids = set()
    scored = []
    for candidate in candidates:
        _require(set(candidate) == {"id", "feasible_interval", "forecast_interval", "cost"}, "candidate fields drift")
        candidate_id = candidate["id"]
        _require(isinstance(candidate_id, str) and candidate_id and candidate_id not in seen_ids, "candidate ids must be nonempty and unique")
        seen_ids.add(candidate_id)
        feasible = _interval(candidate["feasible_interval"], "feasible interval")
        forecast = _interval(candidate["forecast_interval"], "forecast interval")
        cost = Fraction(candidate["cost"])
        _require(cost > 0, "candidate cost must be positive")
        separation = interval_separation(feasible, forecast)
        scored.append({
            "id": candidate_id,
            "feasible_interval": [str(value) for value in feasible],
            "forecast_interval": [str(value) for value in forecast],
            "cost": str(cost),
            "disjoint": separation > 0,
            "separation_margin": str(separation),
            "margin_per_cost": str(separation / cost),
            "_margin": separation,
            "_utility": separation / cost,
        })
    scored.sort(key=lambda item: (-item["_utility"], -item["_margin"], item["id"]))
    for rank, item in enumerate(scored, start=1):
        item["rank"] = rank
        del item["_margin"]
        del item["_utility"]
    return scored


def frozen_candidates() -> list[dict[str, Any]]:
    return [
        {"id": "morphism-sensitive-row", "feasible_interval": ["0", "0"], "forecast_interval": ["1", "1"], "cost": "4"},
        {"id": "semantic-zero-row", "feasible_interval": ["0", "0"], "forecast_interval": ["1/7", "1/7"], "cost": "1"},
        {"id": "overlapping-endpoint-row", "feasible_interval": ["-1/10", "1/10"], "forecast_interval": ["0", "1/5"], "cost": "1/2"},
    ]


def build_result() -> dict[str, Any]:
    candidates = frozen_candidates()
    ranking = rank_candidates(candidates)
    return {
        "schema": SCHEMA,
        "issue": 370,
        "claim_level": "exact_design_control",
        "dependency_group": "synthetic-exact-certificate-boundary-ranking",
        "utility": "closed-interval separation margin divided by declared positive cost",
        "candidates": candidates,
        "ranking": ranking,
        "verification": {
            "candidate_count": len(ranking),
            "disjoint_count": sum(item["disjoint"] for item in ranking),
            "top_candidate": ranking[0]["id"],
            "status": "exact_boundary_ranking_verified",
        },
        "claim_boundary": {
            "included": "exact deterministic ranking for supplied one-dimensional feasible/forecast intervals and declared rational costs",
            "excluded": "forecast probabilities, variance estimation, multidimensional certificate sensitivity, candidate generation, CPU benchmarking, acquisition authorization, or physical conclusions",
            "parent_issue": "remain open",
        },
    }


def validate_result(result: Mapping[str, Any]) -> Mapping[str, Any]:
    expected = build_result()
    _require(result == expected, "next-row ranking artifact does not exactly reproduce")
    return {
        "schema": result["schema"],
        "status": "valid_exact_certificate_boundary_ranking",
        "candidate_count": result["verification"]["candidate_count"],
        "disjoint_count": result["verification"]["disjoint_count"],
        "top_candidate": result["verification"]["top_candidate"],
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--validate", type=Path)
    args = parser.parse_args(argv)
    if args.validate:
        value = json.loads(args.validate.read_text(encoding="utf-8"))
        print(json.dumps(validate_result(value), indent=2, sort_keys=True))
        return 0
    rendered = json.dumps(build_result(), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

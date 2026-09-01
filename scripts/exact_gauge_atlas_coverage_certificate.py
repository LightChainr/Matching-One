#!/usr/bin/env python3
"""Verify exact coverage of a declared one-dimensional rational gauge atlas."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results/model-certificates/framework/gauge-atlas-1d/latest.json"
SCHEMA = "matching-one/exact-one-dimensional-gauge-atlas/v1"


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _parse_interval(value: Mapping[str, Any], label_required: bool) -> dict[str, Any]:
    expected = {"lower", "upper", "lower_closed", "upper_closed"} | ({"id"} if label_required else set())
    _require(set(value) == expected, "interval fields drift")
    lower = Fraction(value["lower"])
    upper = Fraction(value["upper"])
    _require(lower <= upper, "interval endpoints are reversed")
    _require(isinstance(value["lower_closed"], bool) and isinstance(value["upper_closed"], bool), "interval closure flags must be Boolean")
    if lower == upper:
        _require(value["lower_closed"] and value["upper_closed"], "point interval must be closed")
    result = {"lower": lower, "upper": upper, "lower_closed": value["lower_closed"], "upper_closed": value["upper_closed"]}
    if label_required:
        _require(isinstance(value["id"], str) and value["id"], "chart id is required")
        result["id"] = value["id"]
    return result


def _render_interval(value: Mapping[str, Any]) -> Mapping[str, Any]:
    result = {"lower": str(value["lower"]), "upper": str(value["upper"]), "lower_closed": value["lower_closed"], "upper_closed": value["upper_closed"]}
    if "id" in value:
        result["id"] = value["id"]
    return result


def verify_atlas(descriptor: Mapping[str, Any]) -> Mapping[str, Any]:
    _require(set(descriptor) == {"domain", "charts"}, "descriptor fields drift")
    domain = _parse_interval(descriptor["domain"], False)
    charts = [_parse_interval(chart, True) for chart in descriptor["charts"]]
    _require(charts and len({chart["id"] for chart in charts}) == len(charts), "chart ids must be nonempty and unique")
    for chart in charts:
        _require(domain["lower"] <= chart["lower"] <= chart["upper"] <= domain["upper"], "chart lies outside declared domain")
        if chart["lower"] == domain["lower"] and chart["lower_closed"] and not domain["lower_closed"]:
            raise ValueError("chart includes excluded lower domain endpoint")
        if chart["upper"] == domain["upper"] and chart["upper_closed"] and not domain["upper_closed"]:
            raise ValueError("chart includes excluded upper domain endpoint")
    ordered = sorted(charts, key=lambda item: (item["lower"], not item["lower_closed"], item["upper"], item["id"]))
    merged: list[dict[str, Any]] = []
    for chart in ordered:
        current = {key: chart[key] for key in ("lower", "upper", "lower_closed", "upper_closed")}
        if not merged:
            merged.append(current)
            continue
        previous = merged[-1]
        joins = current["lower"] < previous["upper"] or (
            current["lower"] == previous["upper"] and (current["lower_closed"] or previous["upper_closed"])
        )
        if not joins:
            merged.append(current)
            continue
        if current["upper"] > previous["upper"]:
            previous["upper"] = current["upper"]
            previous["upper_closed"] = current["upper_closed"]
        elif current["upper"] == previous["upper"]:
            previous["upper_closed"] = previous["upper_closed"] or current["upper_closed"]
    gaps = []
    first = merged[0]
    if first["lower"] > domain["lower"] or (first["lower"] == domain["lower"] and domain["lower_closed"] and not first["lower_closed"]):
        gaps.append({"lower": domain["lower"], "upper": first["lower"], "lower_closed": domain["lower_closed"], "upper_closed": not first["lower_closed"]})
    for left, right in zip(merged, merged[1:]):
        gaps.append({"lower": left["upper"], "upper": right["lower"], "lower_closed": not left["upper_closed"], "upper_closed": not right["lower_closed"]})
    last = merged[-1]
    if last["upper"] < domain["upper"] or (last["upper"] == domain["upper"] and domain["upper_closed"] and not last["upper_closed"]):
        gaps.append({"lower": last["upper"], "upper": domain["upper"], "lower_closed": not last["upper_closed"], "upper_closed": domain["upper_closed"]})
    gaps = [gap for gap in gaps if gap["lower"] < gap["upper"] or (gap["lower"] == gap["upper"] and gap["lower_closed"] and gap["upper_closed"])]
    return {
        "chart_count": len(charts),
        "merged_components": [_render_interval(value) for value in merged],
        "uncovered_components": [_render_interval(value) for value in gaps],
        "complete_for_declared_domain": not gaps,
        "status": "exact_one_dimensional_gauge_atlas_checked",
    }


def frozen_descriptor() -> dict[str, Any]:
    return {
        "domain": {"lower": "-1", "upper": "1", "lower_closed": True, "upper_closed": True},
        "charts": [
            {"id": "negative-minor", "lower": "-1", "upper": "0", "lower_closed": True, "upper_closed": False},
            {"id": "boundary", "lower": "0", "upper": "0", "lower_closed": True, "upper_closed": True},
            {"id": "positive-minor", "lower": "0", "upper": "1", "lower_closed": False, "upper_closed": True},
        ],
    }


def build_result() -> dict[str, Any]:
    descriptor = frozen_descriptor()
    verification = verify_atlas(descriptor)
    _require(verification["complete_for_declared_domain"], "frozen atlas lost complete coverage")
    return {
        "schema": SCHEMA,
        "issue": 370,
        "claim_level": "exact_control",
        "dependency_group": "synthetic-exact-gauge-atlas-1d",
        "descriptor": descriptor,
        "verification": verification,
        "claim_boundary": {
            "included": "exact union coverage and uncovered-component calculation for supplied rational one-dimensional chart intervals with endpoint semantics",
            "excluded": "matrix gauge construction, multivariate semialgebraic coverage, algebraic boundaries, physical parameter domains, or atlas discovery",
            "parent_issue": "remain open",
        },
    }


def validate_result(result: Mapping[str, Any]) -> Mapping[str, Any]:
    expected = build_result()
    _require(result == expected, "gauge-atlas certificate does not exactly reproduce")
    return {"schema": result["schema"], "status": "valid_exact_one_dimensional_gauge_atlas", "chart_count": result["verification"]["chart_count"]}


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--validate", type=Path)
    args = parser.parse_args(argv)
    if args.validate:
        print(json.dumps(validate_result(json.loads(args.validate.read_text(encoding="utf-8"))), indent=2, sort_keys=True))
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

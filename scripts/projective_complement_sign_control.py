#!/usr/bin/env python3
"""Exact complement/Alexander sign control for projective birth streams."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

from projective_birth_crosswalk import BirthRow, exact_fraction, fraction_text, reconstruct_at_threshold


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "analysis" / "projective_complement_sign_certificate.json"
SCHEMA = "matching-one/projective-complement-sign/v1"


def validate_line_duality(line_duality: Mapping[str, str], required_lines: set[str]) -> None:
    if not isinstance(line_duality, Mapping):
        raise ValueError("line_duality must be a mapping")
    if any(
        not isinstance(source, str)
        or not source.strip()
        or not isinstance(target, str)
        or not target.strip()
        for source, target in line_duality.items()
    ):
        raise ValueError("line-duality labels must be nonempty strings")
    if set(line_duality) != set(line_duality.values()):
        raise ValueError("line_duality must be a bijection on its declared labels")
    if any(line_duality[line_duality[line]] != line for line in line_duality):
        raise ValueError("line_duality must be involutive")
    if not required_lines <= set(line_duality):
        raise ValueError("line_duality does not cover every plateau line")


def row_record(row: BirthRow) -> dict[str, Any]:
    return {
        "tau1": fraction_text(row.tau1),
        "tau2": fraction_text(row.tau2),
        "kind": row.kind,
        "line": row.line,
    }


def dualize_rows(
    records: Sequence[Mapping[str, Any]], line_duality: Mapping[str, str]
) -> list[dict[str, Any]]:
    rows = tuple(BirthRow.from_record(record) for record in records)
    required_lines = {row.line for row in rows if row.kind == "plateau"}
    validate_line_duality(line_duality, required_lines)
    dual = []
    for row in rows:
        line = None if row.kind == "direct_rank2" else line_duality[row.line]
        dual.append(
            row_record(
                BirthRow(
                    tau1=Fraction(1) - row.tau2,
                    tau2=Fraction(1) - row.tau1,
                    kind=row.kind,
                    line=line,
                )
            )
        )
    return dual


def canonical_rows(records: Sequence[Mapping[str, Any]]) -> list[tuple[str, str, str, Optional[str]]]:
    return sorted(
        (fraction_text(row.tau1), fraction_text(row.tau2), row.kind, row.line)
        for row in (BirthRow.from_record(record) for record in records)
    )


def certify_complement(
    records: Sequence[Mapping[str, Any]], threshold_source: Any, line_duality: Mapping[str, str]
) -> dict[str, Any]:
    threshold = exact_fraction(threshold_source, field="threshold")
    if not Fraction(0) <= threshold <= Fraction(1):
        raise ValueError("threshold must lie in [0,1]")
    rows = tuple(BirthRow.from_record(record) for record in records)
    if any(threshold in (row.tau1, row.tau2) for row in rows):
        raise ValueError("complement control requires a tie-free threshold")
    dual_records = dualize_rows(records, line_duality)
    original = reconstruct_at_threshold(records, threshold)
    dual_threshold = Fraction(1) - threshold
    dual = reconstruct_at_threshold(dual_records, dual_threshold)
    original_rank = original["rank_counts"]
    dual_rank = dual["rank_counts"]
    if dual_rank != {
        "P0": original_rank["P2"],
        "P1": original_rank["P1"],
        "P2": original_rank["P0"],
    }:
        raise ArithmeticError("complement rank sectors did not reverse")
    original_m = exact_fraction(original["probabilities"]["M"], field="M")
    dual_m = exact_fraction(dual["probabilities"]["M"], field="dual M")
    if dual_m != -original_m:
        raise ArithmeticError("complement did not reverse Matching-One sign")
    twice_dual = dualize_rows(dual_records, line_duality)
    if canonical_rows(twice_dual) != canonical_rows(records):
        raise ArithmeticError("row duality is not involutive")
    return {
        "threshold": fraction_text(threshold),
        "dual_threshold": fraction_text(dual_threshold),
        "original_rank_counts": original_rank,
        "dual_rank_counts": dual_rank,
        "original_M": fraction_text(original_m),
        "dual_M": fraction_text(dual_m),
        "exact_checks": {
            "rank_sectors_reverse": True,
            "M_changes_sign": True,
            "row_duality_involutive": True,
        },
    }


def analyze_orientation_pair(record: Mapping[str, Any]) -> dict[str, Any]:
    if set(record) != {"threshold", "line_duality", "orientations"}:
        raise ValueError("control fields must be exactly threshold,line_duality,orientations")
    threshold = exact_fraction(record["threshold"], field="threshold")
    orientations = record["orientations"]
    if not isinstance(orientations, Sequence) or isinstance(orientations, (str, bytes)) or len(orientations) != 2:
        raise ValueError("exactly two ordered orientations are required")
    parsed = []
    names = set()
    for orientation in orientations:
        if not isinstance(orientation, Mapping) or set(orientation) != {"name", "covector", "rows"}:
            raise ValueError("orientation fields must be exactly name,covector,rows")
        name = orientation["name"]
        if not isinstance(name, str) or not name.strip() or name in names:
            raise ValueError("orientation names must be unique nonempty strings")
        names.add(name)
        covector = exact_fraction(orientation["covector"], field=f"{name}.covector")
        rows = orientation["rows"]
        if not isinstance(rows, Sequence) or isinstance(rows, (str, bytes)):
            raise ValueError("orientation rows must be a sequence")
        certificate = certify_complement(rows, threshold, record["line_duality"])
        parsed.append((name, covector, rows, certificate))
    if parsed[0][1] == parsed[1][1]:
        raise ValueError("orientation covectors must be distinct")
    if len(parsed[0][2]) != len(parsed[1][2]):
        raise ValueError("orientation sample counts must agree")
    gap = parsed[1][1] - parsed[0][1]
    original_m = [exact_fraction(item[3]["original_M"], field="original M") for item in parsed]
    dual_m = [exact_fraction(item[3]["dual_M"], field="dual M") for item in parsed]
    original_contrast = (original_m[1] - original_m[0]) / gap
    dual_contrast = (dual_m[1] - dual_m[0]) / gap
    if dual_contrast != -original_contrast:
        raise ArithmeticError("ordered orientation contrast did not reverse sign")

    self_matching_rows = [
        list(item[2]) + dualize_rows(item[2], record["line_duality"]) for item in parsed
    ]
    midpoint_m = []
    for rows in self_matching_rows:
        midpoint = reconstruct_at_threshold(rows, "1/2")
        midpoint_m.append(exact_fraction(midpoint["probabilities"]["M"], field="midpoint M"))
    midpoint_contrast = (midpoint_m[1] - midpoint_m[0]) / gap
    if any(midpoint_m) or midpoint_contrast:
        raise ArithmeticError("complement-closed self-matching control has a nonzero odd coordinate")

    return {
        "threshold": fraction_text(threshold),
        "dual_threshold": fraction_text(Fraction(1) - threshold),
        "orientation_order": [item[0] for item in parsed],
        "covectors": {item[0]: fraction_text(item[1]) for item in parsed},
        "covector_gap": fraction_text(gap),
        "orientation_certificates": {item[0]: item[3] for item in parsed},
        "H4_contrasts": {
            "original": fraction_text(original_contrast),
            "dual": fraction_text(dual_contrast),
        },
        "self_matching_midpoint": {
            parsed[0][0]: fraction_text(midpoint_m[0]),
            parsed[1][0]: fraction_text(midpoint_m[1]),
            "odd_H4_contrast": fraction_text(midpoint_contrast),
        },
        "exact_checks": {
            "orientation_descriptor_preserved": True,
            "dual_H4_equals_negative_original": True,
            "self_matching_odd_coordinate_zero": True,
        },
    }


def synthetic_control() -> dict[str, Any]:
    rows = [
        {"tau1": "1/4", "tau2": "1/4", "kind": "direct_rank2", "line": None},
        {"tau1": "3/4", "tau2": "3/4", "kind": "direct_rank2", "line": None},
        {"tau1": "1/4", "tau2": "3/4", "kind": "plateau", "line": "L0"},
        {"tau1": "0", "tau2": "1/4", "kind": "plateau", "line": "L2"},
        {"tau1": "0", "tau2": "3/4", "kind": "plateau", "line": "L0"},
        {"tau1": "1/4", "tau2": "1", "kind": "plateau", "line": "L1"},
    ]
    line_duality = {"L0": "L1", "L1": "L0", "L2": "L2"}
    return {
        "threshold": "3/8",
        "line_duality": line_duality,
        "orientations": [
            {"name": "first", "covector": "-1", "rows": rows},
            {"name": "second", "covector": "1", "rows": dualize_rows(rows, line_duality)},
        ],
    }


def build_artifact() -> dict[str, Any]:
    report = analyze_orientation_pair(synthetic_control())
    if report["H4_contrasts"] != {"original": "-1/6", "dual": "1/6"}:
        raise AssertionError("synthetic complement contrast changed")
    return {
        "schema": SCHEMA,
        "issue": 439,
        "status": "exact_complement_orientation_sign_control",
        "transform": "(tau1,tau2)->(1-tau2,1-tau1)",
        "tie_policy": "fail closed when p equals a birth time",
        "report": report,
        "claim_boundary": {
            "included": "exact rank reversal, M sign, orientation sign, and self-matching zero",
            "excluded": "production archives, empirical amplitudes, transfer fits, or physics claims",
            "parent_issue": "remain open",
        },
    }


def validate_artifact(artifact: Mapping[str, Any]) -> dict[str, Any]:
    expected = build_artifact()
    if artifact != expected:
        raise ValueError("complement sign certificate does not exactly reproduce")
    return {
        "schema": SCHEMA,
        "status": "valid_exact_complement_orientation_sign_control",
        "orientation_count": len(expected["report"]["orientation_order"]),
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--validate", type=Path)
    args = parser.parse_args(argv)
    if args.validate:
        artifact = json.loads(args.validate.read_text(encoding="utf-8"))
        print(json.dumps(validate_artifact(artifact), indent=2, sort_keys=True))
        return 0
    rendered = json.dumps(build_artifact(), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

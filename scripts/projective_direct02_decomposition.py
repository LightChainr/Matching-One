#!/usr/bin/env python3
"""Exact direct-rank2 versus plateau decomposition for Issue 439."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

from projective_birth_crosswalk import BirthRow, exact_fraction, fraction_text


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "analysis" / "projective_direct02_decomposition_certificate.json"
SCHEMA = "matching-one/projective-direct02-decomposition/v1"


def sector_statistics(rows: Sequence[BirthRow], threshold: Fraction) -> dict[str, Any]:
    if not rows:
        raise ValueError("orientation requires at least one birth row")
    direct = tuple(row for row in rows if row.kind == "direct_rank2")
    plateau = tuple(row for row in rows if row.kind == "plateau")
    if not direct or not plateau:
        raise ValueError("both DIRECT_RANK2 and plateau sectors must be present")

    def rank_counts(selected: Sequence[BirthRow]) -> tuple[int, int, int]:
        return tuple(sum(row.rank_at(threshold) == rank for row in selected) for rank in range(3))

    total_counts = rank_counts(rows)
    direct_counts = rank_counts(direct)
    plateau_counts = rank_counts(plateau)
    if any(
        total_counts[rank] != direct_counts[rank] + plateau_counts[rank]
        for rank in range(3)
    ):
        raise ArithmeticError("rank-sector decomposition does not close")

    sample_count = len(rows)
    direct_numerator = direct_counts[2] - direct_counts[0]
    plateau_numerator = plateau_counts[2] - plateau_counts[0]
    total_numerator = total_counts[2] - total_counts[0]
    if total_numerator != direct_numerator + plateau_numerator:
        raise ArithmeticError("Matching-One numerator decomposition does not close")

    direct_births = sum(row.tau2 <= threshold for row in direct)
    return {
        "sample_count": sample_count,
        "sector_sizes": {"DIRECT_RANK2": len(direct), "plateau": len(plateau)},
        "rank_counts": {
            "total": {"P0": total_counts[0], "P1": total_counts[1], "P2": total_counts[2]},
            "DIRECT_RANK2": {
                "P0": direct_counts[0],
                "P1": direct_counts[1],
                "P2": direct_counts[2],
            },
            "plateau": {
                "P0": plateau_counts[0],
                "P1": plateau_counts[1],
                "P2": plateau_counts[2],
            },
        },
        "P_direct02": fraction_text(Fraction(direct_births, sample_count)),
        "M_with_direct02": fraction_text(Fraction(total_numerator, sample_count)),
        "M_without_direct02": fraction_text(Fraction(plateau_numerator, len(plateau))),
        "additive_M_contributions": {
            "DIRECT_RANK2": fraction_text(Fraction(direct_numerator, sample_count)),
            "plateau": fraction_text(Fraction(plateau_numerator, sample_count)),
        },
        "exact_checks": {
            "rank_counts_close": True,
            "total_M_equals_direct_plus_plateau": True,
            "conditioned_plateau_uses_its_own_denominator": True,
        },
    }


def analyze_batch(record: Mapping[str, Any]) -> dict[str, Any]:
    if set(record) != {"batch", "threshold", "orientations"}:
        raise ValueError("batch fields must be exactly batch,threshold,orientations")
    batch = record["batch"]
    if not isinstance(batch, str) or not batch.strip():
        raise ValueError("batch must be a nonempty string")
    threshold = exact_fraction(record["threshold"], field="threshold")
    if not Fraction(0) <= threshold <= Fraction(1):
        raise ValueError("threshold must lie in [0,1]")
    orientations = record["orientations"]
    if not isinstance(orientations, Sequence) or isinstance(orientations, (str, bytes)):
        raise ValueError("orientations must be a two-entry sequence")
    if len(orientations) != 2:
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
        rows_source = orientation["rows"]
        if not isinstance(rows_source, Sequence) or isinstance(rows_source, (str, bytes)):
            raise ValueError("orientation rows must be a sequence")
        rows = tuple(BirthRow.from_record(row) for row in rows_source)
        parsed.append((name, covector, rows, sector_statistics(rows, threshold)))

    if len(parsed[0][2]) != len(parsed[1][2]):
        raise ValueError("ordered orientations must have equal sample counts")
    covector_gap = parsed[1][1] - parsed[0][1]
    if covector_gap == 0:
        raise ValueError("orientation covectors must be distinct")

    def statistic(index: int, key: str, sector: Optional[str] = None) -> Fraction:
        source = parsed[index][3][key]
        if sector is not None:
            source = source[sector]
        return exact_fraction(source, field=key)

    total_contrast = (
        statistic(1, "M_with_direct02") - statistic(0, "M_with_direct02")
    ) / covector_gap
    direct_contrast = (
        statistic(1, "additive_M_contributions", "DIRECT_RANK2")
        - statistic(0, "additive_M_contributions", "DIRECT_RANK2")
    ) / covector_gap
    plateau_contrast = (
        statistic(1, "additive_M_contributions", "plateau")
        - statistic(0, "additive_M_contributions", "plateau")
    ) / covector_gap
    if total_contrast != direct_contrast + plateau_contrast:
        raise ArithmeticError("orientation contrast decomposition does not close")

    return {
        "batch": batch,
        "threshold": fraction_text(threshold),
        "orientation_order": [parsed[0][0], parsed[1][0]],
        "covectors": {
            parsed[0][0]: fraction_text(parsed[0][1]),
            parsed[1][0]: fraction_text(parsed[1][1]),
        },
        "covector_gap": fraction_text(covector_gap),
        "orientations": {parsed[index][0]: parsed[index][3] for index in range(2)},
        "H4_contrasts": {
            "A_M_total": fraction_text(total_contrast),
            "A_M_DIRECT_RANK2_contribution": fraction_text(direct_contrast),
            "A_M_plateau_contribution": fraction_text(plateau_contrast),
        },
        "exact_checks": {
            "equal_orientation_sample_counts": True,
            "A_M_total_equals_direct_plus_plateau": True,
            "M_without_direct02_is_descriptive_conditioning_only": True,
        },
    }


def synthetic_batch() -> dict[str, Any]:
    first_rows = [
        {"tau1": "1/4", "tau2": "1/4", "kind": "direct_rank2", "line": None},
        {"tau1": "1/2", "tau2": "1/2", "kind": "direct_rank2", "line": None},
        {"tau1": "0", "tau2": "1/2", "kind": "plateau", "line": "L0"},
        {"tau1": "1/4", "tau2": "3/4", "kind": "plateau", "line": "L0"},
        {"tau1": "1/2", "tau2": "1", "kind": "plateau", "line": "L1"},
        {"tau1": "3/4", "tau2": "1", "kind": "plateau", "line": "L2"},
    ]
    second_rows = [
        {"tau1": "1/2", "tau2": "1/2", "kind": "direct_rank2", "line": None},
        {"tau1": "3/4", "tau2": "3/4", "kind": "direct_rank2", "line": None},
        {"tau1": "0", "tau2": "1/2", "kind": "plateau", "line": "L0"},
        {"tau1": "1/4", "tau2": "1/2", "kind": "plateau", "line": "L0"},
        {"tau1": "1/2", "tau2": "3/4", "kind": "plateau", "line": "L1"},
        {"tau1": "3/4", "tau2": "1", "kind": "plateau", "line": "L2"},
    ]
    return {
        "batch": "synthetic-b0",
        "threshold": "1/2",
        "orientations": [
            {"name": "first", "covector": "-1", "rows": first_rows},
            {"name": "second", "covector": "1", "rows": second_rows},
        ],
    }


def build_artifact() -> dict[str, Any]:
    report = analyze_batch(synthetic_batch())
    if report["H4_contrasts"] != {
        "A_M_total": "-1/12",
        "A_M_DIRECT_RANK2_contribution": "-1/6",
        "A_M_plateau_contribution": "1/12",
    }:
        raise AssertionError("synthetic orientation contrast changed")
    return {
        "schema": SCHEMA,
        "issue": 439,
        "status": "exact_direct02_mechanism_decomposition_control",
        "report": report,
        "claim_boundary": {
            "included": "exact additive direct-rank2/plateau decomposition and orientation contrast",
            "excluded": "raw archives, empirical rates, covariance, continuum attribution, or physics claims",
            "conditioned_quantity": "M_without_direct02 is descriptive only, not a new Matching-One observable",
            "parent_issue": "remain open",
        },
    }


def validate_artifact(artifact: Mapping[str, Any]) -> dict[str, Any]:
    expected = build_artifact()
    if artifact != expected:
        raise ValueError("direct02 decomposition certificate does not exactly reproduce")
    return {
        "schema": SCHEMA,
        "status": "valid_exact_direct02_mechanism_decomposition_control",
        "orientation_count": len(expected["report"]["orientations"]),
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

#!/usr/bin/env python3
"""Exact same-stream crosswalk from projective birth rows to rank observables."""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
import json
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, Sequence, Union


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "analysis" / "projective_birth_crosswalk_certificate.json"
SCHEMA = "matching-one/projective-birth-crosswalk/v1"
ExactInput = Union[int, str, Fraction]


def exact_fraction(value: ExactInput, *, field: str) -> Fraction:
    if isinstance(value, bool) or isinstance(value, float):
        raise TypeError(f"{field} must be an exact integer, fraction string, or Fraction")
    try:
        result = value if isinstance(value, Fraction) else Fraction(value)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError(f"invalid exact value for {field}") from exc
    return result


def fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}" if value.denominator != 1 else str(value.numerator)


@dataclass(frozen=True)
class BirthRow:
    """One permutation row with exact first/second ambient-rank birth times."""

    tau1: Fraction
    tau2: Fraction
    kind: str
    line: Optional[str]

    def __post_init__(self) -> None:
        if not Fraction(0) <= self.tau1 <= self.tau2 <= Fraction(1):
            raise ValueError("require 0 <= tau1 <= tau2 <= 1")
        if self.kind == "direct_rank2":
            if self.tau1 != self.tau2:
                raise ValueError("DIRECT_RANK2 requires tau1 == tau2")
            if self.line is not None:
                raise ValueError("DIRECT_RANK2 cannot carry a rank-one plateau line")
        elif self.kind == "plateau":
            if not self.tau1 < self.tau2:
                raise ValueError("plateau rows require tau1 < tau2")
            if not isinstance(self.line, str) or not self.line.strip():
                raise ValueError("plateau rows require a nonempty line label")
        else:
            raise ValueError("kind must be 'direct_rank2' or 'plateau'")

    @classmethod
    def from_record(cls, record: Mapping[str, Any]) -> "BirthRow":
        required = {"tau1", "tau2", "kind", "line"}
        if set(record) != required:
            raise ValueError("birth record fields must be exactly tau1,tau2,kind,line")
        return cls(
            tau1=exact_fraction(record["tau1"], field="tau1"),
            tau2=exact_fraction(record["tau2"], field="tau2"),
            kind=record["kind"],
            line=record["line"],
        )

    def rank_at(self, threshold: Fraction) -> int:
        if threshold < self.tau1:
            return 0
        if threshold < self.tau2:
            return 1
        return 2


def reconstruct_at_threshold(
    records: Iterable[Mapping[str, Any]], threshold_source: ExactInput
) -> dict[str, Any]:
    rows = tuple(BirthRow.from_record(record) for record in records)
    if not rows:
        raise ValueError("crosswalk requires at least one birth row")
    threshold = exact_fraction(threshold_source, field="threshold")
    if not Fraction(0) <= threshold <= Fraction(1):
        raise ValueError("threshold must lie in [0,1]")

    ranks = Counter(row.rank_at(threshold) for row in rows)
    line_counts = Counter(
        row.line for row in rows if row.rank_at(threshold) == 1
    )
    kind_counts = Counter(row.kind for row in rows)
    sample_count = len(rows)
    p0, p1, p2 = (ranks[index] for index in range(3))
    f1 = sum(row.tau1 <= threshold for row in rows)
    f2 = sum(row.tau2 <= threshold for row in rows)

    if p0 + p1 + p2 != sample_count:
        raise ArithmeticError("rank counts do not close")
    if f1 != p1 + p2 or f2 != p2:
        raise ArithmeticError("birth CDF and rank counts disagree")
    if sum(line_counts.values()) != p1:
        raise ArithmeticError("rank-one plateau line counts do not close")
    if kind_counts["direct_rank2"] + kind_counts["plateau"] != sample_count:
        raise ArithmeticError("input kind counts do not close")

    matching_numerator = p2 - p0
    cdf_numerator = f1 + f2 - sample_count
    if matching_numerator != cdf_numerator:
        raise ArithmeticError("M=P2-P0 and M=F1+F2-1 disagree")

    return {
        "threshold": fraction_text(threshold),
        "sample_count": sample_count,
        "rank_counts": {"P0": p0, "P1": p1, "P2": p2},
        "birth_cdf_counts": {"F1": f1, "F2": f2},
        "plateau_line_counts": {
            str(line): count for line, count in sorted(line_counts.items())
        },
        "input_kind_counts": {
            "DIRECT_RANK2": kind_counts["direct_rank2"],
            "plateau": kind_counts["plateau"],
        },
        "probabilities": {
            "P0": fraction_text(Fraction(p0, sample_count)),
            "P1": fraction_text(Fraction(p1, sample_count)),
            "P2": fraction_text(Fraction(p2, sample_count)),
            "F1": fraction_text(Fraction(f1, sample_count)),
            "F2": fraction_text(Fraction(f2, sample_count)),
            "M": fraction_text(Fraction(matching_numerator, sample_count)),
        },
        "exact_identities": {
            "rank_partition": True,
            "P1_equals_line_sum": True,
            "M_equals_P2_minus_P0": True,
            "M_equals_F1_plus_F2_minus_1": True,
        },
    }


def synthetic_rows() -> list[dict[str, Any]]:
    return [
        {"tau1": "1/4", "tau2": "3/4", "kind": "plateau", "line": "L0"},
        {"tau1": "1/2", "tau2": "1", "kind": "plateau", "line": "L1"},
        {"tau1": "1/2", "tau2": "1/2", "kind": "direct_rank2", "line": None},
        {"tau1": "3/4", "tau2": "3/4", "kind": "direct_rank2", "line": None},
        {"tau1": "0", "tau2": "1/2", "kind": "plateau", "line": "L0"},
        {"tau1": "3/4", "tau2": "1", "kind": "plateau", "line": "L2"},
    ]


def build_artifact() -> dict[str, Any]:
    rows = synthetic_rows()
    thresholds = ("0", "1/4", "1/2", "3/4", "1")
    return {
        "schema": SCHEMA,
        "issue": 439,
        "status": "exact_same_stream_birth_crosswalk_control",
        "boundary_convention": "birth at tau is included when tau <= threshold",
        "synthetic_rows": rows,
        "threshold_reconstructions": [
            reconstruct_at_threshold(rows, threshold) for threshold in thresholds
        ],
        "claim_boundary": {
            "included": "typed exact row-level P0/P1/P2/F1/F2/M reconstruction",
            "excluded": "raw archive import, covariance, wedge/common-ray/transfer scoring, production or physical claims",
            "parent_issue": "remain open",
        },
    }


def validate_artifact(artifact: Mapping[str, Any]) -> dict[str, Any]:
    regenerated = build_artifact()
    if artifact != regenerated:
        raise ValueError("projective birth crosswalk artifact does not exactly reproduce")
    return {
        "schema": SCHEMA,
        "status": "valid_exact_same_stream_birth_crosswalk_control",
        "row_count": len(regenerated["synthetic_rows"]),
        "threshold_count": len(regenerated["threshold_reconstructions"]),
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

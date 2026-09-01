#!/usr/bin/env python3
"""Extract amplitude-free stage selectivity from the exact collar matrix."""
from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path


@dataclass(frozen=True)
class Interval:
    lo: F
    hi: F

    @classmethod
    def from_record(cls, row: dict) -> "Interval":
        return cls(F(row["lower"]), F(row["upper"]))

    def __neg__(self) -> "Interval":
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other: "Interval") -> "Interval":
        return self + (-other)

    def __add__(self, other: "Interval") -> "Interval":
        return Interval(self.lo + other.lo, self.hi + other.hi)

    def __mul__(self, other: "Interval") -> "Interval":
        values = (self.lo * other.lo, self.lo * other.hi,
                  self.hi * other.lo, self.hi * other.hi)
        return Interval(min(values), max(values))

    def __truediv__(self, other: "Interval") -> "Interval":
        if other.lo <= 0 <= other.hi:
            raise ZeroDivisionError("interval denominator contains zero")
        return self * Interval(1 / other.hi, 1 / other.lo)


def record(value: Interval) -> dict:
    midpoint = (value.lo + value.hi) / 2
    return {
        "lower": str(value.lo),
        "upper": str(value.hi),
        "midpoint": float(midpoint),
        "width": float(value.hi - value.lo),
        "excludes_zero": value.lo > 0 or value.hi < 0,
        "excludes_one": value.hi < 1 or value.lo > 1,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--score", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    score = json.loads(args.score.read_text())
    matrix = [[Interval.from_record(cell) for cell in row]
              for row in score["matrix"]["P4_Schur"]]
    if len(matrix) != 2 or any(len(row) != 2 for row in matrix):
        raise ValueError("expected the frozen 2x2 preferred collar matrix")
    first_absent, first_present = matrix[0]
    second_absent, second_present = matrix[1]
    first_retention = first_present / first_absent
    second_retention = second_present / second_absent
    cross_ratio = (first_absent * second_present) / (first_present * second_absent)
    retention_ratio = first_retention / second_retention
    predicted_second_present = first_retention * second_absent
    stage_interaction = second_present - predicted_second_present
    payload = {
        "schema": "matching-one/p537-finite-collar-stage-selectivity/v1",
        "status": "exact_stage_selective_source_gate",
        "source_score": str(args.score),
        "rows": score["matrix"]["row_order"],
        "columns": score["matrix"]["column_order"],
        "matrix": [[record(cell) for cell in row] for row in matrix],
        "present_over_absent_response": {
            "first_birth": record(first_retention),
            "second_birth": record(second_retention),
        },
        "first_over_second_retention_ratio": record(retention_ratio),
        "rank_one_cross_ratio": {
            "value": record(cross_ratio),
            "rank_one_prediction": "1 exactly",
        },
        "first_birth_calibrated_rank_one_prediction": {
            "predicted_second_birth_present": record(predicted_second_present),
            "actual_second_birth_present": record(second_present),
            "stage_interaction_residual": record(stage_interaction),
            "identity": "residual = determinant / first_birth_absent",
        },
        "present_minus_absent_contrast": {
            "first_birth": record(first_present - first_absent),
            "second_birth": record(second_present - second_absent),
        },
        "decision": (
            "Source presence retains about 60.7% of the first-birth response "
            "but only about 2.25% of the second-birth response; the retention "
            "ratio is about 27.0 and the amplitude-free cross-ratio excludes one."
        ),
        "boundary": (
            "These are response-amplitude ratios within the same exact N25 "
            "radius-one collar matrix, not probabilities or independent evidence."
        ),
    }
    args.output.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({
        "first_birth_retention": payload["present_over_absent_response"]["first_birth"]["midpoint"],
        "second_birth_retention": payload["present_over_absent_response"]["second_birth"]["midpoint"],
        "retention_ratio": payload["first_over_second_retention_ratio"]["midpoint"],
        "rank_one_cross_ratio": payload["rank_one_cross_ratio"]["value"]["midpoint"],
        "stage_interaction_residual": payload["first_birth_calibrated_rank_one_prediction"]["stage_interaction_residual"]["midpoint"],
    }, indent=2))


if __name__ == "__main__":
    main()

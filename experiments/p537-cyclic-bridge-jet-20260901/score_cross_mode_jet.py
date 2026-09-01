#!/usr/bin/env python3
"""Certify the cross-mode thermal jet in the P537 provisional bridge basis."""
from __future__ import annotations

import json
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path


HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "p537-landing-matrix-preflight-20260901" / "result.json"


@dataclass(frozen=True)
class Interval:
    lo: F
    hi: F

    @classmethod
    def from_record(cls, record: dict) -> "Interval":
        return cls(F(record["lower"]), F(record["upper"]))

    def __add__(self, other: "Interval") -> "Interval":
        return Interval(self.lo + other.lo, self.hi + other.hi)

    def __neg__(self) -> "Interval":
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other: "Interval") -> "Interval":
        return self + (-other)

    def __mul__(self, other: "Interval") -> "Interval":
        products = (
            self.lo * other.lo,
            self.lo * other.hi,
            self.hi * other.lo,
            self.hi * other.hi,
        )
        return Interval(min(products), max(products))

    def __truediv__(self, other: "Interval") -> "Interval":
        if other.lo <= 0 <= other.hi:
            raise ZeroDivisionError("interval denominator contains zero")
        reciprocal = Interval(1 / other.hi, 1 / other.lo)
        return self * reciprocal


def record(value: Interval) -> dict:
    return {
        "lower": str(value.lo),
        "upper": str(value.hi),
        "midpoint": float((value.lo + value.hi) / 2),
        "width": float(value.hi - value.lo),
        "excludes_zero": value.lo > 0 or value.hi < 0,
    }


def main() -> None:
    source = json.loads(SOURCE.read_text())
    modes = {row["mode"]: row for row in source["modes"]}

    def jet(mode: str) -> tuple[Interval, Interval]:
        row = modes[mode]
        fixed_m = Interval.from_record(
            row["P4_root_projection"]["source_after_fixed_M_Schur_elimination"]
        )
        thermal = Interval.from_record(row["root_conditioned_mixed_hessian_Tp_over_Mp"])
        return fixed_m, thermal

    same_f, same_t = jet("clean_same")
    reversed_f, reversed_t = jet("clean_reversed")
    determinant = same_f * reversed_t - reversed_f * same_t
    result = {
        "schema": "matching-one/p537-cyclic-bridge-cross-mode-jet/v1",
        "status": "exact_nonzero_cross_mode_thermal_jet_under_provisional_contract",
        "source_artifact": str(SOURCE.relative_to(HERE.parents[1])),
        "basis": {
            "columns": ["clean_same", "clean_reversed"],
            "rows": [
                "F = source response of P4 Y after fixed-M Schur elimination",
                "H = T_p/M_p = normalized thermal derivative of F",
            ],
        },
        "matrix": [
            [record(same_f), record(reversed_f)],
            [record(same_t), record(reversed_t)],
        ],
        "determinant": record(determinant),
        "normalized_thermal_slopes_H_over_F": {
            "clean_same": record(same_t / same_f),
            "clean_reversed": record(reversed_t / reversed_f),
        },
        "decision": (
            "The two cyclic bridge-order modes are linearly independent in the "
            "(fixed-M response, root-conditioned thermal derivative) plane."
        ),
        "boundary": (
            "This reuses the exact N25 populations and the explicit provisional "
            "Bell-port contract of the source artifact; it is not an independent "
            "data block and not the canonical site-flip/no-extra-branch gate."
        ),
    }
    (HERE / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({
        "determinant_midpoint": result["determinant"]["midpoint"],
        "determinant_excludes_zero": result["determinant"]["excludes_zero"],
        "same_H_over_F": result["normalized_thermal_slopes_H_over_F"]["clean_same"]["midpoint"],
        "reversed_H_over_F": result["normalized_thermal_slopes_H_over_F"]["clean_reversed"]["midpoint"],
    }, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Extract and verify the single frozen radius-one collar minor."""
from __future__ import annotations

import argparse
import json
from fractions import Fraction as F
from pathlib import Path


ROWS = ("collar_r1_birth:[0,1]", "collar_r1_birth:[1,2]")
COLS = ("axial2:absent", "axial2:present")


def interval(record: dict[str, object]) -> tuple[F, F]:
    return F(str(record["lower"])), F(str(record["upper"]))


def multiply(x: tuple[F, F], y: tuple[F, F]) -> tuple[F, F]:
    values = (x[0]*y[0], x[0]*y[1], x[1]*y[0], x[1]*y[1])
    return min(values), max(values)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--score", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    score = json.loads(args.score.read_text())
    if score["matrix"]["row_order"] != list(ROWS):
        raise ValueError("unexpected preferred row order")
    if score["matrix"]["column_order"] != list(COLS):
        raise ValueError("unexpected preferred column order")
    matrix = score["matrix"]["P4_Schur"]
    a, b = interval(matrix[0][0]), interval(matrix[0][1])
    c, d = interval(matrix[1][0]), interval(matrix[1][1])
    ad, bc = multiply(a, d), multiply(b, c)
    determinant = ad[0]-bc[1], ad[1]-bc[0]
    reported = score["minors"][0]["determinant"]
    reported_interval = interval(reported)
    if not (reported_interval[0] <= determinant[0] <= determinant[1] <= reported_interval[1]):
        raise ValueError("reported determinant does not enclose the rounded-cell recomputation")
    if determinant[0] <= 0 <= determinant[1]:
        raise SystemExit("preferred collar minor contains zero")
    payload = {
        "schema": "matching-one/p537-finite-collar-preferred-minor/v1",
        "status": "exact_interval_nonzero",
        "collar": "B_inf(z,1)",
        "landing": "four labelled arms before outer reconnection; diagonal corner masks summed",
        "outer_attachment": {"J_B": 1, "J_W": 1},
        "local_source_contact_mask": 0,
        "rows": list(ROWS),
        "columns": list(COLS),
        "determinant": reported,
        "minor_search": score["minor_search"],
        "source_score": str(args.score),
        "boundary": "exact finite N25 radius-one collar rank-two certificate; no asymptotic arm claim",
    }
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Extract the frozen branch-free near-block witness from a full score."""
from __future__ import annotations

import argparse
import json
from fractions import Fraction as F
from pathlib import Path


ROWS = (
    "near_block:[0,1,5,0,0,0,0,0]",
    "near_block:[1,2,5,0,0,0,0,0]",
)
COLS = ("axial2:absent", "axial2:present")


def interval(record):
    return F(record["lower"]), F(record["upper"])


def multiply(x, y):
    values = (x[0]*y[0], x[0]*y[1], x[1]*y[0], x[1]*y[1])
    return min(values), max(values)


def subtract(x, y):
    return x[0]-y[1], x[1]-y[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--score", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    score = json.loads(args.score.read_text())
    row_index = {name: i for i, name in enumerate(score["matrix"]["row_order"])}
    col_index = {name: i for i, name in enumerate(score["matrix"]["column_order"])}
    matrix = score["matrix"]["P4_Schur"]
    a = interval(matrix[row_index[ROWS[0]]][col_index[COLS[0]]])
    b = interval(matrix[row_index[ROWS[0]]][col_index[COLS[1]]])
    c = interval(matrix[row_index[ROWS[1]]][col_index[COLS[0]]])
    d = interval(matrix[row_index[ROWS[1]]][col_index[COLS[1]]])
    determinant = subtract(multiply(a,d), multiply(b,c))
    payload = {
        "schema": "matching-one/p537-siteflip-preferred-minor/v1",
        "status": "exact_interval_nonzero" if determinant[1] < 0 or determinant[0] > 0 else "unresolved",
        "scope": "relaxed near-block; occupied component degree-branch=0, source-port extra-contact=0, globally occupied/vacant landing IDs both merged",
        "rows": list(ROWS),
        "columns": list(COLS),
        "determinant": {
            "lower": str(determinant[0]),
            "upper": str(determinant[1]),
            "midpoint": float(sum(determinant, F(0))/2),
            "excludes_zero": determinant[1] < 0 or determinant[0] > 0,
        },
        "boundary": "this proves finite rank two in the branch-free global-double-merged local-alternation block; the strict finite-collar ordinary four-arm row remains unconstructed",
        "source_score": str(args.score),
    }
    if not payload["determinant"]["excludes_zero"]:
        raise SystemExit("preferred minor interval contains zero")
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Fast point-estimate norm-5 thermal-jet diagnostic.

This script intentionally omits delete-one covariance. It uses the exact
finite-N Krawtchouk/Hermite coordinate map and the observed paired-rank-gap
width to answer one exploratory question quickly:

    Does measured w_can align the matching-odd thermal jet across N->5N?

Formal evidence remains the covariance-aware scorer. The output of this script
is a mechanism diagnostic only and must not be added as an independent primary
evidence block.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Mapping, Sequence

import mpmath as mp

from analyze_matching_parity_derivatives_fast import combine, read
from analyze_rank_gap_thermal_window import pooled_statistics, read_run
from hermite_krawtchouk_scaling_jet import (
    canonical_dimensionless_width,
    cocycle_residual,
    scaling_derivative_jet,
    width_normalized_jet,
)
from threshold_score_modes import project


SIZES = (65, 85, 130, 170, 325, 425)
LINKS = ((65, 325), (85, 425))
TRIPLES = ((65, 130, 325), (85, 170, 425))
ALPHA = mp.mpf(13) / 8
Q2_C = mp.mpf(8) / 5
JORDAN_C = mp.log(5) / mp.log(2)


@dataclass(frozen=True)
class RunSpec:
    n: int
    histogram: Path
    moments: Path
    metadata: Path


def parse_run(text: str) -> RunSpec:
    fields = text.split(":", 3)
    if len(fields) != 4:
        raise argparse.ArgumentTypeError("run must be N:HISTOGRAM:MOMENTS:METADATA")
    return RunSpec(int(fields[0]), Path(fields[1]), Path(fields[2]), Path(fields[3]))


def orientation_totals(histogram: Path, n: int):
    data = read(histogram)
    first = combine([data[key] for key in sorted(data) if key[0] == n and key[1] == "first"])
    second = combine([data[key] for key in sorted(data) if key[0] == n and key[1] == "second"])
    return first, second


def thermal_coefficients(point: Mapping[str, object], max_order: int) -> list[mp.mpf]:
    s = point["P4_S_modes"]
    d = point["P4_D_modes"]
    return [mp.mpf(d[r] if r % 2 == 0 else s[r]) for r in range(max_order + 1)]


def normalized_cross_residual(left: Sequence[mp.mpf], right: Sequence[mp.mpf]) -> list[mp.mpf]:
    if len(left) != len(right):
        raise ValueError("vectors must have equal length")
    if not left:
        return []
    out = []
    for r in range(1, len(left)):
        a = right[r] * left[0]
        b = left[r] * right[0]
        scale = abs(a) + abs(b)
        out.append(mp.mpf(0) if scale == 0 else (a - b) / scale)
    return out


def cosine(left: Sequence[mp.mpf], right: Sequence[mp.mpf]) -> mp.mpf:
    numerator = mp.fsum(a * b for a, b in zip(left, right))
    denominator = mp.sqrt(mp.fsum(a * a for a in left) * mp.fsum(b * b for b in right))
    return numerator / denominator if denominator else mp.nan


def nstr_list(values: Sequence[mp.mpf], digits: int = 16) -> list[str]:
    return [mp.nstr(value, digits) for value in values]


def analyze_run(run: RunSpec, max_order: int) -> dict[str, object]:
    first, second = orientation_totals(run.histogram, run.n)
    projected = project(first, second, max_order)
    coeff = thermal_coefficients(projected, max_order)
    jet = scaling_derivative_jet(coeff, run.n, mp.mpf(projected["p0"]), ALPHA)
    gap_run = read_run(run.n, run.moments, run.metadata)
    gap = pooled_statistics(gap_run)["gap_mean"]
    width = canonical_dimensionless_width(run.n, gap)
    width_jet = width_normalized_jet(jet, width)
    return {
        "N": run.n,
        "p0": mp.nstr(projected["p0"], 20),
        "gap_mean": mp.nstr(gap, 20),
        "w_can": mp.nstr(width, 20),
        "thermal_coefficients": nstr_list(coeff),
        "finite_jet": nstr_list(jet),
        "width_normalized_jet": nstr_list(width_jet),
        "thermal_signs": [0 if value == 0 else (1 if value > 0 else -1) for value in jet],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="append", type=parse_run, required=True)
    parser.add_argument("--max-order", type=int, default=6)
    parser.add_argument("--dps", type=int, default=40)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    mp.mp.dps = args.dps
    by_n = {run.n: run for run in args.run}
    if tuple(sorted(by_n)) != SIZES:
        raise SystemExit(f"exactly sizes {SIZES} are required")

    analyzed = {n: analyze_run(by_n[n], args.max_order) for n in SIZES}
    width_jets = {
        n: [mp.mpf(value) for value in analyzed[n]["width_normalized_jet"]]
        for n in SIZES
    }

    links = []
    for parent, child in LINKS:
        left = width_jets[parent]
        right = width_jets[child]
        residual = normalized_cross_residual(left, right)
        links.append(
            {
                "parent": parent,
                "child": child,
                "orders": list(range(1, args.max_order + 1)),
                "normalized_width_cross_residual": nstr_list(residual),
                "max_abs_residual_r2_plus": mp.nstr(max(abs(value) for value in residual[1:]), 16),
                "cosine_r2_plus": mp.nstr(cosine(left[2:], right[2:]), 16),
            }
        )

    triples = []
    for n, n2, n5 in TRIPLES:
        parent = width_jets[n][2:]
        norm2 = width_jets[n2][2:]
        norm5 = width_jets[n5][2:]
        q2 = cocycle_residual(parent, norm2, norm5, Q2_C)
        jordan = cocycle_residual(parent, norm2, norm5, JORDAN_C)
        scale = [abs(parent[i]) + abs(norm2[i]) + abs(norm5[i]) for i in range(len(parent))]
        q2_rel = [mp.mpf(0) if scale[i] == 0 else q2[i] / scale[i] for i in range(len(q2))]
        jordan_rel = [mp.mpf(0) if scale[i] == 0 else jordan[i] / scale[i] for i in range(len(jordan))]
        triples.append(
            {
                "sizes": [n, n2, n5],
                "orders": list(range(2, args.max_order + 1)),
                "q2_relative_cocycle": nstr_list(q2_rel),
                "jordan_relative_cocycle": nstr_list(jordan_rel),
                "q2_rms_relative": mp.nstr(mp.sqrt(mp.fsum(v * v for v in q2_rel) / len(q2_rel)), 16),
                "jordan_rms_relative": mp.nstr(mp.sqrt(mp.fsum(v * v for v in jordan_rel) / len(jordan_rel)), 16),
            }
        )

    payload = {
        "schema": "matching-one/norm5-fast-thermal-jet/v1",
        "status": "exploratory_point_estimate_only",
        "evidence_guard": "No covariance; do not add as an independent primary evidence block.",
        "by_N": {str(n): analyzed[n] for n in SIZES},
        "width_collapse_links": links,
        "width_corrected_cocycles": triples,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

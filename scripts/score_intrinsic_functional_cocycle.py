#!/usr/bin/env python3
"""Score the frozen intrinsic full-curve q=2/Jordan cocycle.

The input is the long threshold-rank histogram format emitted by
``threshold_rank_orientation_mc``.  Every intrinsic coordinate, projector,
and residual is recomputed inside each delete-one-batch replicate.

Runs that share random counters belong in one covariance group.  Independent
runs must be placed in separate groups; their jackknife covariance
contributions are then added as independent blocks.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Dict, Mapping, Sequence, Tuple

import mpmath as mp

from analyze_p48_retrospective import (
    Histogram,
    add_histograms,
    cos4,
    covariance_of_mean,
    pseudovalues,
    read_histograms,
    tail_and_derivative,
)


SIZES = (65, 85, 130, 170, 325, 425)
LINEAGES = ((65, 130, 325), (85, 170, 425))
LEVELS = (0.0, 0.025, 0.05)
MODELS = (
    ("q2_analytic", 8.0 / 5.0),
    ("jordan_rank2", math.log(5.0) / math.log(2.0)),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def merge_inputs(paths: Sequence[Path]) -> Dict[Tuple[int, str, int], Histogram]:
    merged: Dict[Tuple[int, str, int], Histogram] = {}
    for path in paths:
        rows = read_histograms(path)
        overlap = set(merged) & set(rows)
        if overlap:
            raise ValueError(f"duplicate histogram keys: {sorted(overlap)[:3]}")
        merged.update(rows)
    sizes = tuple(sorted({key[0] for key in merged}))
    if sizes != SIZES:
        raise ValueError(f"expected sizes {SIZES}, got {sizes}")
    return merged


def parse_groups(raw: Sequence[str]) -> Tuple[Tuple[int, ...], ...]:
    groups = tuple(tuple(int(value) for value in item.split(",")) for item in raw)
    if any(not group for group in groups):
        raise ValueError("empty covariance group")
    flat = [n for group in groups for n in group]
    if sorted(flat) != list(SIZES) or len(set(flat)) != len(flat):
        raise ValueError(f"covariance groups must partition {SIZES}")
    return groups


def grouped_rows(records: Mapping[Tuple[int, str, int], Histogram]):
    return {
        n: {
            orientation: sorted(
                (row for key, row in records.items() if key[:2] == (n, orientation)),
                key=lambda row: row.batch,
            )
            for orientation in ("first", "second")
        }
        for n in SIZES
    }


def validate_groups(rows, groups: Sequence[Sequence[int]]) -> None:
    for group in groups:
        signature = None
        for n in group:
            first = rows[n]["first"]
            second = rows[n]["second"]
            if len(first) != len(second) or len(first) < 2:
                raise ValueError(f"N={n}: invalid orientation batch pairing")
            current = (
                tuple(row.batch for row in first),
                tuple(row.samples for row in first),
                tuple(row.batch for row in second),
                tuple(row.samples for row in second),
            )
            if signature is None:
                signature = current
            elif current != signature:
                raise ValueError(f"covariance group {tuple(group)} is not batch aligned")


def aggregate(rows: Sequence[Histogram], omitted: int = -1):
    included = [row for row in rows if row.batch != omitted]
    return {
        "a": rows[0].a,
        "b": rows[0].b,
        "samples": sum(row.samples for row in included),
        "minus": add_histograms(rows, "minus", omitted),
        "plus": add_histograms(rows, "plus", omitted),
    }


def orientation_observables(n: int, row, p: float):
    minus, minus_prime = tail_and_derivative(row["minus"], row["samples"], p)
    plus, plus_prime = tail_and_derivative(row["plus"], row["samples"], p)
    r_primal = plus
    r_matching = 1.0 - minus
    return {
        "S": (r_primal + r_matching) / 2.0,
        "D": (r_primal - r_matching) / 2.0,
        "S_prime": (plus_prime - minus_prime) / 2.0,
        "D_prime": (plus_prime + minus_prime) / 2.0,
        "M": r_primal - r_matching,
        "M_prime": plus_prime + minus_prime,
    }


def solve_target(function, target: float) -> float:
    # The forward binomial recurrence starts at (1-p)^N.  At N=425 the
    # seemingly harmless p=0.9 endpoint underflows in binary64, while the
    # entire percolation transition and frozen u-grid are well inside this
    # narrower interval.
    lower, upper = 0.2, 0.8
    if function(lower) > target or function(upper) < target:
        raise ValueError(f"intrinsic target {target} is not bracketed")
    for _ in range(56):
        midpoint = (lower + upper) / 2.0
        if function(midpoint) < target:
            lower = midpoint
        else:
            upper = midpoint
    return (lower + upper) / 2.0


def size_statistics(by_orientation, omitted: int = -1):
    n = by_orientation["first"][0].n
    data = {
        orientation: aggregate(by_orientation[orientation], omitted)
        for orientation in ("first", "second")
    }

    def obs(orientation: str, p: float):
        return orientation_observables(n, data[orientation], p)

    def mean_matching(p: float) -> float:
        return (obs("first", p)["M"] + obs("second", p)["M"]) / 2.0

    delta_cos4 = cos4(data["first"]["a"], data["first"]["b"]) - cos4(
        data["second"]["a"], data["second"]["b"]
    )
    if delta_cos4 == 0.0:
        raise ValueError(f"N={n}: zero Delta cos(4 theta)")

    def projected(p: float):
        first = obs("first", p)
        second = obs("second", p)
        return {
            "P4_S": (first["S"] - second["S"]) / delta_cos4,
            "P4_D": (first["D"] - second["D"]) / delta_cos4,
            "P4_S_prime": (first["S_prime"] - second["S_prime"]) / delta_cos4,
            "P4_D_prime": (first["D_prime"] - second["D_prime"]) / delta_cos4,
            "mean_M_prime": (first["M_prime"] + second["M_prime"]) / 2.0,
        }

    levels = {}
    for level in LEVELS:
        if level == 0.0:
            p_minus = p_plus = solve_target(mean_matching, 0.0)
        else:
            p_minus = solve_target(mean_matching, -level)
            p_plus = solve_target(mean_matching, level)
        minus = projected(p_minus)
        plus = projected(p_plus)
        d_even = (plus["P4_D"] + minus["P4_D"]) / 2.0
        s_odd = (plus["P4_S"] - minus["P4_S"]) / 2.0
        levels[str(level)] = {
            "p_minus": p_minus,
            "p_plus": p_plus,
            "P4_D_even": d_even,
            "P4_S_odd": s_odd,
            "T": d_even + s_odd,
        }

    center = projected(levels["0.0"]["p_minus"])
    j_value = center["P4_S_prime"] / center["mean_M_prime"]
    xi_value = j_value / center["P4_D"]
    return {
        "N": n,
        "first_rep": [data["first"]["a"], data["first"]["b"]],
        "second_rep": [data["second"]["a"], data["second"]["b"]],
        "samples": data["first"]["samples"],
        "delta_cos4": delta_cos4,
        "center": {
            **center,
            "J": j_value,
            "J_scaled_N13_8": n ** (13.0 / 8.0) * j_value,
            "Xi": xi_value,
        },
        "levels": levels,
    }


def build_samples(rows, groups):
    full = {n: size_statistics(rows[n]) for n in SIZES}
    deleted = {}
    for group in groups:
        batch_ids = [row.batch for row in rows[group[0]]["first"]]
        group_samples = []
        for omitted in batch_ids:
            sample = dict(full)
            for n in group:
                sample[n] = size_statistics(rows[n], omitted)
            group_samples.append(sample)
        deleted[tuple(group)] = group_samples
    return full, deleted


def z_value(sample, n: int, level: float) -> float:
    return n ** (13.0 / 8.0) * sample[n]["levels"][str(level)]["T"]


def residual_function(parent: int, norm2: int, norm5: int, level: float, c: float):
    return lambda sample: (
        z_value(sample, norm5, level)
        - c * z_value(sample, norm2, level)
        + (c - 1.0) * z_value(sample, parent, level)
    )


def grouped_covariance(full, deleted, functions):
    point = [function(full) for function in functions]
    covariance = [[0.0 for _ in functions] for _ in functions]
    for samples in deleted.values():
        deleted_values = [[function(sample) for function in functions] for sample in samples]
        pseudo = [
            [
                pseudovalues(point[column], [row[column] for row in deleted_values])[batch]
                for column in range(len(functions))
            ]
            for batch in range(len(samples))
        ]
        contribution = covariance_of_mean(pseudo)
        covariance = [
            [covariance[i][j] + contribution[i][j] for j in range(len(functions))]
            for i in range(len(functions))
        ]
    return point, covariance


def spectral_quadratic(vector: Sequence[float], covariance: Sequence[Sequence[float]]):
    matrix = mp.matrix(covariance)
    eigenvalues, eigenvectors = mp.eigsy(matrix)
    values = [float(eigenvalues[i]) for i in range(len(vector))]
    scale = max(values)
    if scale <= 0.0:
        raise ValueError("covariance has no positive eigenvalue")
    tolerance = scale * 1e-10
    if min(values) < -tolerance:
        raise ValueError("covariance has a materially negative eigenvalue")
    kept = [i for i, value in enumerate(values) if value > tolerance]
    transformed = [
        math.fsum(float(eigenvectors[row, column]) * vector[row] for row in range(len(vector)))
        for column in kept
    ]
    chi_square = math.fsum(value * value / values[column] for value, column in zip(transformed, kept))
    return {
        "chi_square": chi_square,
        "df_effective": len(kept),
        "covariance_rank": len(kept),
        "covariance_condition": max(values[i] for i in kept) / min(values[i] for i in kept),
        "eigenvalue_cutoff": tolerance,
    }


def calculate(records, groups):
    rows = grouped_rows(records)
    validate_groups(rows, groups)
    full, deleted = build_samples(rows, groups)

    scores = {}
    for model, c_value in MODELS:
        functions = [
            residual_function(parent, norm2, norm5, level, c_value)
            for parent, norm2, norm5 in LINEAGES
            for level in LEVELS
        ]
        residual, covariance = grouped_covariance(full, deleted, functions)
        scores[model] = {
            "c": c_value,
            "order": [
                {"lineage": [parent, norm2, norm5], "u": level}
                for parent, norm2, norm5 in LINEAGES
                for level in LEVELS
            ],
            "residual": residual,
            "covariance": covariance,
            "marginal_z": [
                residual[i] / math.sqrt(covariance[i][i]) if covariance[i][i] > 0.0 else None
                for i in range(len(residual))
            ],
            **spectral_quadratic(residual, covariance),
        }

    center_functions = []
    center_order = []
    for n in SIZES:
        for field in ("J_scaled_N13_8", "Xi"):
            center_functions.append(lambda sample, n=n, field=field: sample[n]["center"][field])
            center_order.append({"N": n, "field": field})
    center_point, center_covariance = grouped_covariance(full, deleted, center_functions)

    return {
        "format_version": 1,
        "classification": "frozen pre-target functional-cocycle scorer",
        "levels": list(LEVELS),
        "covariance_groups": [list(group) for group in groups],
        "model_order": [name for name, _ in MODELS],
        "sizes": full,
        "center_metric": {
            "order": center_order,
            "point": center_point,
            "covariance": center_covariance,
        },
        "models_in_frozen_order": scores,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--histograms", nargs=6, required=True, type=Path)
    parser.add_argument(
        "--covariance-groups", nargs="+", required=True,
        help="comma-separated size groups, e.g. 65,85,130,170 325 425",
    )
    parser.add_argument("--json", required=True, type=Path)
    args = parser.parse_args()
    try:
        groups = parse_groups(args.covariance_groups)
        result = calculate(merge_inputs(args.histograms), groups)
    except (ArithmeticError, ValueError) as exc:
        raise SystemExit(str(exc))
    result["provenance"] = {str(path): sha256(path) for path in args.histograms}
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(args.json)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

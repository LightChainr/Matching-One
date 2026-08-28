#!/usr/bin/env python3
"""Score an H4-null scalar root projector from same-N orientation full curves.

For one same-N pair with c_i=cos(4 theta_i) and individual matching roots p_i,
define

    p_scalar(N) = (c1*p2 - c2*p1)/(c1-c2).

If

    p_i(N) = p_c + a0*N^-beta0 + a4*c_i*N^-2 + ...,

then p_scalar cancels the displayed H4 root term exactly and retains
p_c+a0*N^-beta0.  The V_<1,4> scalar mechanism predicts beta0=7/2
(equivalently L^-7).

For a Gaussian doubling lineage N->2N, q=2^-beta0 and

    p_c_hat(N) = [p_scalar(2N)-q*p_scalar(N)]/(1-q)

contains neither the scalar amplitude nor an externally supplied p_c.  Two or
more lineages therefore give a parameter-free consistency test.  All roots and
projectors are recomputed inside synchronized delete-one-batch replicates.

This scorer is deliberately agnostic about the matching-parity proof.  A pass
supports the finite-size scalar mechanism; it does not establish the OPE/RG
parity assumption.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Dict, Mapping, Sequence, Tuple

from analyze_p48_retrospective import (
    Histogram,
    add_histograms,
    cos4,
    covariance_of_mean,
    inverse,
    matmul,
    pseudovalues,
    quadratic,
    read_histograms,
    transpose,
)
from score_p49_fullcurve_doubling import tail_and_derivative


DEFAULT_LINEAGES = ((65, 130), (85, 170))


def merge_inputs(paths: Sequence[Path]) -> Dict[Tuple[int, str, int], Histogram]:
    merged: Dict[Tuple[int, str, int], Histogram] = {}
    for path in paths:
        records = read_histograms(path)
        overlap = set(merged) & set(records)
        if overlap:
            raise ValueError(f"duplicate histogram keys: {sorted(overlap)[:3]}")
        merged.update(records)
    sizes = sorted({key[0] for key in merged})
    if len(sizes) < 4:
        raise ValueError("need at least four sizes (two doubling lineages)")
    signature = None
    for n in sizes:
        selected = sorted(
            row for key, row in merged.items()
            if key[0] == n and key[1] == "first"
        )
        current = (
            tuple(row.batch for row in selected),
            tuple(row.samples for row in selected),
        )
        if signature is None:
            signature = current
        elif current != signature:
            raise ValueError("cross-size batch/sample alignment is absent")
    return merged


def grouped(records: Mapping[Tuple[int, str, int], Histogram]):
    sizes = sorted({key[0] for key in records})
    return {
        n: {
            orientation: sorted(
                (row for key, row in records.items() if key[:2] == (n, orientation)),
                key=lambda row: row.batch,
            )
            for orientation in ("first", "second")
        }
        for n in sizes
    }


def aggregate(rows: Sequence[Histogram], omitted: int = -1):
    included = [row for row in rows if row.batch != omitted]
    if not included:
        raise ValueError("cannot omit all batches")
    return {
        "a": rows[0].a,
        "b": rows[0].b,
        "samples": sum(row.samples for row in included),
        "minus": add_histograms(rows, "minus", omitted),
        "plus": add_histograms(rows, "plus", omitted),
    }


def matching_value(row, p: float) -> Tuple[float, float]:
    minus, d_minus = tail_and_derivative(row["minus"], row["samples"], p)
    plus, d_plus = tail_and_derivative(row["plus"], row["samples"], p)
    return minus + plus - 1.0, d_minus + d_plus


def solve_root(row) -> Tuple[float, float]:
    lower, upper = 0.55, 0.63
    f_lower = matching_value(row, lower)[0]
    f_upper = matching_value(row, upper)[0]
    if not f_lower <= 0.0 <= f_upper:
        raise ValueError("matching root is not bracketed")
    for _ in range(58):
        midpoint = (lower + upper) / 2.0
        if matching_value(row, midpoint)[0] < 0.0:
            lower = midpoint
        else:
            upper = midpoint
    root = (lower + upper) / 2.0
    slope = matching_value(row, root)[1]
    return root, slope


def size_stat(by_orientation, omitted: int = -1) -> dict:
    first = aggregate(by_orientation["first"], omitted)
    second = aggregate(by_orientation["second"], omitted)
    p1, b1 = solve_root(first)
    p2, b2 = solve_root(second)
    c1 = cos4(first["a"], first["b"])
    c2 = cos4(second["a"], second["b"])
    dc = c1 - c2
    if abs(dc) < 1e-12:
        raise ValueError("orientation pair has vanishing H4 leverage")
    p_scalar = (c1 * p2 - c2 * p1) / dc
    p_h4 = (p1 - p2) / dc
    return {
        "N": by_orientation["first"][0].n,
        "first": [first["a"], first["b"]],
        "second": [second["a"], second["b"]],
        "cos4_first": c1,
        "cos4_second": c2,
        "delta_cos4": dc,
        "root_first": p1,
        "root_second": p2,
        "slope_first": b1,
        "slope_second": b2,
        "p_scalar": p_scalar,
        "p_h4_coefficient": p_h4,
    }


def covariance_and_full(full: Sequence[float], deleted: Sequence[Sequence[float]]):
    pseudo_columns = [
        pseudovalues(full[j], [row[j] for row in deleted])
        for j in range(len(full))
    ]
    pseudo_rows = [
        [pseudo_columns[j][batch] for j in range(len(full))]
        for batch in range(len(deleted))
    ]
    return covariance_of_mean(pseudo_rows), pseudo_rows


def gls_fixed_beta(ns: Sequence[int], values: Sequence[float], covariance, beta: float):
    x = [[1.0, n ** (-beta)] for n in ns]
    inv_c = inverse(covariance)
    xt = transpose(x)
    normal = matmul(xt, matmul(inv_c, x))
    parameter_covariance = inverse(normal)
    rhs = matmul(xt, matmul(inv_c, [[v] for v in values]))
    parameters = matmul(parameter_covariance, rhs)
    fitted = [parameters[0][0] + parameters[1][0] * n ** (-beta) for n in ns]
    residual = [v - f for v, f in zip(values, fitted)]
    return {
        "beta_in_N": beta,
        "pc": parameters[0][0],
        "amplitude": parameters[1][0],
        "parameter_covariance": parameter_covariance,
        "fitted": fitted,
        "residual": residual,
        "chi_square": quadratic(residual, covariance),
        "df": len(ns) - 2,
    }


def lineage_stat(sample: Mapping[int, dict], parent: int, child: int, beta: float):
    if child != 2 * parent:
        raise ValueError("this scorer currently expects exact norm-2 lineages")
    q = 2.0 ** (-beta)
    p_parent = sample[parent]["p_scalar"]
    p_child = sample[child]["p_scalar"]
    pc_hat = (p_child - q * p_parent) / (1.0 - q)
    amplitude = (p_child - p_parent) * (parent ** beta) / (q - 1.0)
    return {
        "parent_N": parent,
        "child_N": child,
        "q": q,
        "p_scalar_parent": p_parent,
        "p_scalar_child": p_child,
        "pc_hat": pc_hat,
        "amplitude": amplitude,
    }


def jackknife_scalar(full_value: float, deleted_values: Sequence[float]) -> Tuple[float, float]:
    pseudo = pseudovalues(full_value, deleted_values)
    mean = math.fsum(pseudo) / len(pseudo)
    variance = math.fsum((x - mean) ** 2 for x in pseudo) / (len(pseudo) * (len(pseudo) - 1))
    return mean, math.sqrt(max(variance, 0.0))


def calculate(records, beta: float = 3.5, lineages=DEFAULT_LINEAGES) -> dict:
    groups = grouped(records)
    sizes = sorted(groups)
    required = {n for pair in lineages for n in pair}
    if not required.issubset(sizes):
        raise ValueError(f"missing lineage sizes: {sorted(required - set(sizes))}")
    full = {n: size_stat(groups[n]) for n in sizes}
    batches = len(groups[sizes[0]]["first"])
    deleted = [
        {n: size_stat(groups[n], omitted=batch) for n in sizes}
        for batch in range(batches)
    ]

    scalar_full = [full[n]["p_scalar"] for n in sizes]
    scalar_deleted = [[sample[n]["p_scalar"] for n in sizes] for sample in deleted]
    scalar_covariance, _ = covariance_and_full(scalar_full, scalar_deleted)
    fixed_beta_fit = gls_fixed_beta(sizes, scalar_full, scalar_covariance, beta)

    lineage_rows = {}
    pc_full = []
    pc_deleted_matrix = []
    amp_full = []
    amp_deleted_matrix = []
    for parent, child in lineages:
        key = f"{parent}->{child}"
        row = lineage_stat(full, parent, child, beta)
        deleted_rows = [lineage_stat(sample, parent, child, beta) for sample in deleted]
        pc_mean, pc_se = jackknife_scalar(row["pc_hat"], [x["pc_hat"] for x in deleted_rows])
        amp_mean, amp_se = jackknife_scalar(row["amplitude"], [x["amplitude"] for x in deleted_rows])
        row.update({
            "pc_hat_jackknife_mean": pc_mean,
            "pc_hat_se": pc_se,
            "amplitude_jackknife_mean": amp_mean,
            "amplitude_se": amp_se,
        })
        lineage_rows[key] = row
        pc_full.append(row["pc_hat"])
        amp_full.append(row["amplitude"])
    for sample in deleted:
        pc_deleted_matrix.append([
            lineage_stat(sample, parent, child, beta)["pc_hat"]
            for parent, child in lineages
        ])
        amp_deleted_matrix.append([
            lineage_stat(sample, parent, child, beta)["amplitude"]
            for parent, child in lineages
        ])

    pc_covariance, _ = covariance_and_full(pc_full, pc_deleted_matrix)
    amp_covariance, _ = covariance_and_full(amp_full, amp_deleted_matrix)
    if len(lineages) == 2:
        pc_diff = pc_full[0] - pc_full[1]
        pc_diff_var = pc_covariance[0][0] + pc_covariance[1][1] - 2.0 * pc_covariance[0][1]
        amp_diff = amp_full[0] - amp_full[1]
        amp_diff_var = amp_covariance[0][0] + amp_covariance[1][1] - 2.0 * amp_covariance[0][1]
        consistency = {
            "pc_difference": pc_diff,
            "pc_difference_se": math.sqrt(max(pc_diff_var, 0.0)),
            "pc_difference_z": pc_diff / math.sqrt(pc_diff_var) if pc_diff_var > 0 else None,
            "amplitude_difference": amp_diff,
            "amplitude_difference_se": math.sqrt(max(amp_diff_var, 0.0)),
            "amplitude_difference_z": amp_diff / math.sqrt(amp_diff_var) if amp_diff_var > 0 else None,
        }
    else:
        consistency = None

    return {
        "format_version": 1,
        "hypothesis": {
            "operator": "V_<1,4> scalar",
            "x": "33/4",
            "central_M_power_in_N": "25/8",
            "root_bias_power_in_N": "7/2",
            "root_bias_power_in_L": "7",
            "beta_in_N": beta,
            "doubling_q": 2.0 ** (-beta),
        },
        "batches": batches,
        "sizes": full,
        "p_scalar_covariance": scalar_covariance,
        "fixed_beta_gls": fixed_beta_fit,
        "lineages": lineage_rows,
        "lineage_pc_covariance": pc_covariance,
        "lineage_amplitude_covariance": amp_covariance,
        "two_lineage_consistency": consistency,
        "limitations": [
            "Two-angle scalar projection cancels H4 exactly but not arbitrary H8/H12/... contamination.",
            "Orientation-dependent slopes and nonlinear root conversion enter as subleading corrections.",
            "Matching-odd parity of V_<1,4> is a separate RG/OPE hypothesis.",
        ],
    }


def write_csv(path: Path, payload: dict) -> None:
    fields = [
        "N", "first", "second", "cos4_first", "cos4_second", "delta_cos4",
        "root_first", "root_second", "p_scalar", "p_h4_coefficient",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for n in sorted(payload["sizes"]):
            row = payload["sizes"][n]
            writer.writerow({field: row[field] for field in fields})


def report(payload: dict) -> str:
    fit = payload["fixed_beta_gls"]
    lines = [
        "# V_<1,4> scalar-root projector score",
        "",
        "Primary fixed hypothesis: `p_scalar(N)=p_c+a*N^(-7/2)`.",
        "The projector cancels the displayed H4 root term without using an external p_c.",
        "",
        "## Size projectors",
        "",
        "| N | p_scalar | H4 root coefficient |",
        "|---:|---:|---:|",
    ]
    for n in sorted(payload["sizes"]):
        row = payload["sizes"][n]
        lines.append(f"| {n} | {row['p_scalar']:.15g} | {row['p_h4_coefficient']:.8g} |")
    lines.extend([
        "",
        "## Fixed beta=7/2 GLS",
        "",
        f"- pc = `{fit['pc']:.15g}`",
        f"- scalar amplitude = `{fit['amplitude']:.8g}`",
        f"- chi-square = `{fit['chi_square']:.6g} / {fit['df']}`",
        "",
        "## Parameter-free doubling-lineage reconstruction",
        "",
        "| lineage | pc_hat | SE | scalar amplitude | SE |",
        "|---|---:|---:|---:|---:|",
    ])
    for key, row in payload["lineages"].items():
        lines.append(
            f"| {key} | {row['pc_hat']:.15g} | {row['pc_hat_se']:.3g} | "
            f"{row['amplitude']:.7g} | {row['amplitude_se']:.3g} |"
        )
    if payload["two_lineage_consistency"]:
        c = payload["two_lineage_consistency"]
        lines.extend([
            "",
            "Two-lineage consistency:",
            f"- pc-hat difference z = `{c['pc_difference_z']:.4g}`" if c["pc_difference_z"] is not None else "- pc-hat difference z unavailable",
            f"- amplitude difference z = `{c['amplitude_difference_z']:.4g}`" if c["amplitude_difference_z"] is not None else "- amplitude difference z unavailable",
        ])
    lines.extend([
        "",
        "## Boundary",
        "",
        "This is a scalar-sector diagnostic. H12/H8 contamination, slope anisotropy, and the conditional matching-parity assignment remain explicit alternatives.",
        "",
    ])
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--histograms", nargs="+", required=True, type=Path)
    parser.add_argument("--beta", type=float, default=3.5)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    records = merge_inputs(args.histograms)
    payload = calculate(records, beta=args.beta)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.csv.parent.mkdir(parents=True, exist_ok=True)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_csv(args.csv, payload)
    args.report.write_text(report(payload), encoding="utf-8")
    print(args.json)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

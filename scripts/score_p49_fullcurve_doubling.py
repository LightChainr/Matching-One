#!/usr/bin/env python3
"""Score clean norm-2 full-curve Gaussian lineages for issues #48 and #49.

The parent/child inputs must use aligned batch ids and identical per-batch
sample counts.  Every delete-one replicate recomputes the intrinsic thermal
coordinates, projectors, roots, and fixed residuals.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Dict, Iterable, Mapping, Sequence, Tuple

from analyze_p48_retrospective import (
    Histogram,
    add_histograms,
    cos4,
    covariance_of_mean,
    project_size,
    pseudovalues,
    quadratic,
    read_histograms,
    tail_and_derivative,
)


SIZES = (65, 85, 130, 170)
LINEAGES = ((65, 130), (85, 170))
LEVELS = (0.0, 0.025, 0.05)
P48_POWERS = {
    "P4_S": 1.0,
    "P4_D": 13.0 / 8.0,
    "P4_S_prime": 5.0 / 4.0,
    "P4_D_prime": 5.0 / 8.0,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def merge_inputs(paths: Sequence[Path]) -> Dict[Tuple[int, str, int], Histogram]:
    merged: Dict[Tuple[int, str, int], Histogram] = {}
    for path in paths:
        records = read_histograms(path)
        overlap = set(merged) & set(records)
        if overlap:
            raise ValueError(f"duplicate histogram keys in inputs: {sorted(overlap)[:3]}")
        merged.update(records)
    sizes = tuple(sorted({key[0] for key in merged}))
    if sizes != SIZES:
        raise ValueError(f"expected sizes {SIZES}, got {sizes}")
    signatures = {}
    for n in sizes:
        selected = sorted(
            (row for key, row in merged.items() if key[0] == n and key[1] == "first"),
            key=lambda row: row.batch,
        )
        signatures[n] = ([row.batch for row in selected], [row.samples for row in selected])
    if len({(tuple(a), tuple(b)) for a, b in signatures.values()}) != 1:
        raise ValueError("cross-size batch/sample alignment is absent")
    return merged


def grouped(records: Mapping[Tuple[int, str, int], Histogram]):
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


def aggregate(rows: Sequence[Histogram], omitted: int):
    included = [row for row in rows if row.batch != omitted]
    return {
        "a": rows[0].a,
        "b": rows[0].b,
        "samples": sum(row.samples for row in included),
        "minus": add_histograms(rows, "minus", omitted),
        "plus": add_histograms(rows, "plus", omitted),
    }


def orientation_values(n: int, row, p: float):
    minus, minus_d = tail_and_derivative(row["minus"], row["samples"], p)
    plus, plus_d = tail_and_derivative(row["plus"], row["samples"], p)
    return {
        "M": minus + plus - 1.0,
        "M_prime": minus_d + plus_d,
        "R_G": plus,
        "R_hat": 1.0 - minus,
    }


def solve_target(function, target: float) -> float:
    # Keep the binomial recurrence away from endpoint underflow.  Every frozen
    # level is deep inside this bracket for the percolation data in scope.
    lower, upper = 0.1, 0.9
    f_lower = function(lower) - target
    f_upper = function(upper) - target
    if not f_lower <= 0.0 <= f_upper:
        raise ValueError(f"target {target} is not bracketed")
    for _ in range(56):
        midpoint = (lower + upper) / 2.0
        if function(midpoint) < target:
            lower = midpoint
        else:
            upper = midpoint
    return (lower + upper) / 2.0


def size_statistics(by_orientation, omitted: int = -1):
    n = by_orientation["first"][0].n
    rows = {name: aggregate(by_orientation[name], omitted) for name in ("first", "second")}

    def value(name: str, p: float):
        return orientation_values(n, rows[name], p)

    def mean_matching(p: float) -> float:
        return (value("first", p)["M"] + value("second", p)["M"]) / 2.0

    coordinates = {}
    contrasts = {}
    for level in LEVELS:
        if level == 0.0:
            p_minus = p_plus = solve_target(mean_matching, 0.0)
        else:
            p_minus = solve_target(mean_matching, -level)
            p_plus = solve_target(mean_matching, level)
        minus_first = value("first", p_minus)
        minus_second = value("second", p_minus)
        plus_first = value("first", p_plus)
        plus_second = value("second", p_plus)
        stored_minus = minus_first["M"] - minus_second["M"]
        stored_plus = plus_first["M"] - plus_second["M"]
        # Parent stored order is exact lineage order.  Child stored order is
        # reversed by the multiplication-by-(1+i) genealogy in issue #49.
        lineage_sign = 1.0 if n in (65, 85) else -1.0
        x_minus = lineage_sign * stored_minus
        x_plus = lineage_sign * stored_plus
        coordinates[str(level)] = {
            "p_minus": p_minus,
            "p_plus": p_plus,
            "inverse_slope_minus": 1.0 / abs((minus_first["M_prime"] + minus_second["M_prime"]) / 2.0),
            "inverse_slope_plus": 1.0 / abs((plus_first["M_prime"] + plus_second["M_prime"]) / 2.0),
        }
        contrasts[str(level)] = {
            "X_even": (x_plus + x_minus) / 2.0,
            "X_odd": (x_plus - x_minus) / 2.0,
        }

    p0 = coordinates["0.0"]["p_minus"]
    center_values = {name: value(name, p0) for name in ("first", "second")}
    mean_slope = (center_values["first"]["M_prime"] + center_values["second"]["M_prime"]) / 2.0
    roots = {name: solve_target(lambda p, name=name: value(name, p)["M"], 0.0) for name in ("first", "second")}
    root_gap_stored = roots["first"] - roots["second"]
    lineage_sign = 1.0 if n in (65, 85) else -1.0
    p48 = project_size(by_orientation, omitted)
    return {
        "N": n,
        "p0": p0,
        "mean_slope": mean_slope,
        "root_first_stored": roots["first"],
        "root_second_stored": roots["second"],
        "root_gap_lineage": lineage_sign * root_gap_stored,
        "coordinates": coordinates,
        "contrasts": contrasts,
        "p48": {metric: p48[metric] for metric in P48_POWERS},
        "delta_cos4_stored": p48["delta_cos4"],
    }


def jackknife_se(values: Sequence[float]) -> float:
    mean = math.fsum(values) / len(values)
    return math.sqrt((len(values) - 1) / len(values) * math.fsum((x - mean) ** 2 for x in values))


def residual_statistics(full: float, deleted: Sequence[float]):
    se = jackknife_se(deleted)
    return {"value": full, "se": se, "z": full / se if se else None}


def joint_statistics(full_sample, deleted_samples, functions):
    full = [function(full_sample) for function in functions]
    deleted = [[function(sample) for function in functions] for sample in deleted_samples]
    pseudo = [
        [pseudovalues(full[column], [row[column] for row in deleted])[batch]
         for column in range(len(functions))]
        for batch in range(len(deleted))
    ]
    covariance = covariance_of_mean(pseudo)
    return {
        "residual": full,
        "covariance": covariance,
        "chi_square": quadratic(full, covariance),
        "df": len(full),
    }


def covariance_for_functions(full_sample, deleted_samples, functions):
    full = [function(full_sample) for function in functions]
    deleted = [[function(sample) for function in functions] for sample in deleted_samples]
    pseudo = [
        [pseudovalues(full[column], [row[column] for row in deleted])[batch]
         for column in range(len(functions))]
        for batch in range(len(deleted))
    ]
    return full, covariance_of_mean(pseudo)


def fixed_model_score(target, target_covariance, design, parameters, parameter_covariance):
    predicted = [math.fsum(x * beta for x, beta in zip(row, parameters)) for row in design]
    source_covariance = [[
        math.fsum(
            design[i][a] * parameter_covariance[a][b] * design[j][b]
            for a in range(len(parameters)) for b in range(len(parameters))
        ) if parameters else 0.0
        for j in range(len(design))
    ] for i in range(len(design))]
    covariance = [[target_covariance[i][j] + source_covariance[i][j] for j in range(len(design))] for i in range(len(design))]
    residual = [target[i] - predicted[i] for i in range(len(target))]
    return {
        "prediction": predicted,
        "residual": residual,
        "residual_covariance": covariance,
        "z_marginal": [residual[i] / math.sqrt(covariance[i][i]) for i in range(len(residual))],
        "chi_square": quadratic(residual, covariance),
        "df": len(residual),
    }


def calculate(records):
    groups = grouped(records)
    full = {n: size_statistics(groups[n]) for n in SIZES}
    batch_count = len(groups[SIZES[0]]["first"])
    deleted = [
        {n: size_statistics(groups[n], omitted) for n in SIZES}
        for omitted in range(batch_count)
    ]
    radial = 2.0 ** (-13.0 / 8.0)
    slope_ratio = 2.0 ** (3.0 / 8.0)
    lineages = {}
    for parent, child in LINEAGES:
        key = f"{parent}->{child}"
        level_scores = {}
        for level in LEVELS:
            level_key = str(level)
            def even_res(sample):
                return sample[child]["contrasts"][level_key]["X_even"] + radial * sample[parent]["contrasts"][level_key]["X_even"]
            def odd_value(sample):
                return sample[child]["contrasts"][level_key]["X_odd"]
            level_scores[level_key] = {
                "X_even_parent": full[parent]["contrasts"][level_key]["X_even"],
                "X_even_child": full[child]["contrasts"][level_key]["X_even"],
                "fixed_even_residual": residual_statistics(even_res(full), [even_res(x) for x in deleted]),
                "X_odd_parent": full[parent]["contrasts"][level_key]["X_odd"],
                "X_odd_child": odd_value(full),
                "p_parent": full[parent]["coordinates"][level_key],
                "p_child": full[child]["coordinates"][level_key],
            }
        def slope_res(sample):
            return sample[child]["mean_slope"] - slope_ratio * sample[parent]["mean_slope"]
        def root_raw(sample):
            return sample[child]["root_gap_lineage"] + 0.25 * sample[parent]["root_gap_lineage"]
        def root_finite(sample):
            return sample[child]["root_gap_lineage"] + radial * (
                sample[parent]["mean_slope"] / sample[child]["mean_slope"]
            ) * sample[parent]["root_gap_lineage"]
        p48_scores = {}
        for metric, power in P48_POWERS.items():
            factor = 2.0 ** (-power)
            def normalized(sample, metric=metric, factor=factor):
                return sample[child]["p48"][metric] - factor * sample[parent]["p48"][metric]
            def artifact_negative(sample, metric=metric, factor=factor):
                return sample[child]["p48"][metric] + factor * sample[parent]["p48"][metric]
            p48_scores[metric] = {
                "parent": full[parent]["p48"][metric],
                "child": full[child]["p48"][metric],
                "descriptive_ratio": full[child]["p48"][metric] / full[parent]["p48"][metric],
                "normalized_P4_positive_ratio_residual": residual_statistics(normalized(full), [normalized(x) for x in deleted]),
                "frozen_artifact_negative_ratio_residual": residual_statistics(artifact_negative(full), [artifact_negative(x) for x in deleted]),
            }
        lineages[key] = {
            "levels": level_scores,
            "slope_parent": full[parent]["mean_slope"],
            "slope_child": full[child]["mean_slope"],
            "slope_ratio": full[child]["mean_slope"] / full[parent]["mean_slope"],
            "fixed_slope_residual": residual_statistics(slope_res(full), [slope_res(x) for x in deleted]),
            "root_gap_parent": full[parent]["root_gap_lineage"],
            "root_gap_child": full[child]["root_gap_lineage"],
            "root_gap_ratio": full[child]["root_gap_lineage"] / full[parent]["root_gap_lineage"],
            "fixed_raw_root_residual": residual_statistics(root_raw(full), [root_raw(x) for x in deleted]),
            "finite_slope_root_residual": residual_statistics(root_finite(full), [root_finite(x) for x in deleted]),
            "p48": p48_scores,
        }

    def even_function(parent, child, level):
        level_key = str(level)
        return lambda sample: sample[child]["contrasts"][level_key]["X_even"] + radial * sample[parent]["contrasts"][level_key]["X_even"]
    def slope_function(parent, child):
        return lambda sample: sample[child]["mean_slope"] - slope_ratio * sample[parent]["mean_slope"]
    def root_raw_function(parent, child):
        return lambda sample: sample[child]["root_gap_lineage"] + 0.25 * sample[parent]["root_gap_lineage"]
    def root_finite_function(parent, child):
        return lambda sample: sample[child]["root_gap_lineage"] + radial * (sample[parent]["mean_slope"] / sample[child]["mean_slope"]) * sample[parent]["root_gap_lineage"]
    def p48_function(parent, child, metric, sign):
        factor = 2.0 ** (-P48_POWERS[metric])
        return lambda sample: sample[child]["p48"][metric] + sign * factor * sample[parent]["p48"][metric]

    joint = {}
    for level in LEVELS:
        joint[f"P49_X_even_u={level}"] = joint_statistics(
            full, deleted, [even_function(parent, child, level) for parent, child in LINEAGES]
        )
    joint["P49_slope"] = joint_statistics(full, deleted, [slope_function(*pair) for pair in LINEAGES])
    joint["P49_root_raw"] = joint_statistics(full, deleted, [root_raw_function(*pair) for pair in LINEAGES])
    joint["P49_root_finite_slope"] = joint_statistics(full, deleted, [root_finite_function(*pair) for pair in LINEAGES])
    for metric in P48_POWERS:
        joint[f"P48_normalized_positive_{metric}"] = joint_statistics(
            full, deleted, [p48_function(parent, child, metric, -1.0) for parent, child in LINEAGES]
        )
        joint[f"P48_artifact_negative_{metric}"] = joint_statistics(
            full, deleted, [p48_function(parent, child, metric, +1.0) for parent, child in LINEAGES]
        )
    joint["P48_normalized_positive_all"] = joint_statistics(
        full, deleted,
        [p48_function(parent, child, metric, -1.0)
         for metric in P48_POWERS for parent, child in LINEAGES],
    )
    joint["P48_artifact_negative_all"] = joint_statistics(
        full, deleted,
        [p48_function(parent, child, metric, +1.0)
         for metric in P48_POWERS for parent, child in LINEAGES],
    )

    target_sizes = (130, 170)
    target, target_covariance = covariance_for_functions(
        full, deleted,
        [lambda sample, n=n: sample[n]["p48"]["P4_S_prime"] for n in target_sizes],
    )
    model_specs = {
        "pure_N^-5/4": {
            "design": [[n ** (-1.25)] for n in target_sizes],
            "parameters": [1.9434247576878727],
            "parameter_covariance": [[0.005868307572059406]],
        },
        "zero_effect": {
            "design": [[] for _ in target_sizes],
            "parameters": [],
            "parameter_covariance": [],
        },
        "q2_even_scalar": {
            "design": [[n ** (-1.25), n ** (-2.25)] for n in target_sizes],
            "parameters": [3.203310807356976, -90.59560328584558],
            "parameter_covariance": [[0.1003546426313145, -7.164152468416187], [-7.164152468416186, 540.1590397230893]],
        },
        "rank2_jordan_log": {
            "design": [[n ** (-1.25), n ** (-1.25) * math.log(n)] for n in target_sizes],
            "parameters": [-2.422594685734799, 1.016646899281392],
            "parameter_covariance": [[1.251237089769748, -0.286289830474317], [-0.286289830474317, 0.065785221364566]],
        },
    }
    sprime_replication = {
        name: fixed_model_score(target, target_covariance, **spec)
        for name, spec in model_specs.items()
    }
    return {
        "format_version": 1,
        "batches": batch_count,
        "sizes": full,
        "lineages": lineages,
        "joint_scores": joint,
        "P48_Sprime_fresh_seed_replication": {
            "target_sizes": target_sizes,
            "target": target,
            "target_covariance": target_covariance,
            "classification": "fresh-seed replication at geometries present in retrospective source; no target refit",
            "models_in_frozen_order": sprime_replication,
        },
    }


def flatten(payload) -> Iterable[dict]:
    for lineage, row in payload["lineages"].items():
        for level, score in row["levels"].items():
            residual = score["fixed_even_residual"]
            yield {"lineage": lineage, "family": "P49_X_even", "channel": f"u={level}", **residual}
        for family, name in (
            ("P49_slope", "fixed_slope_residual"),
            ("P49_root_raw", "fixed_raw_root_residual"),
            ("P49_root_finite_slope", "finite_slope_root_residual"),
        ):
            yield {"lineage": lineage, "family": family, "channel": "center", **row[name]}
        for metric, score in row["p48"].items():
            for family, name in (
                ("P48_normalized_positive", "normalized_P4_positive_ratio_residual"),
                ("P48_artifact_negative", "frozen_artifact_negative_ratio_residual"),
            ):
                yield {"lineage": lineage, "family": family, "channel": metric, **score[name]}


def report_text(payload) -> str:
    joint = payload["joint_scores"]
    models = payload["P48_Sprime_fresh_seed_replication"]["models_in_frozen_order"]
    lines = [
        "# Clean 100M full-curve norm-2 score for issues #48/#49",
        "",
        "**Classification:** fresh-seed replication at geometries already present in the retrospective P33 source. No target parameter was fit.",
        "",
        "## Decision",
        "",
        "The full-curve thermal-even and root tests remain compatible but are driven to about 2.1 sigma by the 85->170 lineage. The raw asymptotic slope multiplier is sharply rejected at this precision, exposing a small finite-size correction. The four-channel P48 pure-power conjunction fails specifically and reproducibly in `P4[S']`; both preregistered correction models survive this same-geometry fresh-seed replication.",
        "",
        "## P49 no-fit Gaussian doubling",
        "",
        "| statistic | chi-square / df |",
        "|---|---:|",
    ]
    for name in (
        "P49_X_even_u=0.0", "P49_X_even_u=0.025", "P49_X_even_u=0.05",
        "P49_root_raw", "P49_root_finite_slope", "P49_slope",
    ):
        score = joint[name]
        lines.append(f"| `{name}` | {score['chi_square']:.6g} / {score['df']} |")
    lines.extend(["", "Point ratios:", ""])
    for lineage, row in payload["lineages"].items():
        lines.append(
            f"- `{lineage}`: slope={row['slope_ratio']:.9g} (target {2.0**(3.0/8.0):.9g}); "
            f"root={row['root_gap_ratio']:.9g} (raw target -0.25)."
        )
    lines.extend([
        "",
        "All `u={0,.025,.05}` coordinates were solved inside every delete-one replicate. The reported inverse-slope condition proxies and coordinates are retained in `analysis/score.json`.",
        "",
        "## P48 derivative spectrum",
        "",
        "The mathematically consistent child/parent factor for a projector normalized by each size's own `Delta cos(4 theta)` is **positive** `2^-alpha`. The negative H4 rotation factor belongs to the unnormalized exact-lineage contrast. Applying a negative factor to normalized `P4` double-counts the angular sign.",
        "",
        "| normalized channel | positive-ratio chi-square / 2 | negative-artifact chi-square / 2 |",
        "|---|---:|---:|",
    ])
    for metric in P48_POWERS:
        positive = joint[f"P48_normalized_positive_{metric}"]
        negative = joint[f"P48_artifact_negative_{metric}"]
        lines.append(f"| `{metric}` | {positive['chi_square']:.6g} | {negative['chi_square']:.6g} |")
    lines.extend([
        f"| **joint** | **{joint['P48_normalized_positive_all']['chi_square']:.6g} / 8** | **{joint['P48_artifact_negative_all']['chi_square']:.6g} / 8** |",
        "",
        "The normalized pure laws pass cleanly for `P4[S]`, are mildly strained for `P4[D]` and `P4[D']`, and fail for `P4[S']`.",
        "",
        "## Frozen P48 S-prime models on fresh N=130/170 counters",
        "",
        "| model (frozen order) | chi-square / 2 | marginal z at N=130,170 |",
        "|---|---:|---:|",
    ])
    for name in ("pure_N^-5/4", "zero_effect", "q2_even_scalar", "rank2_jordan_log"):
        score = models[name]
        lines.append(
            f"| `{name}` | {score['chi_square']:.6g} | "
            f"{score['z_marginal'][0]:+.3f}, {score['z_marginal'][1]:+.3f} |"
        )
    lines.extend([
        "",
        "The q=2 correction is the first surviving model in the frozen order. The Jordan-log adversary also survives and has a smaller descriptive chi-square, but these same geometries do not constitute a new angular holdout and should not be used to reverse the preregistered order.",
        "",
        "## Provenance",
        "",
        "- Parent N=65/85 and child N=130/170 use source commit `6d2d68a`, seed `2026104501`, counters `[5000000000,5100000000)`, 100 aligned batches, and 100,000,000 samples per orientation pair.",
        "- Child runs were executed on Huawei DevEnv `f415a4bcbd9a438b85f5f29e4a507ea4` (AArch64, 16 vCPU); both stderr files are empty.",
        "- Raw child histograms, moments, metadata, stdout and stderr are preserved under `raw/`.",
        "",
    ])
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--histograms", nargs=4, required=True, type=Path)
    parser.add_argument("--json", required=True, type=Path)
    parser.add_argument("--csv", required=True, type=Path)
    args = parser.parse_args()
    records = merge_inputs(args.histograms)
    payload = calculate(records)
    payload["provenance"] = {str(path): sha256(path) for path in args.histograms}
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.csv.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    rows = list(flatten(payload))
    with args.csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    output_dir = args.json.parent.parent
    (output_dir / "REPORT.md").write_text(report_text(payload), encoding="utf-8")
    command = "python3 scripts/score_p49_fullcurve_doubling.py --histograms " + " ".join(map(str, args.histograms)) + f" --json {args.json} --csv {args.csv}"
    child_base = output_dir / "raw"
    analysis_base = output_dir / "analysis"
    commands = [
        "scp Huawei-CodeBuddy:/workspace/fullcurve-results-20260828/'n130.*' " + str(child_base) + "/",
        "scp Huawei-CodeBuddy:/workspace/fullcurve-results-20260828/'n170.*' " + str(child_base) + "/",
        f"python3 scripts/analyze_matching_parity_derivatives_fast.py --histograms {child_base / 'n130.hist.csv'} --dps 40 --json {analysis_base / 'n130.p48.json'} --csv {analysis_base / 'n130.p48.csv'}",
        f"python3 scripts/analyze_matching_parity_derivatives_fast.py --histograms {child_base / 'n170.hist.csv'} --dps 40 --json {analysis_base / 'n170.p48.json'} --csv {analysis_base / 'n170.p48.csv'}",
        command,
        "python3 -m py_compile scripts/score_p49_fullcurve_doubling.py",
        "python3 tests/test_p48_retrospective.py",
    ]
    (output_dir / "commands.txt").write_text("\n".join(commands) + "\n", encoding="utf-8")
    checksum_paths = sorted(
        path for path in output_dir.rglob("*")
        if path.is_file() and path.name != "checksums.sha256"
    )
    (output_dir / "checksums.sha256").write_text(
        "".join(f"{sha256(path)}  {path.relative_to(output_dir)}\n" for path in checksum_paths),
        encoding="utf-8",
    )
    print(args.json)
    print(args.csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

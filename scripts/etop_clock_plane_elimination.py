#!/usr/bin/env python3
"""Profile the production E_top ray against the measured C/W clock plane.

The fit is an errors-in-variables distance to a linear hyperplane.  Each row
uses its complete same-batch covariance, so uncertainty in A_top, C and W is
not treated as fixed design.  The P205 quotient prism is reconstructed from
its raw aligned histogram/moment batches and used as an independent geometry
block rather than as another fitted lineage.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import tempfile
from typing import Any, Iterable, Mapping, Sequence

import mpmath as mp

from etop_rank1_elimination import _chi2_survival
from p205_quotient_etop_ray_test import SIZES, materialize_sources
from rank_plane_crosswalk import (
    _clock_projection,
    _projected_p_derivative,
    _projected_state,
    combine,
    cos4,
    endpoint_observables,
    intrinsic_center,
    read_histograms,
    read_moments,
)


ROOT = Path(__file__).resolve().parents[1]
CROSSWALK = ROOT / "results/rank-plane-crosswalk/latest.json"
P205_EXISTING = ROOT / "results/etop-p205-quotient-ray-test/latest.json"
PRODUCTION_GROUPS = ("P49", "P43", "P50", "P57")
METRICS = ("P4_A_top", "P4_E_top", "P4_C", "P4_W")
MODELS = {
    "A_ray": ("P4_A_top",),
    "A_plus_C": ("P4_A_top", "P4_C"),
    "A_plus_W": ("P4_A_top", "P4_W"),
    "A_plus_C_plus_W": ("P4_A_top", "P4_C", "P4_W"),
}
ALPHA = mp.mpf("0.01")


def number(value: Any, digits: int = 16) -> float:
    return float(mp.nstr(mp.mpf(str(value)), digits))


def covariance_of_mean(rows: Sequence[Sequence[mp.mpf]]) -> list[list[mp.mpf]]:
    count = len(rows)
    return [
        [
            mp.fsum(row[i] * row[j] for row in rows) / (count * (count - 1))
            for j in range(len(METRICS))
        ]
        for i in range(len(METRICS))
    ]


def production_rows(covariance_key: str) -> list[dict[str, Any]]:
    payload = json.loads(CROSSWALK.read_text())
    output = []
    for dataset in payload["datasets"]:
        group = dataset["id"].split("-", 1)[0]
        if group not in PRODUCTION_GROUPS:
            continue
        order = dataset["covariance_metric_order"]
        positions = [order.index(metric) for metric in METRICS]
        covariance = dataset[covariance_key]
        output.append({
            "id": dataset["id"],
            "dependency_group": group,
            "N": dataset["N"],
            "estimate": [mp.mpf(str(dataset["point"][metric])) for metric in METRICS],
            "covariance": [
                [mp.mpf(str(covariance[i][j])) for j in positions]
                for i in positions
            ],
        })
    if len(output) != 8:
        raise ValueError(f"expected eight high-statistics production rows, got {len(output)}")
    return output


def reconstruct_p205(covariance_kind: str) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    existing = json.loads(P205_EXISTING.read_text())
    expected = {row["N"]: row for row in existing["pairs"]}
    output = []
    with tempfile.TemporaryDirectory(prefix="p205-clock-plane-") as directory:
        paths, hashes = materialize_sources(Path(directory))
        for n in SIZES:
            histograms = read_histograms(paths[n]["histogram"])
            moments = read_moments(paths[n]["moments"])
            batches = sorted({key[2] for key in histograms})
            by_orientation = {
                orientation: [histograms[(n, orientation, batch)] for batch in batches]
                for orientation in ("first", "second")
            }
            total_first = combine(by_orientation["first"])
            total_second = combine(by_orientation["second"])
            delta_cos4 = cos4(total_first.a, total_first.b) - cos4(
                total_second.a, total_second.b
            )
            p0 = intrinsic_center(total_first, total_second)
            first = endpoint_observables(total_first, p0)
            second = endpoint_observables(total_second, p0)
            state = _projected_state(first, second, delta_cos4)
            state.update({
                metric: mp.fsum(
                    _clock_projection(
                        moments[(n, "first", batch)],
                        moments[(n, "second", batch)],
                        delta_cos4,
                    )[metric]
                    for batch in batches
                ) / len(batches)
                for metric in ("P4_C", "P4_W")
            })
            derivative = _projected_p_derivative(first, second, delta_cos4)
            center_slope = (first["S_birth"] + second["S_birth"]) / 2
            fixed, intrinsic = [], []
            for batch in batches:
                left = endpoint_observables(by_orientation["first"][batch], p0)
                right = endpoint_observables(by_orientation["second"][batch], p0)
                batch_state = _projected_state(left, right, delta_cos4)
                batch_state.update(_clock_projection(
                    moments[(n, "first", batch)],
                    moments[(n, "second", batch)],
                    delta_cos4,
                ))
                fixed_row = [batch_state[metric] - state[metric] for metric in METRICS]
                center_shift = -(left["A_top"] + right["A_top"]) / (2 * center_slope)
                fixed.append(fixed_row)
                intrinsic.append([
                    fixed_row[index] + derivative[metric] * center_shift
                    for index, metric in enumerate(METRICS)
                ])
            reconstructed_ae = [state["P4_A_top"], state["P4_E_top"]]
            reference_ae = [mp.mpf(str(value)) for value in expected[n]["intrinsic"]["estimate"]]
            replay_error = max(abs(a - b) for a, b in zip(reconstructed_ae, reference_ae))
            if replay_error > mp.mpf("1e-14"):
                raise ValueError(f"P205 N{n} A/E replay changed by {replay_error}")
            output.append({
                "id": f"P205-prism-N{n}",
                "dependency_group": "P205-prism",
                "N": n,
                "estimate": [state[metric] for metric in METRICS],
                "covariance": covariance_of_mean(
                    intrinsic if covariance_kind == "intrinsic" else fixed
                ),
                "replay_AE_max_abs_error": replay_error,
            })
    return output, {"raw_sha256": hashes, "existing_result_sha256": hashlib.sha256(
        P205_EXISTING.read_bytes()).hexdigest()}


def quadratic_form(row: Mapping[str, Any], predictors: Sequence[str], beta: Sequence[float]) -> float:
    metric_index = {metric: index for index, metric in enumerate(METRICS)}
    coefficient = [0.0] * len(METRICS)
    coefficient[metric_index["P4_E_top"]] = 1.0
    for value, metric in zip(beta, predictors):
        coefficient[metric_index[metric]] -= value
    estimate = [float(value) for value in row["estimate"]]
    covariance = [[float(value) for value in line] for line in row["covariance"]]
    residual = sum(left * right for left, right in zip(coefficient, estimate))
    variance = sum(
        coefficient[i] * covariance[i][j] * coefficient[j]
        for i in range(len(METRICS)) for j in range(len(METRICS))
    )
    if variance <= 0:
        return math.inf
    return residual * residual / variance


def objective(rows: Sequence[Mapping[str, Any]], predictors: Sequence[str], beta: Sequence[float]) -> float:
    return sum(quadratic_form(row, predictors, beta) for row in rows)


def nelder_mead(function, start: Sequence[float], step: float = 1.0,
                max_iter: int = 900, tolerance: float = 1e-11):
    dimension = len(start)
    simplex = [list(start)]
    for index in range(dimension):
        vertex = list(start)
        vertex[index] += step
        simplex.append(vertex)
    values = [function(point) for point in simplex]
    for _ in range(max_iter):
        ranked = sorted(zip(values, simplex), key=lambda item: item[0])
        values = [item[0] for item in ranked]
        simplex = [item[1] for item in ranked]
        if max(abs(value - values[0]) for value in values) < tolerance:
            break
        centroid = [sum(point[j] for point in simplex[:-1]) / dimension
                    for j in range(dimension)]
        reflected = [centroid[j] + (centroid[j] - simplex[-1][j])
                     for j in range(dimension)]
        reflected_value = function(reflected)
        if values[0] <= reflected_value < values[-2]:
            simplex[-1], values[-1] = reflected, reflected_value
            continue
        if reflected_value < values[0]:
            expanded = [centroid[j] + 2 * (reflected[j] - centroid[j])
                        for j in range(dimension)]
            expanded_value = function(expanded)
            if expanded_value < reflected_value:
                simplex[-1], values[-1] = expanded, expanded_value
            else:
                simplex[-1], values[-1] = reflected, reflected_value
            continue
        contracted = [centroid[j] + 0.5 * (simplex[-1][j] - centroid[j])
                      for j in range(dimension)]
        contracted_value = function(contracted)
        if contracted_value < values[-1]:
            simplex[-1], values[-1] = contracted, contracted_value
            continue
        best = simplex[0]
        simplex = [best] + [
            [best[j] + 0.5 * (point[j] - best[j]) for j in range(dimension)]
            for point in simplex[1:]
        ]
        values = [function(point) for point in simplex]
    ranked = sorted(zip(values, simplex), key=lambda item: item[0])
    return ranked[0][0], ranked[0][1]


def fit(rows: Sequence[Mapping[str, Any]], predictors: Sequence[str]) -> dict[str, Any]:
    function = lambda beta: objective(rows, predictors, beta)
    starts = itertools.product((-8.0, -2.0, 0.0, 2.0, 8.0), repeat=len(predictors))
    candidates = []
    for start in starts:
        value, beta = nelder_mead(function, start, step=1.25)
        if all(math.isfinite(item) and abs(item) < 1e4 for item in beta):
            candidates.append((value, beta))
    candidates.sort(key=lambda item: item[0])
    if not candidates:
        raise ValueError("optimizer found no finite candidate")
    value, beta = candidates[0]
    return {
        "chi2": value,
        "beta": {metric: coefficient for metric, coefficient in zip(predictors, beta)},
        "optimizer": {
            "method": "deterministic 5^k multistart Nelder-Mead",
            "starts": 5 ** len(predictors),
            "best_second_gap": candidates[1][0] - value if len(candidates) > 1 else None,
            "objective_replay_error": abs(function(beta) - value),
        },
    }


def score_block(covariance_kind: str) -> dict[str, Any]:
    production_key = (
        "covariance_intrinsic_center_first_order_influence"
        if covariance_kind == "intrinsic"
        else "covariance_fixed_center_exact_batch_estimator"
    )
    production = production_rows(production_key)
    p205, provenance = reconstruct_p205(covariance_kind)
    joint = production + p205
    models = {}
    for name, predictors in MODELS.items():
        production_fit = fit(production, predictors)
        joint_fit = fit(joint, predictors)
        production_df = len(production) - len(predictors)
        joint_df = len(joint) - len(predictors)
        predictive = max(0.0, joint_fit["chi2"] - production_fit["chi2"])
        fixed_prediction = objective(
            p205, predictors,
            [production_fit["beta"][metric] for metric in predictors],
        )
        models[name] = {
            "predictors": list(predictors),
            "production": {**production_fit, "df": production_df,
                           "p": number(_chi2_survival(production_fit["chi2"], production_df))},
            "joint": {**joint_fit, "df": joint_df,
                      "p": number(_chi2_survival(joint_fit["chi2"], joint_df))},
            "P205_profiled_prediction": {
                "chi2": predictive, "df": len(p205),
                "p": number(_chi2_survival(predictive, len(p205))),
            },
            "P205_fixed_production_coefficients": {
                "chi2": fixed_prediction, "df": len(p205),
                "p": number(_chi2_survival(fixed_prediction, len(p205))),
                "role": "stress sensitivity; ignores coefficient uncertainty",
            },
        }
    ray = models["A_ray"]
    for name in ("A_plus_C", "A_plus_W"):
        for block in ("production", "joint"):
            delta = max(0.0, ray[block]["chi2"] - models[name][block]["chi2"])
            models[name][block]["improvement_over_A_ray"] = {
                "delta_chi2": delta, "df": 1,
                "p": number(_chi2_survival(delta, 1)),
            }
    for baseline in ("A_plus_C", "A_plus_W"):
        for block in ("production", "joint"):
            delta = max(0.0, models[baseline][block]["chi2"] -
                        models["A_plus_C_plus_W"][block]["chi2"])
            models["A_plus_C_plus_W"][block][f"improvement_over_{baseline}"] = {
                "delta_chi2": delta, "df": 1,
                "p": number(_chi2_survival(delta, 1)),
            }
    return {
        "covariance_semantics": covariance_kind,
        "models": models,
        "P205_rows": [{
            "id": row["id"], "N": row["N"],
            "metric_order": list(METRICS),
            "estimate": [number(value) for value in row["estimate"]],
            "covariance": [[number(value) for value in line] for line in row["covariance"]],
            "replay_AE_max_abs_error": number(row["replay_AE_max_abs_error"]),
        } for row in p205],
        "provenance": provenance,
    }


def build_report() -> dict[str, Any]:
    mp.mp.dps = 60
    primary = score_block("intrinsic")
    sensitivity = score_block("fixed")
    models = primary["models"]
    return {
        "schema": "matching-one/etop-clock-plane-elimination/v1",
        "status": "zero_new_sample_production_plus_independent_archive_analysis",
        "alpha": number(ALPHA),
        "question": "does a measured C/W clock-plane coordinate explain the cross-geometry A_top/E_top ray failure, and does that explanation transport to the independent P205 prism?",
        "source": {
            "production_crosswalk": str(CROSSWALK.relative_to(ROOT)),
            "production_crosswalk_sha256": hashlib.sha256(CROSSWALK.read_bytes()).hexdigest(),
            "P205_existing_result": str(P205_EXISTING.relative_to(ROOT)),
        },
        "primary": primary,
        "fixed_center_sensitivity": sensitivity,
        "decision": {
            "A_ray": "eliminated jointly",
            "A_plus_C": "a common plane survives jointly and under coefficient-uncertainty-profiled P205 prediction",
            "A_plus_W": "survives at alpha .01 but is substantially weaker on P205 transport",
            "A_plus_C_plus_W": "does not improve significantly over A_plus_C",
            "mechanism_update": "one measured clock-plane direction absorbs the apparent E_top ray rotation; the independent prism favors a C-like translation-clock plane over a W-only lifetime plane, without yet fixing universal numerical coefficients",
        },
        "claim_boundary": "This is a covariance-aware linear source decomposition of the measured topology/clock coordinates. The fixed production best-fit coefficients do not directly predict the much more precise P205 rows; survival is after profiling their training uncertainty. The result does not identify C as a continuum primary, equate it with an ordinary local Potts energy field, or prove a unique microscopic mixing law.",
    }


def markdown(report: Mapping[str, Any]) -> str:
    p = report["primary"]["models"]
    lines = [
        "# E_top clock-plane production elimination", "",
        "## Answer", "",
        "The one-coordinate `E_top = beta_A A_top` ray is eliminated after the",
        "independent P205 prism is added.  A second measured clock coordinate",
        "changes the model ranking without adding new Monte Carlo samples.", "",
        "| model | production chi2/df (p) | joint chi2/df (p) | profiled P205 p |", "|---|---:|---:|---:|",
    ]
    labels = {"A_ray": "A only", "A_plus_C": "A + C", "A_plus_W": "A + W", "A_plus_C_plus_W": "A + C + W"}
    for name in ("A_ray", "A_plus_C", "A_plus_W", "A_plus_C_plus_W"):
        row = p[name]
        lines.append(
            f"| {labels[name]} | {row['production']['chi2']:.6g}/{row['production']['df']} ({row['production']['p']:.5g}) "
            f"| {row['joint']['chi2']:.6g}/{row['joint']['df']} ({row['joint']['p']:.5g}) "
            f"| {row['P205_profiled_prediction']['p']:.5g} |"
        )
    c = p["A_plus_C"]
    w = p["A_plus_W"]
    cw = p["A_plus_C_plus_W"]
    lines += [
        "", "The production-only coefficients are:", "",
        f"- `A+C`: `{c['production']['beta']}`",
        f"- `A+W`: `{w['production']['beta']}`",
        "",
        f"Adding `C` to the ray improves the joint profile by `Delta chi2={c['joint']['improvement_over_A_ray']['delta_chi2']:.6g}` on one df.  Adding `W` also rescues the ray at the declared alpha, but its coefficient-uncertainty-profiled P205 score is near the boundary (`p={w['P205_profiled_prediction']['p']:.5g}`) whereas `A+C` remains comfortable (`p={c['P205_profiled_prediction']['p']:.5g}`).  This is common-plane compatibility, not a claim that the production point estimates of the coefficients predict P205: the fixed-coefficient stress rows reject both because the production block leaves those coefficients broad.",
        "",
        f"Adding `W` after `A+C` gains only `Delta chi2={cw['joint']['improvement_over_A_plus_C']['delta_chi2']:.6g}` on one df (`p={cw['joint']['improvement_over_A_plus_C']['p']:.5g}`), so the current data do not require a third coordinate.",
        "", "## Scientific card", "",
        "- **Mechanism space changed:** universal topology ray versus a topology-plus-measured-clock plane.",
        "- **Result:** one measured clock coordinate absorbs the four-lineage ray failure; after profiling coefficient uncertainty, a common plane remains compatible with the independent quotient prism.  The `C` coordinate is the parsimonious current survivor, while `W` alone is not formally excluded at alpha .01.",
        "- **Not proved:** `C` is not thereby identified with a continuum primary or an independently marked local Potts energy field.",
        "- **Observer / sector / source:** `P4(A_top,E_top,C,W)` / threshold-rank and event-clock plane / aligned production batches.",
        "- **Dependency:** P49/P43/P50/P57 train the coefficients; P205 is an independent raw archive and enters only as a profiled external block.",
        "- **Next discriminator:** record an independently marked same-batch local singlet/energy row and ask whether it replaces `C` and pins the plane coefficients, rather than collecting another untyped E_top ray.",
        "", "## Reproduction", "", "```bash",
        "python3 scripts/etop_clock_plane_elimination.py --json results/etop-clock-plane-elimination/latest.json --markdown results/etop-clock-plane-elimination/REPORT.md",
        "python3 -m unittest discover -s tests -p 'test_etop_clock_plane_elimination.py'", "```", "",
    ]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args()
    report = build_report()
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.markdown.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    args.markdown.write_text(markdown(report))
    print(json.dumps(report["decision"], indent=2))


if __name__ == "__main__":
    main()

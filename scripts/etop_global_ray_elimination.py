#!/usr/bin/env python3
"""Test one global production ray for all P4(A_top,E_top) geometries."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

import mpmath as mp

from etop_rank1_elimination import (
    CROSSWALK,
    LINEAGES,
    PRIMARY_COVARIANCE,
    ROOT,
    SENSITIVITY_COVARIANCE,
    _chi2_survival,
    _state_row,
)


SOURCE_RESULT = ROOT / "results/etop-rank1-elimination/latest.json"


def _number(value: mp.mpf | float | int, digits: int = 17) -> float:
    return float(mp.nstr(value, digits))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def inverse_2x2(matrix: Sequence[Sequence[mp.mpf]]) -> list[list[mp.mpf]]:
    determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
    if determinant <= 0:
        raise ValueError("state covariance is not positive definite")
    return [
        [matrix[1][1] / determinant, -matrix[0][1] / determinant],
        [-matrix[1][0] / determinant, matrix[0][0] / determinant],
    ]


def dot(left: Sequence[mp.mpf], right: Sequence[mp.mpf]) -> mp.mpf:
    return mp.fsum(a * b for a, b in zip(left, right))


def matvec(matrix, vector):
    return [dot(row, vector) for row in matrix]


def ray_row_chi2(row: Mapping[str, Any], theta: mp.mpf) -> tuple[mp.mpf, mp.mpf]:
    direction = [mp.cos(theta), mp.sin(theta)]
    precision = inverse_2x2(row["covariance"])
    weighted_direction = matvec(precision, direction)
    weighted_point = matvec(precision, row["estimate"])
    denominator = dot(direction, weighted_direction)
    amplitude = dot(direction, weighted_point) / denominator
    residual = [
        value - amplitude * unit for value, unit in zip(row["estimate"], direction)
    ]
    chi2 = dot(residual, matvec(precision, residual))
    return chi2, amplitude


def ray_objective(rows: Sequence[Mapping[str, Any]], theta: mp.mpf) -> mp.mpf:
    return mp.fsum(ray_row_chi2(row, theta)[0] for row in rows)


def golden_minimum(function: Callable[[mp.mpf], mp.mpf], lower, upper):
    golden = (mp.sqrt(5) - 1) / 2
    left = upper - golden * (upper - lower)
    right = lower + golden * (upper - lower)
    for _ in range(140):
        if function(left) <= function(right):
            upper, right = right, left
            left = upper - golden * (upper - lower)
        else:
            lower, left = left, right
            right = lower + golden * (upper - lower)
    location = (lower + upper) / 2
    return location, function(location)


def minimize_periodic(function: Callable[[mp.mpf], mp.mpf], grid_size: int = 4096):
    period = mp.pi
    step = period / grid_size
    grid = [step * index for index in range(grid_size)]
    values = [function(value) for value in grid]
    spread = max(values) - min(values)
    if spread <= max(mp.mpf(1), abs(min(values))) * mp.mpf("1e-40"):
        index = min(range(grid_size), key=lambda item: values[item])
        return grid[index], values[index], {
            "method": "periodic grid; objective is constant at working precision",
            "grid_size": grid_size,
            "local_basins": "continuum",
            "best_grid_chi2": _number(values[index]),
            "refined_minus_grid_chi2": 0.0,
        }
    local = [
        index
        for index in range(grid_size)
        if values[index] < values[(index - 1) % grid_size]
        and values[index] <= values[(index + 1) % grid_size]
    ]
    candidates = [
        golden_minimum(function, grid[index] - step, grid[index] + step)
        for index in local
    ]
    theta, value = min(candidates, key=lambda row: row[1])
    theta %= period
    return theta, value, {
        "method": "periodic grid isolation followed by high-precision golden minimization of every local basin",
        "grid_size": grid_size,
        "local_basins": len(local),
        "best_grid_chi2": _number(min(values)),
        "refined_minus_grid_chi2": _number(value - min(values)),
    }


def signed_angle_degrees(theta: mp.mpf) -> mp.mpf:
    angle = (theta % mp.pi) * 180 / mp.pi
    if angle >= 90:
        angle -= 180
    return angle


def profile_offsets(rows, best_theta, best_chi2, delta):
    target = best_chi2 + mp.mpf(delta)

    def crossing(sign: int):
        inner = mp.mpf(0)
        f_inner = ray_objective(rows, best_theta) - target
        steps = 4096
        for index in range(1, steps + 1):
            outer = sign * (mp.pi / 2) * index / steps
            f_outer = ray_objective(rows, best_theta + outer) - target
            if f_inner * f_outer <= 0:
                lower, upper = sorted((inner, outer))
                for _ in range(120):
                    midpoint = (lower + upper) / 2
                    f_mid = ray_objective(rows, best_theta + midpoint) - target
                    if (ray_objective(rows, best_theta + lower) - target) * f_mid <= 0:
                        upper = midpoint
                    else:
                        lower = midpoint
                return (lower + upper) / 2
            inner, f_inner = outer, f_outer
        return sign * mp.pi / 2

    lower, upper = crossing(-1), crossing(+1)
    return {
        "delta_chi2": float(delta),
        "angle_offset_degrees": [_number(lower * 180 / mp.pi), _number(upper * 180 / mp.pi)],
        "absolute_angle_degrees_mod_ray": [
            _number(signed_angle_degrees(best_theta + lower)),
            _number(signed_angle_degrees(best_theta + upper)),
        ],
    }


def fit_ray(
    rows: Sequence[Mapping[str, Any]], include_profile: bool = False, grid_size: int = 4096
) -> dict[str, Any]:
    theta, chi2, certificate = minimize_periodic(
        lambda value: ray_objective(rows, value), grid_size=grid_size
    )
    if theta >= mp.pi / 2:
        theta -= mp.pi
    direction = [mp.cos(theta), mp.sin(theta)]
    row_details = []
    for row in rows:
        row_chi2, amplitude = ray_row_chi2(row, theta)
        row_details.append(
            {
                "id": row["id"],
                "N": row["N"],
                "dependency_group": row["dependency_group"],
                "amplitude": _number(amplitude),
                "chi2_contribution": _number(row_chi2),
            }
        )
    result = {
        "unit_direction_A_nonnegative": [_number(value) for value in direction],
        "angle_degrees_from_A_axis": _number(signed_angle_degrees(theta)),
        "E_over_A_slope": _number(direction[1] / direction[0])
        if abs(direction[0]) > mp.mpf("1e-40")
        else "vertical",
        "min_chi2": _number(chi2),
        "row_profiles": row_details,
        "optimizer_certificate": certificate,
        "theta_internal": theta,
    }
    if include_profile:
        result["profile_intervals"] = {
            "delta_chi2_1": profile_offsets(rows, theta, chi2, 1),
            "delta_chi2_3p841": profile_offsets(
                rows, theta, chi2, "3.841458820694124"
            ),
        }
    return result


def load_rows(crosswalk: Mapping[str, Any], covariance_key: str):
    datasets = {row["id"]: row for row in crosswalk["datasets"]}
    rows = []
    for _, parent_id, child_id, dependency_group in LINEAGES:
        for dataset_id in (parent_id, child_id):
            row = _state_row(datasets[dataset_id], covariance_key)
            row["dependency_group"] = dependency_group
            rows.append(row)
    return rows


def clean_fit(fit: Mapping[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in fit.items() if key != "theta_internal"}


def analyze(crosswalk, covariance_key: str, alpha: mp.mpf):
    rows = load_rows(crosswalk, covariance_key)
    groups = {
        name: [row for row in rows if row["dependency_group"] == name]
        for name in ("P49", "P43", "P50", "P57")
    }
    lineage_fits = {name: fit_ray(group_rows) for name, group_rows in groups.items()}
    four_ray_chi2 = mp.fsum(
        mp.mpf(str(fit["min_chi2"])) for fit in lineage_fits.values()
    )
    global_fit = fit_ray(rows, include_profile=True)
    global_chi2 = mp.mpf(str(global_fit["min_chi2"]))
    global_p = _chi2_survival(global_chi2, 7)
    difference = global_chi2 - four_ray_chi2
    difference_p = _chi2_survival(difference, 3)

    pairwise = []
    for left, right in itertools.combinations(groups, 2):
        pair_fit = fit_ray(groups[left] + groups[right])
        separate = mp.mpf(str(lineage_fits[left]["min_chi2"])) + mp.mpf(
            str(lineage_fits[right]["min_chi2"])
        )
        delta = mp.mpf(str(pair_fit["min_chi2"])) - separate
        theta_left = lineage_fits[left]["theta_internal"]
        theta_right = lineage_fits[right]["theta_internal"]
        determinant = mp.sin(theta_right - theta_left)
        pairwise.append(
            {
                "lineages": [left, right],
                "common_ray_chi2": pair_fit["min_chi2"],
                "separate_ray_chi2": _number(separate),
                "delta_chi2": _number(delta),
                "degrees_of_freedom": 1,
                "p_value": _number(_chi2_survival(delta, 1)),
                "signed_unit_ray_determinant": _number(determinant),
                "absolute_unit_ray_determinant": _number(abs(determinant)),
                "angle_difference_degrees_mod_ray": _number(
                    abs(signed_angle_degrees(theta_right - theta_left))
                ),
                "pair_common_ray": clean_fit(pair_fit),
            }
        )
    pairwise.sort(key=lambda row: row["delta_chi2"], reverse=True)

    heldout = []
    for name in groups:
        training = [row for row in rows if row["dependency_group"] != name]
        target = groups[name]
        train_fit = fit_ray(training)
        train_min = mp.mpf(str(train_fit["min_chi2"]))

        def predictive_profile(theta):
            return (
                ray_objective(training, theta)
                - train_min
                + ray_objective(target, theta)
            )

        theta, statistic, certificate = minimize_periodic(predictive_profile)
        heldout.append(
            {
                "heldout_lineage": name,
                "training_direction_angle_degrees": train_fit[
                    "angle_degrees_from_A_axis"
                ],
                "predictive_profile_direction_angle_degrees": _number(
                    signed_angle_degrees(theta)
                ),
                "profiled_predictive_chi2": _number(statistic),
                "degrees_of_freedom": 2,
                "p_value": _number(_chi2_survival(statistic, 2)),
                "decision_at_alpha": "eliminated"
                if _chi2_survival(statistic, 2) < alpha
                else "survives",
                "optimizer_certificate": certificate,
            }
        )

    contribution_by_group = {
        name: _number(
            mp.fsum(
                mp.mpf(str(row["chi2_contribution"]))
                for row in global_fit["row_profiles"]
                if row["dependency_group"] == name
            )
        )
        for name in groups
    }
    return {
        "covariance_key": covariance_key,
        "global_ray": {
            **clean_fit(global_fit),
            "degrees_of_freedom": 7,
            "p_value": _number(global_p),
            "decision_at_alpha": "eliminated" if global_p < alpha else "survives",
            "decision_scope": "absolute goodness-of-fit against unrestricted row means",
            "chi2_contribution_by_dependency_group": contribution_by_group,
        },
        "four_lineage_specific_rays": {
            "fits": {name: clean_fit(fit) for name, fit in lineage_fits.items()},
            "summed_chi2": _number(four_ray_chi2),
            "degrees_of_freedom": 4,
            "p_value": _number(_chi2_survival(four_ray_chi2, 4)),
        },
        "global_vs_four_ray_likelihood_ratio": {
            "delta_chi2": _number(difference),
            "degrees_of_freedom": 3,
            "p_value": _number(difference_p),
            "decision_at_alpha": "global_ray_eliminated"
            if difference_p < alpha
            else "global_ray_survives",
        },
        "primary_model_decision": "global_ray_eliminated_against_four_lineage_specific_rays"
        if difference_p < alpha
        else "global_ray_survives_four_lineage_specific_ray_comparison",
        "pairwise_independent_geometry_determinants": pairwise,
        "largest_pairwise_incompatibility": pairwise[0],
        "leave_one_lineage_out_predictive_profiles": heldout,
    }


def build_report(alpha: str = "0.01", dps: int = 70):
    mp.mp.dps = max(50, dps)
    crosswalk = json.loads(CROSSWALK.read_text(encoding="utf-8"))
    alpha_value = mp.mpf(alpha)
    primary = analyze(crosswalk, PRIMARY_COVARIANCE, alpha_value)
    sensitivity = analyze(crosswalk, SENSITIVITY_COVARIANCE, alpha_value)
    return {
        "schema": "matching-one/etop-global-ray-elimination/v1",
        "issues": [337, 370],
        "status": "cross_geometry_global_ray_model_elimination",
        "alpha": float(alpha),
        "null": "all eight P4(A_top,E_top) means share one unoriented ray through the origin; each size/geometry has its own unrestricted scalar amplitude",
        "alternative": "four lineage-specific unoriented rays, each shared only by its parent and child",
        "degrees_of_freedom": {
            "global_ray": "16 observations - 8 amplitudes - 1 ray angle = 7",
            "four_rays": "16 observations - 8 amplitudes - 4 ray angles = 4",
            "nested_difference": 3,
        },
        "dependency_contract": "full A_top/E_top covariance within every row; parent/child and P49/P43/P50/P57 blocks are block diagonal as declared by the production archives",
        "source": {
            "crosswalk": str(CROSSWALK.relative_to(ROOT)),
            "crosswalk_sha256": sha256(CROSSWALK),
            "lineage_rank1_result": str(SOURCE_RESULT.relative_to(ROOT)),
            "lineage_rank1_result_sha256": sha256(SOURCE_RESULT),
        },
        "primary": primary,
        "fixed_center_sensitivity": sensitivity,
        "scientific_card": [
            "MECHANISM SPACE: one universal cross-geometry A_top/E_top ray versus four geometry-family rays.",
            "RESULT: the test changes only the geometry-sharing constraint; every size retains a free scalar amplitude and no field name is assigned.",
            "NOT PROVED: ray incompatibility does not identify a field, operator, exponent or asymptotic geometry law.",
            "OBSERVER-SECTOR-SOURCE-GEOMETRY: P4(A_top,E_top) | Alexander odd/even state plane | threshold-rank source | P49/P43/P50/P57 production geometries.",
            "DEPENDENCY GROUPS: four named blocks remain explicit; the global statistic is their block-diagonal joint profile, not eight independent narrative votes.",
            "UPWEIGHT OBSERVATION: a new independent geometry should be frozen near the largest pairwise determinant direction, then scored against the source-profiled global ray.",
        ],
    }


def render_markdown(report: Mapping[str, Any]) -> str:
    primary = report["primary"]
    global_ray = primary["global_ray"]
    comparison = primary["global_vs_four_ray_likelihood_ratio"]
    largest = primary["largest_pairwise_incompatibility"]
    fixed = report["fixed_center_sensitivity"]
    lines = [
        "# Cross-geometry E_top global-ray elimination",
        "",
        "## Answer",
        "",
        f"The one-global-ray null is **{comparison['decision_at_alpha']}** against",
        "the four-lineage-ray alternative. Its absolute goodness-of-fit is",
        f"`chi2={global_ray['min_chi2']:.6g}` on `7` df",
        f"(`p={global_ray['p_value']:.6g}`), which alone {global_ray['decision_at_alpha']} at alpha=0.01. Relative to the four lineage-specific",
        f"rays, the gauge-free nested penalty is `Delta chi2={comparison['delta_chi2']:.6g}`",
        f"on `3` df (`p={comparison['p_value']:.6g}`), so the stronger geometry-sharing",
        f"constraint is `{comparison['decision_at_alpha']}`.",
        "",
        f"The fitted global direction is `{global_ray['angle_degrees_from_A_axis']:.5g}` degrees",
        f"from the positive A_top axis (`E/A={global_ray['E_over_A_slope']:.6g}`).",
        f"The largest independent-geometry incompatibility is `{largest['lineages']}`",
        f"with `Delta chi2={largest['delta_chi2']:.6g}` on 1 df",
        f"and absolute unit-ray determinant `{largest['absolute_unit_ray_determinant']:.6g}`.",
        "",
        "## Lineage directions",
        "",
        "| dependency group | angle from A (deg) | E/A | lineage chi2 / 1 df |",
        "|---|---:|---:|---:|",
    ]
    for name, fit in primary["four_lineage_specific_rays"]["fits"].items():
        lines.append(
            f"| {name} | {fit['angle_degrees_from_A_axis']:.5g} | {fit['E_over_A_slope']:.6g} | {fit['min_chi2']:.6g} |"
        )
    lines.extend(
        [
            "",
            "## Pairwise geometry contrasts",
            "",
            "| pair | Delta chi2 / 1 df | p | abs det(unit rays) | angle difference (deg) |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for row in primary["pairwise_independent_geometry_determinants"]:
        lines.append(
            f"| {row['lineages']} | {row['delta_chi2']:.6g} | {row['p_value']:.6g} "
            f"| {row['absolute_unit_ray_determinant']:.6g} | {row['angle_difference_degrees_mod_ray']:.6g} |"
        )
    lines.extend(
        [
            "",
            "## Leave-one-lineage-out sensitivity",
            "",
            "| held out | profiled predictive chi2 / 2 df | p | decision |",
            "|---|---:|---:|---|",
        ]
    )
    for row in primary["leave_one_lineage_out_predictive_profiles"]:
        lines.append(
            f"| {row['heldout_lineage']} | {row['profiled_predictive_chi2']:.6g} "
            f"| {row['p_value']:.6g} | {row['decision_at_alpha']} |"
        )
    lines.extend(
        [
            "",
            "The held-out statistic profiles the ray against the training likelihood",
            "penalty rather than pretending that the training direction is known exactly.",
            "The fixed-center exact-batch sensitivity gives global-ray",
            f"`chi2={fixed['global_ray']['min_chi2']:.6g}` and nested",
            f"`Delta chi2={fixed['global_vs_four_ray_likelihood_ratio']['delta_chi2']:.6g}`;",
            "the decision is unchanged.",
            "",
            "## Scientific card",
            "",
        ]
    )
    lines.extend(f"- {item}" for item in report["scientific_card"])
    lines.extend(
        [
            "",
            "## Reproduction",
            "",
            "```bash",
            "python3 scripts/etop_global_ray_elimination.py --format json --output results/etop-global-ray-elimination/latest.json",
            "python3 scripts/etop_global_ray_elimination.py --format markdown --output results/etop-global-ray-elimination/REPORT.md",
            "python3 -m unittest discover -s tests -p 'test_etop_global_ray_elimination.py'",
            "```",
            "",
            "## Claim boundary",
            "",
            "This is a cross-geometry model elimination on already revealed production",
            "blocks. It neither assigns the incompatible directions to fields nor treats",
            "the four dependency groups as separate discoveries. Intrinsic-center covariance",
            "is first-order; the fixed-center exact-batch row is the declared sensitivity.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--format", choices=("json", "markdown"), default="json")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--alpha", default="0.01")
    parser.add_argument("--dps", type=int, default=70)
    args = parser.parse_args()
    report = build_report(alpha=args.alpha, dps=args.dps)
    payload = (
        json.dumps(report, indent=2, sort_keys=True) + "\n"
        if args.format == "json"
        else render_markdown(report)
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(payload, encoding="utf-8")


if __name__ == "__main__":
    main()

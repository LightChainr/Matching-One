#!/usr/bin/env python3
"""Freeze multi-u thermal-response templates from the #101 coordinate map.

Issue #119 asks for the shape of a correction sector across the already
frozen intrinsic levels rather than a new derivative at the centre.  The
nonlinear bare-to-thermal reparametrization of issue #101 predicts, at
leading order and without a fit,

    w_u / u           ~  B N^{-3/8}
    (c_u - c_0) / u^2 ~  A N^{-3/4}

with the frozen grid `u ∈ {0, 0.025, 0.05}`.  Ordinary analytic q=2
corrections add higher odd/even powers of u; a rank-2 Jordan factor
multiplies the same scaling functions by log N and therefore changes the
N-power of the even piece.

This module freezes those monomials before any P43 / N185 / N265 or
Issue #57 coordinate is read.  P49 N=130/170 reconstructions are
development-only diagnostics of the template, not a held-out score.
N=10 Beta(3,3) is the exact odd-around-1/2 oracle: every midpoint is
1/2, so the even u-shape vanishes identically.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from pathlib import Path
from typing import Dict, Mapping, Sequence, Tuple

ROOT = Path(__file__).resolve().parents[1]
P49_RAW = ROOT / "results" / "server-20260828" / "P49-fullcurve-doubling-100m" / "raw"

FROZEN_U: Tuple[float, ...] = (0.0, 0.025, 0.05)
EVEN_U_POWER = 2
ODD_U_POWER = 1
EVEN_N_EXPONENT = Fraction(-3, 4)
ODD_N_EXPONENT = Fraction(-3, 8)
DOUBLING_EVEN = Fraction(1, 2) ** Fraction(3, 4)
DOUBLING_ODD = Fraction(1, 2) ** Fraction(3, 8)


def n10_matching(p: float) -> float:
    return 12.0 * p**5 - 30.0 * p**4 + 20.0 * p**3 - 1.0


def solve_matching(function, target: float, lower: float = 0.1, upper: float = 0.9) -> float:
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


def center_width(function, u: float) -> Dict[str, float]:
    if u == 0.0:
        root = solve_matching(function, 0.0)
        return {"u": 0.0, "p_minus": root, "p_plus": root, "c": root, "w": 0.0}
    p_minus = solve_matching(function, -u)
    p_plus = solve_matching(function, u)
    return {
        "u": u,
        "p_minus": p_minus,
        "p_plus": p_plus,
        "c": (p_plus + p_minus) / 2.0,
        "w": (p_plus - p_minus) / 2.0,
    }


def shape_row(levels: Mapping[float, Mapping[str, float]]) -> Dict[str, float]:
    c0 = levels[0.0]["c"]
    odd_ratios = []
    even_ratios = []
    for u in FROZEN_U:
        if u == 0.0:
            continue
        odd_ratios.append(levels[u]["w"] / (u ** ODD_U_POWER))
        even_ratios.append((levels[u]["c"] - c0) / (u ** EVEN_U_POWER))
    return {
        "w_over_u_0.025": odd_ratios[0],
        "w_over_u_0.05": odd_ratios[1],
        "c_shift_over_u2_0.025": even_ratios[0],
        "c_shift_over_u2_0.05": even_ratios[1],
        "odd_u_ratio": odd_ratios[0] / odd_ratios[1] if odd_ratios[1] else float("nan"),
        "even_u_ratio": (
            even_ratios[0] / even_ratios[1] if abs(even_ratios[1]) > 1e-18 else 1.0
        ),
    }


def n10_oracle() -> Dict[str, object]:
    levels = {u: center_width(n10_matching, u) for u in FROZEN_U}
    shape = shape_row(levels)
    return {
        "N": 10,
        "oracle": "Beta(3,3)",
        "levels": levels,
        "shape": shape,
        "even_shape_vanishes": all(abs(levels[u]["c"] - 0.5) < 1e-15 for u in FROZEN_U),
        "Q": levels[0.05]["c"] - levels[0.025]["c"],
    }


def p49_descriptive(n: int) -> Dict[str, object]:
    from analyze_p48_retrospective import read_histograms
    from score_p49_fullcurve_doubling import aggregate, orientation_values, solve_target

    path = P49_RAW / f"n{n}.hist.csv"
    records = read_histograms(path)
    by_orientation = {
        name: sorted((row for key, row in records.items() if key[1] == name), key=lambda row: row.batch)
        for name in ("first", "second")
    }
    rows = {name: aggregate(by_orientation[name], -1) for name in ("first", "second")}

    def mean_matching(p: float) -> float:
        return (
            orientation_values(n, rows["first"], p)["M"]
            + orientation_values(n, rows["second"], p)["M"]
        ) / 2.0

    levels = {}
    for u in FROZEN_U:
        if u == 0.0:
            root = solve_target(mean_matching, 0.0)
            levels[u] = {"u": 0.0, "p_minus": root, "p_plus": root, "c": root, "w": 0.0}
        else:
            p_minus = solve_target(mean_matching, -u)
            p_plus = solve_target(mean_matching, u)
            levels[u] = {
                "u": u,
                "p_minus": p_minus,
                "p_plus": p_plus,
                "c": (p_plus + p_minus) / 2.0,
                "w": (p_plus - p_minus) / 2.0,
            }
    shape = shape_row(levels)
    odd_scaled = {
        str(u): levels[u]["w"] * (n ** float(-ODD_N_EXPONENT)) / u
        for u in FROZEN_U
        if u
    }
    even_scaled = {
        str(u): (levels[u]["c"] - levels[0.0]["c"]) * (n ** float(-EVEN_N_EXPONENT)) / (u * u)
        for u in FROZEN_U
        if u
    }
    return {
        "N": n,
        "role": "development_descriptive_only",
        "not_a_target_for_p43_or_issue_57": True,
        "not_a_doubling_pair_with_the_other_p49_child": True,
        "levels": levels,
        "shape": shape,
        "odd_scaled_w_over_u": odd_scaled,
        "even_scaled_c_shift_over_u2": even_scaled,
        "frozen_u": list(FROZEN_U),
    }


def frozen_templates() -> Dict[str, object]:
    return {
        "u_grid": list(FROZEN_U),
        "do_not_add_levels_after_looking": True,
        "coordinate_nonlinearity": {
            "odd": "w_u / u ~ N^{-3/8}",
            "even": "(c_u - c_0) / u^2 ~ N^{-3/4}",
            "doubling_odd": float(DOUBLING_ODD),
            "doubling_even": float(DOUBLING_EVEN),
            "source": "issue_101_bare_to_thermal_reparametrization",
        },
        "q2_analytic": {
            "odd": "polynomial in odd powers of u",
            "even": "polynomial in even powers of u",
            "leading_agrees_with_coordinate_nonlinearity": True,
        },
        "jordan_log": {
            "multiplies_the_same_u_shape_by_log_N": True,
            "changes_the_N_power_of_the_even_piece": True,
        },
        "joint_covariance": "recompute_p_plus_minus_inside_each_delete_one_replicate",
        "not_a_target_for_p43_or_issue_57": True,
    }


def run_suite() -> Dict[str, object]:
    n10 = n10_oracle()
    p49 = [p49_descriptive(130), p49_descriptive(170)]
    return {
        "schema": "multi-u thermal response templates v1",
        "templates": frozen_templates(),
        "n10_oracle": n10,
        "p49_descriptive": p49,
        "not_a_target_for_p43_or_issue_57": True,
    }


def render_report(payload: Dict[str, object]) -> str:
    n10 = payload["n10_oracle"]
    lines = [
        "# Multi-u thermal-response templates",
        "",
        "Source: `scripts/multi_u_thermal_response.py`.",
        "Claim level: C0 template freeze, C1 N=10 oracle. P49 rows are development-only.",
        "Not a P43 / Issue #57 target. Frozen `u={0, 0.025, 0.05}`; do not add levels.",
        "",
        "## Frozen monomials",
        "",
        "```text",
        "w_u / u           ~ B N^{-3/8}",
        "(c_u - c_0) / u^2 ~ A N^{-3/4}",
        "```",
        "",
        "These are the leading #101 coordinate-nonlinearity shapes. Ordinary q=2",
        "analytic corrections add higher powers of u. A Jordan log multiplies the",
        "same u-shape by log N. Templates must be used jointly with delete-one",
        "covariance; they are not three independent tests.",
        "",
        "## N=10 Beta(3,3) oracle",
        "",
        f"Even midpoints vanish: `Q_10={n10['Q']}`, every `c_u=1/2`.",
        "",
        "## Descriptive P49 N=130/170",
        "",
        "These sizes are children of different doubling lineages, not a doubling pair.",
        "",
        "| N | `w/u` at 0.025 | `w/u` at 0.05 | `(c-c0)/u^2` at 0.025 | `(c-c0)/u^2` at 0.05 |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in payload["p49_descriptive"]:
        shape = row["shape"]
        lines.append(
            "| {N} | {w025:.8f} | {w05:.8f} | {e025:.6f} | {e05:.6f} |".format(
                N=row["N"],
                w025=shape["w_over_u_0.025"],
                w05=shape["w_over_u_0.05"],
                e025=shape["c_shift_over_u2_0.025"],
                e05=shape["c_shift_over_u2_0.05"],
            )
        )
    lines.extend(
        [
            "",
            "On both sizes `w_{0.025}/w_{0.05}≈1/2` and the two even ratios agree to",
            "relative O(10^{-3}). That is compatible with the frozen monomials. It is",
            "not a prospective Jordan discriminator and is not a P43 score.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    payload = run_suite()
    n10 = payload["n10_oracle"]
    print("frozen u", FROZEN_U)
    print("N=10 even shape vanishes", n10["even_shape_vanishes"], "Q", n10["Q"])
    for row in payload["p49_descriptive"]:
        shape = row["shape"]
        print(
            "P49 N={N} w/u={w:.6f}/{w2:.6f} even/u^2={e:.6f}/{e2:.6f} labeled {role}".format(
                N=row["N"],
                w=shape["w_over_u_0.025"],
                w2=shape["w_over_u_0.05"],
                e=shape["c_shift_over_u2_0.025"],
                e2=shape["c_shift_over_u2_0.05"],
                role=row["role"],
            )
        )
    if args.report is not None:
        args.report.write_text(render_report(payload), encoding="utf-8")
        print("wrote " + str(args.report))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

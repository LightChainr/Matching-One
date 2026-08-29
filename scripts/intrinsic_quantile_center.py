#!/usr/bin/env python3
"""Intrinsic quantile-center spectroscopy of the nonlinear thermal field.

Issue #101. Frozen before any P43 coordinate is read as a target:

    Mbar_N(p_-^u) = -u,   Mbar_N(p_+^u) = +u
    c_u = (p_+^u + p_-^u)/2
    w_u = (p_+^u - p_-^u)/2
    Q_N = c_{0.05} - c_{0.025}  ~  C N^{-3/4}

The level set is exactly ``u = {0.025, 0.05}``. Do not add levels after looking.
The no-fit doubling prediction on a true doubling lineage is

    Q_{2N} / Q_N = 2^{-3/4}.

N=10 C4 self-matching is the exact oracle: M(p)=2 I_p(3,3)-1 is odd about
1/2, so every midpoint is 1/2 and Q_10=0. Clean P49 N=130/170 reconstructions
are development-only and are not a P43 or Issue #57 target.
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Callable, Mapping, Sequence


ROOT = Path(__file__).resolve().parents[1]
P49_RAW = ROOT / "results" / "server-20260828" / "P49-fullcurve-doubling-100m" / "raw"
P49_HISTOGRAMS = {
    130: P49_RAW / "n130.hist.csv",
    170: P49_RAW / "n170.hist.csv",
}

FROZEN_U = (0.025, 0.05)
Q_LEVELS = (0.025, 0.05)
Q_EXPONENT = Fraction(-3, 4)
DOUBLING_RATIO = 2 ** float(Q_EXPONENT)
WIDTH_EXPONENT = Fraction(3, 8)
N10_DEGREE = 10


def n10_matching(p: float) -> float:
    """Exact N=10 C4 self-matching polynomial M(p)=2 I_p(3,3)-1."""

    return (((12.0 * p - 30.0) * p + 20.0) * p * p * p) - 1.0


def solve_level(function: Callable[[float], float], target: float) -> float:
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


@dataclass(frozen=True)
class QuantileLevel:
    u: float
    p_minus: float
    p_plus: float

    @property
    def c(self) -> float:
        return 0.5 * (self.p_plus + self.p_minus)

    @property
    def w(self) -> float:
        return 0.5 * (self.p_plus - self.p_minus)


def quantile_levels(
    matching: Callable[[float], float], levels: Sequence[float] = FROZEN_U
) -> dict[float, QuantileLevel]:
    if tuple(levels) != FROZEN_U:
        raise ValueError("issue #101 freezes u={0.025, 0.05}; do not add levels")
    output = {}
    for u in levels:
        output[u] = QuantileLevel(
            u=u,
            p_minus=solve_level(matching, -u),
            p_plus=solve_level(matching, u),
        )
    return output


def midpoint_difference(levels: Mapping[float, QuantileLevel]) -> float:
    return levels[0.05].c - levels[0.025].c


def n10_oracle() -> dict[str, object]:
    levels = quantile_levels(n10_matching)
    q_n = midpoint_difference(levels)
    return {
        "N": N10_DEGREE,
        "oracle": "M(p)=2 I_p(3,3)-1",
        "self_matching": True,
        "levels": {
            str(u): {
                "p_minus": level.p_minus,
                "p_plus": level.p_plus,
                "c": level.c,
                "w": level.w,
            }
            for u, level in levels.items()
        },
        "Q": q_n,
        "Q_vanishes_by_oddness": abs(q_n) < 1e-15,
        "c_u_equals_half": all(abs(level.c - 0.5) < 1e-15 for level in levels.values()),
        "not_a_numerical_target_for_p43_or_issue_57": True,
    }


def frozen_monomials() -> dict[str, object]:
    return {
        "u": list(FROZEN_U),
        "Q_N": "c_0.05 - c_0.025",
        "Q_N_n_exponent": str(Q_EXPONENT),
        "doubling_ratio": "2^{-3/4}",
        "doubling_ratio_float": DOUBLING_RATIO,
        "w_u_n_exponent": str(WIDTH_EXPONENT),
        "do_not_add_levels_after_looking": True,
        "not_a_numerical_target_for_p43_or_issue_57": True,
        "true_doubling_lineages": ["65->130", "85->170", "145->290"],
        "p49_130_170_are_not_a_doubling_pair": True,
    }


def reconstruct_histogram(path: Path) -> dict[str, object]:
    from analyze_p48_retrospective import read_histograms
    from score_p49_fullcurve_doubling import aggregate, orientation_values

    records = read_histograms(path)
    sizes = {key[0] for key in records}
    if len(sizes) != 1:
        raise ValueError(f"{path} must contain a single N, got {sorted(sizes)}")
    n = next(iter(sizes))
    by_orientation = {
        orientation: sorted(
            (row for key, row in records.items() if key[1] == orientation),
            key=lambda row: row.batch,
        )
        for orientation in ("first", "second")
    }
    rows = {
        name: aggregate(by_orientation[name], -1) for name in ("first", "second")
    }

    def mean_matching(p: float) -> float:
        first = orientation_values(n, rows["first"], p)["M"]
        second = orientation_values(n, rows["second"], p)["M"]
        return 0.5 * (first + second)

    levels = quantile_levels(mean_matching)
    q_n = midpoint_difference(levels)
    first = by_orientation["first"][0]
    second = by_orientation["second"][0]
    return {
        "N": n,
        "path": str(path),
        "first": [first.a, first.b],
        "second": [second.a, second.b],
        "batches": len(by_orientation["first"]),
        "samples": sum(row.samples for row in by_orientation["first"]),
        "channel": "rank-2 cross",
        "levels": {
            str(u): {
                "p_minus": level.p_minus,
                "p_plus": level.p_plus,
                "c": level.c,
                "w": level.w,
                "M_minus": mean_matching(level.p_minus),
                "M_plus": mean_matching(level.p_plus),
            }
            for u, level in levels.items()
        },
        "Q": q_n,
        "Q_times_N_to_3_over_4": q_n * (n ** 0.75),
        "w_0.025_times_N_to_3_over_8": levels[0.025].w * (n ** 0.375),
        "w_0.05_times_N_to_3_over_8": levels[0.05].w * (n ** 0.375),
        "role": "development_descriptive_only",
        "not_a_numerical_target_for_p43_or_issue_57": True,
        "not_a_doubling_pair_with_the_other_p49_child": True,
    }


def descriptive_p49() -> dict[str, object]:
    rows = {
        str(n): reconstruct_histogram(path) for n, path in P49_HISTOGRAMS.items()
    }
    return {
        "schema": "issue-101 intrinsic quantile center v1",
        "claim_level": "C0 definition freeze / C1 N=10 oracle",
        "frozen": frozen_monomials(),
        "n10_oracle": n10_oracle(),
        "p49_descriptive": rows,
        "p49_role": "development_descriptive_only",
        "not_a_numerical_target_for_p43_or_issue_57": True,
        "p43_was_not_read_as_a_target": True,
    }


def render_report(data: dict[str, object]) -> str:
    n10 = data["n10_oracle"]
    frozen = data["frozen"]
    lines = [
        "# Intrinsic quantile-center spectroscopy",
        "",
        "Source: `scripts/intrinsic_quantile_center.py`.",
        "Claim level: C0 definition freeze, C1 N=10 oracle. Not a P43/#57 target.",
        "",
        "Frozen levels `u={0.025, 0.05}` and",
        "",
        "```text",
        "Q_N = c_{0.05} - c_{0.025}  ~  C N^{-3/4}",
        "Q_{2N}/Q_N = 2^{-3/4}     on a true doubling lineage",
        "```",
        "",
        f"Numeric doubling ratio `{frozen['doubling_ratio_float']:.16f}`.",
        "Do not add quantile levels after looking at outcomes.",
        "",
        "## N=10 Beta(3,3) oracle",
        "",
        "The C4 self-matching control has `M(p)=2 I_p(3,3)-1`, which is odd about",
        "`p=1/2`. Every intrinsic midpoint is therefore `1/2` and `Q_10=0` exactly.",
        "",
        "```text",
        f"Q_10 = {n10['Q']}",
        "c_u = 1/2 for both frozen u",
        "```",
        "",
        "## Descriptive P49 N=130/170",
        "",
        "These two sizes are children of different doubling lineages, not a",
        "doubling pair. Scaled `Q_N N^{3/4}` is reported as a development",
        "diagnostic only. Do not score P43 against these numbers.",
        "",
        "| N | `c_{0.025}` | `c_{0.05}` | `Q_N` | `Q_N N^{3/4}` | `w_{0.025} N^{3/8}` | `w_{0.05} N^{3/8}` |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for n in ("130", "170"):
        row = data["p49_descriptive"][n]
        lines.append(
            "| {N} | {c025:.10f} | {c05:.10f} | {Q:.6e} | {scaled:.6e} | {w025:.6f} | {w05:.6f} |".format(
                N=row["N"],
                c025=row["levels"]["0.025"]["c"],
                c05=row["levels"]["0.05"]["c"],
                Q=row["Q"],
                scaled=row["Q_times_N_to_3_over_4"],
                w025=row["w_0.025_times_N_to_3_over_8"],
                w05=row["w_0.05_times_N_to_3_over_8"],
            )
        )
    lines.extend(
        [
            "",
            "## Boundary",
            "",
            "- Retrospective P49 numbers are development only.",
            "- A claim-bearing score must recompute `p_±^u` inside each delete-one",
            "  replicate and must be frozen before the target coordinates are read.",
            "- P43 N=185/265 and Issue #57 norm-5 are not targets of this freeze.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path)
    parser.add_argument("--report", type=Path)
    parser.add_argument("--skip-p49", action="store_true")
    args = parser.parse_args()
    if args.skip_p49:
        data = {
            "frozen": frozen_monomials(),
            "n10_oracle": n10_oracle(),
            "p49_descriptive": {},
            "not_a_numerical_target_for_p43_or_issue_57": True,
        }
    else:
        data = descriptive_p49()
    print("frozen u = {0.025, 0.05}")
    print("Q_N ~ N^{-3/4}; doubling ratio 2^{-3/4} = " + f"{DOUBLING_RATIO:.16f}")
    oracle = data["n10_oracle"]
    print("N=10 oracle Q = " + str(oracle["Q"]) + " (vanishes by oddness)")
    for n, row in data.get("p49_descriptive", {}).items():
        print(
            "P49 descriptive N={N}: Q={Q:.6e} Q N^{{3/4}}={scaled:.6e}".format(
                N=row["N"], Q=row["Q"], scaled=row["Q_times_N_to_3_over_4"]
            )
        )
        print("  labeled development-only; not a P43/#57 target")
    if args.json is not None:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
        print("wrote " + str(args.json))
    if args.report is not None:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(render_report(data), encoding="utf-8")
        print("wrote " + str(args.report))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

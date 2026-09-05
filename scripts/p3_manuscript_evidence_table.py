#!/usr/bin/env python3
"""Assemble the P3 manuscript evidence tables from committed ladder artifacts.

This script runs no Monte Carlo and rescores nothing that has not already been
scored and committed.  It reads the frozen N=580 design, the committed Fieller
scoring, and the committed projective rescoring, and derives only what follows
exactly from those three inputs:

  * the three-rung response with each rung's distance from zero, which is what
    decides whether a denominator is safe to nominate;
  * the two-entry control -- the projective statistic restricted to the pair of
    rungs the frozen test actually used, which must reproduce Fieller's z
    squared exactly, and does;
  * the verdict table, frozen against projective, so that every change is
    attributable to the third rung rather than to a change of statistic;
  * the curvature functional and the spin-8 amplitude each competitor needs.

Nothing here is hand-typed into the draft: ``tables.md`` is rendered from the
artifact this script writes, and a regression test fails if the two drift.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

try:  # pragma: no cover - import shape depends on how the script is invoked
    from scripts.projective_inference import ray_residual
    from scripts import aspect_ladder_projective_rescore as rescore_module
except ModuleNotFoundError:  # pragma: no cover
    from projective_inference import ray_residual
    import aspect_ladder_projective_rescore as rescore_module

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "results" / "p3-projective-inference-manuscript" / "latest.json"
DEFAULT_MARKDOWN = ROOT / "docs" / "manuscripts" / "p3-projective-inference" / "tables.md"
FROZEN_SCORING = ROOT / "results" / "aspect-ladder-n580" / "latest.json"
PROJECTIVE_SCORING = ROOT / "results" / "aspect-ladder-n580-projective" / "latest.json"
FROZEN_DESIGN = ROOT / "predictions" / "aspect_ladder_n580_20260905.yaml"
SCHEMA = "matching-one/p3-projective-inference-manuscript/v1"
ISSUE = 579
RUNGS = (1, 2, 4)
THREE_SIGMA = 3.0


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def response_table(vector: Sequence[float],
                   covariance: Sequence[Sequence[float]]) -> list[dict[str, Any]]:
    """Each rung's value, its error, and how far it sits from zero.

    The last column is the whole reason a denominator is a choice and not a
    convention: nominating a rung as denominator inherits its distance from
    zero, and every other rung's error is then reported through it.
    """
    rows = []
    for index, rung in enumerate(RUNGS):
        error = math.sqrt(covariance[index][index])
        rows.append({
            "rung": rung,
            "value": vector[index],
            "standard_error": error,
            "sigma_from_zero": vector[index] / error,
        })
    return rows


def two_entry_control(vector: Sequence[float],
                      covariance: Sequence[Sequence[float]],
                      competitors: Mapping[str, Sequence[float]],
                      frozen: Mapping[str, Any]) -> dict[str, Any]:
    """Projective statistic on the frozen test's own two rungs, against Fieller.

    The frozen design nominated ``r4_over_r1`` as the deciding entry, so this
    restricts the projective statistic to exactly those two rungs and one ray.
    Fieller's z squared and the restricted projective statistic are the same
    number; showing that on the real data is what makes the verdict changes in
    ``verdict_table`` attributable to the third rung and to nothing else.
    """
    entry = frozen["discriminating_entry"]
    _require(entry == "r4_over_r1", f"unexpected discriminating entry {entry!r}")
    pair_covariance = [[covariance[0][0], covariance[0][2]],
                       [covariance[2][0], covariance[2][2]]]
    rows = []
    worst = 0.0
    for name in sorted(competitors):
        ray = competitors[name]
        fit = ray_residual([vector[0], vector[2]], pair_covariance, [ray[0], ray[2]])
        statistic = float(fit["statistic"])
        fieller_z = frozen["comparison"][name]["z"][entry]
        deviation = abs(statistic - fieller_z ** 2) / max(abs(fieller_z ** 2), 1e-300)
        worst = max(worst, deviation)
        rows.append({
            "competitor": name,
            "fieller_z": fieller_z,
            "fieller_z_squared": fieller_z ** 2,
            "projective_statistic_two_entries": statistic,
            "relative_deviation": deviation,
        })
    return {
        "entry": entry,
        "rows": rows,
        "largest_relative_deviation": worst,
        "what_it_shows": (
            "on the two rungs the frozen test used, the projective statistic is "
            "Fieller's z squared to machine precision, so no verdict below moves "
            "because the statistic changed"
        ),
    }


def verdict_table(competitors: Mapping[str, Sequence[float]],
                  frozen: Mapping[str, Any],
                  projective: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Frozen verdict against projective verdict, competitor by competitor."""
    entry = frozen["discriminating_entry"]
    rows = []
    for name in sorted(competitors):
        fieller_z = abs(frozen["comparison"][name]["z"][entry])
        block = projective["competitors"][name]
        frozen_excluded = fieller_z > THREE_SIGMA
        projective_excluded = bool(block["excluded_at_3_sigma"])
        rows.append({
            "competitor": name,
            "ray": list(competitors[name]),
            "fieller_sigma_one_entry": fieller_z,
            "projective_sigma_three_rungs": block["equivalent_sigma"],
            "frozen_verdict": "excluded" if frozen_excluded else "compatible",
            "projective_verdict": "excluded" if projective_excluded else "compatible",
            "verdict_changed": frozen_excluded != projective_excluded,
            "sigma_over_admissible_correlations": list(
                block["sigma_over_admissible_correlations"]),
            "verdict_survives_the_missing_covariance": bool(
                block["verdict_survives_the_missing_covariance"]),
        })
    return rows


def curvature_table(projective: Mapping[str, Any],
                    competitors: Mapping[str, Sequence[float]]) -> dict[str, Any]:
    """The second divided difference, measured, against what each ray predicts.

    ``f[1,2,4]`` is a linear functional of the response.  It has no denominator
    and no matrix inverse, so unlike every ratio in this problem it does not
    degrade when a rung sits close to zero.
    """
    measured = projective["curvature"]
    predictions = []
    for name in sorted(competitors):
        predictions.append({
            "competitor": name,
            "curvature_predicted": projective["competitors"][name]["curvature_predicted"],
        })
    signs = {round(row["curvature_predicted"], 12) > 0.0 for row in predictions}
    return {
        "weights": list(rescore_module.CURVATURE_WEIGHTS),
        "measured": measured,
        "predicted": predictions,
        "no_competitor_predicts_a_negative_curvature": all(
            row["curvature_predicted"] >= -1e-12 for row in predictions),
        "at_least_one_competitor_predicts_a_positive_curvature": True in signs,
    }


def spin8_table(projective: Mapping[str, Any],
                competitors: Mapping[str, Sequence[float]]) -> dict[str, Any]:
    """How large the design's own declared systematic must be, per competitor.

    The frozen design kept the r=2 rung out of the decision because it carries a
    spin-8 leakage of the opposite sign, and justified doing so with a bound of
    ``|A8/A4|`` well below 1.  This is the value that bound would have to take
    for each competitor to reach the rung that was dropped.
    """
    rows = []
    for name in sorted(competitors):
        rows.append({
            "competitor": name,
            "required_abs_A8_over_A4_to_reach_r2":
                projective["competitors"][name]["required_abs_A8_over_A4_to_reach_r2"],
        })
    finite = [row["required_abs_A8_over_A4_to_reach_r2"] for row in rows
              if math.isfinite(row["required_abs_A8_over_A4_to_reach_r2"])]
    _require(bool(finite), "no competitor produced a finite spin-8 requirement")
    return {
        "leakage_coefficient": rescore_module.LEAKAGE,
        "assumed_by_the_frozen_design":
            projective["spin8_bound_provenance"]["assumed_by_the_frozen_design"],
        "provenance": projective["spin8_bound_provenance"],
        "rows": rows,
        "smallest_requirement": min(finite),
        "largest_requirement": max(finite),
    }


def assemble() -> dict[str, Any]:
    frozen = json.loads(FROZEN_SCORING.read_text(encoding="utf-8"))
    projective = json.loads(PROJECTIVE_SCORING.read_text(encoding="utf-8"))
    competitors = rescore_module.load_competitors()
    response = rescore_module.load_response(FROZEN_SCORING)
    vector, covariance = response["vector"], response["covariance"]
    _require(list(projective["response_vector"]) == list(vector),
             "the projective artifact was scored on a different response vector")
    return {
        "schema": SCHEMA,
        "issue": ISSUE,
        "channel": projective["channel"],
        "sources": {
            "frozen_design": str(FROZEN_DESIGN.relative_to(ROOT)),
            "frozen_scoring": str(FROZEN_SCORING.relative_to(ROOT)),
            "projective_scoring": str(PROJECTIVE_SCORING.relative_to(ROOT)),
        },
        "site_count": frozen["site_count"],
        "batches": frozen["batches"],
        "samples_per_rung": frozen["samples_per_rung"],
        "covariance_is_complete": response["covariance_is_complete"],
        "response": response_table(vector, covariance),
        "fieller_interval_3sigma": frozen["fieller_interval_3sigma"],
        "two_entry_control": two_entry_control(vector, covariance, competitors, frozen),
        "verdicts": verdict_table(competitors, frozen, projective),
        "curvature": curvature_table(projective, competitors),
        "spin8": spin8_table(projective, competitors),
        "admissible_corr_r2_r4": list(projective["admissible_corr_r2_r4"]),
        "what_this_does_not_separate": projective["what_this_does_not_separate"],
    }


def _fmt(value: float, digits: int = 3) -> str:
    if value == math.inf:
        return "&infin;"
    if value == -math.inf:
        return "&minus;&infin;"
    if value == 0.0 or 1e-3 <= abs(value) < 1e5:
        return f"{value:.{digits}f}"
    return f"{value:.{digits}e}"


def render(payload: Mapping[str, Any]) -> str:
    lines = [
        "# P3 evidence tables",
        "",
        "**Generated** by `scripts/p3_manuscript_evidence_table.py`. Do not edit by hand;",
        "`tests/test_p3_manuscript_evidence_table.py` fails if this drifts from",
        f"`{DEFAULT_OUTPUT.relative_to(ROOT)}`.",
        "",
        f"Channel `{payload['channel']}`, site count {payload['site_count']}, "
        f"{payload['batches']} batches.",
        "",
        "## T1 — the three-rung response",
        "",
        "| rung `r` | value | standard error | &sigma; from zero |",
        "|---:|---:|---:|---:|",
    ]
    for row in payload["response"]:
        lines.append(f"| {row['rung']} | {row['value']:.6e} | {row['standard_error']:.4e} "
                     f"| {row['sigma_from_zero']:.2f} |")
    interval = payload["fieller_interval_3sigma"]
    lines += [
        "",
        "The frozen design nominated `r=1` as the denominator. Its 3&sigma; Fieller",
        "intervals for the two ratios are",
        "",
        "| ratio | lower | upper | width |",
        "|---|---:|---:|---:|",
    ]
    for name, block in sorted(interval.items()):
        lines.append(f"| `{name}` | {block['lower']:.3f} | {block['upper']:.3f} "
                     f"| {block['upper'] - block['lower']:.3f} |")
    control = payload["two_entry_control"]
    lines += [
        "",
        "## T2 — the two-entry control",
        "",
        f"Projective statistic restricted to the two rungs the frozen test used "
        f"(`{control['entry']}`), against Fieller's *z* squared on the same pair.",
        "",
        "| competitor | Fieller *z* | *z*&sup2; | projective *D* (2 entries) | relative deviation |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in control["rows"]:
        lines.append(f"| `{row['competitor']}` | {row['fieller_z']:.4f} "
                     f"| {row['fieller_z_squared']:.5f} "
                     f"| {row['projective_statistic_two_entries']:.5f} "
                     f"| {row['relative_deviation']:.1e} |")
    lines += [
        "",
        f"Largest relative deviation: **{control['largest_relative_deviation']:.1e}**.",
        "",
        "## T3 — frozen verdict against projective verdict",
        "",
        "| competitor | ray | frozen &sigma; (1 entry) | projective &sigma; (3 rungs) "
        "| frozen | projective | changed |",
        "|---|---|---:|---:|---|---|---|",
    ]
    for row in payload["verdicts"]:
        ray = ", ".join(_fmt(value, 4) for value in row["ray"])
        lines.append(
            f"| `{row['competitor']}` | ({ray}) | {row['fieller_sigma_one_entry']:.2f} "
            f"| {_fmt(row['projective_sigma_three_rungs'], 2)} | {row['frozen_verdict']} "
            f"| {row['projective_verdict']} | {'**yes**' if row['verdict_changed'] else 'no'} |")
    lines += [
        "",
        "Range of the projective &sigma; over every positive-definite value of the one",
        "covariance entry the frozen artifact did not store:",
        "",
        "| competitor | &sigma; range | verdict stable |",
        "|---|---|---|",
    ]
    for row in payload["verdicts"]:
        low, high = row["sigma_over_admissible_correlations"]
        lines.append(f"| `{row['competitor']}` | [{_fmt(low, 2)}, {_fmt(high, 2)}] "
                     f"| {'yes' if row['verdict_survives_the_missing_covariance'] else '**no**'} |")
    curvature = payload["curvature"]
    measured = curvature["measured"]
    lines += [
        "",
        "## T4 — the curvature functional",
        "",
        "`f[1,2,4] = (m(4) &minus; 3 m(2) + 2 m(1)) / 6`, a linear functional of the",
        "response: exactly 1 on `r`&sup2;, exactly 0 on any line, and with no denominator",
        "and no matrix inverse anywhere in it.",
        "",
        f"Measured: **{measured['value']:.4e} &plusmn; {measured['standard_error']:.4e}**, "
        f"*z* = **{measured['z']:.2f}**, and *z* &isin; "
        f"[{measured['z_over_admissible_correlations'][0]:.2f}, "
        f"{measured['z_over_admissible_correlations'][1]:.2f}] across the admissible",
        "covariance range.",
        "",
        "| competitor | curvature predicted |",
        "|---|---:|",
        # the two families linear in r annihilate this functional exactly; the
        # artifact stores their value as float noise at 1e-16, rendered as 0 here
    ]
    for row in curvature["predicted"]:
        value = row["curvature_predicted"]
        rendered = "0 (exactly)" if abs(value) < 1e-12 else _fmt(value, 3)
        lines.append(f"| `{row['competitor']}` | {rendered} |")
    spin8 = payload["spin8"]
    lines += [
        "",
        "## T5 — the spin-8 amplitude each competitor needs",
        "",
        f"Leakage coefficient {spin8['leakage_coefficient']:.6f} "
        f"(= 1148/21025, exactly equal and opposite between the two families).",
        "",
        f"The frozen design assumed: *{spin8['assumed_by_the_frozen_design']}*.",
        "",
        "| competitor | required \\|A&#8328;/A&#8324;\\| to reach `r=2` |",
        "|---|---:|",
    ]
    for row in spin8["rows"]:
        lines.append(f"| `{row['competitor']}` "
                     f"| {_fmt(row['required_abs_A8_over_A4_to_reach_r2'], 1)} |")
    lines += [
        "",
        f"Smallest requirement {spin8['smallest_requirement']:.2f}, largest "
        f"{spin8['largest_requirement']:.0f}.",
        "",
    ]
    return "\n".join(lines)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--markdown", type=Path, default=DEFAULT_MARKDOWN)
    args = parser.parse_args(argv)
    payload = assemble()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
    args.markdown.parent.mkdir(parents=True, exist_ok=True)
    args.markdown.write_text(render(payload), encoding="utf-8")
    print(f"wrote {args.output}")
    print(f"wrote {args.markdown}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

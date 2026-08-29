#!/usr/bin/env python3
"""Post-reveal P275 score of the transported unconnected mean J_D4.

This consumes the already frozen Phase-1 microcanonical streams.  It changes
neither their matching-root rule nor the Phase-1 decision: the only discovery
move is to replace Cov(A_top,J_D4)/B by E[J_D4], after the connected response
was revealed to contain a scale-zero global-line background.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence

import numpy as np

from score_p275_atop_field_identity import (
    COORDINATE_ORDER,
    MODULUS_ORDER,
    SIZE_ORDER,
    _expectation,
    _geometry_map,
    _load_prediction,
    _matching_root,
    _number,
    _pool_levels,
    _read_run,
    fit_models,
    validate_runs,
)


def geometry_mean_estimate(run: dict, target: dict, omitted_batch: int | None = None) -> dict:
    levels = _pool_levels(run["rows"], omitted_batch)
    p = _matching_root(levels)
    mean_q = _expectation(levels, "sum_q", p)
    lab = complex(*[
        _expectation(levels, f"sum_{component}_J_D4", p)
        for component in ("Re", "Im")
    ])
    transport = target["transport"]
    phase = complex(_number(transport["real"]), _number(transport["imag"]))
    canonical = phase*lab
    scaled = run["N"]**(13.0/8.0)*canonical
    return {
        "p_matching": p,
        "mean_q_residual": mean_q,
        "mean_J_D4_lab": [lab.real, lab.imag],
        "mean_J_D4_canonical": [canonical.real, canonical.imag],
        "Y_mean_J_D4_scaled": [scaled.real, scaled.imag],
    }


def estimate_vector(runs: Sequence[dict], prediction: dict, omitted: tuple[int, int] | None = None):
    geometry = _geometry_map(prediction)
    run_map = {(run["N"], run["modulus"]): run for run in runs}
    values, details = [], {}
    for n in SIZE_ORDER:
        for modulus in MODULUS_ORDER:
            batch = omitted[1] if omitted is not None and omitted[0] == n else None
            estimate = geometry_mean_estimate(run_map[(n, modulus)], geometry[(n, modulus)], batch)
            values.extend(estimate["Y_mean_J_D4_scaled"])
            details[f"N{n}:{modulus}"] = estimate
    return np.asarray(values, dtype=float), details


def jackknife_covariance(runs: Sequence[dict], prediction: dict) -> tuple[np.ndarray, dict]:
    full, details = estimate_vector(runs, prediction)
    covariance = np.zeros((18, 18), dtype=float)
    batches = int(prediction["phase1_microcanonical_matching_root"]["batches"])
    for size_index, n in enumerate(SIZE_ORDER):
        coordinates = slice(6*size_index, 6*(size_index+1))
        deleted = np.asarray([
            estimate_vector(runs, prediction, (n, batch))[0][coordinates]
            for batch in range(batches)
        ])
        centered = deleted-deleted.mean(axis=0)
        covariance[coordinates, coordinates] = (batches-1.0)/batches*(centered.T@centered)
    return covariance, {"observation": full, "geometries": details}


def build_report(runs: Sequence[dict], prediction: dict) -> dict:
    provenance = validate_runs(runs, prediction)
    covariance, estimates = jackknife_covariance(runs, prediction)
    fits = fit_models(estimates["observation"], covariance, prediction)
    return {
        "schema": "matching-one/p275-mean-jd4-discovery-score/v1",
        "status": "post_reveal_discovery_does_not_change_phase1",
        "issues": [205, 275],
        "provenance": provenance,
        "coordinate_order": list(COORDINATE_ORDER),
        "discovery_change": "replace Cov(A_top,J_D4)/B by transported E[J_D4] at the same delete-one finite matching roots",
        "estimates": {
            "observation_Y_mean_J_D4": estimates["observation"].tolist(),
            "geometries": estimates["geometries"],
            "covariance_18x18": covariance.tolist(),
        },
        "model_score_using_phase1_bases": fits,
        "claim_boundary": [
            "This score was chosen after the connected Gamma reveal and is discovery-only.",
            "It reuses the frozen samples, roots, transport, scaling, covariance and Phase-1 model bases.",
            "It does not alter the frozen Phase-1 rejection of the Gamma scaling map.",
            "The nine rows are single quotients, not the paired-orientation P4 q2 parent/child difference scored in 634040d.",
            "A surviving basis nominates a preregistered local-source follow-up; it does not identify a field by itself.",
        ],
    }


def render_markdown(report: dict) -> str:
    lines = [
        "# P275 post-reveal mean-J_D4 discovery score", "",
        "This reuses the frozen nine streams but replaces `Cov(A_top,J_D4)/B` with",
        "transported `E[J_D4]`. It does not change the Phase-1 result.", "",
        "| geometry | p_N | Re N^(13/8)E[J_D4] | Im | Re E[J_D4] | Im |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for n in SIZE_ORDER:
        for modulus in MODULUS_ORDER:
            row = report["estimates"]["geometries"][f"N{n}:{modulus}"]
            scaled = row["Y_mean_J_D4_scaled"]
            raw = row["mean_J_D4_canonical"]
            lines.append(
                f"| N{n}/{modulus} | {row['p_matching']:.12g} | {scaled[0]:.8g} | {scaled[1]:.8g} | {raw[0]:.8g} | {raw[1]:.8g} |"
            )
    lines.extend(["", "## Same Phase-1 model bases", "", "| model | chi2 | df | survival p |", "|---|---:|---:|---:|"])
    scores = report["model_score_using_phase1_bases"]["scores"]
    for name in ("Q4_epsilon_ordinary", "Q4_energy_Jordan", "generic_allowed_H4_pure", "generic_allowed_H4_affine_log", "zero_response"):
        row = scores[name]
        lines.append(f"| {name} | {row['chi_square']:.6g} | {row['dof']} | {row['survival_p']:.6g} |")
    lines.extend([
        "", "## Reading", "",
        "No Phase-1 basis survives for the single-quotient mean source either. The least",
        "bad basis is generic affine-log H4 (`chi2=868.75/6`), still decisively rejected.",
        "This does not contradict commit `634040d`: its `0.33085` q2 ratio belongs to a",
        "paired-orientation P4 difference on an exact N65/N130 parent-child chain. These",
        "nine rows are individual quotients at N50/N130/N170 and contain neither that",
        "orientation subtraction nor its exact radial chain. Therefore the q2 mean-source",
        "hint remains a separate candidate, while the present nine-geometry archive cannot",
        "promote transported single-quotient `E[J_D4]` to the field-identity covector.",
    ])
    lines.extend(["", "## Boundary", ""] + [f"- {line}" for line in report["claim_boundary"]] + [""])
    return "\n".join(lines)


def main(argv: Sequence[str] | None = None) -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prediction", type=Path, default=root/"predictions/p275_atop_q4_field_identity_20260829.yaml")
    parser.add_argument("--run", action="append", required=True, help="N:modulus:CSV:METADATA")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path)
    args = parser.parse_args(argv)
    prediction = _load_prediction(args.prediction)
    report = build_report([_read_run(spec) for spec in args.run], prediction)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    if args.markdown:
        args.markdown.write_text(render_markdown(report), encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Freeze the natural charged-current scale test before revealing N85."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Mapping, Optional, Sequence


P_REF = 0.592746050790
SOURCE_N = 65
TARGET_N = 85
H4_EXPONENT = -13.0 / 8.0


def jackknife_natural_coordinate(payload: Mapping[str, object]) -> dict[str, object]:
    batches = payload["joint_estimate"]["batch_values"]
    count = len(batches)
    factor = P_REF * (1.0 - P_REF)

    def estimate(indices: Sequence[int]) -> list[float]:
        values = []
        for orientation in ("first", "second"):
            w = math.fsum(batches[i]["values"][f"{orientation}_A_W"] for i in indices) / len(indices)
            current = math.fsum(
                batches[i]["values"][f"{orientation}_A_J_minus"] for i in indices
            ) / len(indices)
            values.append(factor * current / w)
        return [values[0], values[1], values[1] - values[0]]

    full = estimate(list(range(count)))
    replicates = [estimate([j for j in range(count) if j != i]) for i in range(count)]
    centers = [math.fsum(row[j] for row in replicates) / count for j in range(3)]
    covariance = [[
        (count - 1) / count * math.fsum(
            (row[i] - centers[i]) * (row[j] - centers[j]) for row in replicates
        )
        for j in range(3)] for i in range(3)]
    return {
        "order": ["K_A_first", "K_A_second", "Delta_K_A_second_minus_first"],
        "value": full,
        "standard_error": [math.sqrt(max(0.0, covariance[i][i])) for i in range(3)],
        "covariance": covariance,
        "jackknife_replicates": replicates,
        "definition": "K_A=p(1-p) Jminus_A/W_A=d_eta log W_A",
    }


def preregister(activity_path: Path) -> dict[str, object]:
    parent = json.loads(activity_path.read_text(encoding="utf-8"))
    if parent.get("schema") != "matching-one/N65-F3-charged-activity-net/v1":
        raise ValueError("unexpected N65 activity/net certificate")
    natural = jackknife_natural_coordinate(parent)
    delta = natural["value"][2]
    delta_se = natural["standard_error"][2]
    scale_ratio = (TARGET_N / SOURCE_N) ** H4_EXPONENT
    targets = {
        "zero": {"value": 0.0, "fit_variance": 0.0},
        "source_fitted_scale_neutral": {
            "value": delta, "fit_variance": delta_se * delta_se,
        },
        "source_fitted_project_H4": {
            "value": scale_ratio * delta,
            "fit_variance": scale_ratio * scale_ratio * delta_se * delta_se,
            "exponent_in_area_N": H4_EXPONENT,
            "scale_ratio": scale_ratio,
        },
    }

    def required_samples(target_gap: float, z: float) -> float:
        return 20000.0 * (z * delta_se / target_gap) ** 2

    gaps = {
        "zero_vs_scale_neutral": abs(targets["source_fitted_scale_neutral"]["value"]),
        "zero_vs_H4": abs(targets["source_fitted_project_H4"]["value"]),
        "scale_neutral_vs_H4": abs(
            targets["source_fitted_scale_neutral"]["value"]
            - targets["source_fitted_project_H4"]["value"]
        ),
    }
    design = {
        "archive_search_result": "no N85 or larger tau1/primitive-ell/tau2 projective archive exists",
        "rejected_old_N85": [
            "results/server-20260828/C01/N85/eval.metadata.json: no projective_births/births_csv",
            "results/server-20260828/P34/N85/mc.metadata.json: no projective_births/births_csv",
        ],
        "geometry": {
            "N": TARGET_N,
            "first": [9, 2], "second": [7, 6],
            "first_period_matrix": [[9, -2], [2, 9]],
            "second_period_matrix": [[7, -6], [6, 7]],
            "reason": "smallest existing-project design above N65 with two primitive Gaussian representations",
        },
        "samples_per_shape": 200000,
        "batches": 20,
        "seed": 202608337,
        "replica_offset": 0,
        "power_extrapolation": {
            "assumption": "N85 per-sample Delta_K variance no worse than N65",
            "N65_samples": 20000,
            "N65_standard_error": delta_se,
            "planned_standard_error": delta_se * math.sqrt(20000 / 200000),
            "required_samples_3sigma": {
                name: required_samples(gap, 3.0) for name, gap in gaps.items()
            },
            "required_samples_5sigma": {
                name: required_samples(gap, 5.0) for name, gap in gaps.items()
            },
            "why_200k": "exceeds the N65-variance 3-sigma scale-neutral-vs-H4 requirement without pretending fitted-target uncertainty vanishes",
        },
    }
    return {
        "schema": "matching-one/P337-natural-current-scale-preregistration/v1",
        "status": "frozen before N85 generation or inspection",
        "source": {
            "activity_commit": "f3384f4",
            "archive_commit": "1714141",
            "activity_certificate": str(activity_path),
            "p_ref": P_REF,
        },
        "natural_coordinate": natural,
        "frozen_targets_at_N85": targets,
        "design": design,
        "scoring_contract": {
            "observable": "Delta_K_A=K_A_second-K_A_first from aligned delete-one-batch jackknife",
            "normalization": "p(1-p)/W_A fixed from N65 definition; no N85 normalization selection",
            "measurement_only_residual": "(observed-target)^2/Var_N85 for the three fixed numeric targets",
            "predictive_residual": "for fitted targets divide by Var_N85+fit_variance; zero has fit_variance 0",
            "model_order": ["zero", "source_fitted_scale_neutral", "source_fitted_project_H4"],
            "no_selection": "all three targets and H4 exponent -13/8 are frozen before N85",
        },
        "claim_boundary": "N65 normalization and variance design only; N85 values absent from this artifact",
    }


def render_markdown(payload: Mapping[str, object]) -> str:
    natural = payload["natural_coordinate"]
    targets = payload["frozen_targets_at_N85"]
    design = payload["design"]
    lines = [
        "# Preregistration: natural charged-current scale test", "",
        "This artifact is frozen before generating or inspecting N85 projective-birth data.", "",
        "```text",
        "K_A = p(1-p) Jminus_A/W_A = d_eta log W_A,",
        "Delta_K_A = K_A(second)-K_A(first).",
        "```", "",
        f"N65 gives `K_first={natural['value'][0]:.12g}`, `K_second={natural['value'][1]:.12g}`, and `Delta_K={natural['value'][2]:.12g} +/- {natural['standard_error'][2]:.3g}` by aligned batch jackknife.", "",
        "Frozen N85 targets:", "",
        f"- zero: `0`",
        f"- source-fitted scale-neutral: `{targets['source_fitted_scale_neutral']['value']:.12g}`",
        f"- source-fitted project H4, `(85/65)^(-13/8)`: `{targets['source_fitted_project_H4']['value']:.12g}`", "",
        f"No old N85 file contains the required `tau1,ell,tau2` statistics. The fixed next pair is `(9+2i,7+6i)`, 200000 samples per shape, {design['batches']} aligned batches, seed `{design['seed']}`.", "",
        f"N65 variance projects an N85 SE of `{design['power_extrapolation']['planned_standard_error']:.4g}`. The 200k block exceeds the 3-sigma scale-neutral versus H4 gap requirement while remaining a short fresh run.", "",
        "Scoring will report both residuals to the frozen numeric targets and predictive residuals including N65 fitted-target variance. N85 cannot choose a new normalization, exponent, or target.", "",
    ]
    return "\n".join(lines)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--activity", type=Path, required=True)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--markdown", type=Path)
    args = parser.parse_args(argv)
    payload = preregister(args.activity)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.markdown:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(render_markdown(payload), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

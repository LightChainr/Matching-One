#!/usr/bin/env python3
"""Frozen closure-first scorer for the N325 Z5 charged three-point stream."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mpmath as mp

from norm5_chiral_hecke_phase import gaussian_ratio_power
from score_norm5_chiral_phase import inverse, matvec, quadratic
from z5_charged_threepoint_mc import PRIMARY_REAL_ORDER, PRODUCTION_ID


SPINS = (4, 8, 12)


def chi_square_survival(chi_square: float, dof: int) -> float:
    value = mp.gammainc(mp.mpf(dof) / 2, mp.mpf(chi_square) / 2, mp.inf) / mp.gamma(mp.mpf(dof) / 2)
    return float(value)


def gls_model(mean, covariance, q: complex) -> dict:
    precision = inverse(covariance)
    design = [
        [q.real, -q.imag, 0.0, 0.0],
        [q.imag, q.real, 0.0, 0.0],
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, q.real, -q.imag],
        [0.0, 0.0, q.imag, q.real],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ]
    precision_design = [
        [sum(precision[i][k] * design[k][j] for k in range(8)) for j in range(4)]
        for i in range(8)
    ]
    normal = [
        [sum(design[k][i] * precision_design[k][j] for k in range(8)) for j in range(4)]
        for i in range(4)
    ]
    rhs_precision = matvec(precision, mean)
    rhs = [sum(design[k][i] * rhs_precision[k] for k in range(8)) for i in range(4)]
    beta = matvec(inverse(normal), rhs)
    fitted = matvec(design, beta)
    residual = [mean[i] - fitted[i] for i in range(8)]
    chi_square = quadratic(residual, precision)
    return {
        "fitted_amplitudes_A_re_im_B_re_im": beta,
        "fitted_vector": fitted,
        "residual": residual,
        "chi_square": chi_square,
        "degrees_of_freedom": 4,
        "survival_p": chi_square_survival(chi_square, 4),
    }


def zero_score(point, covariance) -> dict:
    chi_square = quadratic(point, inverse(covariance))
    return {
        "chi_square": chi_square,
        "degrees_of_freedom": len(point),
        "survival_p": chi_square_survival(chi_square, len(point)),
    }


def validate(payload: dict, manifest: dict | None) -> None:
    if payload.get("schema") != "matching-one/z5-charged-threepoint-response/v1":
        raise ValueError("unexpected response schema")
    if tuple(payload["analysis"]["primary_order"]) != PRIMARY_REAL_ORDER:
        raise ValueError("primary order changed")
    if not payload["mapping_gate"].get("passed"):
        raise ValueError("mapping gate failed")
    if payload["analysis"]["conjugacy_max_abs"] > 1e-12:
        raise ValueError("charged DFT conjugacy failed")
    if manifest is not None:
        if not manifest.get("production_authorized") or manifest.get("production_id") != PRODUCTION_ID:
            raise ValueError("manifest is not authorized")
        observed = payload["run"]
        for key in ("samples", "batches", "workers", "p", "seed", "radius", "replica_offset"):
            if observed[key] != manifest["run"][key]:
                raise ValueError(f"run differs from manifest for {key}")
        if payload.get("manifest_runner_commit") != manifest.get("runner_commit"):
            raise ValueError("runner commit mismatch")


def score(payload: dict, manifest: dict | None = None) -> dict:
    mp.mp.dps = 60
    validate(payload, manifest)
    analysis = payload["analysis"]
    mean = [float(value) for value in analysis["primary_point"]]
    covariance = [[float(value) for value in row] for row in analysis["primary_covariance_of_mean"]]
    closure = analysis["closure"]
    closure_score = zero_score(
        [float(value) for value in closure["point_re_im"]],
        [[float(value) for value in row] for row in closure["delete_one_covariance_re_im"]],
    )
    models = {}
    for spin in SPINS:
        real, imag, denominator = gaussian_ratio_power(3 * spin)
        q = complex(real / denominator, imag / denominator)
        row = gls_model(mean, covariance, q)
        row.update({
            "spin_each": spin,
            "cubic_spin": 3 * spin,
            "exact_q_re": f"{real}/{denominator}",
            "exact_q_im": f"{imag}/{denominator}",
        })
        models[f"H{spin}"] = row
    ranking = sorted(models, key=lambda name: models[name]["chi_square"])
    nonneutral = analysis["nonneutral_controls"]
    nonneutral_score = zero_score(
        [float(value) for value in nonneutral["point"]],
        [[float(value) for value in row] for row in nonneutral["covariance_of_mean"]],
    )
    return {
        "schema": "matching-one/z5-charged-threepoint-score/v1",
        "status": "frozen_production_reveal" if manifest else "variance_smoke_only",
        "score_order": ["cross_product_closure", "joint_H4_H8_H12_GLS"],
        "primary_order": list(PRIMARY_REAL_ORDER),
        "primary_point": mean,
        "primary_covariance_8x8": covariance,
        "cross_product_closure": {
            **closure_score,
            "point_re_im": closure["point_re_im"],
            "covariance_2x2": closure["delete_one_covariance_re_im"],
            "relation": closure["relation"],
        },
        "joint_models": models,
        "ranking": ranking,
        "nonneutral_joint_zero_control": {
            **nonneutral_score,
            "order": nonneutral["order"],
            "point": nonneutral["point"],
        },
        "DFT_conjugacy_max_abs": analysis["conjugacy_max_abs"],
        "claim_boundary": [
            "Closure is scored before the conditional single-spin models.",
            "The three GLS targets are frozen q_(3s) phases with two complex channel amplitudes; no per-channel hand phase is fitted.",
            "This does not claim a universal raw cubic phase or a specific continuum primary.",
        ],
    }


def render(result: dict) -> str:
    closure = result["cross_product_closure"]
    lines = [
        "# Z5 charged three-point closure-first score", "",
        "## Cross-product closure", "",
        f"- point: `{closure['point_re_im']}`",
        f"- chi-square: `{closure['chi_square']}/{closure['degrees_of_freedom']}`, p `{closure['survival_p']}`",
        "", "## Frozen joint GLS", "",
        "| model | chi-square/df | p |", "|---|---:|---:|",
    ]
    for name in result["ranking"]:
        row = result["joint_models"][name]
        lines.append(f"| {name} | {row['chi_square']}/{row['degrees_of_freedom']} | {row['survival_p']} |")
    control = result["nonneutral_joint_zero_control"]
    lines += [
        "", "## Controls", "",
        f"- nonneutral joint zero: `{control['chi_square']}/{control['degrees_of_freedom']}`, p `{control['survival_p']}`",
        f"- DFT conjugacy maximum residual: `{result['DFT_conjugacy_max_abs']}`", "",
        "## Boundary", "",
        *[f"- {line}" for line in result["claim_boundary"]], "",
    ]
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("response", type=Path)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text()) if args.manifest else None
    result = score(json.loads(args.response.read_text()), manifest)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    args.markdown.write_text(render(result))
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

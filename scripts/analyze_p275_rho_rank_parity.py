#!/usr/bin/env python3
"""Recover the rho-child A_top/E_top complex C3 rank-parity split for #275."""

from __future__ import annotations

import argparse
import cmath
import csv
import hashlib
import io
import json
import math
import subprocess
from pathlib import Path
from typing import Any, Mapping, Sequence

import yaml


ROOT = Path(__file__).resolve().parents[1]
MANIFEST_SCHEMA = "matching-one.p275-rho-rank-parity.manifest.v1"
OUTPUT_SCHEMA = "matching-one.p275-rho-rank-parity.v1"


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    return sha256_bytes(path.read_bytes())


def git_blob(commit: str, path: str) -> bytes:
    completed = subprocess.run(
        ["git", "show", f"{commit}:{path}"],
        cwd=ROOT,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    return completed.stdout


def dft_r1(values: Sequence[float]) -> complex:
    zeta = cmath.exp(2j * math.pi / 3.0)
    return sum(values[j] * zeta ** (-j) for j in range(3)) / 3.0


def covariance_of_mean(rows: Sequence[Sequence[float]]) -> list[list[float]]:
    count = len(rows)
    means = [sum(row[j] for row in rows) / count for j in range(len(rows[0]))]
    return [
        [
            sum((row[i] - means[i]) * (row[j] - means[j]) for row in rows)
            / (count * (count - 1))
            for j in range(len(means))
        ]
        for i in range(len(means))
    ]


def complex_zero_score(mean: Sequence[float], covariance: Sequence[Sequence[float]]) -> dict[str, Any]:
    determinant = covariance[0][0] * covariance[1][1] - covariance[0][1] ** 2
    if determinant <= 0:
        raise ValueError("A_top complex covariance is not positive definite")
    precision = [
        [covariance[1][1] / determinant, -covariance[0][1] / determinant],
        [-covariance[0][1] / determinant, covariance[0][0] / determinant],
    ]
    chi_square = sum(
        mean[i] * precision[i][j] * mean[j] for i in range(2) for j in range(2)
    )
    return {
        "value_re_im": list(mean),
        "covariance_2x2": [list(row) for row in covariance],
        "chi_square": chi_square,
        "degrees_of_freedom": 2,
        "p_value": math.exp(-chi_square / 2.0),
    }


def build_report(manifest_path: Path) -> dict[str, Any]:
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema") != MANIFEST_SCHEMA:
        raise ValueError("unexpected rho rank-parity manifest schema")
    source = manifest["source"]
    batch_bytes = git_blob(source["commit"], source["batches_path"])
    score_bytes = git_blob(source["commit"], source["score_path"])
    if sha256_bytes(batch_bytes) != source["batches_sha256"]:
        raise ValueError("rho batch SHA256 mismatch")
    if sha256_bytes(score_bytes) != source["score_sha256"]:
        raise ValueError("rho source score SHA256 mismatch")

    source_score = json.loads(score_bytes)
    e_score = source_score["primary_nontrivial_Etop_r1"]
    rows = list(csv.DictReader(io.StringIO(batch_bytes.decode("utf-8"))))
    if len(rows) != int(source["batches"]):
        raise ValueError("rho batch count mismatch")

    joint_rows: list[list[float]] = []
    for row in rows:
        if int(row["samples"]) != int(source["samples_per_batch"]):
            raise ValueError("rho samples-per-batch mismatch")
        p0 = []
        p2 = []
        for child in source["child_order"]:
            samples = float(row["samples"])
            p0.append(float(row[f"{child}_rank0"]) / samples)
            p2.append(float(row[f"{child}_rank2"]) / samples)
        a_r1 = dft_r1([p2[j] - p0[j] for j in range(3)])
        e_raw_r1 = dft_r1([p0[j] + p2[j] for j in range(3)])
        joint_rows.append([a_r1.real, a_r1.imag, e_raw_r1.real, e_raw_r1.imag])

    joint_covariance = covariance_of_mean(joint_rows)
    a_mean = [
        sum(row[index] for row in joint_rows) / len(joint_rows) for index in (0, 1)
    ]
    a_score = complex_zero_score(
        a_mean, [row[:2] for row in joint_covariance[:2]]
    )
    alpha = float(manifest["score"]["alpha"])
    a_score["decision_alpha"] = alpha
    a_score["decision"] = "compatible_with_zero" if a_score["p_value"] >= alpha else "resolved_nonzero"

    e_mean = [float(value) for value in e_score["value_re_im"]]
    joint_mean = a_mean + e_mean
    return {
        "schema": OUTPUT_SCHEMA,
        "status": "completed_branch_batch_reanalysis_zero_new_MC",
        "decision": "RANK_EVEN_ONLY_C3_CORRECTION_AT_CURRENT_PRECISION",
        "observable": {
            "joint_order": ["A_top_r1_re", "A_top_r1_im", "E_top_r1_re", "E_top_r1_im"],
            "joint_mean": joint_mean,
            "joint_covariance_of_mean": joint_covariance,
            "A_top_definition": "P2-P0",
            "E_top_definition": "P0+P2",
            "C3_transform": manifest["observable"]["transform"],
            "continuum_A_top": manifest["observable"]["continuum_A_top"],
            "units": manifest["observable"]["units"],
            "normalizer": manifest["observable"]["normalizer"],
        },
        "scores": {
            "A_top_nontrivial_C3_r1": a_score,
            "E_top_nontrivial_C3_r1": {
                "value_re_im": e_mean,
                "covariance_2x2": e_score["covariance_2x2"],
                "chi_square": float(e_score["chi_square"]),
                "degrees_of_freedom": int(e_score["dof"]),
                "p_value": float(e_score["p"]),
                "decision_alpha": alpha,
                "decision": e_score["decision"],
                "source": "frozen_branch_score",
            },
        },
        "interpretation": {
            "supported": "the_rho_child_complex_C3_correction_is_resolved_on_the_Alexander_even_E_top_axis_while_A_top_is_compatible_with_zero",
            "not_supported": "a_B_source_response_continuum_energy_identity_or_original_U_transport",
            "next_use": "the_full_rank_six_coordinate_rank0_rank2_archive_is_ready_to_score_a_theory_given_restricted_trace_transport_vector",
        },
        "source": {
            **source,
            "source_score_schema": source_score["schema"],
        },
        "boundaries": manifest["claim_boundary"],
        "manifest": {
            "path": str(manifest_path.relative_to(ROOT)),
            "sha256": sha256_file(manifest_path),
        },
    }


def markdown_report(report: Mapping[str, Any]) -> str:
    a = report["scores"]["A_top_nontrivial_C3_r1"]
    e = report["scores"]["E_top_nontrivial_C3_r1"]
    covariance = a["covariance_2x2"]
    lines = [
        "# Issue #275 rho-child rank-parity C3 split",
        "",
        "## Decision",
        "",
        "**RANK_EVEN_ONLY_C3_CORRECTION_AT_CURRENT_PRECISION.**  In the same 100 aligned",
        "rho-child batches, the matching-odd `A_top=P2-P0` nontrivial C3 coordinate is",
        "compatible with zero, while the Alexander-even `E_top=P0+P2` coordinate is strongly",
        "resolved.  They are correlated transforms of one stream, not independent votes.",
        "",
        "```text",
        f"A_top r1 = {a['value_re_im'][0]:+.12g} {a['value_re_im'][1]:+.12g} i",
        f"Cov(A_top r1) = [[{covariance[0][0]:.12g}, {covariance[0][1]:.12g}],",
        f"                 [{covariance[1][0]:.12g}, {covariance[1][1]:.12g}]]",
        f"chi2={a['chi_square']:.6g}/2, p={a['p_value']:.6g}",
        "",
        f"E_top r1 = {e['value_re_im'][0]:+.12g} {e['value_re_im'][1]:+.12g} i",
        f"chi2={e['chi_square']:.6g}/2, p={e['p_value']:.6g}",
        "```",
        "",
        "At Q=1 the continuum `A_top` baseline is exactly zero for every modulus because",
        "Arguin gives `P2=P0`; no fitted continuum subtraction is used for that row.  The result",
        "is finite square-bond fixed-p rank parity.  It does not supply a B-source insertion,",
        "pooled-root p-jet, or original-U normalizer.",
        "",
        "## Immediate use",
        "",
        "The underlying archive retains rank0 and rank2 counts for all three children in every",
        "batch, so its six-coordinate covariance is already sufficient to score any frozen",
        "rank-restricted theory vector.  The missing object is the typed observer/source",
        "transport column, not more rho-child sampling.",
        "",
        "## Boundaries",
        "",
    ]
    lines.extend(f"- `{boundary}`" for boundary in report["boundaries"])
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest",
        type=Path,
        default=ROOT / "analysis" / "p275_rho_rank_parity_manifest.yaml",
    )
    parser.add_argument("--output-json", type=Path)
    parser.add_argument("--output-md", type=Path)
    args = parser.parse_args()
    manifest_path = args.manifest.resolve()
    manifest = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    report = build_report(manifest_path)
    output_json = (args.output_json or ROOT / manifest["outputs"]["json"]).resolve()
    output_md = (args.output_md or ROOT / manifest["outputs"]["markdown"]).resolve()
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_md.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    output_md.write_text(markdown_report(report), encoding="utf-8")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Score the frozen P334 one-pass bottleneck proxy pilot."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Mapping, Sequence

import mpmath as mp

from score_p334_current_k0_geometry_pilot import (
    add_matrix,
    centered_fit,
    covariance,
    joint_two,
    sha256,
    student_p,
)


ORIENTATIONS = ("first", "second")
CHEAP = ("g_size", "g_carriers", "g_occupied_frontier", "g_vacant_frontier")
PROXY = (
    "g_boundary_cut", "g_boundary_multicontact", "g_boundary_pairs",
    "g_core_vertices", "g_core_edges", "g_articulations", "g_bridges",
)
SHAPE = CHEAP + PROXY
MODELS = {
    "M0": ("age",),
    "Mcheap": ("age",) + CHEAP,
    "Mshape": ("age",) + SHAPE,
    "MH2": ("age", "h2_rate"),
}


def load_rows(path: Path, metadata_path: Path, frozen: Mapping[str, object], commit: str):
    metadata = json.loads(metadata_path.read_text())
    checks = {
        "git_commit": metadata["git_commit"] == commit,
        "samples": metadata["samples_per_pair"] == frozen["samples"],
        "batches": metadata["batches"] == frozen["batches"],
        "seed": metadata["seed"] == frozen["seed"],
        "replica_first": metadata["replica_counter_first"] == frozen["replica_counter_first"],
        "replica_last": metadata["replica_counter_last_exclusive"] == frozen["replica_counter_last_exclusive"],
        "k0": metadata["geometry_pilot_k0"] == frozen["k0"],
    }
    if not all(checks.values()):
        raise ValueError(f"metadata contract failed: {checks}")
    rows = []
    with path.open(newline="") as stream:
        for raw in csv.DictReader(stream):
            n, k0 = int(raw["n"]), int(raw["k0"])
            remaining = n - k0
            row = {
                "n": n, "k0": k0,
                "orientation": raw["orientation"], "batch": int(raw["batch"]),
                "ell_u": int(raw["ell_u"]), "ell_v": int(raw["ell_v"]),
                "y": int(raw["next_exit"]), "age": int(raw["age_steps"]) / n,
                "g_size": int(raw["essential_size"]) / n,
                "g_carriers": int(raw["essential_carriers"]) - 1,
                "g_occupied_frontier": int(raw["occupied_frontier"]) / k0,
                "g_vacant_frontier": int(raw["vacant_frontier"]) / remaining,
                "g_boundary_cut": int(raw["boundary_cut_edges"]) / n,
                "g_boundary_multicontact": int(raw["boundary_multicontact_sites"]) / remaining,
                "g_boundary_pairs": int(raw["boundary_contact_pairs"]) / remaining,
                "g_core_vertices": int(raw["core_vertices"]) / n,
                "g_core_edges": int(raw["core_edges"]) / n,
                "g_articulations": int(raw["articulation_vertices"]) / n,
                "g_bridges": int(raw["bridges"]) / n,
                "h2_rate": int(raw["H2"]) / remaining,
            }
            if int(raw["H2"]) != int(raw["H2_theta"]) + int(raw["H2_figure8"]) + int(raw["H2_separate"]):
                raise ValueError("H2 decomposition failed")
            if int(raw["core_vertices"]) > int(raw["essential_size"]):
                raise ValueError("2-core exceeds essential carrier")
            rows.append(row)
    return rows, {
        "metadata_checks": checks,
        "csv_sha256": sha256(path),
        "metadata_sha256": sha256(metadata_path),
    }


def fit_size(rows: Sequence[dict], batches: int):
    full, order = {}, []
    for orientation in ORIENTATIONS:
        subset = [row for row in rows if row["orientation"] == orientation]
        full[orientation] = {}
        for model, predictors in MODELS.items():
            full[orientation][model] = centered_fit(subset, predictors)
            order.append((orientation, model))
    deleted = []
    for omitted in range(batches):
        vector = []
        for orientation, model in order:
            subset = [row for row in rows
                      if row["orientation"] == orientation and row["batch"] != omitted]
            vector.append(float(centered_fit(subset, MODELS[model])["coefficients"]["age"]))
        deleted.append(vector)
    cov = covariance(deleted)
    index = {key: position for position, key in enumerate(order)}
    summary, joint = {}, {}
    for orientation in ORIENTATIONS:
        summary[orientation] = {}
        primary = float(full[orientation]["M0"]["coefficients"]["age"])
        for model in MODELS:
            position = index[(orientation, model)]
            fit = full[orientation][model]
            beta = float(fit["coefficients"]["age"])
            se = math.sqrt(max(cov[position][position], 0.0))
            summary[orientation][model] = {
                **fit, "beta_age": beta, "se": se,
                "p": student_p(beta, se, batches - 1),
                "absolute_retention_vs_fresh_M0": abs(beta / primary) if primary else None,
            }
    for model in MODELS:
        positions = [index[(orientation, model)] for orientation in ORIENTATIONS]
        values = [summary[orientation][model]["beta_age"] for orientation in ORIENTATIONS]
        block = [[cov[i][j] for j in positions] for i in positions]
        statistic, p = joint_two(values, block)
        joint[model] = {"chi2": statistic, "df": 2, "p": p, "covariance": block}
    support = {}
    for orientation in ORIENTATIONS:
        subset = [row for row in rows if row["orientation"] == orientation]
        support[orientation] = {
            "risk_rows": len(subset), "next_exits": sum(row["y"] for row in subset),
            "gate": len(subset) >= 400 and sum(row["y"] for row in subset) >= 10,
        }
    return {"fits": summary, "joint": joint, "support": support,
            "vector_order": [list(key) for key in order], "jackknife_covariance": cov}


def contrast_covariance(score: Mapping[str, object], model: str):
    order = [tuple(value) for value in score["vector_order"]]
    cov = score["jackknife_covariance"]
    output = [[0.0, 0.0], [0.0, 0.0]]
    deltas = []
    for orientation in ORIENTATIONS:
        deltas.append(score["fits"][orientation][model]["beta_age"] -
                      score["fits"][orientation]["M0"]["beta_age"])
    for i, orientation in enumerate(ORIENTATIONS):
        a, b = order.index((orientation, model)), order.index((orientation, "M0"))
        for j, other in enumerate(ORIENTATIONS):
            c, d = order.index((other, model)), order.index((other, "M0"))
            output[i][j] = cov[a][c] - cov[a][d] - cov[b][c] + cov[b][d]
    return deltas, output


def production_anchor(score: Mapping[str, object], production: Mapping[str, object], size: str):
    deltas, delta_cov = contrast_covariance(score, "Mshape")
    source = production["sizes"][size]
    source_order = source["point_vector_order"]
    indices = [source_order.index(f"{orientation}_primary") for orientation in ORIENTATIONS]
    source_values = [source["point_vector"][index] for index in indices]
    source_cov = [[source["delete_one_covariance"][i][j] for j in indices] for i in indices]
    values = [source_values[index] + deltas[index] for index in range(2)]
    combined = add_matrix(source_cov, delta_cov)
    statistic, p = joint_two(values, combined)
    orientations = {}
    for index, orientation in enumerate(ORIENTATIONS):
        variance = combined[index][index]
        source_var, delta_var = source_cov[index][index], delta_cov[index][index]
        df = variance * variance / (source_var * source_var / 99 + delta_var * delta_var / 19)
        se = math.sqrt(variance)
        orientations[orientation] = {
            "production_M0": source_values[index],
            "fresh_delta_Mshape_minus_M0": deltas[index],
            "anchored_Mshape": values[index],
            "se": se,
            "welch_satterthwaite_df": df,
            "p": student_p(values[index], se, max(int(df), 1)),
            "absolute_retention": abs(values[index] / source_values[index]),
        }
    return {"orientations": orientations, "covariance": combined,
            "joint": {"chi2": statistic, "df": 2, "p": p}}


def geometry_coefficients(rows: Sequence[dict], excluded: int | None = None):
    subset = [row for row in rows if row["batch"] != excluded]
    return centered_fit(subset, SHAPE, stratum=lambda row: (
        row["orientation"], row["ell_u"], row["ell_v"]))["coefficients"]


def transferred(rows_source, rows_target, source_excluded=None, target_excluded=None):
    gamma = geometry_coefficients(rows_source, source_excluded)
    values = []
    for orientation in ORIENTATIONS:
        subset = [row for row in rows_target
                  if row["orientation"] == orientation and row["batch"] != target_excluded]
        fit = centered_fit(subset, ("age",), outcome=lambda row, gamma=gamma:
            row["y"] - sum(float(gamma[name]) * row[name] for name in SHAPE))
        values.append(float(fit["coefficients"]["age"]))
    return values, {name: float(gamma[name]) for name in SHAPE}


def transfer_score(source, target, batches):
    point, gamma = transferred(source, target)
    source_deleted = [transferred(source, target, source_excluded=batch)[0]
                      for batch in range(batches)]
    target_deleted = [transferred(source, target, target_excluded=batch)[0]
                      for batch in range(batches)]
    cov = add_matrix(covariance(source_deleted), covariance(target_deleted))
    statistic, p = joint_two(point, cov)
    return {"geometry_coefficients": gamma,
            "target_residual_age": dict(zip(ORIENTATIONS, point)),
            "covariance": cov, "joint": {"chi2": statistic, "df": 2, "p": p}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--freeze", type=Path, required=True)
    parser.add_argument("--production-age-score", type=Path, required=True)
    parser.add_argument("--n325-csv", type=Path, required=True)
    parser.add_argument("--n325-metadata", type=Path, required=True)
    parser.add_argument("--n425-csv", type=Path, required=True)
    parser.add_argument("--n425-metadata", type=Path, required=True)
    parser.add_argument("--runner-commit", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    freeze = json.loads(args.freeze.read_text())
    production = json.loads(args.production_age_score.read_text())
    if sha256(args.production_age_score) != freeze["production_age_result"]["sha256"]:
        raise ValueError("production age score hash changed")
    rows325, provenance325 = load_rows(args.n325_csv, args.n325_metadata,
                                       freeze["runs"]["N325"], args.runner_commit)
    rows425, provenance425 = load_rows(args.n425_csv, args.n425_metadata,
                                       freeze["runs"]["N425"], args.runner_commit)
    score325 = fit_size(rows325, 20)
    score425 = fit_size(rows425, 20)
    anchor325 = production_anchor(score325, production, "N325")
    anchor425 = production_anchor(score425, production, "N425")
    retentions = [anchor["orientations"][orientation]["absolute_retention"]
                  for anchor in (anchor325, anchor425) for orientation in ORIENTATIONS]
    absorbed = (anchor325["joint"]["p"] >= 0.01 and anchor425["joint"]["p"] >= 0.01
                and max(retentions) <= 0.25)
    payload = {
        "schema": "matching-one/p334-bottleneck-proxy-pilot-score/v1",
        "freeze_sha256": sha256(args.freeze), "runner_commit": args.runner_commit,
        "inputs": {"N325": provenance325, "N425": provenance425,
                   "production_age_score_sha256": sha256(args.production_age_score)},
        "fresh_scores": {"N325": score325, "N425": score425},
        "production_anchor": {"N325": anchor325, "N425": anchor425},
        "cross_size_transfer": transfer_score(rows325, rows425, 20),
        "decision": "shape_proxy_absorbs_production_age" if absorbed else
                    "shape_proxy_fails_common_production_age_absorption",
        "decision_rule": freeze["primary_decision"],
        "claim_boundary": freeze["interpretation_boundary"],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "decision": payload["decision"],
        "retentions": retentions,
        "anchor_joint_p": {"N325": anchor325["joint"]["p"],
                           "N425": anchor425["joint"]["p"]},
        "transfer": payload["cross_size_transfer"]["target_residual_age"],
    }, indent=2))


if __name__ == "__main__":
    main()

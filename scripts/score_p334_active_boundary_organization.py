#!/usr/bin/env python3
"""Score the frozen P334 active-boundary organization pilot."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess

import numpy as np

from score_p334_current_k0_geometry_pilot import add_matrix, covariance, sha256
from score_p334_morphology_state_transport import (
    FEATURES, ORIENTATIONS, SIZES, chi_square, centered_arrays, load_rows,
    train_subspace,
)


ORG = ("axis_anisotropy", "corner_balance", "frontier_components_L",
       "largest_frontier_arc_L", "frontier_concentration")


def load_fresh(path: Path):
    rows = []
    with path.open(newline="") as stream:
        for raw in csv.DictReader(stream):
            n = int(raw["n"])
            length = math.sqrt(n)
            cut = int(raw["boundary_cut_edges"])
            pairs = int(raw["boundary_contact_pairs"])
            frontier = int(raw["vacant_frontier"])
            rows.append({
                "size": f"N{n}", "orientation": raw["orientation"],
                "batch": int(raw["batch"]),
                "line": (int(raw["ell_u"]), int(raw["ell_v"])),
                "age": int(raw["age_steps"]) / n, "y": float(raw["next_exit"]),
                "k1": int(raw["k1"]), "k2": int(raw["k2"]),
                "area_size": int(raw["essential_size"]) / n,
                "carriers": float(int(raw["essential_carriers"]) - 1),
                "occupied_frontier_L": int(raw["occupied_frontier"]) / length,
                "vacant_frontier_L": frontier / length,
                "boundary_cut_L": cut / length,
                "boundary_multicontact_L": int(raw["boundary_multicontact_sites"]) / length,
                "boundary_pairs_L": pairs / length,
                "core_vertices_area": int(raw["core_vertices"]) / n,
                "core_edges_area": int(raw["core_edges"]) / n,
                "articulations_L": int(raw["articulation_vertices"]) / length,
                "bridges_L": int(raw["bridges"]) / length,
                "axis_anisotropy": abs(int(raw["boundary_axis_imbalance"])) / max(cut, 1),
                "corner_balance": int(raw["boundary_corner_balance"]) / max(pairs, 1),
                "frontier_components_L": int(raw["frontier_components"]) / length,
                "largest_frontier_arc_L": int(raw["largest_frontier_component"]) / length,
                "frontier_concentration": int(raw["frontier_component_sumsq"]) /
                                            max(frontier * frontier, 1),
            })
    return rows


def centered(rows, outcome):
    groups = {}
    for row in rows:
        groups.setdefault(row["line"], []).append(row)
    age, y, physical, organization = [], [], [], []
    for members in groups.values():
        raw_y = [float(outcome(row)) for row in members]
        means = {
            "age": sum(float(row["age"]) for row in members) / len(members),
            "y": sum(raw_y) / len(members),
        }
        physical_means = [sum(float(row[name]) for row in members) / len(members)
                          for name in FEATURES]
        org_means = [sum(float(row[name]) for row in members) / len(members)
                     for name in ORG]
        for index, row in enumerate(members):
            age.append(float(row["age"]) - means["age"])
            y.append(raw_y[index] - means["y"])
            physical.append([float(row[name]) - physical_means[j]
                             for j, name in enumerate(FEATURES)])
            organization.append([float(row[name]) - org_means[j]
                                 for j, name in enumerate(ORG)])
    return tuple(np.asarray(value) for value in (age, y, physical, organization))


def cross_size_subspaces(old_rows):
    arrays = {key: centered_arrays(value) for key, value in old_rows.items()}
    output = {}
    for target, source in (("N325", "N425"), ("N425", "N325")):
        output[target] = train_subspace(
            arrays, [(source, orientation) for orientation in ORIENTATIONS], 1)
    return output


def fit_nested(rows, subspace, outcome):
    age, y, physical, organization = centered(rows, outcome)
    canonical = (physical / subspace["scales"]) @ subspace["canonical_basis"]
    base_score = canonical @ subspace["directions"].T
    base_design = np.column_stack((age, base_score))
    full_design = np.column_stack((age, base_score, organization))
    base_beta, *_ = np.linalg.lstsq(base_design, y, rcond=None)
    full_beta, *_ = np.linalg.lstsq(full_design, y, rcond=None)
    base_residual = y - base_design @ base_beta
    full_residual = y - full_design @ full_beta
    total = float(y @ y)
    return {
        "baseline_beta_age": float(base_beta[0]),
        "organization_beta_age": float(full_beta[0]),
        "delta_beta_age": float(full_beta[0] - base_beta[0]),
        "baseline_R2": 1.0 - float(base_residual @ base_residual) / total,
        "organization_R2": 1.0 - float(full_residual @ full_residual) / total,
        "incremental_R2": float((base_residual @ base_residual -
                                  full_residual @ full_residual) / total),
        "organization_coefficients": dict(zip(ORG, full_beta[2:].tolist())),
        "rows": len(rows),
    }


def old_anchor(old_transport):
    order = []
    values = []
    matrix = np.zeros((4, 4))
    for size, regime in (("N325", "N425_to_N325"), ("N425", "N325_to_N425")):
        block = old_transport["transport_scores"][regime]["rank1"]["production_anchor"]
        offset = len(order)
        for orientation, value in zip(ORIENTATIONS, block["age_slopes"]):
            order.append((size, orientation))
            values.append(float(value))
        matrix[offset:offset + 2, offset:offset + 2] = np.asarray(block["covariance"])
    return order, values, matrix.tolist()


def age_score(fresh, subspaces, old_transport, batches=20):
    order, anchor_values, anchor_cov = old_anchor(old_transport)

    def vector(omitted_size=None, omitted_batch=None):
        output = []
        fits = {}
        for key in order:
            rows = [row for row in fresh[key]
                    if not (key[0] == omitted_size and row["batch"] == omitted_batch)]
            fit = fit_nested(rows, subspaces[key[0]], lambda row: row["y"])
            output.append(fit["delta_beta_age"])
            fits[f"{key[0]}_{key[1]}"] = fit
        return output, fits

    delta, fits = vector()
    delta_cov = np.zeros((4, 4))
    for size in SIZES:
        deleted = [vector(size, batch)[0] for batch in range(batches)]
        delta_cov += np.asarray(covariance(deleted))
    adjusted = [anchor_values[index] + delta[index] for index in range(4)]
    combined = add_matrix(anchor_cov, delta_cov.tolist())
    by_size = {}
    for index, size in enumerate(SIZES):
        positions = [2 * index, 2 * index + 1]
        block = [[combined[i][j] for j in positions] for i in positions]
        values = [adjusted[i] for i in positions]
        statistic, df, p = chi_square(values, block)
        by_size[size] = {"age_slopes": values, "covariance": block,
                         "chi2": statistic, "df": df, "p": p}
    retention = [abs(adjusted[index] / anchor_values[index]) for index in range(4)]
    absorbed = all(by_size[size]["p"] >= 0.01 for size in SIZES) and max(retention) <= 0.25
    return {
        "vector_order": [list(key) for key in order], "fresh_fits": fits,
        "old_rank1_anchor_age_slopes": anchor_values,
        "old_rank1_anchor_covariance": anchor_cov,
        "fresh_paired_delta": delta,
        "fresh_paired_delta_covariance": delta_cov.tolist(),
        "organization_anchor_age_slopes": adjusted,
        "combined_covariance": combined,
        "absolute_retention_vs_old_rank1": retention,
        "size_joint": by_size,
        "decision": "organization_absorbs_remaining_rank1_age" if absorbed else
                    "organization_fails_common_absorption_of_remaining_rank1_age",
    }


def upper_matrix(vector):
    dimension = int((math.sqrt(8 * len(vector) + 1) - 1) / 2)
    output = np.zeros((dimension, dimension))
    index = 0
    for i in range(dimension):
        for j in range(i, dimension):
            output[i, j] = output[j, i] = vector[index]
            index += 1
    return output


def temporal_modes(temporal):
    output = {}
    for size in temporal["sizes"]:
        key = f"N{size['N']}"
        dimension = len(size["layers"])
        upper_count = dimension * (dimension + 1) // 2
        output[key] = {"layers": size["layers"], "vectors": {}}
        for oi, orientation in enumerate(ORIENTATIONS):
            kernel = upper_matrix(size["kernel_vector"][oi * upper_count:(oi + 1) * upper_count])
            values, vectors = np.linalg.eigh(kernel)
            vectors = vectors[:, np.argsort(values)[::-1]]
            output[key]["vectors"][orientation] = {
                "mode2": vectors[:, 1], "mode3": vectors[:, 2]}
    return output


def temporal_score(fresh, subspaces, temporal, batches=20):
    modes = temporal_modes(temporal)
    order = [(size, orientation, mode) for size in SIZES
             for orientation in ORIENTATIONS for mode in ("mode2", "mode3")]

    def score_vector(omitted_size=None, omitted_batch=None):
        values, fits = [], {}
        for size, orientation, mode in order:
            rows = [row for row in fresh[(size, orientation)]
                    if not (size == omitted_size and row["batch"] == omitted_batch)]
            layers = modes[size]["layers"]
            vector = modes[size]["vectors"][orientation][mode]
            fit = fit_nested(rows, subspaces[size], lambda row, layers=layers, vector=vector:
                float(np.asarray([int(row["k1"] <= layer) + int(row["k2"] <= layer)
                                  for layer in layers]) @ vector))
            values.append(fit["incremental_R2"])
            fits[f"{size}_{orientation}_{mode}"] = fit
        return values, fits

    values, fits = score_vector()
    cov = np.zeros((len(order), len(order)))
    for size in SIZES:
        deleted = [score_vector(size, batch)[0] for batch in range(batches)]
        cov += np.asarray(covariance(deleted))
    return {"vector_order": [list(key) for key in order],
            "incremental_R2": values, "fits": fits,
            "delete_one_covariance": cov.tolist(),
            "standard_errors": np.sqrt(np.maximum(np.diag(cov), 0.0)).tolist()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--freeze", type=Path, required=True)
    parser.add_argument("--old-transport-score", type=Path, required=True)
    parser.add_argument("--old-raw-lock", type=Path, required=True)
    parser.add_argument("--n325-csv", type=Path, required=True)
    parser.add_argument("--n325-metadata", type=Path, required=True)
    parser.add_argument("--n425-csv", type=Path, required=True)
    parser.add_argument("--n425-metadata", type=Path, required=True)
    parser.add_argument("--runner-commit", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    freeze = json.loads(args.freeze.read_text())
    if sha256(args.old_transport_score) != freeze["old_transport_score"]["sha256"]:
        raise ValueError("old transport score hash changed")
    if sha256(args.old_raw_lock) != freeze["old_raw_lock"]["sha256"]:
        raise ValueError("old raw lock hash changed")
    old_transport = json.loads(args.old_transport_score.read_text())
    old_lock = json.loads(args.old_raw_lock.read_text())
    old_rows = {}
    for size in SIZES:
        loaded = load_rows(Path(old_lock["runs"][size]["csv"]))
        for orientation in ORIENTATIONS:
            old_rows[(size, orientation)] = [row for row in loaded
                                             if row["orientation"] == orientation]
    subspaces = cross_size_subspaces(old_rows)

    fresh = {}
    inputs = {}
    for size, csv_path, metadata_path in (
            ("N325", args.n325_csv, args.n325_metadata),
            ("N425", args.n425_csv, args.n425_metadata)):
        metadata = json.loads(metadata_path.read_text())
        contract = freeze["runs"][size]
        checks = {
            "commit": metadata["git_commit"] == args.runner_commit,
            "samples": metadata["samples_per_pair"] == contract["samples"],
            "batches": metadata["batches"] == contract["batches"],
            "seed": metadata["seed"] == contract["seed"],
            "replica_first": metadata["replica_counter_first"] ==
                             contract["replica_counter_first"],
            "replica_last": metadata["replica_counter_last_exclusive"] ==
                            contract["replica_counter_last_exclusive"],
            "k0": metadata["geometry_pilot_k0"] == contract["k0"],
        }
        if not all(checks.values()):
            raise ValueError(f"{size} metadata contract failed: {checks}")
        loaded = load_fresh(csv_path)
        for orientation in ORIENTATIONS:
            fresh[(size, orientation)] = [row for row in loaded
                                          if row["orientation"] == orientation]
        inputs[size] = {"csv_sha256": sha256(csv_path),
                        "metadata_sha256": sha256(metadata_path),
                        "rows": len(loaded), "checks": checks}

    temporal_bytes = subprocess.check_output(
        ["git", "show", freeze["temporal_score"]["git_spec"]])
    if hashlib.sha256(temporal_bytes).hexdigest() != freeze["temporal_score"]["sha256"]:
        raise ValueError("temporal score hash changed")
    temporal = json.loads(temporal_bytes)
    payload = {
        "schema": "matching-one/p334-active-boundary-organization-score/v1",
        "freeze_sha256": sha256(args.freeze), "runner_commit": args.runner_commit,
        "inputs": inputs,
        "age_primary": age_score(fresh, subspaces, old_transport),
        "temporal_secondary": temporal_score(fresh, subspaces, temporal),
        "claim_boundary": freeze["claim_boundary"],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "decision": payload["age_primary"]["decision"],
        "retention": payload["age_primary"]["absolute_retention_vs_old_rank1"],
        "joint_p": {size: payload["age_primary"]["size_joint"][size]["p"]
                    for size in SIZES},
        "temporal_incremental_R2": payload["temporal_secondary"]["incremental_R2"],
    }, indent=2))


if __name__ == "__main__":
    main()

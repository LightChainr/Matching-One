#!/usr/bin/env python3
"""Score the frozen P334 local active-boundary motif pilot."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

from score_p334_current_k0_geometry_pilot import add_matrix, covariance, sha256
from score_p334_morphology_state_transport import (
    FEATURES, ORIENTATIONS, SIZES, chi_square, centered_arrays, load_rows,
    train_subspace,
)

BASE = ("occupied", "essential", "same_side_pairs", "opposite_side_pairs")
MOTIFS = tuple(f"{parity}_{name}" for name in BASE for parity in ("S", "D"))


def load_fresh(path: Path):
    rows = []
    with path.open(newline="") as stream:
        for raw in csv.DictReader(stream):
            n = int(raw["n"])
            length = math.sqrt(n)
            row = {
                "size": f"N{n}", "orientation": raw["orientation"],
                "batch": int(raw["batch"]),
                "line": (int(raw["ell_u"]), int(raw["ell_v"])),
                "age": int(raw["age_steps"]) / n, "y": float(raw["next_exit"]),
                "area_size": int(raw["essential_size"]) / n,
                "carriers": float(int(raw["essential_carriers"]) - 1),
                "occupied_frontier_L": int(raw["occupied_frontier"]) / length,
                "vacant_frontier_L": int(raw["vacant_frontier"]) / length,
                "boundary_cut_L": int(raw["boundary_cut_edges"]) / length,
                "boundary_multicontact_L":
                    int(raw["boundary_multicontact_sites"]) / length,
                "boundary_pairs_L": int(raw["boundary_contact_pairs"]) / length,
                "core_vertices_area": int(raw["core_vertices"]) / n,
                "core_edges_area": int(raw["core_edges"]) / n,
                "articulations_L": int(raw["articulation_vertices"]) / length,
                "bridges_L": int(raw["bridges"]) / length,
            }
            for name in BASE:
                birth = float(raw[f"birth_r1_{name}"])
                exit_ = float(raw[f"exit_r1_{name}"])
                row[f"S_{name}"] = 0.5 * (birth + exit_)
                row[f"D_{name}"] = 0.5 * (exit_ - birth)
            rows.append(row)
    return rows


def centered(rows):
    groups = {}
    for row in rows:
        groups.setdefault(row["line"], []).append(row)
    age, y, physical, motifs = [], [], [], []
    for members in groups.values():
        means = {name: sum(float(row[name]) for row in members) / len(members)
                 for name in ("age", "y", *FEATURES, *MOTIFS)}
        for row in members:
            age.append(float(row["age"]) - means["age"])
            y.append(float(row["y"]) - means["y"])
            physical.append([float(row[name]) - means[name] for name in FEATURES])
            motifs.append([float(row[name]) - means[name] for name in MOTIFS])
    return tuple(np.asarray(value) for value in (age, y, physical, motifs))


def cross_size_subspaces(old_rows):
    arrays = {key: centered_arrays(value) for key, value in old_rows.items()}
    return {
        "N325": train_subspace(arrays, [("N425", o) for o in ORIENTATIONS], 1),
        "N425": train_subspace(arrays, [("N325", o) for o in ORIENTATIONS], 1),
    }


def fit_nested(rows, subspace):
    age, y, physical, motifs = centered(rows)
    canonical = (physical / subspace["scales"]) @ subspace["canonical_basis"]
    rank1 = canonical @ subspace["directions"].T
    base = np.column_stack((age, rank1))
    full = np.column_stack((age, rank1, motifs))
    beta0, *_ = np.linalg.lstsq(base, y, rcond=None)
    beta1, *_ = np.linalg.lstsq(full, y, rcond=None)
    r0, r1 = y - base @ beta0, y - full @ beta1
    total = float(y @ y)
    return {
        "baseline_beta_age": float(beta0[0]),
        "motif_beta_age": float(beta1[0]),
        "delta_beta_age": float(beta1[0] - beta0[0]),
        "baseline_R2": 1.0 - float(r0 @ r0) / total,
        "motif_R2": 1.0 - float(r1 @ r1) / total,
        "incremental_R2": float((r0 @ r0 - r1 @ r1) / total),
        "motif_coefficients": dict(zip(MOTIFS, beta1[2:].tolist())),
        "rows": len(rows),
    }


def old_anchor(old_transport):
    order, values = [], []
    matrix = np.zeros((4, 4))
    for si, (size, regime) in enumerate(
            (("N325", "N425_to_N325"), ("N425", "N325_to_N425"))):
        block = old_transport["transport_scores"][regime]["rank1"]["production_anchor"]
        for orientation, value in zip(ORIENTATIONS, block["age_slopes"]):
            order.append((size, orientation)); values.append(float(value))
        matrix[2 * si:2 * si + 2, 2 * si:2 * si + 2] = block["covariance"]
    return order, values, matrix.tolist()


def score(fresh, subspaces, old_transport, batches=20):
    order, anchor, anchor_cov = old_anchor(old_transport)

    def vector(omit_size=None, omit_batch=None):
        values, fits = [], {}
        for size, orientation in order:
            rows = [row for row in fresh[(size, orientation)]
                    if not (size == omit_size and row["batch"] == omit_batch)]
            fit = fit_nested(rows, subspaces[size])
            values.append(fit["delta_beta_age"])
            fits[f"{size}_{orientation}"] = fit
        return values, fits

    delta, fits = vector()
    delta_cov = np.zeros((4, 4))
    for size in SIZES:
        delta_cov += np.asarray(covariance(
            [vector(size, batch)[0] for batch in range(batches)]))
    adjusted = [anchor[i] + delta[i] for i in range(4)]
    combined = add_matrix(anchor_cov, delta_cov.tolist())
    by_size = {}
    for si, size in enumerate(SIZES):
        pos = [2 * si, 2 * si + 1]
        cov = [[combined[i][j] for j in pos] for i in pos]
        values = [adjusted[i] for i in pos]
        stat, df, p = chi_square(values, cov)
        by_size[size] = {"age_slopes": values, "covariance": cov,
                         "chi2": stat, "df": df, "p": p}
    retention = [abs(adjusted[i] / anchor[i]) for i in range(4)]
    absorbed = all(by_size[s]["p"] >= 0.01 for s in SIZES) and max(retention) <= 0.25
    coherent_hint = max(retention) < 0.50
    return {
        "vector_order": [list(key) for key in order], "fresh_fits": fits,
        "old_rank1_anchor_age_slopes": anchor,
        "old_rank1_anchor_covariance": anchor_cov,
        "fresh_paired_delta": delta,
        "fresh_paired_delta_covariance": delta_cov.tolist(),
        "local_motif_anchor_age_slopes": adjusted,
        "combined_covariance": combined,
        "absolute_retention_vs_old_rank1": retention,
        "size_joint": by_size,
        "decision": "local_motif_absorbs_remaining_rank1_age" if absorbed else
                    "local_motif_fails_common_absorption_of_remaining_rank1_age",
        "100k_extension_precondition_partial": coherent_hint,
    }


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
    if sha256(Path(__file__).parents[1] / "src/threshold_rank_integer_period_mc.cpp") != \
            freeze["runner_source_sha256"]:
        raise ValueError("runner source hash differs from frozen source")

    old_transport = json.loads(args.old_transport_score.read_text())
    old_lock = json.loads(args.old_raw_lock.read_text())
    old_rows = {}
    for size in SIZES:
        loaded = load_rows(Path(old_lock["runs"][size]["csv"]))
        for orientation in ORIENTATIONS:
            old_rows[(size, orientation)] = [r for r in loaded
                                             if r["orientation"] == orientation]
    subspaces = cross_size_subspaces(old_rows)

    fresh, inputs = {}, {}
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
            fresh[(size, orientation)] = [r for r in loaded
                                          if r["orientation"] == orientation]
        inputs[size] = {"csv_sha256": sha256(csv_path),
                        "metadata_sha256": sha256(metadata_path),
                        "rows": len(loaded), "checks": checks}

    primary = score(fresh, subspaces, old_transport)
    payload = {
        "schema": "matching-one/p334-local-active-boundary-motif-score/v1",
        "freeze_sha256": sha256(args.freeze), "runner_commit": args.runner_commit,
        "inputs": inputs, "typed_coordinates": list(MOTIFS),
        "age_primary": primary, "claim_boundary": freeze["claim_boundary"],
        "extension_rule": freeze["extension_rule"],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"decision": primary["decision"],
                      "retention": primary["absolute_retention_vs_old_rank1"],
                      "joint_p": {s: primary["size_joint"][s]["p"] for s in SIZES}},
                     indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Post-reveal cross-scale transport of the P334 morphology state."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Mapping, Sequence

import mpmath as mp
import numpy as np

from score_p334_current_k0_geometry_pilot import add_matrix, covariance, sha256


ORIENTATIONS = ("first", "second")
SIZES = ("N325", "N425")
FEATURES = (
    "area_size", "carriers", "occupied_frontier_L", "vacant_frontier_L",
    "boundary_cut_L", "boundary_multicontact_L", "boundary_pairs_L",
    "core_vertices_area", "core_edges_area", "articulations_L", "bridges_L",
)
GROUPS = {
    "mass_core": ("area_size", "core_vertices_area", "core_edges_area"),
    "boundary": ("occupied_frontier_L", "vacant_frontier_L", "boundary_cut_L",
                 "boundary_multicontact_L", "boundary_pairs_L"),
    "bottleneck": ("carriers", "articulations_L", "bridges_L"),
}


def chi_square(values: Sequence[float], cov: Sequence[Sequence[float]]):
    vector = np.asarray(values, dtype=float)
    matrix = np.asarray(cov, dtype=float)
    statistic = float(vector @ np.linalg.pinv(matrix, rcond=1e-12) @ vector)
    rank = int(np.linalg.matrix_rank(matrix, tol=max(matrix.shape) *
                                     np.linalg.svd(matrix, compute_uv=False)[0] * 1e-12))
    p = float(mp.gammainc(mp.mpf(rank) / 2, mp.mpf(statistic) / 2,
                          mp.inf, regularized=True))
    return statistic, rank, p


def load_rows(path: Path):
    rows = []
    with path.open(newline="") as stream:
        for raw in csv.DictReader(stream):
            n = int(raw["n"])
            length = math.sqrt(n)
            rows.append({
                "size": f"N{n}", "orientation": raw["orientation"],
                "batch": int(raw["batch"]), "line": (int(raw["ell_u"]), int(raw["ell_v"])),
                "age": int(raw["age_steps"]) / n, "y": float(raw["next_exit"]),
                "area_size": int(raw["essential_size"]) / n,
                "carriers": float(int(raw["essential_carriers"]) - 1),
                "occupied_frontier_L": int(raw["occupied_frontier"]) / length,
                "vacant_frontier_L": int(raw["vacant_frontier"]) / length,
                "boundary_cut_L": int(raw["boundary_cut_edges"]) / length,
                "boundary_multicontact_L": int(raw["boundary_multicontact_sites"]) / length,
                "boundary_pairs_L": int(raw["boundary_contact_pairs"]) / length,
                "core_vertices_area": int(raw["core_vertices"]) / n,
                "core_edges_area": int(raw["core_edges"]) / n,
                "articulations_L": int(raw["articulation_vertices"]) / length,
                "bridges_L": int(raw["bridges"]) / length,
            })
    return rows


def centered_arrays(rows: Sequence[Mapping[str, object]], omitted_batch: int | None = None):
    kept = [row for row in rows if row["batch"] != omitted_batch]
    groups: dict[tuple[int, int], list[Mapping[str, object]]] = {}
    for row in kept:
        groups.setdefault(row["line"], []).append(row)
    age, y, x = [], [], []
    for members in groups.values():
        age_mean = sum(float(row["age"]) for row in members) / len(members)
        y_mean = sum(float(row["y"]) for row in members) / len(members)
        feature_means = [sum(float(row[name]) for row in members) / len(members)
                         for name in FEATURES]
        for row in members:
            age.append(float(row["age"]) - age_mean)
            y.append(float(row["y"]) - y_mean)
            x.append([float(row[name]) - feature_means[index]
                      for index, name in enumerate(FEATURES)])
    return np.asarray(age), np.asarray(y), np.asarray(x)


def fit_age_shape(age, y, scores=None):
    design = age[:, None] if scores is None else np.column_stack((age, scores))
    beta, _, _, singular = np.linalg.lstsq(design, y, rcond=None)
    residual = y - design @ beta
    return {
        "beta_age": float(beta[0]),
        "shape_amplitudes": [] if scores is None else beta[1:].tolist(),
        "rss": float(residual @ residual),
        "rank": int(np.linalg.matrix_rank(design)),
        "singular_values": singular.tolist(),
    }


def train_subspace(arrays, train_keys, rank):
    pooled = np.vstack([arrays[key][2] for key in train_keys])
    scales = np.sqrt(np.mean(pooled * pooled, axis=0))
    if np.any(scales <= 1e-14):
        raise ValueError("zero-information physical morphology column")
    standardized = pooled / scales
    _, singular, vh = np.linalg.svd(standardized, full_matrices=False)
    active = singular > singular[0] * max(standardized.shape) * np.finfo(float).eps
    basis = vh[active].T
    coefficients = []
    for key in train_keys:
        age, y, x = arrays[key]
        z = (x / scales) @ basis
        fit = fit_age_shape(age, y, z)
        coefficients.append(fit["shape_amplitudes"])
    coefficient_matrix = np.asarray(coefficients)
    _, predictive_singular, predictive_vh = np.linalg.svd(coefficient_matrix,
                                                           full_matrices=False)
    if rank > predictive_vh.shape[0]:
        raise ValueError("requested predictive rank exceeds training environments")
    return {
        "scales": scales, "canonical_basis": basis,
        "directions": predictive_vh[:rank],
        "morphology_rank": int(active.sum()),
        "morphology_singular_values": singular.tolist(),
        "predictive_singular_values": predictive_singular.tolist(),
        "rank_energy": float(np.sum(predictive_singular[:rank] ** 2) /
                             np.sum(predictive_singular ** 2)),
    }


def regime_point(all_rows, train_keys, target_keys, rank, omitted=None):
    omitted = omitted or {}
    arrays = {}
    for key in set(train_keys) | set(target_keys):
        arrays[key] = centered_arrays(all_rows[key], omitted.get(key[0]))
    subspace = train_subspace(arrays, train_keys, rank)
    targets = []
    for key in target_keys:
        age, y, x = arrays[key]
        z = (x / subspace["scales"]) @ subspace["canonical_basis"]
        scores = z @ subspace["directions"].T
        base = fit_age_shape(age, y)
        transported = fit_age_shape(age, y, scores)
        targets.append({
            "environment": f"{key[0]}_{key[1]}",
            "baseline_beta_age": base["beta_age"],
            "transport_beta_age": transported["beta_age"],
            "delta_beta_age": transported["beta_age"] - base["beta_age"],
            "baseline_rss": base["rss"], "transport_rss": transported["rss"],
            "rss_gain_fraction": (base["rss"] - transported["rss"]) / base["rss"],
            "shape_amplitudes": transported["shape_amplitudes"],
        })
    return targets, subspace


def production_block(production, target_keys):
    values, cov = [], np.zeros((len(target_keys), len(target_keys)))
    positions = []
    for size, orientation in target_keys:
        block = production["sizes"][size]
        position = block["point_vector_order"].index(f"{orientation}_primary")
        values.append(float(block["point_vector"][position]))
        positions.append((size, position))
    for i, (size_i, index_i) in enumerate(positions):
        for j, (size_j, index_j) in enumerate(positions):
            if size_i == size_j:
                cov[i, j] = production["sizes"][size_i]["delete_one_covariance"][index_i][index_j]
    return values, cov.tolist()


def score_regime(all_rows, production, train_keys, target_keys, rank, batches=20):
    point, subspace = regime_point(all_rows, train_keys, target_keys, rank)
    delta = [entry["delta_beta_age"] for entry in point]
    delta_cov = np.zeros((len(target_keys), len(target_keys)))
    for size in SIZES:
        deleted = []
        for batch in range(batches):
            estimate, _ = regime_point(all_rows, train_keys, target_keys, rank,
                                       omitted={size: batch})
            deleted.append([entry["delta_beta_age"] for entry in estimate])
        delta_cov += np.asarray(covariance(deleted))
    production_values, production_cov = production_block(production, target_keys)
    anchored = [production_values[index] + delta[index]
                for index in range(len(target_keys))]
    combined = add_matrix(production_cov, delta_cov.tolist())
    statistic, df, p = chi_square(anchored, combined)
    base_statistic, base_df, base_p = chi_square(production_values, production_cov)
    return {
        "train": [f"{size}_{orientation}" for size, orientation in train_keys],
        "target": [f"{size}_{orientation}" for size, orientation in target_keys],
        "rank": rank,
        "training_subspace": {
            "physical_morphology_rank": subspace["morphology_rank"],
            "predictive_singular_values": subspace["predictive_singular_values"],
            "rank_energy": subspace["rank_energy"],
        },
        "fresh_target_fits": point,
        "paired_delta_covariance": delta_cov.tolist(),
        "production_baseline": {
            "age_slopes": production_values, "covariance": production_cov,
            "chi2": base_statistic, "df": base_df, "p": base_p,
        },
        "production_anchor": {
            "age_slopes": anchored, "covariance": combined,
            "absolute_retention": [abs(anchored[index] / production_values[index])
                                   for index in range(len(target_keys))],
            "chi2": statistic, "df": df, "p": p,
        },
    }


def global_rotation(all_rows, omitted=None):
    omitted = omitted or {}
    arrays = {key: centered_arrays(rows, omitted.get(key[0]))
              for key, rows in all_rows.items()}
    keys = [(size, orientation) for size in SIZES for orientation in ORIENTATIONS]
    pooled = np.vstack([arrays[key][2] for key in keys])
    scales = np.sqrt(np.mean(pooled * pooled, axis=0))
    standardized = pooled / scales
    _, singular, vh = np.linalg.svd(standardized, full_matrices=False)
    active = singular > singular[0] * max(standardized.shape) * np.finfo(float).eps
    basis = vh[active].T
    beta = {}
    for key in keys:
        age, y, x = arrays[key]
        z = (x / scales) @ basis
        beta[key] = np.asarray(fit_age_shape(age, y, z)["shape_amplitudes"])
    matrix = np.vstack([beta[key] for key in keys])
    predictive_singular = np.linalg.svd(matrix, compute_uv=False)

    def cosine(left, right):
        return float(left @ right / math.sqrt((left @ left) * (right @ right)))

    cosines = {
        "scale_first": cosine(beta[("N325", "first")], beta[("N425", "first")]),
        "scale_second": cosine(beta[("N325", "second")], beta[("N425", "second")]),
        "orientation_N325": cosine(beta[("N325", "first")], beta[("N325", "second")]),
        "orientation_N425": cosine(beta[("N425", "first")], beta[("N425", "second")]),
    }

    def group_energy(pairs):
        original_differences = []
        for left, right in pairs:
            original_differences.append(basis @ (beta[left] - beta[right]))
        total = sum(float(value @ value) for value in original_differences)
        shares = {}
        for group, names in GROUPS.items():
            indices = [FEATURES.index(name) for name in names]
            shares[group] = sum(float(np.sum(value[indices] ** 2))
                                for value in original_differences) / total
        return shares

    return {
        "predictive_singular_values": predictive_singular.tolist(),
        "rank1_energy": float(predictive_singular[0] ** 2 /
                              np.sum(predictive_singular ** 2)),
        "rank2_energy": float(np.sum(predictive_singular[:2] ** 2) /
                              np.sum(predictive_singular ** 2)),
        "coefficient_cosines": cosines,
        "scale_rotation_group_energy": group_energy([
            (("N325", "first"), ("N425", "first")),
            (("N325", "second"), ("N425", "second")),
        ]),
        "orientation_rotation_group_energy": group_energy([
            (("N325", "first"), ("N325", "second")),
            (("N425", "first"), ("N425", "second")),
        ]),
    }


def rotation_with_uncertainty(all_rows, batches=20):
    point = global_rotation(all_rows)
    names = ["rank1_energy", "rank2_energy"] + [
        f"cos_{name}" for name in point["coefficient_cosines"]
    ] + [f"scale_{name}" for name in GROUPS] + [f"orientation_{name}" for name in GROUPS]

    def vector(value):
        return ([value["rank1_energy"], value["rank2_energy"]] +
                [value["coefficient_cosines"][name]
                 for name in point["coefficient_cosines"]] +
                [value["scale_rotation_group_energy"][name] for name in GROUPS] +
                [value["orientation_rotation_group_energy"][name] for name in GROUPS])

    cov = np.zeros((len(names), len(names)))
    for size in SIZES:
        deleted = [vector(global_rotation(all_rows, {size: batch}))
                   for batch in range(batches)]
        cov += np.asarray(covariance(deleted))
    point["diagnostic_vector_order"] = names
    point["diagnostic_vector"] = vector(point)
    point["delete_one_covariance"] = cov.tolist()
    point["standard_errors"] = dict(zip(names, np.sqrt(np.maximum(np.diag(cov), 0.0))))
    return point


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-lock", type=Path, required=True)
    parser.add_argument("--production-age-score", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    lock = json.loads(args.raw_lock.read_text())
    production = json.loads(args.production_age_score.read_text())
    if sha256(args.production_age_score) != "43901ffdc597328207cc31a4e1307a2f801011652d443ed24abed8e7bd6ad6eb":
        raise ValueError("production score hash changed")
    all_rows = {}
    inputs = {}
    for size in SIZES:
        path = Path(lock["runs"][size]["csv"])
        if sha256(path) != lock["runs"][size]["csv_sha256"]:
            raise ValueError(f"{size} raw hash changed")
        loaded = load_rows(path)
        for orientation in ORIENTATIONS:
            all_rows[(size, orientation)] = [row for row in loaded
                                             if row["orientation"] == orientation]
        inputs[size] = {"path": str(path), "sha256": sha256(path), "rows": len(loaded)}

    regimes = {
        "N325_to_N425": ([('N325', 'first'), ('N325', 'second')],
                          [('N425', 'first'), ('N425', 'second')]),
        "N425_to_N325": ([('N425', 'first'), ('N425', 'second')],
                          [('N325', 'first'), ('N325', 'second')]),
        "first_to_second": ([('N325', 'first'), ('N425', 'first')],
                             [('N325', 'second'), ('N425', 'second')]),
        "second_to_first": ([('N325', 'second'), ('N425', 'second')],
                             [('N325', 'first'), ('N425', 'first')]),
    }
    scores = {}
    for name, (train, target) in regimes.items():
        scores[name] = {f"rank{rank}": score_regime(all_rows, production, train, target, rank)
                        for rank in (1, 2)}
    payload = {
        "schema": "matching-one/p334-morphology-state-transport/v1",
        "status": "post_reveal_zero_new_sample_mechanism_exploration",
        "physical_coordinates": {
            "area_like_divisor": "N",
            "boundary_or_bottleneck_divisor": "L=sqrt(N)",
            "dimensionless": ["essential_carriers-1"],
            "feature_order": list(FEATURES),
        },
        "inputs": inputs,
        "production_age_score_sha256": sha256(args.production_age_score),
        "transport_scores": scores,
        "rotation_diagnostic": rotation_with_uncertainty(all_rows),
        "interpretation_boundary": "Held-out means the morphology subspace is learned without target outcomes; target amplitudes are then fitted inside that fixed rank-one or rank-two subspace. This tests subspace transport, not a point forecast. The analysis is post-reveal and does not add an independent evidence block.",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    compact = {name: {rank: value["production_anchor"]
                      for rank, value in ranks.items()}
               for name, ranks in scores.items()}
    print(json.dumps({"transport": compact,
                      "rotation": payload["rotation_diagnostic"]}, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Post-reveal P334 production anchor and temporal-mode geometry description."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.stats import chi2, t


ORIENTATIONS = ("first", "second")
CHEAP = ("g_size", "g_carriers", "g_occupied_frontier", "g_vacant_frontier")


def sha256(path: Path) -> str:
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    return digest


def load_rows(path: Path) -> list[dict[str, float]]:
    rows = []
    with path.open(newline="") as stream:
        for raw in csv.DictReader(stream):
            n, k0 = int(raw["n"]), int(raw["k0"])
            rows.append({
                "orientation": raw["orientation"],
                "batch": int(raw["batch"]),
                "ell_u": int(raw["ell_u"]),
                "ell_v": int(raw["ell_v"]),
                "k1": int(raw["k1"]),
                "k2": int(raw["k2"]),
                "g_size": int(raw["essential_size"]) / n,
                "g_carriers": int(raw["essential_carriers"]) - 1,
                "g_occupied_frontier": int(raw["occupied_frontier"]) / k0,
                "g_vacant_frontier": int(raw["vacant_frontier"]) / (n - k0),
            })
    return rows


def upper_matrix(vector: list[float]) -> np.ndarray:
    dimension = int((np.sqrt(8 * len(vector) + 1) - 1) / 2)
    if dimension * (dimension + 1) // 2 != len(vector):
        raise ValueError("upper-triangle vector has nontriangular length")
    output = np.zeros((dimension, dimension))
    index = 0
    for i in range(dimension):
        for j in range(i, dimension):
            output[i, j] = output[j, i] = vector[index]
            index += 1
    return output


def line_center(rows: list[dict[str, float]], values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    geometry = np.asarray([[row[name] for name in CHEAP] for row in rows], dtype=float)
    centered_y = values.copy()
    centered_x = geometry.copy()
    groups: dict[tuple[int, int], list[int]] = {}
    for index, row in enumerate(rows):
        groups.setdefault((int(row["ell_u"]), int(row["ell_v"])), []).append(index)
    for indices in groups.values():
        centered_y[indices] -= np.mean(centered_y[indices])
        centered_x[indices] -= np.mean(centered_x[indices], axis=0)
    return centered_x, centered_y


def explain(rows: list[dict[str, float]], values: np.ndarray) -> dict[str, object]:
    x, y = line_center(rows, values)
    active = np.sum(x * x, axis=0) > 1e-14
    coefficients, *_ = np.linalg.lstsq(x[:, active], y, rcond=None)
    fitted = x[:, active] @ coefficients
    total = float(y @ y)
    residual = float((y - fitted) @ (y - fitted))
    covariance = (x.T @ y / len(y)).tolist()
    standard = np.sqrt(np.sum(x * x, axis=0) * total)
    correlations = np.divide(x.T @ y, standard, out=np.zeros(len(CHEAP)), where=standard > 0)
    return {
        "rows": len(rows),
        "line_centered_mode_variance": total / len(y),
        "cheap_geometry_R2": 1.0 - residual / total,
        "cheap_coefficients": {name: (float(coefficients[list(np.flatnonzero(active)).index(i)])
                                              if active[i] else 0.0)
                               for i, name in enumerate(CHEAP)},
        "mode_geometry_covariance": dict(zip(CHEAP, covariance)),
        "mode_geometry_correlation": dict(zip(CHEAP, correlations.tolist())),
        "dropped_zero_information": [name for name, keep in zip(CHEAP, active) if not keep],
    }


def temporal_description(rows: list[dict[str, float]], temporal_size: dict[str, object]) -> dict[str, object]:
    dimension = len(temporal_size["layers"])
    upper_count = dimension * (dimension + 1) // 2
    kernel_vector = temporal_size["kernel_vector"]
    output = {}
    r2_order = []
    deleted = []
    mode_vectors = {}
    for o_index, orientation in enumerate(ORIENTATIONS):
        kernel = upper_matrix(kernel_vector[o_index * upper_count:(o_index + 1) * upper_count])
        eigenvalues, eigenvectors = np.linalg.eigh(kernel)
        order = np.argsort(eigenvalues)[::-1]
        eigenvalues, eigenvectors = eigenvalues[order], eigenvectors[:, order]
        mode_vectors[orientation] = eigenvectors
        subset = [row for row in rows if row["orientation"] == orientation]
        ranks = np.asarray([[int(row["k1"] <= layer) + int(row["k2"] <= layer)
                             for layer in temporal_size["layers"]] for row in subset], dtype=float)
        output[orientation] = {"production_kernel_eigenvalues": eigenvalues.tolist()}
        for mode in (1, 2):
            score = ranks @ eigenvectors[:, mode]
            item = explain(subset, score)
            item["production_mode_index_one_based"] = mode + 1
            output[orientation][f"mode{mode + 1}"] = item
            r2_order.append((orientation, mode))
    for omitted in range(50):
        values = []
        for orientation, mode in r2_order:
            subset = [row for row in rows
                      if row["orientation"] == orientation and row["batch"] != omitted]
            ranks = np.asarray([[int(row["k1"] <= layer) + int(row["k2"] <= layer)
                                 for layer in temporal_size["layers"]] for row in subset], dtype=float)
            values.append(explain(subset, ranks @ mode_vectors[orientation][:, mode])["cheap_geometry_R2"])
        deleted.append(values)
    deleted = np.asarray(deleted)
    centered = deleted - np.mean(deleted, axis=0)
    covariance = 49 / 50 * centered.T @ centered
    for index, (orientation, mode) in enumerate(r2_order):
        output[orientation][f"mode{mode + 1}"]["R2_jackknife_se"] = float(
            np.sqrt(covariance[index, index]))
    output["R2_vector_order"] = [[orientation, mode + 1] for orientation, mode in r2_order]
    output["R2_jackknife_covariance"] = covariance.tolist()
    return output


def production_anchor(pilot: dict[str, object], production: dict[str, object]) -> dict[str, object]:
    output = {}
    for size in ("N325", "N425"):
        pilot_size = pilot["scores"][size]
        order = [tuple(value) for value in pilot_size["vector_order"]]
        pilot_cov = np.asarray(pilot_size["jackknife_covariance"], dtype=float)
        production_size = production["sizes"][size]
        production_order = production_size["point_vector_order"]
        production_covariance = np.asarray(production_size["delete_one_covariance"], dtype=float)
        values, delta_values = [], []
        delta_covariance = np.zeros((2, 2))
        for oi, orientation in enumerate(ORIENTATIONS):
            primary = order.index((orientation, "M0"))
            cheap = order.index((orientation, "Mcheap"))
            delta_values.append(
                pilot_size["fits"][orientation]["Mcheap"]["beta_age"] -
                pilot_size["fits"][orientation]["M0"]["beta_age"])
            for oj, other in enumerate(ORIENTATIONS):
                other_primary = order.index((other, "M0"))
                other_cheap = order.index((other, "Mcheap"))
                delta_covariance[oi, oj] = (
                    pilot_cov[cheap, other_cheap] - pilot_cov[cheap, other_primary] -
                    pilot_cov[primary, other_cheap] + pilot_cov[primary, other_primary])
        production_indices = [production_order.index(f"{orientation}_primary")
                              for orientation in ORIENTATIONS]
        production_values = np.asarray([production_size["point_vector"][index]
                                        for index in production_indices])
        production_block = production_covariance[np.ix_(production_indices, production_indices)]
        adjusted = production_values + np.asarray(delta_values)
        combined_covariance = production_block + delta_covariance
        inverse = np.linalg.inv(combined_covariance)
        joint_chi2 = float(adjusted @ inverse @ adjusted)
        for index, orientation in enumerate(ORIENTATIONS):
            prod_var = production_block[index, index]
            delta_var = delta_covariance[index, index]
            variance = prod_var + delta_var
            df = variance ** 2 / (prod_var ** 2 / 99 + delta_var ** 2 / 49)
            se = np.sqrt(variance)
            values.append({
                "orientation": orientation,
                "production_beta_M0": float(production_values[index]),
                "pilot_paired_delta_Mcheap_minus_M0": float(delta_values[index]),
                "production_anchored_beta_Mcheap": float(adjusted[index]),
                "se": float(se),
                "welch_satterthwaite_df": float(df),
                "two_sided_p": float(2 * t.sf(abs(adjusted[index] / se), df)),
                "absolute_retention_vs_production_M0": float(
                    abs(adjusted[index] / production_values[index])),
            })
        output[size] = {
            "orientations": values,
            "combined_covariance": combined_covariance.tolist(),
            "joint_chi_square": joint_chi2,
            "joint_df": 2,
            "joint_p": float(chi2.sf(joint_chi2, 2)),
            "boundary": "Post-reveal independent production anchor: it replaces the noisy 50k M0 point by the disjoint 2M production M0 and adds only the paired 50k control shift. It is not the frozen primary score."
        }
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pilot-score", type=Path, required=True)
    parser.add_argument("--production-age-score", type=Path, required=True)
    parser.add_argument("--temporal-score", type=Path, required=True)
    parser.add_argument("--n325-csv", type=Path, required=True)
    parser.add_argument("--n425-csv", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    pilot = json.loads(args.pilot_score.read_text())
    production = json.loads(args.production_age_score.read_text())
    temporal = json.loads(args.temporal_score.read_text())
    temporal_by_size = {f"N{value['N']}": value for value in temporal["sizes"]}
    payload = {
        "schema": "matching-one/p334-current-geometry-post-reveal-description/v1",
        "inputs": {
            "pilot_score_sha256": sha256(args.pilot_score),
            "production_age_score_sha256": sha256(args.production_age_score),
            "temporal_score_sha256": sha256(args.temporal_score),
            "temporal_result_commit": "5a7f2d9",
        },
        "production_anchored_age_attenuation": production_anchor(pilot, production),
        "temporal_modes": {
            "N325": temporal_description(load_rows(args.n325_csv), temporal_by_size["N325"]),
            "N425": temporal_description(load_rows(args.n425_csv), temporal_by_size["N425"]),
        },
        "boundary": "Post-reveal descriptive reuse of the fresh pilot. It neither changes the frozen age decision nor turns temporal modes into exact states or memory variables."
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "production_anchor": {size: {
            "retentions": [row["absolute_retention_vs_production_M0"]
                           for row in payload["production_anchored_age_attenuation"][size]["orientations"]],
            "joint_p": payload["production_anchored_age_attenuation"][size]["joint_p"],
        } for size in ("N325", "N425")},
        "temporal_R2": {size: {orientation: {
            mode: payload["temporal_modes"][size][orientation][mode]["cheap_geometry_R2"]
            for mode in ("mode2", "mode3")}
            for orientation in ORIENTATIONS} for size in ("N325", "N425")},
    }, indent=2))


if __name__ == "__main__":
    main()

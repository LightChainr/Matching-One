#!/usr/bin/env python3
"""Score the frozen P334 current-k0 geometry pilot.

The scorer uses only small cross-products after exact line fixed-effect
centering.  It deliberately keeps H2 separate from the cheap geometry model:
H2/(N-k0) is the conditional one-step hazard by construction.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Mapping, Sequence, Tuple

import mpmath as mp


CHEAP = ("g_size", "g_carriers", "g_occupied_frontier", "g_vacant_frontier")
TYPED = ("h2_theta_rate", "h2_figure8_rate", "h2_separate_rate")
MODEL_PREDICTORS = {
    "M0": ("age",),
    "Mcheap": ("age",) + CHEAP,
    "Mtyped": ("age",) + CHEAP + TYPED,
    "MH2": ("age", "h2_rate"),
}
ORIENTATIONS = ("first", "second")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def solve(matrix: Sequence[Sequence[float]], vector: Sequence[float]) -> List[float]:
    n = len(vector)
    work = [list(matrix[row]) + [float(vector[row])] for row in range(n)]
    scale = max((abs(value) for row in matrix for value in row), default=1.0)
    for column in range(n):
        pivot = max(range(column, n), key=lambda row: abs(work[row][column]))
        if abs(work[pivot][column]) <= 1e-12 * max(scale, 1.0):
            raise ValueError("singular centered regression system")
        work[column], work[pivot] = work[pivot], work[column]
        divisor = work[column][column]
        for index in range(column, n + 1):
            work[column][index] /= divisor
        for row in range(n):
            if row == column:
                continue
            factor = work[row][column]
            for index in range(column, n + 1):
                work[row][index] -= factor * work[column][index]
    return [work[row][n] for row in range(n)]


def centered_fit(
    rows: Sequence[Mapping[str, float]],
    predictors: Sequence[str],
    outcome: Callable[[Mapping[str, float]], float] = lambda row: row["y"],
    stratum: Callable[[Mapping[str, float]], Tuple[object, ...]] =
        lambda row: (row["ell_u"], row["ell_v"]),
) -> Dict[str, object]:
    names = list(predictors)
    groups: Dict[Tuple[object, ...], List[Mapping[str, float]]] = {}
    for row in rows:
        groups.setdefault(stratum(row), []).append(row)
    centered: List[Tuple[List[float], float]] = []
    for members in groups.values():
        means = [sum(float(row[name]) for row in members) / len(members) for name in names]
        ybar = sum(outcome(row) for row in members) / len(members)
        for row in members:
            centered.append((
                [float(row[name]) - means[index] for index, name in enumerate(names)],
                outcome(row) - ybar,
            ))
    diagonal = [sum(x[index] * x[index] for x, _ in centered)
                for index in range(len(names))]
    active = [index for index, value in enumerate(diagonal) if value > 1e-12]
    if 0 not in active:
        raise ValueError("age is unidentifiable after line centering")
    matrix = [[sum(x[i] * x[j] for x, _ in centered) for j in active] for i in active]
    vector = [sum(x[i] * y for x, y in centered) for i in active]
    coefficients = solve(matrix, vector)
    by_name = {name: 0.0 for name in names}
    for index, value in zip(active, coefficients):
        by_name[names[index]] = value
    residual_ss = 0.0
    for x, y in centered:
        fitted = sum(by_name[name] * x[index] for index, name in enumerate(names))
        residual_ss += (y - fitted) ** 2
    return {
        "coefficients": by_name,
        "dropped_zero_information": [names[index] for index in range(len(names))
                                     if index not in active],
        "rows": len(rows),
        "strata": len(groups),
        "residual_ss": residual_ss,
        "age_information": diagonal[0],
    }


def covariance(deleted: Sequence[Sequence[float]]) -> List[List[float]]:
    count = len(deleted)
    width = len(deleted[0])
    means = [sum(row[column] for row in deleted) / count for column in range(width)]
    factor = (count - 1) / count
    return [[factor * sum((row[i] - means[i]) * (row[j] - means[j])
                          for row in deleted)
             for j in range(width)] for i in range(width)]


def add_matrix(first: Sequence[Sequence[float]], second: Sequence[Sequence[float]]) -> List[List[float]]:
    return [[first[i][j] + second[i][j] for j in range(len(first))]
            for i in range(len(first))]


def student_p(value: float, se: float, df: int) -> float:
    if se <= 0:
        return 0.0 if value else 1.0
    t2 = (value / se) ** 2
    return float(mp.betainc(df / 2, mp.mpf("0.5"), 0, df / (df + t2), regularized=True))


def chi2_p(value: float, df: int) -> float:
    return float(mp.gammainc(df / 2, value / 2, mp.inf, regularized=True))


def joint_two(values: Sequence[float], cov: Sequence[Sequence[float]]) -> Tuple[float, float]:
    determinant = cov[0][0] * cov[1][1] - cov[0][1] * cov[1][0]
    if determinant <= 0:
        raise ValueError("non-positive two-orientation covariance")
    chi2 = (values[0] ** 2 * cov[1][1] - 2 * values[0] * values[1] * cov[0][1] +
            values[1] ** 2 * cov[0][0]) / determinant
    return chi2, chi2_p(chi2, 2)


def load_rows(path: Path, metadata_path: Path, frozen: Mapping[str, object], commit: str) -> Tuple[List[dict], dict]:
    metadata = json.loads(metadata_path.read_text())
    expected = frozen
    checks = {
        "git_commit": metadata["git_commit"] == commit,
        "samples": metadata["samples_per_pair"] == expected["samples"],
        "batches": metadata["batches"] == expected["batches"],
        "seed": metadata["seed"] == expected["seed"],
        "replica_first": metadata["replica_counter_first"] == expected["replica_counter_first"],
        "replica_last": metadata["replica_counter_last_exclusive"] == expected["replica_counter_last_exclusive"],
        "k0": metadata["geometry_pilot_k0"] == expected["k0"],
    }
    if not all(checks.values()):
        raise ValueError(f"metadata contract failed: {checks}")
    rows: List[dict] = []
    with path.open(newline="") as stream:
        for raw in csv.DictReader(stream):
            n = int(raw["n"])
            k0 = int(raw["k0"])
            remaining = n - k0
            row = {key: float(value) for key, value in raw.items()
                   if key not in ("orientation",)}
            row["orientation"] = raw["orientation"]
            row["batch"] = int(raw["batch"])
            row["ell_u"] = int(raw["ell_u"])
            row["ell_v"] = int(raw["ell_v"])
            row["y"] = int(raw["next_exit"])
            row["age"] = int(raw["age_steps"]) / n
            row["g_size"] = int(raw["essential_size"]) / n
            row["g_carriers"] = int(raw["essential_carriers"]) - 1
            row["g_occupied_frontier"] = int(raw["occupied_frontier"]) / k0
            row["g_vacant_frontier"] = int(raw["vacant_frontier"]) / remaining
            row["h2_rate"] = int(raw["H2"]) / remaining
            row["h2_theta_rate"] = int(raw["H2_theta"]) / remaining
            row["h2_figure8_rate"] = int(raw["H2_figure8"]) / remaining
            row["h2_separate_rate"] = int(raw["H2_separate"]) / remaining
            if int(raw["H2"]) != (int(raw["H2_theta"]) + int(raw["H2_figure8"]) +
                                  int(raw["H2_separate"])):
                raise ValueError("H2 trigger decomposition failed")
            if int(raw["H2"]) != (int(raw["H2_direction_positive"]) +
                                  int(raw["H2_direction_negative"]) +
                                  int(raw["H2_direction_mixed"])):
                raise ValueError("H2 direction decomposition failed")
            rows.append(row)
    return rows, {"metadata_checks": checks, "raw_sha256": sha256(path),
                  "metadata_sha256": sha256(metadata_path)}


def fit_size(rows: Sequence[dict], batches: int) -> Dict[str, object]:
    full: Dict[str, Dict[str, object]] = {}
    order: List[Tuple[str, str]] = []
    for orientation in ORIENTATIONS:
        subset = [row for row in rows if row["orientation"] == orientation]
        full[orientation] = {}
        for model, predictors in MODEL_PREDICTORS.items():
            full[orientation][model] = centered_fit(subset, predictors)
            order.append((orientation, model))
    deleted = []
    for batch in range(batches):
        values = []
        for orientation, model in order:
            subset = [row for row in rows
                      if row["orientation"] == orientation and row["batch"] != batch]
            fit = centered_fit(subset, MODEL_PREDICTORS[model])
            values.append(float(fit["coefficients"]["age"]))
        deleted.append(values)
    cov = covariance(deleted)
    index = {key: position for position, key in enumerate(order)}
    summary: Dict[str, object] = {}
    for orientation in ORIENTATIONS:
        summary[orientation] = {}
        for model in MODEL_PREDICTORS:
            position = index[(orientation, model)]
            fit = full[orientation][model]
            beta = float(fit["coefficients"]["age"])
            se = math.sqrt(max(cov[position][position], 0.0))
            summary[orientation][model] = {
                **fit,
                "beta_age": beta,
                "se": se,
                "p": student_p(beta, se, batches - 1),
            }
        primary = summary[orientation]["M0"]["beta_age"]
        for model in ("Mcheap", "Mtyped", "MH2"):
            summary[orientation][model]["absolute_slope_retention"] = (
                abs(summary[orientation][model]["beta_age"]) / abs(primary)
                if primary else None)
    joint = {}
    for model in MODEL_PREDICTORS:
        positions = [index[(orientation, model)] for orientation in ORIENTATIONS]
        values = [summary[orientation][model]["beta_age"] for orientation in ORIENTATIONS]
        block = [[cov[i][j] for j in positions] for i in positions]
        chi2, p = joint_two(values, block)
        joint[model] = {"chi2": chi2, "df": 2, "p": p, "covariance": block}
    support = {
        orientation: {
            "risk_rows": sum(row["orientation"] == orientation for row in rows),
            "next_exits": sum(row["y"] for row in rows if row["orientation"] == orientation),
            "mean_h2_rate": sum(row["h2_rate"] for row in rows if row["orientation"] == orientation) /
                            sum(row["orientation"] == orientation for row in rows),
            "mean_next_exit": sum(row["y"] for row in rows if row["orientation"] == orientation) /
                              sum(row["orientation"] == orientation for row in rows),
            "essential_carrier_support": sorted({int(row["g_carriers"] + 1) for row in rows
                                                  if row["orientation"] == orientation}),
            "figure8_total": int(sum(row["h2_figure8_rate"] for row in rows
                                     if row["orientation"] == orientation) *
                                 (int(rows[0]["n"]) - int(rows[0]["k0"]))),
            "separate_total": int(sum(row["h2_separate_rate"] for row in rows
                                      if row["orientation"] == orientation) *
                                  (int(rows[0]["n"]) - int(rows[0]["k0"]))),
        } for orientation in ORIENTATIONS
    }
    for value in support.values():
        value["support_gate"] = value["risk_rows"] >= 1000 and value["next_exits"] >= 25
    return {"fits": summary, "joint": joint, "support": support,
            "vector_order": [list(value) for value in order], "jackknife_covariance": cov}


def geometry_gamma(rows: Sequence[dict]) -> Dict[str, object]:
    return centered_fit(rows, CHEAP, stratum=lambda row: (
        row["orientation"], row["ell_u"], row["ell_v"]))


def transferred_slopes(source: Sequence[dict], target: Sequence[dict],
                       source_exclude: int | None = None,
                       target_exclude: int | None = None) -> Tuple[List[float], Dict[str, float]]:
    source_rows = [row for row in source if row["batch"] != source_exclude]
    target_rows = [row for row in target if row["batch"] != target_exclude]
    gamma = geometry_gamma(source_rows)["coefficients"]
    values = []
    for orientation in ORIENTATIONS:
        subset = [row for row in target_rows if row["orientation"] == orientation]
        fit = centered_fit(
            subset, ("age",),
            outcome=lambda row, gamma=gamma: row["y"] -
                sum(float(gamma[name]) * row[name] for name in CHEAP))
        values.append(float(fit["coefficients"]["age"]))
    return values, {name: float(gamma[name]) for name in CHEAP}


def transfer_score(source: Sequence[dict], target: Sequence[dict], batches: int) -> Dict[str, object]:
    point, gamma = transferred_slopes(source, target)
    source_deleted = [transferred_slopes(source, target, source_exclude=batch)[0]
                      for batch in range(batches)]
    target_deleted = [transferred_slopes(source, target, target_exclude=batch)[0]
                      for batch in range(batches)]
    cov = add_matrix(covariance(source_deleted), covariance(target_deleted))
    chi2, p = joint_two(point, cov)
    return {
        "source": "N325",
        "heldout_target": "N425",
        "geometry_coefficients": gamma,
        "target_residual_age_slopes": dict(zip(ORIENTATIONS, point)),
        "covariance_source_plus_target": cov,
        "joint": {"chi2": chi2, "df": 2, "p": p},
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--freeze", type=Path, required=True)
    parser.add_argument("--n325-csv", type=Path, required=True)
    parser.add_argument("--n325-metadata", type=Path, required=True)
    parser.add_argument("--n425-csv", type=Path, required=True)
    parser.add_argument("--n425-metadata", type=Path, required=True)
    parser.add_argument("--runner-commit", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    freeze = json.loads(args.freeze.read_text())
    rows325, provenance325 = load_rows(
        args.n325_csv, args.n325_metadata, freeze["runs"]["N325"], args.runner_commit)
    rows425, provenance425 = load_rows(
        args.n425_csv, args.n425_metadata, freeze["runs"]["N425"], args.runner_commit)
    score325 = fit_size(rows325, int(freeze["runs"]["N325"]["batches"]))
    score425 = fit_size(rows425, int(freeze["runs"]["N425"]["batches"]))
    transfer = transfer_score(rows325, rows425, 50)
    payload = {
        "schema": "matching-one/p334-current-k0-geometry-pilot-score/v1",
        "freeze_sha256": sha256(args.freeze),
        "runner_commit": args.runner_commit,
        "inputs": {"N325": provenance325, "N425": provenance425},
        "scores": {"N325": score325, "N425": score425},
        "cross_size_transfer": transfer,
        "claim_boundary": freeze["interpretation_boundary"],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "N325_Mcheap_retention": [score325["fits"][o]["Mcheap"]["absolute_slope_retention"]
                                  for o in ORIENTATIONS],
        "N425_Mcheap_retention": [score425["fits"][o]["Mcheap"]["absolute_slope_retention"]
                                  for o in ORIENTATIONS],
        "N325_MH2_joint_p": score325["joint"]["MH2"]["p"],
        "N425_MH2_joint_p": score425["joint"]["MH2"]["p"],
        "transfer_joint_p": transfer["joint"]["p"],
        "output": str(args.output),
    }, indent=2))


if __name__ == "__main__":
    main()

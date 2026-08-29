#!/usr/bin/env python3
"""Score the frozen P275 A_top field-identity microcanonical stream.

The finite matching root is solved independently for every geometry and every
delete-one replicate before gamma, frame transport, scaling, or GLS.  The
18-real covariance is block diagonal across independent N seeds and complete
within each same-field three-modulus block.
"""

from __future__ import annotations

import argparse
import csv
from fractions import Fraction
import json
import math
from pathlib import Path
from typing import Iterable, Sequence

import mpmath as mp
import numpy as np
import yaml


MODULUS_ORDER = ("i", "2i", "5i_over_2")
SIZE_ORDER = (50, 130, 170)
COORDINATE_ORDER = tuple(
    f"N{n}:{modulus}:{component}"
    for n in SIZE_ORDER
    for modulus in MODULUS_ORDER
    for component in ("Re", "Im")
)
REQUIRED_COLUMNS = {
    "n", "modulus", "batch", "k", "counter_first", "counter_last_exclusive",
    "samples", "sum_q", "sum_q2", "sum_I01", "sum_I12", "sum_I02",
    "sum_Re_J_S4", "sum_Im_J_S4", "sum_Re_J_D4", "sum_Im_J_D4",
    "sum_q_Re_J_S4", "sum_q_Im_J_S4", "sum_q_Re_J_D4",
    "sum_q_Im_J_D4", "sum_birth_mass", "priority_field_digest",
}
NUMERIC_FIELDS = tuple(sorted(REQUIRED_COLUMNS - {
    "modulus", "priority_field_digest",
}))


def _number(value) -> float:
    return float(Fraction(str(value)))


def _load_prediction(path: Path) -> dict:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not payload.get("production_authorized"):
        raise ValueError("prediction does not authorize production scoring")
    phase = payload.get("phase1_microcanonical_matching_root", {})
    if phase.get("status") != "frozen_prereveal_production_authorized":
        raise ValueError("Phase 1 is not a frozen authorized prereveal")
    if phase.get("evaluation_p") != "finite_matching_root_inside_each_delete_one":
        raise ValueError("evaluation-p contract changed")
    return payload


def _read_run(spec: str) -> dict:
    parts = spec.split(":", 3)
    if len(parts) != 4:
        raise ValueError("--run must be N:modulus:CSV:METADATA")
    n, modulus, csv_path, metadata_path = int(parts[0]), parts[1], Path(parts[2]), Path(parts[3])
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    with csv_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not REQUIRED_COLUMNS.issubset(reader.fieldnames or ()):  # pragma: no cover
            raise ValueError(f"missing columns in {csv_path}")
        rows = []
        for raw in reader:
            row = {key: float(raw[key]) for key in NUMERIC_FIELDS}
            row["modulus"] = raw["modulus"]
            row["priority_field_digest"] = raw["priority_field_digest"]
            rows.append(row)
    return {
        "N": n, "modulus": modulus, "csv": csv_path, "metadata_path": metadata_path,
        "metadata": metadata, "rows": rows,
    }


def _geometry_map(prediction: dict) -> dict[tuple[int, str], dict]:
    result = {}
    for n_name, block in prediction["geometries"].items():
        n = int(n_name.removeprefix("N"))
        for modulus in MODULUS_ORDER:
            result[(n, modulus)] = block[modulus]
    return result


def validate_runs(runs: Sequence[dict], prediction: dict) -> dict:
    expected = {(n, modulus) for n in SIZE_ORDER for modulus in MODULUS_ORDER}
    incoming = {(run["N"], run["modulus"]) for run in runs}
    if incoming != expected or len(runs) != len(expected):
        raise ValueError(f"run set mismatch: missing={expected-incoming}, extra={incoming-expected}")
    geometry = _geometry_map(prediction)
    phase = prediction["phase1_microcanonical_matching_root"]
    samples = int(phase["samples_per_geometry"])
    batches = int(phase["batches"])
    counter = tuple(int(v) for v in phase["replica_counter"])
    binary_hashes = set()
    commit_hashes = set()
    common_fields = {}
    for run in runs:
        n, modulus, metadata = run["N"], run["modulus"], run["metadata"]
        target = geometry[(n, modulus)]
        checks = {
            "schema": metadata.get("schema") == "matching-one/p275-atop-field-identity-microcanonical/v1",
            "phase": metadata.get("phase") == phase["id"],
            "N": metadata.get("N") == n,
            "modulus": metadata.get("modulus") == modulus,
            "matrix": metadata.get("period_matrix") == target["matrix"],
            "smith": metadata.get("smith_invariants") == target["smith"],
            "z": metadata.get("z") == target["z"],
            "samples": metadata.get("samples") == samples,
            "batches": metadata.get("batches") == batches,
            "seed": metadata.get("seed") == prediction["geometries"][f"N{n}"]["seed"],
            "counter_first": metadata.get("replica_counter_first") == counter[0],
            "counter_last": metadata.get("replica_counter_last_exclusive") == counter[1],
            "root_estimator": metadata.get("root_estimator") == "k*last_active_mark+(N-k)*next_inactive_mark",
        }
        if not all(checks.values()):
            raise ValueError(f"metadata contract failed for N{n}/{modulus}: {checks}")
        binary_hashes.add(metadata["binary_sha256"])
        commit_hashes.add(metadata["git_commit"])
        by_batch: dict[int, set[str]] = {}
        row_counts: dict[tuple[int, int], int] = {}
        for row in run["rows"]:
            batch, k = int(row["batch"]), int(row["k"])
            if int(row["n"]) != n or row["modulus"] != modulus or not 0 <= k <= n:
                raise ValueError(f"row identity failed in {run['csv']}")
            by_batch.setdefault(batch, set()).add(row["priority_field_digest"])
            row_counts[(batch, k)] = row_counts.get((batch, k), 0) + 1
        if set(by_batch) != set(range(batches)) or any(len(v) != 1 for v in by_batch.values()):
            raise ValueError(f"batch/digest structure failed in {run['csv']}")
        if set(row_counts) != {(b, k) for b in range(batches) for k in range(n+1)} or any(
            count != 1 for count in row_counts.values()
        ):
            raise ValueError(f"microcanonical grid failed in {run['csv']}")
        common_fields[(n, modulus)] = {b: next(iter(v)) for b, v in by_batch.items()}
    if len(binary_hashes) != 1 or len(commit_hashes) != 1:
        raise ValueError("all nine geometries must use one runner commit and binary")
    for n in SIZE_ORDER:
        reference = common_fields[(n, "i")]
        for modulus in MODULUS_ORDER[1:]:
            if common_fields[(n, modulus)] != reference:
                raise ValueError(f"N{n} priority field is not byte-identical across moduli")
    return {
        "runner_commit": next(iter(commit_hashes)),
        "binary_sha256": next(iter(binary_hashes)),
        "metadata_contract": "pass",
        "same_N_priority_field_digest": "pass",
    }


def _pool_levels(rows: Sequence[dict], omitted_batch: int | None) -> list[dict[str, float]]:
    n = int(rows[0]["n"])
    fields = [field for field in NUMERIC_FIELDS if field not in {
        "n", "batch", "k", "counter_first", "counter_last_exclusive",
    }]
    pooled = [{field: 0.0 for field in fields} for _ in range(n+1)]
    for row in rows:
        if omitted_batch is not None and int(row["batch"]) == omitted_batch:
            continue
        target = pooled[int(row["k"])]
        for field in fields:
            target[field] += row[field]
    return pooled


def _binomial_weights(n: int, p: float) -> np.ndarray:
    values = np.array([
        math.comb(n, k) * p**k * (1.0-p)**(n-k)
        for k in range(n+1)
    ], dtype=float)
    values /= values.sum()
    return values


def _expectation(levels: Sequence[dict[str, float]], field: str, p: float) -> float:
    weights = _binomial_weights(len(levels)-1, p)
    means = np.array([row[field]/row["samples"] for row in levels], dtype=float)
    return float(weights @ means)


def _matching_root(levels: Sequence[dict[str, float]]) -> float:
    lo, hi = 0.0, 1.0
    flo, fhi = _expectation(levels, "sum_q", lo), _expectation(levels, "sum_q", hi)
    if not flo < 0 < fhi:
        raise ValueError(f"matching root is not bracketed: {flo}, {fhi}")
    for _ in range(90):
        mid = (lo+hi)/2.0
        value = _expectation(levels, "sum_q", mid)
        if value > 0:
            hi = mid
        else:
            lo = mid
    return (lo+hi)/2.0


def geometry_estimate(run: dict, target: dict, omitted_batch: int | None = None) -> dict:
    levels = _pool_levels(run["rows"], omitted_batch)
    p = _matching_root(levels)
    mean_q = _expectation(levels, "sum_q", p)
    birth = _expectation(levels, "sum_birth_mass", p)
    if not birth > 0:
        raise ValueError("birth mass is not positive")
    gamma = []
    s_control = []
    for component in ("Re", "Im"):
        mean_j = _expectation(levels, f"sum_{component}_J_D4", p)
        mean_qj = _expectation(levels, f"sum_q_{component}_J_D4", p)
        gamma.append((mean_qj-mean_q*mean_j)/birth)
        mean_s = _expectation(levels, f"sum_{component}_J_S4", p)
        mean_qs = _expectation(levels, f"sum_q_{component}_J_S4", p)
        s_control.append((mean_qs-mean_q*mean_s)/birth)
    transport = target["transport"]
    phase = complex(_number(transport["real"]), _number(transport["imag"]))
    canonical = phase*complex(*gamma)
    scaled = run["N"]**(13.0/8.0)*canonical
    return {
        "p_matching": p,
        "mean_q_residual": mean_q,
        "birth_mass": birth,
        "gamma_lab": [gamma[0], gamma[1]],
        "gamma_S4_control_lab": s_control,
        "Gamma_canonical": [canonical.real, canonical.imag],
        "Y_scaled": [scaled.real, scaled.imag],
    }


def estimate_vector(runs: Sequence[dict], prediction: dict, omitted: tuple[int, int] | None = None):
    geometry = _geometry_map(prediction)
    run_map = {(run["N"], run["modulus"]): run for run in runs}
    values, details = [], {}
    for n in SIZE_ORDER:
        for modulus in MODULUS_ORDER:
            batch = omitted[1] if omitted is not None and omitted[0] == n else None
            estimate = geometry_estimate(run_map[(n, modulus)], geometry[(n, modulus)], batch)
            values.extend(estimate["Y_scaled"])
            details[f"N{n}:{modulus}"] = estimate
    return np.asarray(values, dtype=float), details


def jackknife_covariance(runs: Sequence[dict], prediction: dict) -> tuple[np.ndarray, dict]:
    full, details = estimate_vector(runs, prediction)
    covariance = np.zeros((18,18), dtype=float)
    phase = prediction["phase1_microcanonical_matching_root"]
    batches = int(phase["batches"])
    for size_index, n in enumerate(SIZE_ORDER):
        coordinates = slice(6*size_index, 6*(size_index+1))
        deleted = np.asarray([
            estimate_vector(runs, prediction, (n, batch))[0][coordinates]
            for batch in range(batches)
        ])
        centered = deleted-deleted.mean(axis=0)
        covariance[coordinates, coordinates] = (batches-1.0)/batches*(centered.T@centered)
    return covariance, {"observation": full, "geometries": details}


def _design_matrix(model: str, prediction: dict) -> np.ndarray:
    modulus_values = {
        row["id"]: (_number(row["E4hat_over_i"]), _number(row["eta_coordinate"]))
        for row in prediction["moduli"]
    }
    complex_features = []
    for n in SIZE_ORDER:
        for modulus in MODULUS_ORDER:
            f, eta = modulus_values[modulus]
            if model == "Q4_epsilon_ordinary":
                features = [f]
            elif model == "Q4_energy_Jordan":
                features = [f, f*math.log(n), f*eta]
            elif model == "generic_allowed_H4_pure":
                features = [float(modulus == m) for m in MODULUS_ORDER]
            elif model == "generic_allowed_H4_affine_log":
                features = []
                for m in MODULUS_ORDER:
                    features.extend([float(modulus == m), float(modulus == m)*math.log(n)])
            elif model == "zero_response":
                features = []
            else:  # pragma: no cover
                raise ValueError(model)
            complex_features.append(features)
    p = len(complex_features[0]) if complex_features else 0
    design = np.zeros((18,2*p), dtype=float)
    for index, features in enumerate(complex_features):
        design[2*index, :p] = features
        design[2*index+1, p:] = features
    return design


def _survival(chi_square: float, dof: int) -> float | None:
    if dof <= 0:
        return None
    return float(mp.gammainc(mp.mpf(dof)/2, mp.mpf(chi_square)/2, mp.inf, regularized=True))


def fit_models(observation: np.ndarray, covariance: np.ndarray, prediction: dict) -> dict:
    rcond = float(prediction["phase1_microcanonical_matching_root"]["gls_rcond"])
    weight = np.linalg.pinv(covariance, rcond=rcond, hermitian=True)
    scores = {}
    for model in prediction["model_order"]:
        design = _design_matrix(model, prediction)
        if design.shape[1]:
            normal = design.T@weight@design
            normal_inv = np.linalg.pinv(normal, rcond=rcond, hermitian=True)
            beta = normal_inv@design.T@weight@observation
            residual = observation-design@beta
            beta_cov = normal_inv
        else:
            beta = np.zeros(0); beta_cov = np.zeros((0,0)); residual = observation
        chi_square = float(residual@weight@residual)
        dof = len(observation)-design.shape[1]
        scores[model] = {
            "chi_square": chi_square,
            "dof": dof,
            "survival_p": _survival(chi_square, dof),
            "coefficients_real_then_imag": beta.tolist(),
            "coefficient_covariance": beta_cov.tolist(),
        }
    ordinary = scores["Q4_epsilon_ordinary"]
    jordan = scores["Q4_energy_Jordan"]
    delta = max(0.0, ordinary["chi_square"]-jordan["chi_square"])
    scores["ordinary_to_Jordan_gain"] = {
        "delta_chi_square": delta,
        "delta_dof": 4,
        "survival_p": _survival(delta, 4),
    }
    alpha = float(prediction["phase1_microcanonical_matching_root"]["decision_alpha"])
    supported = [
        name for name in prediction["model_order"]
        if scores[name]["survival_p"] is not None and scores[name]["survival_p"] >= alpha
    ]
    selected = supported[0] if supported else "none"
    return {"scores": scores, "supported_in_frozen_order": supported, "selected": selected}


def build_report(runs: Sequence[dict], prediction: dict) -> dict:
    provenance = validate_runs(runs, prediction)
    covariance, estimates = jackknife_covariance(runs, prediction)
    fits = fit_models(estimates["observation"], covariance, prediction)
    return {
        "schema": "matching-one/p275-atop-field-identity-score/v1",
        "status": "frozen_phase1_reveal",
        "issues": [205, 275],
        "provenance": provenance,
        "coordinate_order": list(COORDINATE_ORDER),
        "estimates": {
            "observation_Y": estimates["observation"].tolist(),
            "geometries": estimates["geometries"],
            "covariance_18x18": covariance.tolist(),
        },
        "model_score": fits,
        "scientific_card": [
            "Question: which field completes the already-selected global H4 channel?",
            "Observable: gamma=Cov(A_top,J_D4)/B at each delete-one finite matching root.",
            "Selector: three moduli times three cyclic sizes with full same-field covariance.",
            "Q4 fingerprint: E4hat modulus vector; Jordan adds fixed log-N and eta directions.",
            "Stop rule: failure promotes another H4 completion; H8/H12 are not rescored.",
        ],
    }


def render_markdown(report: dict) -> str:
    lines = ["# P275 A_top field-identity Phase 1", "", f"Selected: `{report['model_score']['selected']}`", "",
             "| geometry | p_N | Re Y | Im Y | B |", "|---|---:|---:|---:|---:|"]
    values = report["estimates"]["observation_Y"]
    index = 0
    for n in SIZE_ORDER:
        for modulus in MODULUS_ORDER:
            row = report["estimates"]["geometries"][f"N{n}:{modulus}"]
            lines.append(f"| N{n}/{modulus} | {row['p_matching']:.12g} | {values[index]:.8g} | {values[index+1]:.8g} | {row['birth_mass']:.8g} |")
            index += 2
    lines.extend(["", "## Frozen model scores", "", "| model | chi2 | dof | survival p |", "|---|---:|---:|---:|"])
    for model in ("Q4_epsilon_ordinary", "Q4_energy_Jordan", "generic_allowed_H4_pure", "generic_allowed_H4_affine_log", "zero_response"):
        row = report["model_score"]["scores"][model]
        lines.append(f"| {model} | {row['chi_square']:.6g} | {row['dof']} | {row['survival_p']:.6g} |")
    lines.extend(["", "## Scientific card", ""] + [f"- {line}" for line in report["scientific_card"]] + [""])
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
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(render_markdown(report), encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

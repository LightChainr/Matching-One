#!/usr/bin/env python3
"""Analyze conditional winding-line responses without replaying or finding roots."""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
import platform
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.stats import chi2

import analyze_norm4_source_thermal as old
from analyze_norm4_source_endpoint_1m import load_profile
from norm4_source_line_core import INPUT_COLUMNS, COMPLEX_COMPONENTS, evaluate

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis/norm4_source_line_contract.json"
OUTPUT = ROOT / "results/norm4-source-line"
NS = (65, 85, 130, 170, 260, 340)
SCALAR_FIELDS = ("p0", "rootdot_fugacity", "root_comoving_rank1_fugacity", "rank1_common_E")
COMPONENTS = (("first", "re"), ("first", "im"), ("second", "re"), ("second", "im"))
EIGEN_RCOND = 1e-10


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def display(path):
    path = path.resolve()
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


def input_record(path, role, **extra):
    return {"path": display(path), "role": role, "sha256": digest(path), **extra}


def load_line_profile(path, n, per_batch, run, source_profiles):
    """One input pass, including the two already available rank1 identities."""
    profile = np.zeros((100, 2, n + 1, 6), dtype=float)
    seen = set()
    for row in csv.DictReader(path.open()):
        batch, k = int(row["batch"]), int(row["k"])
        direction = row["orientation"]
        g = ("first", "second").index(direction)
        key = (batch, g, k)
        if (key in seen or int(row["n"]) != n or not 0 <= batch < 100 or not 0 <= k <= n
                or int(row["samples"]) != per_batch
                or [int(row["a"]), int(row["b"])] != run[direction]):
            raise ValueError(f"{path}: incompatible or duplicate conditional line row {key}")
        seen.add(key)
        profile[batch, g, k] = [float(row[field]) for field in INPUT_COLUMNS]
    if len(seen) != 100 * 2 * (n + 1) or not np.isfinite(profile).all():
        raise ValueError(f"{path}: incomplete or nonfinite all-K conditional line profile")
    expected_rank1 = per_batch - source_profiles[..., 1]
    expected_rank1_s = n * (source_profiles[..., 2] - source_profiles[..., 4])
    if not np.array_equal(profile[..., 0], expected_rank1):
        raise ValueError(f"N{n}: new rank1 support differs from the same old source configurations")
    # Only the old read_raw density conversion introduces roundoff here.
    if not np.allclose(profile[..., 1], expected_rank1_s, rtol=1e-13, atol=1e-6):
        raise ValueError(f"N{n}: rank1*s does not match N*(S-ES) in its old paired source batches")
    return profile, {
        "rank1_count_equals_samples_minus_E": True,
        "rank1_s_equals_bulk_S_minus_ES": True,
        "maximum_bulk_source_roundoff": float(np.max(np.abs(profile[..., 1] - expected_rank1_s))),
        "operation": "same input pass; source/line identities on each batch, direction and K",
    }


def point(line_sums, source_sums, samples, n, scalar):
    p0 = float(scalar["p0"])
    direction = {
        name: old.direction_values(old.binomial_moments(source_sums[g], samples, p0, n))
        for g, name in enumerate(("first", "second"))
    }
    return evaluate(line_sums, samples, n, p0, scalar, direction)


def vectorize(points):
    return {f"N{n}.{key}": float(value) for n in NS for key, value in points[n].items()}


def supported_zero(values, covariance, labels=None, maximum_rank=None):
    """Nominal zero comparison on the supported covariance subspace."""
    values = np.asarray(values, dtype=float)
    covariance = np.asarray(covariance, dtype=float)
    covariance = (covariance + covariance.T) / 2
    result = {"values": values.tolist(), "covariance": covariance.tolist(),
              "labels": labels, "eigenvalue_relative_tolerance": EIGEN_RCOND}
    if not np.isfinite(values).all() or not np.isfinite(covariance).all():
        return {**result, "status": "nonfinite_covariance", "nominal_p": None}
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    maximum = float(max(0, eigenvalues[-1]))
    tolerance = maximum * EIGEN_RCOND
    support = eigenvalues > tolerance
    rank = int(np.count_nonzero(support))
    result.update(eigenvalues=eigenvalues.tolist(), supported_rank=rank)
    if eigenvalues[0] < -max(maximum, np.finfo(float).tiny) * 1e-8:
        return {**result, "status": "nonpositive_covariance", "nominal_p": None}
    if maximum_rank is not None and rank > maximum_rank:
        return {**result, "status": "rank_exceeds_null_tangent", "nominal_p": None}
    if rank == 0:
        return {**result, "status": "zero_covariance_support", "nominal_p": None}
    coordinates = eigenvectors[:, support].T @ values
    remainder = values - eigenvectors[:, support] @ coordinates
    outside = float(np.linalg.norm(remainder))
    result["outside_support_norm"] = outside
    reference = max(float(np.linalg.norm(values)), math.sqrt(maximum), np.finfo(float).tiny)
    if outside > 1e-7 * reference:
        return {**result, "status": "residual_outside_covariance_support", "nominal_p": None}
    statistic = float(np.sum(coordinates**2 / eigenvalues[support]))
    return {**result, "status": "computed", "chi_square": statistic, "df": rank,
            "nominal_p": float(chi2.sf(statistic, rank))}


def selected_zero(labels, central, covariance, selected):
    index = [labels.index(label) for label in selected]
    return supported_zero(central[index], covariance[np.ix_(index, index)], selected)


def wedge_null(n, labels, central, covariance, minimum_thermal_snr):
    """Regular-null wedge score using T=mu_p and beta_clock=b*p*(1-p)."""
    c_labels = [f"N{n}.{g}.conditional_cov_s_{part}" for g, part in COMPONENTS]
    t_labels = [f"N{n}.{g}.mu_p_{part}" for g, part in COMPONENTS]
    w_labels = [f"N{n}.clock_wedge.{COMPLEX_COMPONENTS[i]}__{COMPLEX_COMPONENTS[j]}"
                for i, j in itertools.combinations(range(4), 2)]
    indices = [labels.index(label) for label in c_labels + t_labels]
    w_indices = [labels.index(label) for label in w_labels]
    c, t = central[indices[:4]], central[indices[4:]]
    sigma = covariance[np.ix_(indices, indices)]
    cc, ct, tc, tt = sigma[:4, :4], sigma[:4, 4:], sigma[4:, :4], sigma[4:, 4:]
    t_norm2 = float(t @ t)
    t_variance_trace = float(max(0, np.trace(tt)))
    thermal_snr = math.sqrt(t_norm2 / t_variance_trace) if t_variance_trace > 0 else None
    weak = t_norm2 == 0 or (thermal_snr is not None and thermal_snr < minimum_thermal_snr)
    w = central[w_indices]
    result = {
        "vector_order": list(COMPLEX_COMPONENTS), "C": c.tolist(), "T_mu_p": t.tolist(),
        "C_T_labels": c_labels + t_labels, "C_T_covariance": sigma.tolist(),
        "wedge_labels": w_labels, "wedge_values": w.tolist(),
        "original_nonlinear_wedge_covariance": covariance[np.ix_(w_indices, w_indices)].tolist(),
        "within_geometry_wedges": [w_labels[0], w_labels[-1]],
        "cross_geometry_wedges": w_labels[1:-1],
        "thermal_norm": math.sqrt(t_norm2), "thermal_covariance_trace": t_variance_trace,
        "thermal_norm_to_rms_SE": thermal_snr,
        "zero_thermal_variance_with_nonzero_vector": bool(t_variance_trace == 0 and t_norm2 > 0),
        "minimum_thermal_norm_snr": minimum_thermal_snr,
        "weak_thermal_vector": bool(weak),
        "beta_definition": "beta_clock=b*p0*(1-p0), because T=mu_p rather than p0*(1-p0)*mu_p",
        "null": "C=beta_clock*T, common real beta for the four geometry/quadrature components",
        "score_rule": "Omega=L_T Cov(C-beta_clock*T) L_T^T; supported null rank <=3, not six nonlinear wedge eigenvalues",
    }
    if t_norm2 == 0:
        return {**result, "status": "weak_thermal_vector_descriptive_only", "beta_clock": None,
                "null_tangent_covariance": None, "nominal_p": None,
                "interpretation": "The wedge is uninformative at zero T; inspect the saved C vector, not a closure claim"}
    beta = float((t @ c) / t_norm2)
    l_t = np.zeros((6, 4))
    for row, (i, j) in enumerate(itertools.combinations(range(4), 2)):
        l_t[row, i] = t[j]
        l_t[row, j] = -t[i]
    residual_covariance = cc - beta * ct - beta * tc + beta**2 * tt
    omega = l_t @ residual_covariance @ l_t.T
    omega = (omega + omega.T) / 2
    result.update(beta_clock=beta, null_tangent_covariance=omega.tolist(),
                  C_minus_beta_T_covariance=residual_covariance.tolist())
    if weak:
        return {**result, "status": "weak_thermal_vector_descriptive_only", "nominal_p": None,
                "interpretation": "C,T and all wedges are reported; noisy thermal direction prevents a regular-null closure interpretation"}
    score = supported_zero(w, omega, w_labels, maximum_rank=3)
    result.update(null_support_score=score, status=score["status"],
                  supported_rank=score.get("supported_rank"),
                  chi_square=score.get("chi_square"), df=score.get("df"),
                  nominal_p=score.get("nominal_p"))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT)
    args = parser.parse_args()
    output_dir = args.output_dir.resolve()
    destination = output_dir / "latest.json"
    if destination.exists():
        raise ValueError("saved line result exists; refusing overwrite or repeated analysis")
    started = time.perf_counter()
    contract = json.loads(CONTRACT.read_text())
    if tuple(contract["Ns"]) != NS or contract["batches"] != 100:
        raise ValueError("line contract must preserve the six sizes and 100 aligned batches")
    minimum_thermal_snr = float(contract.get("wedge_min_thermal_norm_snr", 3.0))
    if not math.isfinite(minimum_thermal_snr) or minimum_thermal_snr <= 0:
        raise ValueError("positive finite descriptive weak-thermal threshold required")
    source_path = ROOT / contract["source_result"]
    source = json.loads(source_path.read_text())
    source_index = {label: i for i, label in enumerate(source["labels"])}
    if len(source_index) != len(source["labels"]):
        raise ValueError("saved source result labels must be unique")
    manifest_path = ROOT / contract.get("source_manifest", str(old.MANIFEST.relative_to(ROOT)))
    manifest = json.loads(manifest_path.read_text())
    runs = {run["N"]: run for run in manifest["runs"]}
    replay_path = output_dir / "run.json"
    replay = json.loads(replay_path.read_text())
    if replay.get("status") != "completed":
        raise ValueError("conditional line replay receipt is not completed")

    source_profiles, line_profiles, inputs, agreements = {}, {}, [], {}
    per_batch = {n: int(contract["marked_permutations_by_N"][str(n)]) // 100 for n in NS}
    for n in NS:
        expected = 10000 if n in (260, 340) else 1000
        if per_batch[n] != expected:
            raise ValueError(f"N{n}: unexpected marked sample union")
        original_path = old.OUTPUT / "raw" / f"n{n}.csv"
        source_profiles[n] = load_profile(original_path, n, 1000, runs[n])
        inputs.append(input_record(original_path, "original_100k_source", N=n))
        if n in (260, 340):
            increment_path = ROOT / "results/norm4-source-endpoint-1m/increment/raw" / f"n{n}.csv"
            source_profiles[n] += load_profile(increment_path, n, 9000, runs[n])
            inputs.append(input_record(increment_path, "same_batch_900k_source_increment", N=n))
        line_path = output_dir / "raw" / f"n{n}.csv"
        line_profiles[n], agreements[n] = load_line_profile(line_path, n, per_batch[n], runs[n], source_profiles[n])
        inputs.append(input_record(line_path, "conditional_line_source_products", N=n))
    # Freeze the exact source artifacts whose central/LOO roots are reused.
    saved_source_inputs = {item["path"]: item["sha256"] for item in source["inputs"]}
    for item in inputs:
        if item["role"] != "conditional_line_source_products" and saved_source_inputs.get(item["path"]) != item["sha256"]:
            raise ValueError(f"source profile does not match the saved scalar analysis: {item['path']}")

    totals_source = {n: source_profiles[n].sum(axis=0) for n in NS}
    totals_line = {n: line_profiles[n].sum(axis=0) for n in NS}
    central_points, central_diagnostics = {}, {}
    for n in NS:
        scalar = {field: source["by_N"][str(n)]["points"][field] for field in SCALAR_FIELDS}
        central_points[n], central_diagnostics[n] = point(
            totals_line[n], totals_source[n], per_batch[n] * 100, n, scalar)
    central_map = vectorize(central_points)
    labels = list(central_map)
    central = np.asarray(list(central_map.values()))
    covariance = np.zeros((len(labels), len(labels)))
    contributions = {}
    covered = []
    for group in contract["dependency_groups"]:
        group_id, sizes = group["id"], list(group["Ns"])
        covered.extend(sizes)
        saved = source["covariance_contributions"]["source:" + group_id]
        if saved["Ns"] != sizes or saved["delete_one_batch_ids"] != list(range(100)):
            raise ValueError(f"{group_id}: saved source delete-one indices or group membership differ")
        for n in sizes:
            if saved["batch_counts"] != [per_batch[n]] * 100:
                raise ValueError(f"{group_id}: saved source batches differ from conditional line batches")
        saved_vectors = np.asarray(saved["delete_one_vectors"], dtype=float)
        if saved_vectors.shape != (100, len(source["labels"])):
            raise ValueError(f"{group_id}: incomplete saved scalar delete-one vectors")
        vectors, roots, alphas = [], {n: [] for n in sizes}, {n: [] for n in sizes}
        for batch in range(100):
            changed = dict(central_points)
            for n in sizes:
                scalar = {field: float(saved_vectors[batch, source_index[f"N{n}.{field}"]])
                          for field in SCALAR_FIELDS}
                changed[n], _ = point(totals_line[n] - line_profiles[n][batch],
                                      totals_source[n] - source_profiles[n][batch],
                                      per_batch[n] * 99, n, scalar)
                roots[n].append(changed[n]["p0"])
                alphas[n].append(changed[n]["alpha_E"])
            vectors.append(list(vectorize(changed).values()))
        vectors = np.asarray(vectors)
        centered = vectors - vectors.mean(axis=0)
        component = 99 / 100 * centered.T @ centered
        covariance += component
        contributions[group_id] = {
            "Ns": sizes, "stage": "same_source_and_line_marks",
            "batch_counts": saved["batch_counts"], "delete_one_batch_ids": list(range(100)),
            "delete_one_vectors": vectors.tolist(), "covariance": component.tolist(),
            "saved_source_root_vector_key": "source:" + group_id,
            "root_ranges": {n: [min(roots[n]), max(roots[n])] for n in sizes},
            "alpha_E_ranges": {n: [min(alphas[n]), max(alphas[n])] for n in sizes},
            "operation": "omit the identical source/line batch; reuse its saved self-rooted p0 and scalar responses; other groups fixed",
        }
    if sorted(covered) != sorted(NS):
        raise ValueError("every N must belong to exactly one preserved source dependency group")
    errors = np.sqrt(np.maximum(0, np.diag(covariance)))
    estimates = {label: {"value": float(value), "se": float(error),
                          "z": float(value / error) if error > 0 else None}
                 for label, value, error in zip(labels, central, errors)}
    primary_labels = {
        n: [f"N{n}.{g}.resid_E_{part}" for g, part in COMPONENTS] for n in NS
    }
    primary_by_n = {n: selected_zero(labels, central, covariance, primary_labels[n]) for n in NS}
    all_primary = [label for n in NS for label in primary_labels[n]]
    primary_joint = selected_zero(labels, central, covariance, all_primary)
    wedge_by_n = {n: wedge_null(n, labels, central, covariance, minimum_thermal_snr) for n in NS}
    summary = {
        "max_abs_primary_component_z": max(abs(estimates[label]["z"] or 0) for label in all_primary),
        "primary_joint_nominal_p": primary_joint.get("nominal_p"),
        "primary_joint_supported_rank": primary_joint.get("supported_rank"),
        "primary_per_N_nominal_p": {n: primary_by_n[n].get("nominal_p") for n in NS},
        "wedge_per_N": {n: {"status": wedge_by_n[n]["status"],
                              "nominal_p": wedge_by_n[n].get("nominal_p"),
                              "supported_rank": wedge_by_n[n].get("supported_rank"),
                              "thermal_norm_to_rms_SE": wedge_by_n[n]["thermal_norm_to_rms_SE"]}
                         for n in NS},
    }
    code_paths = [Path(__file__).resolve(), ROOT / "scripts/norm4_source_line_core.py",
                  ROOT / "scripts/analyze_norm4_source_thermal.py",
                  ROOT / "scripts/analyze_p40_source_thermal.py",
                  ROOT / "scripts/analyze_norm4_source_endpoint_1m.py"]
    result = {
        "schema": "matching-one.norm4-source-line.v1",
        "status": "computed_conditional_line_source_response",
        "execution_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "command": [sys.executable, *sys.argv], "contract": contract,
        "labels": labels, "estimates": estimates, "covariance": covariance.tolist(),
        "by_N": central_diagnostics, "covariance_contributions": contributions,
        "primary_E_plus_clock": {"by_N": primary_by_n, "joint_all_24_components": primary_joint,
                                  "loading": "alpha_E from the additional pooled-rank1 equation; no extra degree of freedom subtracted from the four-real conditional residual"},
        "common_topology_plus_clock": {"by_N": wedge_by_n,
                                        "dependence": "all six wedges and all N retain joint covariance; per-N p values are not independent evidence votes"},
        "summary": summary,
        "input_identity": agreements, "inputs": inputs,
        "source_result": input_record(source_path, "saved_current_self_rooted_scalars_and_delete_ones"),
        "source_manifest": input_record(manifest_path, "original_source_geometry_and_counter_contract"),
        "line_replay_receipt": {**input_record(replay_path, "new_line_products_on_old_counters"), "content": replay},
        "contract_source": input_record(CONTRACT, "frozen_conditional_observer_contract"),
        "code": [input_record(path, "analysis_code") for path in code_paths],
        "uncertainty": "Three source-group aligned delete-one contributions only; same source/line omissions and saved root transport. No independent complement anchor is used or propagated. Components, paired contrasts, E residuals and wedge views share the same observations.",
        "wedge_null_support": "T=mu_p, beta_clock=(T dot C)/(T dot T); null-tangent rank at most3. Original six-wedge nonlinear covariance is retained but does not set the null degrees of freedom.",
        "weak_thermal_rule": {"minimum_norm_to_rms_SE": minimum_thermal_snr,
                              "effect": "descriptive wedges only when thermal direction is weak; primary E residual remains reported; no closure claim and no additional run"},
        "interpretation": "Finite conditional primitive-line readout, not an energy-field identity, continuum scaling result or proof of pointwise microscopic source equality. Every geometry is conditioned before pairing and both complex quadratures are retained.",
        "root_finder_calls": 0, "configuration_replays_by_this_script": 0, "new_samples": 0,
        "server_actions": 0, "test_suites": [],
        "environment": {"python": platform.python_version(), "numpy": np.__version__,
                        "scipy": scipy.__version__, "machine": platform.machine()},
        "elapsed_seconds": time.perf_counter() - started,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    with destination.open("x") as handle:
        handle.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"elapsed_seconds": result["elapsed_seconds"], **summary}, indent=2))


if __name__ == "__main__":
    main()

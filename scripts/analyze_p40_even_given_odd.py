#!/usr/bin/env python3
"""Join one old-counter E=q^2 supplement to the saved P40 production Gram."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import time
from pathlib import Path

import numpy as np

from analyze_p40_absolute_cluster import cos4, digest, git_bytes, load_input


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "analysis/p40_even_given_odd_replay.json"
COORDS = ["q", "S", "K/N", "KK/[N(N-1)]", "T/(2N)", "chi/N", "E"]
MODES = ("raw", "clock", "full")
BASE_METRICS = ["q_mean", "E_mean", "S_mean", "var_q", "var_E", "cov_E_q", "var_E_given_q"]
METRICS = BASE_METRICS + [f"{field}_{mode}" for mode in MODES
                          for field in ("Jq", "JE", "C", "curvature_c")]
JOINT_FIELDS = [f"{field}_{mode}" for mode in MODES for field in ("C", "JE")]
H4_FIELDS = ["q_mean", "E_mean", "Jq_raw", "Jq_clock", "Jq_full", *JOINT_FIELDS]
LABELS = [f"H4.{field}" for field in H4_FIELDS] + [
    f"{direction}.{field}" for direction in ("first", "second") for field in METRICS
]
REQUIRED_COLUMNS = {
    "n", "a", "b", "orientation", "batch", "samples", "sum_q", "sum_e",
    "sum_s", "sum_k", "sum_kk", "sum_edges", "sum_chi", "sum_e_s",
    "sum_e_k", "sum_e_kk", "sum_e_edges", "sum_e_chi",
}
MATCH_ATOL = 2e-7
MATCH_RTOL = 2e-11


def path_label(path):
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def join_supplement(path, run, counts, old_sums, old_grams):
    """Retain marginal seven-coordinate Grams, never fabricate cross-E cells."""
    payload = path.read_bytes()
    reader = csv.DictReader(payload.decode("utf-8").splitlines())
    if not REQUIRED_COLUMNS.issubset(reader.fieldnames or []):
        raise ValueError(f"{path}: missing supplement columns")
    rows = {}
    for raw in reader:
        row = {key: (value if key == "orientation" else int(value))
               for key, value in raw.items() if key in REQUIRED_COLUMNS}
        key = (row["batch"], row["orientation"])
        if key in rows:
            raise ValueError(f"{path}: duplicate batch/orientation {key}")
        rows[key] = row
    expected = {(batch, direction) for batch in range(100)
                for direction in ("first", "second")}
    if set(rows) != expected:
        raise ValueError(f"{path}: expected 100 aligned two-orientation batches")

    n = run["N"]
    sums = np.zeros((100, 2, 7))
    grams = np.zeros((100, 2, 7, 7))
    maximum_first_error = maximum_e_error = 0.0
    for batch in range(100):
        for index, direction in enumerate(("first", "second")):
            row = rows[(batch, direction)]
            if (row["n"] != n or [row["a"], row["b"]] != run[direction]
                    or row["samples"] != 10000 or row["samples"] != counts[batch]):
                raise ValueError(f"{path}: changed source geometry/count at {batch}/{direction}")
            recovered = np.array([
                row["sum_q"], row["sum_s"] / n, row["sum_k"] / n,
                row["sum_kk"] / (n * (n - 1)), row["sum_edges"] / (2 * n),
                row["sum_chi"] / n,
            ], dtype=float)
            block = slice(6 * index, 6 * (index + 1))
            parent_sum = old_sums[batch, block]
            parent_gram = old_grams[batch, block, block]
            first_error = float(np.max(np.abs(recovered - parent_sum)))
            e_error = abs(float(row["sum_e"]) - float(parent_gram[0, 0]))
            maximum_first_error = max(maximum_first_error, first_error)
            maximum_e_error = max(maximum_e_error, e_error)
            if not np.allclose(recovered, parent_sum, atol=MATCH_ATOL, rtol=MATCH_RTOL):
                raise ValueError(f"{path}: old first moments differ at {batch}/{direction}")
            if not math.isclose(row["sum_e"], parent_gram[0, 0],
                                abs_tol=MATCH_ATOL, rel_tol=MATCH_RTOL):
                raise ValueError(f"{path}: E sum differs from old q squared at {batch}/{direction}")

            sums[batch, index, :6] = parent_sum
            sums[batch, index, 6] = row["sum_e"]
            grams[batch, index, :6, :6] = parent_gram
            e_cross = np.array([
                row["sum_q"], row["sum_e_s"] / n, row["sum_e_k"] / n,
                row["sum_e_kk"] / (n * (n - 1)), row["sum_e_edges"] / (2 * n),
                row["sum_e_chi"] / n,
            ], dtype=float)
            grams[batch, index, 6, :6] = e_cross
            grams[batch, index, :6, 6] = e_cross
            grams[batch, index, 6, 6] = row["sum_e"]  # E^2=E; E*q=q above.

    diagnostic = {
        "aligned_rows": len(rows), "maximum_normalized_first_moment_error": maximum_first_error,
        "maximum_E_vs_parent_q_squared_error": maximum_e_error,
        "absolute_tolerance": MATCH_ATOL, "relative_tolerance": MATCH_RTOL,
        "cross_orientation_E_Gram": "not supplied or imputed; marginal estimands use paired batch jackknife covariance",
    }
    source = {"path": path_label(path), "sha256": digest(payload), "bytes": len(payload)}
    return sums, grams, source, diagnostic


def single(count, sums, gram):
    mean = sums / count
    covariance = (gram - np.outer(sums, sums) / count) / (count - 1)
    covariance = (covariance + covariance.T) / 2
    var_q, var_e, cov_eq = covariance[0, 0], covariance[6, 6], covariance[6, 0]
    if var_q <= 0:
        raise ValueError("q variance is nonpositive")
    var_even_given_odd = var_e - cov_eq * cov_eq / var_q
    if var_even_given_odd <= 0:
        raise ValueError("E given q variance is nonpositive")
    probabilities = np.array([(mean[6] - mean[0]) / 2, 1 - mean[6],
                              (mean[6] + mean[0]) / 2])
    if np.any(probabilities <= 0):
        raise ValueError("all three rank sectors are required for sector curvature")

    vectors = {"raw": np.array([0., 1., 0., 0., 0., 0., 0.])}
    control_diagnostics = {}
    for mode, indices in (("clock", [2, 3]), ("full", [2, 3, 4, 5])):
        control_cov = covariance[np.ix_(indices, indices)]
        scales = np.sqrt(np.diag(control_cov))
        if np.any(scales <= 0):
            raise ValueError(f"{mode}: degenerate control variance")
        correlation = control_cov / np.outer(scales, scales)
        eigenvalues = np.linalg.eigvalsh(correlation)
        if eigenvalues[0] <= 1e-10:
            raise ValueError(f"{mode}: declared control span numerically unresolved")
        beta = np.linalg.solve(correlation, covariance[indices, 1] / scales) / scales
        vector = vectors["raw"].copy()
        vector[indices] -= beta
        vectors[mode] = vector
        control_diagnostics[mode] = {
            "controls": [COORDS[i] for i in indices], "coefficients": beta.tolist(),
            "condition_correlation": float(eigenvalues[-1] / eigenvalues[0]),
            "source_covariance_with_controls": (covariance[indices] @ vector).tolist(),
        }

    values = {
        "q_mean": mean[0], "E_mean": mean[6], "S_mean": mean[1],
        "var_q": var_q, "var_E": var_e, "cov_E_q": cov_eq,
        "var_E_given_q": var_even_given_odd,
    }
    mode_diagnostics = {}
    for mode in MODES:
        vector = vectors[mode]
        jq = float(covariance[0] @ vector)
        je = float(covariance[6] @ vector)
        conditional = je - cov_eq / var_q * jq
        source_mean = float(mean @ vector)
        # Convert unbiased sample covariances back to empirical raw moments.
        qs = jq * (count - 1) / count + mean[0] * source_mean
        es = je * (count - 1) / count + mean[6] * source_mean
        sector_means = np.array([(es - qs) / 2, source_mean - es, (es + qs) / 2]) / probabilities
        curvature = float((sector_means[0] + sector_means[2]) / 2 - sector_means[1])
        values.update({f"Jq_{mode}": jq, f"JE_{mode}": je,
                       f"C_{mode}": conditional, f"curvature_c_{mode}": curvature})
        mode_diagnostics[mode] = {
            "source_vector": vector.tolist(), "sector_source_means_minus_zero_plus": sector_means.tolist(),
            "curvature_c": curvature, "C_minus_c_times_var_E_given_q": float(conditional - curvature * var_even_given_odd),
            "q_compensation_h_prime": float(-jq / var_q),
            "source_variance": float(vector @ covariance @ vector),
        }
    diagnostic = {"controls": control_diagnostics, "modes": mode_diagnostics,
                  "sector_probabilities_minus_zero_plus": probabilities.tolist(),
                  "var_E_given_q": float(var_even_given_odd)}
    return values, diagnostic


def point(count, sums, grams, delta):
    values, diagnostic = {}, {}
    for index, direction in enumerate(("first", "second")):
        values[direction], diagnostic[direction] = single(count, sums[index], grams[index])
    output = [(values["first"][key] - values["second"][key]) / delta for key in H4_FIELDS]
    output += [values[direction][key] for direction in ("first", "second") for key in METRICS]
    return np.asarray(output, dtype=float), diagnostic


def analyze(output_dir):
    start = time.perf_counter()
    replay_text = MANIFEST.read_text(encoding="utf-8")
    replay_contract = json.loads(replay_text)
    source_manifest = Path(replay_contract.get("source_analysis_manifest", "analysis/p40_absolute_cluster_reanalysis.json"))
    if not source_manifest.is_absolute():
        source_manifest = ROOT / source_manifest
    contract = json.loads(source_manifest.read_text(encoding="utf-8"))
    if replay_contract["source_commit"] != contract["input_commit"]:
        raise ValueError("supplement and old Gram name different source commits")
    if sorted(run["N"] for run in contract["runs"]) != [65, 85]:
        raise ValueError("only the two archived P40 million-sample blocks are supported")

    backend_path = "src/gaussian_orientation_mc.cpp"
    backend = git_bytes(contract["input_commit"], backend_path)
    backend_blob = hashlib.sha1(f"blob {len(backend)}\0".encode() + backend).hexdigest()
    if replay_contract.get("source_backend_blob") not in (None, backend_blob):
        raise ValueError("source backend blob differs from replay declaration")
    inputs = [{"commit": contract["input_commit"], "path": backend_path,
               "sha256": digest(backend), "git_blob_sha1": backend_blob}]
    results = {}
    for run in contract["runs"]:
        counts, old_sums, old_grams, source, metadata = load_input(contract, run)
        sums, grams, supplement, agreement = join_supplement(
            output_dir / "raw" / f"n{run['N']}.csv", run, counts, old_sums, old_grams)
        total_count = int(counts.sum())
        total_sum, total_gram = sums.sum(axis=0), grams.sum(axis=0)
        delta_exact = cos4(*run["first"]) - cos4(*run["second"])
        delta = float(delta_exact)
        central, diagnostic = point(total_count, total_sum, total_gram, delta)
        loo = np.array([point(total_count - int(count), total_sum - s, total_gram - g, delta)[0]
                        for count, s, g in zip(counts, sums, grams)])
        deviations = loo - loo.mean(axis=0)
        covariance = 99 / 100 * deviations.T @ deviations
        errors = np.sqrt(np.maximum(0, np.diag(covariance)))
        estimates = {key: {"value": float(value), "se": float(error),
                            "z": float(value / error) if error > 0 else None}
                     for key, value, error in zip(LABELS, central, errors)}
        results[str(run["N"])] = {
            "estimates": estimates, "covariance": covariance.tolist(),
            "delete_one_vectors": loo.tolist(), "diagnostics": diagnostic,
            "supplement_parent_agreement": agreement, "delta_cos4_exact": str(delta_exact),
            "parent_metadata": metadata, "samples": total_count,
            "orientation_order": ["first", "second"],
            "transformed_sample_sum": total_sum.tolist(),
            "transformed_marginal_sample_gram": total_gram.tolist(),
        }
        inputs.extend(source)
        inputs.append(supplement)

    joint = {}
    for field in JOINT_FIELDS:
        scores = [result["estimates"][f"H4.{field}"]["z"] for result in results.values()]
        if any(score is None for score in scores):
            raise ValueError(f"zero estimated uncertainty for predeclared {field}")
        statistic = sum(score * score for score in scores)
        joint[field] = {"chi_square": statistic, "df": 2, "nominal_p": math.exp(-statistic / 2)}

    code_paths = [Path(__file__), ROOT / "scripts/analyze_p40_absolute_cluster.py", MANIFEST, source_manifest]
    runner = ROOT / "src/p40_even_given_odd_replay.cpp"
    if runner.exists():
        code_paths.append(runner)
    return {
        "schema": "matching-one.p40-even-given-odd.v1", "labels": LABELS,
        "coordinate_order_per_orientation": COORDS, "by_N": results,
        "joint_zero": joint, "predeclared_joint_fields": JOINT_FIELDS,
        "inputs": inputs, "replay_contract": replay_contract,
        "replay_contract_text": replay_text, "source_contract": contract,
        "definitions": {
            "primary_source": "S=(c_blackNN+c_whiteMatching)/N, identical raw definition in both geometries",
            "auxiliary_sources": "clock/full empirical covariance projections fitted separately by geometry and refitted in each delete-one; neither span includes q or E",
            "E": "q squared, with q in {-1,0,1}; the motif E_edges is not this observable",
            "C": "Cov(E,S_mode)-Cov(E,q)/Var(q)*Cov(q,S_mode)",
            "physical_tangent": "P proportional to Bernoulli(p) exp(lambda*S_mode+h(lambda)*q), h_prime(0)=-Cov(q,S_mode)/Var(q), holds mean q fixed to first order; rank-one probability derivative is -C",
            "sector_curvature_c": "(E[S_mode|q=-1]+E[S_mode|q=1])/2-E[S_mode|q=0]",
            "identity": "C=sector_curvature_c * (Var(E)-Cov(E,q)^2/Var(q))",
            "direction": "(first-second)/exact Delta cos(4theta); no free exponent",
        },
        "uncertainty": "100 aligned delete-one batches per N; omit same batch in both directions and refit all source projections; full output covariance retained; six correlated nominal two-N zero summaries are not pooled as six independent findings",
        "rng_boundary": "counter_uniform includes splitmix64(N); different N blocks use the ordinary N-domain-separated PRNG independence assumption, while same-N orientations and derived views remain paired",
        "coupling_boundary": "reported C and its H4 contrast are marginal responses invariant to cross-geometry cyclic coupling; that coupling still determines paired sampling covariance; missing cross-orientation E Gram entries are not fabricated",
        "interpretation_boundary": "clock/full compensation and subsequent q holding are sequential; adding the q counterterm need not keep the original control means fixed. Nonzero C distinguishes a pure q tangent only at first order in this finite topological marginal, not an independent continuum field or universal operator",
        "not_scoreable": {"N130_norm4_child": "this source has only N65/N85 parents", "thermal_response_derivative": "the required mixed third source/readout/K moment is not supplied by these seven-coordinate Grams"},
        "new_samples": 0, "oldconfigurationreplays": 2000000,
        "scientific_test_suites": [],
        "code": [{"path": path_label(path), "sha256": digest(path.read_bytes())} for path in code_paths],
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "machine": platform.machine()},
        "elapsed_seconds": time.perf_counter() - start,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results/p40-even-given-odd")
    args = parser.parse_args()
    destination = args.output_dir / "latest.json"
    if destination.exists():
        raise ValueError("refusing to overwrite the saved result")
    result = analyze(args.output_dir)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with destination.open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"joint_zero": result["joint_zero"], "elapsed_seconds": result["elapsed_seconds"]}, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""First source derivative of root-normalized U from archived fixed-p samples."""
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
import scipy
from scipy.optimize import brentq

from analyze_p40_absolute_cluster import cos4, digest, git_bytes, load_input


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "analysis/p40_source_thermal_replay.json"
PREVIOUS_EVEN_COMMIT = "56a6267d6a6826a165f93ed3a64a670ca7088180"
ATOMS = ["q", "E", "S", "qS", "ES"]
DIRECTIONS = ("first", "second")
COMMON_LABELS = ["Udot", "Udot_direct", "Udot_rootmotion", "Udot_slope_source",
                 "Udot_slope_root", "U", "p0", "rootdot", "D", "B", "pooled_q",
                 "pooled_Jq", "pooled_Jq_p", "pooled_q_pp", "P4_JE_p", "P4_E_pp"]
DIRECTION_LABELS = ["q", "E", "S", "q_p", "E_p", "S_p", "q_pp", "E_pp",
                    "S_pp", "Jq", "JE", "Jq_p", "JE_p", "C"]
LABELS = COMMON_LABELS + [f"{direction}.{field}" for direction in DIRECTIONS
                          for field in DIRECTION_LABELS]
ROOT_HALF_WIDTH = 0.01
MINIMUM_ESS = 100000
MATCH_ATOL, MATCH_RTOL = 2e-7, 2e-11


class NotScoreable(Exception):
    def __init__(self, stage, reason, diagnostic):
        super().__init__(reason)
        self.stage, self.reason, self.diagnostic = stage, reason, diagnostic


def path_label(path):
    try:
        return str(path.relative_to(ROOT))
    except ValueError:
        return str(path)


def read_previous_even(commit, n):
    path = f"results/p40-even-given-odd/raw/n{n}.csv"
    payload = git_bytes(commit, path)
    rows = {}
    for row in csv.DictReader(payload.decode("utf-8").splitlines()):
        key = (int(row["batch"]), row["orientation"])
        if key in rows:
            raise ValueError("duplicate row in previous E supplement")
        rows[key] = row
    expected = {(batch, direction) for batch in range(100) for direction in DIRECTIONS}
    if set(rows) != expected:
        raise ValueError("previous E supplement is not the full aligned block")
    return rows, {"commit": commit, "path": path, "sha256": digest(payload), "bytes": len(payload)}


def load_binned(path, run, old_counts, old_sums, old_grams, previous_commit):
    n = run["N"]
    previous, previous_input = read_previous_even(previous_commit, n)
    payload = path.read_bytes()
    required = {"n", "a", "b", "orientation", "batch", "k", "samples",
                "sum_q", "sum_e", "sum_s", "sum_qs", "sum_es"}
    reader = csv.DictReader(payload.decode("utf-8").splitlines())
    if not required.issubset(reader.fieldnames or []):
        raise ValueError(f"{path}: missing K-stratified fields")
    counts = np.zeros((100, 2, n + 1), dtype=np.int64)
    integer_sums = np.zeros((100, 2, n + 1, 5), dtype=np.int64)
    seen = set()
    for raw in reader:
        row = {key: (value if key == "orientation" else int(value))
               for key, value in raw.items() if key in required}
        batch, direction, k = row["batch"], row["orientation"], row["k"]
        key = (batch, direction, k)
        if key in seen or direction not in DIRECTIONS or not 0 <= batch < 100 or not 0 <= k <= n:
            raise ValueError(f"{path}: duplicate or invalid stratified key {key}")
        seen.add(key)
        if row["n"] != n or [row["a"], row["b"]] != run[direction] or row["samples"] < 0:
            raise ValueError(f"{path}: changed geometry or sample count at {key}")
        values = [row[name] for name in ("sum_q", "sum_e", "sum_s", "sum_qs", "sum_es")]
        if row["samples"] == 0 and any(values):
            raise ValueError(f"{path}: nonzero sums in empty K bin {key}")
        index = DIRECTIONS.index(direction)
        counts[batch, index, k] = row["samples"]
        integer_sums[batch, index, k] = values
    if len(seen) != 100 * 2 * (n + 1):
        raise ValueError(f"{path}: every K bin, including zeros, must be present")
    if not np.array_equal(counts[:, 0], counts[:, 1]):
        raise ValueError(f"{path}: directions do not share the archived occupation counts")

    sums = integer_sums.astype(float)
    sums[:, :, :, 2:] /= n
    maximum_error = 0.0
    k_values = np.arange(n + 1, dtype=float)
    for batch in range(100):
        for index, direction in enumerate(DIRECTIONS):
            count = int(counts[batch, index].sum())
            old_even = previous[(batch, direction)]
            if (count != 10000 or count != old_counts[batch]
                    or int(old_even["samples"]) != count or int(old_even["n"]) != n
                    or [int(old_even["a"]), int(old_even["b"])] != run[direction]):
                raise ValueError(f"{path}: changed source row at {batch}/{direction}")
            block = slice(index * 6, (index + 1) * 6)
            old_sum, old_gram = old_sums[batch, block], old_grams[batch, block, block]
            observed = np.concatenate([
                sums[batch, index].sum(axis=0),
                [counts[batch, index] @ k_values / n,
                 counts[batch, index] @ (k_values * (k_values - 1)) / (n * (n - 1))],
            ])
            expected = np.array([old_sum[0], old_gram[0, 0], old_sum[1], old_gram[0, 1],
                                 int(old_even["sum_e_s"]) / n, old_sum[2], old_sum[3]])
            maximum_error = max(maximum_error, float(np.max(np.abs(observed - expected))))
            if not np.allclose(observed, expected, atol=MATCH_ATOL, rtol=MATCH_RTOL):
                raise ValueError(f"{path}: K aggregation differs from archived moments at {batch}/{direction}")
            for column, atom in (("sum_q", 0), ("sum_e", 1), ("sum_s", 2), ("sum_e_s", 4)):
                if int(integer_sums[batch, index, :, atom].sum()) != int(old_even[column]):
                    raise ValueError(f"{path}: differs from exact previous {column} at {batch}/{direction}")
    diagnostic = {"aligned_batch_orientation_rows": 200, "K_rows_including_zeros": len(seen),
                  "checked_coordinates": [*ATOMS, "K/N", "KK/[N(N-1)]"],
                  "maximum_normalized_aggregate_error": maximum_error,
                  "absolute_tolerance": MATCH_ATOL, "relative_tolerance": MATCH_RTOL}
    inputs = [previous_input, {"path": path_label(path), "sha256": digest(payload), "bytes": len(payload)}]
    return counts, sums, inputs, diagnostic


def thermal_moments(counts, sums, p, p_star, n):
    """Ratio derivatives with both normalizer derivatives included analytically."""
    k = np.arange(n + 1, dtype=float)
    log_weights = k * math.log(p / p_star) + (n - k) * math.log((1 - p) / (1 - p_star))
    weights = np.exp(log_weights)
    score = k / p - (n - k) / (1 - p)
    score_prime = -k / (p * p) - (n - k) / ((1 - p) ** 2)
    weight_p = weights * score
    weight_pp = weights * (score * score + score_prime)
    z, zp, zpp = counts @ weights, counts @ weight_p, counts @ weight_pp
    if not np.isfinite([z, zp, zpp]).all() or z <= 0:
        raise NotScoreable("importance_weights", "nonfinite importance normalizer", {"p": p})
    a, ap, app = weights @ sums, weight_p @ sums, weight_pp @ sums
    mean = a / z
    first = (ap - mean * zp) / z
    second = (app - mean * zpp - 2 * first * zp) / z
    observed = counts > 0
    observed_weights = weights[observed]
    sample_count = int(counts.sum())
    ess = float(z * z / (counts @ (weights * weights)))
    diagnostic = {
        "ESS": ess, "ESS_fraction": ess / sample_count, "sample_count": sample_count,
        "weight_min_observed": float(observed_weights.min()),
        "weight_max_observed": float(observed_weights.max()),
        "weight_min_all_K": float(weights.min()), "weight_max_all_K": float(weights.max()),
        "maximum_sample_weight_fraction": float(observed_weights.max() / z),
        "maximum_K_bin_weight_fraction": float(np.max(counts * weights) / z),
        "normalizer_per_sample": float(z / sample_count),
        "normalizer_derivative_ratios": [float(zp / z), float(zpp / z)],
    }
    return mean, first, second, diagnostic


def direction_values(packet):
    mean, first, second, _ = packet
    q, e, s, qs, es = mean
    qp, ep, sp, qsp, esp = first
    jq, je = qs - q * s, es - e * s
    jq_p = qsp - qp * s - q * sp
    je_p = esp - ep * s - e * sp
    var_q = e - q * q
    if var_q <= 0:
        raise NotScoreable("rank_marginal", "nonpositive weighted q variance", {"var_q": float(var_q)})
    c = je - q * (1 - e) / var_q * jq
    return dict(zip(DIRECTION_LABELS, [q, e, s, qp, ep, sp, *second[:3], jq, je, jq_p, je_p, c]))


def point(counts, sums, p_star, n, delta):
    bracket = (p_star - ROOT_HALF_WIDTH, p_star + ROOT_HALF_WIDTH)

    def pooled_q(p):
        return float(sum(thermal_moments(counts[i], sums[i], p, p_star, n)[0][0]
                         for i in range(2)) / 2)

    endpoint_values = [pooled_q(p) for p in bracket]
    root_diagnostic = {"bracket": list(bracket), "pooled_q_at_endpoints": endpoint_values}
    if not np.isfinite(endpoint_values).all() or endpoint_values[0] * endpoint_values[1] > 0:
        raise NotScoreable("pooled_root", "pooled q root is not bracketed in the declared interval", root_diagnostic)
    try:
        p0 = float(brentq(pooled_q, *bracket, xtol=5e-14, rtol=5e-14, maxiter=100))
    except (ValueError, RuntimeError) as error:
        raise NotScoreable("pooled_root", str(error), root_diagnostic) from error
    packets = [thermal_moments(counts[i], sums[i], p0, p_star, n) for i in range(2)]
    root_diagnostic.update({"p0": p0, "importance": {d: packets[i][3] for i, d in enumerate(DIRECTIONS)}})
    if any(packet[3]["ESS"] < MINIMUM_ESS for packet in packets):
        raise NotScoreable("importance_ESS", "root importance ESS is below 100000", root_diagnostic)
    rows = [direction_values(packet) for packet in packets]
    d = float(sum(row["q_p"] for row in rows) / 2)
    root_diagnostic["positive_local_slope_D"] = d
    if not math.isfinite(d) or d <= 0:
        raise NotScoreable("pooled_root_slope", "pooled root has nonpositive or nonfinite slope", root_diagnostic)
    jq_bar = float(sum(row["Jq"] for row in rows) / 2)
    jq_p_bar = float(sum(row["Jq_p"] for row in rows) / 2)
    q_pp_bar = float(sum(row["q_pp"] for row in rows) / 2)
    b = float((rows[0]["E_p"] - rows[1]["E_p"]) / delta)
    p4_je_p = float((rows[0]["JE_p"] - rows[1]["JE_p"]) / delta)
    p4_e_pp = float((rows[0]["E_pp"] - rows[1]["E_pp"]) / delta)
    rootdot = -jq_bar / d
    prefactor = n ** (13 / 8) / 2
    pieces = [prefactor * p4_je_p / d, prefactor * p4_e_pp * rootdot / d,
              -prefactor * b * jq_p_bar / (d * d),
              -prefactor * b * q_pp_bar * rootdot / (d * d)]
    udot = math.fsum(pieces)
    common = [udot, *pieces, prefactor * b / d, p0, rootdot, d, b,
              float(sum(row["q"] for row in rows) / 2), jq_bar, jq_p_bar,
              q_pp_bar, p4_je_p, p4_e_pp]
    vector = np.array(common + [row[field] for row in rows for field in DIRECTION_LABELS], dtype=float)
    if not np.isfinite(vector).all():
        raise NotScoreable("root_response", "nonfinite root response", {"p0": p0})
    root_diagnostic["Udot_four_term_sum_roundoff"] = float(udot - sum(pieces))
    root_diagnostic["atomic_moments"] = {
        direction: {"mean": packets[i][0].tolist(), "p_derivative": packets[i][1].tolist(),
                    "p_second_derivative": packets[i][2].tolist()}
        for i, direction in enumerate(DIRECTIONS)
    }
    return vector, root_diagnostic


def analyze_size(counts, sums, n, p_star, delta):
    total_counts, total_sums = counts.sum(axis=0), sums.sum(axis=0)
    try:
        central, central_diagnostic = point(total_counts, total_sums, p_star, n, delta)
    except NotScoreable as error:
        return {"status": "not_scoreable", "not_scoreable": {field: error.reason for field in LABELS},
                "failure_stage": error.stage, "diagnostics": error.diagnostic,
                "covariance": None, "delete_one_vectors": None}

    loo, loo_diagnostics, failures = [], [], []
    for batch in range(100):
        try:
            vector, diagnostic = point(total_counts - counts[batch], total_sums - sums[batch], p_star, n, delta)
            loo.append(vector.tolist())
            loo_diagnostics.append({"omitted_batch": batch, **diagnostic})
        except NotScoreable as error:
            loo.append(None)
            failures.append({"omitted_batch": batch, "stage": error.stage,
                             "reason": error.reason, "diagnostic": error.diagnostic})
    if failures:
        return {"status": "not_scoreable",
                "not_scoreable": {"full_covariance_and_all_standard_errors": "one or more aligned delete-one roots violate the declared root/ESS boundary", "joint_Udot_zero": "complete scoreable delete-one covariance is unavailable"},
                "central_point_estimates_without_scoreable_uncertainty": dict(zip(LABELS, central.tolist())),
                "diagnostics": central_diagnostic, "delete_one_failures": failures,
                "covariance": None, "delete_one_vectors": loo, "delete_one_diagnostics": loo_diagnostics}

    loo_array = np.asarray(loo)
    deviations = loo_array - loo_array.mean(axis=0)
    covariance = 99 / 100 * deviations.T @ deviations
    errors = np.sqrt(np.maximum(0, np.diag(covariance)))
    estimates = {field: {"value": float(value), "se": float(error),
                         "z": float(value / error) if error > 0 else None}
                 for field, value, error in zip(LABELS, central, errors)}
    roots = [item["p0"] for item in loo_diagnostics]
    return {"status": "scoreable", "estimates": estimates, "covariance": covariance.tolist(),
            "delete_one_vectors": loo, "diagnostics": central_diagnostic,
            "delete_one_diagnostics": loo_diagnostics, "delete_one_root_range": [min(roots), max(roots)],
            "minimum_delete_one_ESS": min(item["importance"][direction]["ESS"]
                                          for item in loo_diagnostics for direction in DIRECTIONS)}


def analyze(output_dir):
    start = time.perf_counter()
    manifest_text = MANIFEST.read_text(encoding="utf-8")
    manifest = json.loads(manifest_text)
    source_manifest = Path(manifest.get("source_analysis_manifest", "analysis/p40_absolute_cluster_reanalysis.json"))
    if not source_manifest.is_absolute():
        source_manifest = ROOT / source_manifest
    contract = json.loads(source_manifest.read_text(encoding="utf-8"))
    if manifest["source_commit"] != contract["input_commit"]:
        raise ValueError("binned replay and old Gram source commits differ")
    if sorted(run["N"] for run in contract["runs"]) != [65, 85]:
        raise ValueError("only original N65/N85 blocks are supported")
    p_star = float(contract["p"])
    previous_commit = manifest.get("previous_even_commit", PREVIOUS_EVEN_COMMIT)
    backend_path = manifest.get("source_backend_path", "src/gaussian_orientation_mc.cpp")
    backend = git_bytes(contract["input_commit"], backend_path)
    backend_blob = hashlib.sha1(f"blob {len(backend)}\0".encode() + backend).hexdigest()
    if manifest.get("source_backend_blob") not in (None, backend_blob):
        raise ValueError("replay source backend blob differs from declaration")
    inputs = [{"commit": contract["input_commit"], "path": backend_path,
               "sha256": digest(backend), "git_blob_sha1": backend_blob}]
    results = {}
    for run in contract["runs"]:
        old_counts, old_sums, old_grams, old_inputs, metadata = load_input(contract, run)
        n = run["N"]
        counts, sums, new_inputs, agreement = load_binned(
            output_dir / "raw" / f"n{n}.csv", run, old_counts, old_sums, old_grams, previous_commit)
        delta_exact = cos4(*run["first"]) - cos4(*run["second"])
        result = analyze_size(counts, sums, n, p_star, float(delta_exact))
        result.update({"delta_cos4_exact": str(delta_exact), "source_agreement": agreement,
                       "parent_metadata": metadata, "samples_per_orientation": int(counts[:, 0].sum())})
        results[str(n)] = result
        inputs.extend(old_inputs + new_inputs)
    if all(result["status"] == "scoreable" and result["estimates"]["Udot"]["z"] is not None
           for result in results.values()):
        statistic = sum(result["estimates"]["Udot"]["z"] ** 2 for result in results.values())
        joint = {"Udot": {"status": "scoreable", "chi_square": statistic, "df": 2,
                            "nominal_p": math.exp(-statistic / 2)}}
    else:
        joint = {"Udot": {"status": "not_scoreable", "reason": "both N blocks require complete scoreable root/ESS/delete-one uncertainty"}}

    code_paths = [Path(__file__), ROOT / "scripts/analyze_p40_absolute_cluster.py", MANIFEST, source_manifest]
    runner = ROOT / "src/p40_source_thermal_replay.cpp"
    if runner.exists():
        code_paths.append(runner)
    return {
        "schema": "matching-one.p40-source-thermal.v1", "labels": LABELS,
        "atomic_coordinate_order": ATOMS, "by_N": results, "joint_zero": joint,
        "inputs": inputs, "manifest": manifest, "manifest_text": manifest_text,
        "source_contract": contract, "previous_even_commit": previous_commit,
        "definitions": {
            "source": "S=(CB+CW)/N is the common raw source; no geometry-adapted clock/full projection",
            "importance_weight": "(p/p_star)^K ((1-p)/(1-p_star))^(N-K)",
            "ratio_derivatives": "mu=A/Z; mu_p=(A_p-mu Z_p)/Z; mu_pp=(A_pp-mu Z_pp-2 mu_p Z_p)/Z",
            "moment_convention": "self-normalized empirical tilted moments; no ad hoc sample-covariance Bessel correction in the ratio derivatives",
            "root": "pooled q=0 within declared p_star +/- .01 with positive local slope; not a full-interval or global-unique-root proof",
            "U": "N^(13/8)/2 * P4[E_p]/mean(q_p), at the locally bracketed pooled matching root",
            "rootdot": "-mean(Jq)/D, Jq=mean(qS)-mean(q)mean(S)",
            "Udot_terms": "direct + rootmotion + slope_source + slope_root; all four terms share one source and full covariance",
            "C": "JE-q*(1-E)/(E-q*q)*Jq in the weighted finite three-sector marginal",
        },
        "thermal_clock_null": "If the same-N source is common S=a+bK in both geometries, its tilt only shifts Bernoulli log-odds; refinding the matching root makes numerator and denominator Jacobians cancel, so Udot=0. The nominal joint zero statistic tests this necessary condition, not all K^2 clocks, geometry-dependent clocks, or CFT mechanisms; compatibility does not prove the affine-K model.",
        "scoreability": {"root_bracket": [p_star - ROOT_HALF_WIDTH, p_star + ROOT_HALF_WIDTH],
                          "minimum_ESS": MINIMUM_ESS, "scope": "central and every aligned delete-one root; report failure without a substituted target"},
        "uncertainty": "100 aligned delete-one batches per N, refind pooled root in each omission, full output covariance and all vectors; N-domain-separated PRNG independence for the single nominal two-N Udot zero statistic; four additive terms are not independent evidence",
        "estimator_boundary": "Same original U definition but a different and noisier fixed-p self-normalized importance estimator, not a replacement or reproduction of threshold-histogram Rao-Blackwell results; N65/N85 alone are not the complete65-130-260 or85-170-340 norm4 lineages",
        "source_deformation_boundary": "K-binned first source moments support lambda=0 first derivatives only; exp(lambda*mean_S_given_K) is not finite-lambda configuration reweighting. No assumed q2/Jordan source-deformation law, fitted exponent, finite difference or energy identity.",
        "new_samples": 0, "oldconfigurationreplays": 2000000, "scientific_test_suites": [],
        "code": [{"path": path_label(path), "sha256": digest(path.read_bytes())} for path in code_paths],
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__, "machine": platform.machine()},
        "elapsed_seconds": time.perf_counter() - start,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results/p40-source-thermal")
    args = parser.parse_args()
    destination = args.output_dir / "latest.json"
    if destination.exists():
        raise ValueError("refusing to overwrite the saved source-thermal result")
    result = analyze(args.output_dir)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with destination.open("x", encoding="utf-8") as handle:
        handle.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"joint_zero": result["joint_zero"], "elapsed_seconds": result["elapsed_seconds"],
                      "by_N_status": {n: row["status"] for n, row in result["by_N"].items()}}, indent=2))


if __name__ == "__main__":
    main()

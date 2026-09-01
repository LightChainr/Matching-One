#!/usr/bin/env python3
"""Evaluate the fixed source-jet split on P40's independent old million streams.

Reuse saved full/delete-one roots and raw U derivatives. No root solving,
counter replay, compensator selection, or additional probability-point sample.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
DIRECTORY = ROOT / "results/p40-source-thermal"
PARENT = DIRECTORY / "latest.json"
DESTINATION = DIRECTORY / "jet-split.json"
DEFINITION = ROOT / "notes/p40-thermal-clock-source-quotient.md"
NS = (65, 85)
DIRECTIONS = ("first", "second")
STATES = np.array([-1.0, 0.0, 1.0])
FIELDS = ("raw_bulk", "within_sector_bulk", "topology_only_bulk")
SAVED_FIELDS = ("p0", "D", "U", "Udot", "first.q", "first.E", "second.q", "second.E")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_profile(path, n, run):
    counts = np.zeros((100, 2, n + 1), dtype=np.int64)
    sums = np.zeros((100, 2, n + 1, 5), dtype=float)
    seen = set()
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            batch, k = int(row["batch"]), int(row["k"])
            g = DIRECTIONS.index(row["orientation"])
            key = (batch, g, k)
            if (key in seen or not 0 <= batch < 100 or not 0 <= k <= n
                    or int(row["n"]) != n
                    or [int(row["a"]), int(row["b"])] != run[DIRECTIONS[g]]):
                raise ValueError(f"invalid archived K row: {path} {key}")
            seen.add(key)
            counts[batch, g, k] = int(row["samples"])
            sums[batch, g, k] = [int(row[field]) for field in
                                 ("sum_q", "sum_e", "sum_s", "sum_qs", "sum_es")]
    if len(seen) != 100 * 2 * (n + 1) or np.any(counts < 0):
        raise ValueError(f"incomplete archived K table: {path}")
    if not np.all(counts.sum(axis=2) == 10000):
        raise ValueError(f"changed P40 batch weights: {path}")
    sums[:, :, :, 2:] /= n
    return counts, sums


def conditional_table(counts, sums, p, p_star, n):
    """Importance-integrated conditional moments using S=s/N, without S^2."""
    k = np.arange(n + 1, dtype=float)
    log_weights = k * math.log(p / p_star) + (n - k) * math.log((1 - p) / (1 - p_star))
    log_weights -= float(log_weights.max())
    weights = np.exp(log_weights)
    normalizer = float(weights @ counts)
    if not math.isfinite(normalizer) or normalizer <= 0:
        raise ValueError(f"N{n}: invalid importance normalization")
    q, e, s, qs, es = sums.T
    sector_counts = np.column_stack(((e - q) / 2, counts - e, (e + q) / 2))
    sector_sums = np.column_stack(((es - qs) / 2, s - es, (es + qs) / 2))
    weighted_counts = weights @ sector_counts
    if np.any(weighted_counts <= 0):
        raise ValueError(f"N{n}: all three sectors are needed for the fixed quotient")
    probabilities = weighted_counts / normalizer
    mean_k = weights @ (sector_counts * k[:, None]) / weighted_counts
    mean_s = weights @ sector_sums / weighted_counts
    dk = k[:, None] - mean_k[None, :]
    var_k = weights @ (sector_counts * dk * dk) / weighted_counts
    cov_sk = weights @ ((sector_sums - sector_counts * mean_s[None, :]) * dk) / weighted_counts
    return {"probabilities": probabilities, "mean_K": mean_k,
            "mean_S_density": mean_s, "var_K": var_k, "cov_S_density_K": cov_sk}


def at_saved_root(counts, sums, state, n, p_star, delta):
    p0, d, u = state["p0"], state["D"], state["U"]
    tables = [conditional_table(counts[g], sums[g], p0, p_star, n) for g in range(2)]
    numerator = math.fsum(float(t["probabilities"] @ t["cov_S_density_K"]) for t in tables)
    denominator = math.fsum(float(t["probabilities"] @ t["var_K"]) for t in tables)
    if denominator < 0 or not math.isfinite(denominator) or d <= 0:
        raise ValueError(f"N{n}: invalid saved slope or conditional clock variance")
    if denominator == 0 and numerator != 0:
        raise ValueError(f"N{n}: zero clock variance with nonzero clock covariance")
    beta = numerator / denominator if denominator > 0 else 0.0
    jq_p, je_p, directions = [], [], {}
    for g, table in zip(DIRECTIONS, tables):
        gamma = table["cov_S_density_K"] - beta * table["var_K"]
        probabilities = table["probabilities"]
        jq_p.append(float(probabilities @ ((STATES - state[f"{g}.q"]) * gamma) / (p0 * (1 - p0))))
        je_p.append(float(probabilities @ ((STATES**2 - state[f"{g}.E"]) * gamma) / (p0 * (1 - p0))))
        directions[g] = {
            "sector_order": [-1, 0, 1],
            **{key: value.tolist() for key, value in table.items()},
            "gamma_density": gamma.tolist(),
            "topology_counterterm_density":
                (table["mean_S_density"] - beta * table["mean_K"]).tolist(),
            "Jq_prime_R_density": jq_p[-1], "JE_prime_R_density": je_p[-1],
        }
    direct = n ** (13 / 8) / 2 * (je_p[0] - je_p[1]) / (delta * d)
    slope = -u * math.fsum(jq_p) / (2 * d)
    within = n * math.fsum((direct, slope))
    raw = n * state["Udot"]
    values = dict(zip(FIELDS, (raw, within, raw - within)))
    if not all(math.isfinite(value) for value in values.values()):
        raise ValueError(f"N{n}: nonfinite split")
    return values, {
        "saved_state": state, "b_star_density": beta, "b_star_bulk": n * beta,
        "common_clock_projection_numerator": numerator,
        "common_clock_projection_denominator": denominator,
        "directions": directions,
        "within_sector_pieces_bulk": {"direct": n * direct, "slope_source": n * slope,
                                      "rootmotion": 0.0, "slope_root": 0.0},
    }


def vectorize(points):
    return {f"N{n}.{field}": points[n][field] for n in NS for field in FIELDS}


def main():
    started = time.perf_counter()
    if DESTINATION.exists():
        raise ValueError("saved P40 jet split exists; do not overwrite or repeat")
    previous = json.loads(PARENT.read_text())
    old_labels = previous["labels"]
    old_indices = {field: old_labels.index(field) for field in SAVED_FIELDS}
    runs = {run["N"]: run for run in previous["manifest"]["runs"]}
    p_star = float(previous["manifest"]["p_ref"])
    input_by_path = {item["path"]: item for item in previous["inputs"]}
    raw, totals, saved_full, saved_loo, root_mapping, inputs = {}, {}, {}, {}, {}, []
    for n in NS:
        row = previous["by_N"][str(n)]
        if row["status"] != "scoreable":
            raise ValueError(f"N{n}: the saved parent root family is not scoreable")
        path = DIRECTORY / "raw" / f"n{n}.csv"
        relative = str(path.relative_to(ROOT))
        digest = sha(path)
        if digest != input_by_path[relative]["sha256"]:
            raise ValueError(f"changed raw input since saved roots: {path}")
        inputs.append({"path": relative, "sha256": digest})
        raw[n] = read_profile(path, n, runs[n])
        totals[n] = tuple(array.sum(axis=0) for array in raw[n])
        saved_full[n] = {field: float(row["estimates"][field]["value"]) for field in SAVED_FIELDS}
        diagnostics = row["delete_one_diagnostics"]
        vectors = row["delete_one_vectors"]
        if len(diagnostics) != 100 or len(vectors) != 100:
            raise ValueError(f"N{n}: incomplete saved delete-one mapping")
        indices = {int(item["omitted_batch"]): i for i, item in enumerate(diagnostics)}
        if set(indices) != set(range(100)):
            raise ValueError(f"N{n}: ambiguous saved batch/root mapping")
        saved_loo[n], root_mapping[n] = {}, []
        for batch in range(100):
            index = indices[batch]
            vector = vectors[index]
            state = {field: float(vector[column]) for field, column in old_indices.items()}
            if state["p0"] != diagnostics[index]["p0"]:
                raise ValueError(f"N{n}: parent vector and root diagnostic disagree")
            saved_loo[n][batch] = state
            root_mapping[n].append({"omitted_batch": batch, "parent_vector_index": index,
                                    "parent_saved_state": state})

    deltas = {n: float(Fraction(previous["by_N"][str(n)]["delta_cos4_exact"])) for n in NS}
    central_by_n, central_diagnostics = {}, {}
    for n in NS:
        central_by_n[n], central_diagnostics[n] = at_saved_root(
            *totals[n], saved_full[n], n, p_star, deltas[n])
    central_map = vectorize(central_by_n)
    labels = list(central_map)
    central = np.array(list(central_map.values()))
    covariance = np.zeros((len(labels), len(labels)))
    groups = {}
    for n in NS:
        vectors, betas = [], []
        for batch in range(100):
            omitted = dict(central_by_n)
            counts, sums = (totals[n][i] - raw[n][i][batch] for i in range(2))
            omitted[n], diagnostic = at_saved_root(
                counts, sums, saved_loo[n][batch], n, p_star, deltas[n])
            betas.append(diagnostic["b_star_density"])
            vectors.append(list(vectorize(omitted).values()))
        vectors = np.asarray(vectors)
        deviations = vectors - vectors.mean(axis=0)
        component = 99 / 100 * deviations.T @ deviations
        covariance += component
        groups[f"P40-N{n}-N-keyed-million"] = {
            "Ns": [n], "delete_one_batch_ids": list(range(100)),
            "delete_one_vectors": vectors.tolist(), "covariance_contribution": component.tolist(),
            "b_star_density_delete_one": betas, "saved_root_mapping": root_mapping[n],
            "operation": "Remove the same batch in both directions, use its saved parent root/U/D/q/E, recompute conditional moments and common b_star, hold the other N at its central value",
        }
    errors = np.sqrt(np.maximum(0, np.diag(covariance)))
    estimates = {label: {"value": float(value), "se": float(error),
                          "z": float(value / error) if error > 0 else None}
                 for label, value, error in zip(labels, central, errors)}
    result = {
        "schema": "matching-one.p40-source-jet-split.v1",
        "status": "computed_independent_existing_stream_contrast",
        "labels": labels, "estimates": estimates, "covariance": covariance.tolist(),
        "by_N": central_diagnostics, "dependency_groups": groups,
        "definitions": {
            "source": "S=(CB+CW)/N; all displayed responses are bulk v_N=N*dU/dlambda_density",
            "within_source": "R_g=S_g-E_g[S_g|q]-b_star*(K-E_g[K|q])",
            "common_b_star": "sum_g E_g[Cov(S,K|q)] / sum_g E_g[Var(K|q)], equal direction weights",
            "within_response": "W_bulk=N*(N^(13/8)*P4[JE_prime(R)]/(2D)-U*mean(Jq_prime(R))/D)",
            "source_identity": "S_g=R_g+b_star*K+f_g(q), with fixed f_g(q)=E_g[S_g|q]-b_star*E_g[K|q]; common K is annihilated by the root-normalized readout",
            "remainder": "topology_only_bulk=v-W; a decomposition of one source, not another independently measured source",
            "derivative_convention": "At every saved root refit conditional functions/b_star, then freeze them for the p/lambda derivative; R has zero q/E first response at that point, not zero thermal mixed response",
        },
        "comparison_role": "The same already-defined W on an independent old P40 random stream after the NZ N85 hint was known. This is an archived-stream contrast, not prospective confirmation of 2.37 SE, fresh Monte Carlo, or an independent vote from each component.",
        "estimator": "Fixed-p Bernoulli samples reweighted with the exact likelihood ratio and a self-normalizing denominator, unlike the complete-prefix NZ Binomial estimator. Saved central/100-delete-one roots and raw U derivatives are reused without solving them again.",
        "uncertainty": "N65 and N85 are separate N-keyed PRNG domains; sum their separate aligned-delete-one contributions, other N held central. Full six-coordinate covariance is redundant by raw=within+topology and is not inverted.",
        "scope": "Each size uses one million old configurations and the same two directions; only two parents, not full norm4 lineages. A nonzero W concerns within-sector cluster/occupancy coupling in this thermal jet, not field identity. No S^2, source variance, efficiency, attribution percentage, or component omnibus test.",
        "parent_result": {"path": str(PARENT.relative_to(ROOT)), "sha256": sha(PARENT),
                          "schema": previous["schema"], "code": previous["code"]},
        "inputs": inputs,
        "code": [{"path": str(path.relative_to(ROOT)), "sha256": sha(path)}
                 for path in (Path(__file__).resolve(), DEFINITION)],
        "environment": {"python": platform.python_version(), "executable": sys.executable,
                        "numpy": np.__version__, "machine": platform.machine()},
        "elapsed_seconds": time.perf_counter() - started,
        "new_samples": 0, "configuration_replays": 0, "root_solves": 0,
        "server_actions": 0, "test_suites": [],
    }
    with DESTINATION.open("x") as handle:
        handle.write(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"output": str(DESTINATION), "elapsed_seconds": result["elapsed_seconds"],
                      "estimates": estimates}, indent=2))


if __name__ == "__main__":
    main()

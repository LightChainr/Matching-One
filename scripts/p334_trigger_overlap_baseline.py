#!/usr/bin/env python3
"""Read real checkpoint sufficient statistics against a fixed-edge graph null.

No new MC: the source is the frozen P334 cooperative pilot. The graph null
forgets physical geometry but preserves each checkpoint's numbers of safe
vertices and minimal triggering pairs. No null graph simulation is needed.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import subprocess

import numpy as np

SOURCE = "e81dd59ff6be69056e504e0e81cfeccf73dc5e97"
PREFIX = "results/local-20260831/P334-cooperative-closure/raw"
HASHES = {325: "266e6c3c225cb507ba11c45acaa55c5eb5e132692a37ca3baccb82702a459c03",
          425: "e8120d7c717a28e950ce9d4deb69398c7c3774c5f34e8ce948ab288a988f55b2"}
METRICS = ("delta_observed", "delta_fixed_edge_null", "delta_excess",
           "degree_variance_observed", "degree_variance_fixed_edge_null",
           "two_stars_observed", "two_stars_fixed_edge_null", "two_stars_excess",
           "safe_vertices", "minimal_pair_triggers")


def load_blob(path):
    return subprocess.check_output(["git", "cat-file", "blob", f"{SOURCE}:{path}"])


def graph_statistics(a, b2_safe, sum_child_degree_sq, vacancies):
    """Complement graph edges are minimal rank-two-triggering pairs."""
    slots = a * (a-1) // 2
    triggers = slots-b2_safe
    assert 0 <= triggers <= slots
    # c_v = (a-1)-t_v, with t_v the trigger-graph degree.
    sum_trigger_sq = a*(a-1)**2 - 4*(a-1)*b2_safe + sum_child_degree_sq
    twice_stars = sum_trigger_sq-2*triggers
    assert twice_stars >= 0 and twice_stars % 2 == 0
    stars = twice_stars//2
    variance_num = a*sum_child_degree_sq-(2*b2_safe)**2
    assert variance_num >= 0
    observed_var = variance_num/(a*a) if a else 0.0
    # A vertex degree is Hypergeometric(C(a,2), a-1, triggers).
    # Its mean is fixed across all graphs; hence this is E[population variance].
    p = triggers/slots if slots else 0.0
    null_var = (a-1)**2*p*(1-p)/(a+1) if a > 1 else 0.0
    null_stars = (a*(a-1)*(a-2)/2 * triggers*(triggers-1)
                  / (slots*(slots-1))) if slots > 1 else 0.0
    weight = a/vacancies/(vacancies-1)**2
    observed_delta = weight*observed_var
    null_delta = weight*null_var
    excess = observed_delta-null_delta
    star_excess = stars-null_stars
    assert abs(excess-2*star_excess/(vacancies*(vacancies-1)**2)) < 2e-14
    return np.array([observed_delta, null_delta, excess, observed_var, null_var,
                     stars, null_stars, star_excess, a, triggers])


def score_size(size):
    path = f"{PREFIX}/N{size}.geometry_pilot.csv"
    payload = load_blob(path)
    digest = hashlib.sha256(payload).hexdigest()
    assert digest == HASHES[size]
    metadata_bytes = load_blob(f"{PREFIX}/N{size}.metadata.json")
    meta = json.loads(metadata_bytes)
    n = meta["samples_per_pair"]
    width = len(METRICS)
    x = np.zeros((n, 2, width))
    eligible = np.zeros((n, 2))
    above = np.zeros(2, dtype=int)
    for row in csv.DictReader(io.StringIO(payload.decode())):
        i = int(row["replica"])-meta["replica_counter_first"]
        o = ("first", "second").index(row["orientation"])
        assert 0 <= i < n and eligible[i, o] == 0
        d = int(row["n"])-int(row["k0"])
        a = int(row["checkpoint_b1_safe_count"])
        assert a == d-int(row["H2"])
        x[i, o] = graph_statistics(a, int(row["checkpoint_b2_safe_pairs"]),
                                   int(row["checkpoint_sum_child_b1_sq"]), d)
        eligible[i, o] = 1
        above[o] += int(x[i, o, 2] > 1e-15)
    counts = eligible.sum(0)
    mean = x.sum(0)/counts[:, None]
    # Retain one cluster per base permutation, with both orientations aligned.
    influence = (x-eligible[:, :, None]*mean)/(counts/n)[None, :, None]
    flat = influence.reshape(n, 2*width)
    cov = flat.T @ flat/(n*(n-1))
    se = np.sqrt(np.diag(cov)).reshape(2, width)
    ratio_jac = np.zeros((2, 2*width))
    ratios = mean[:, 0]/mean[:, 1]
    for o in range(2):
        ratio_jac[o, o*width] = 1/mean[o, 1]
        ratio_jac[o, o*width+1] = -mean[o, 0]/mean[o, 1]**2
    ratio_cov = ratio_jac @ cov @ ratio_jac.T
    excess_indices = [2, width+2]
    excess_cov = cov[np.ix_(excess_indices, excess_indices)]
    excess = mean[:, 2]
    stat = float(excess @ np.linalg.solve(excess_cov, excess))
    rows = {}
    for o, name in enumerate(("first", "second")):
        rows[name] = {
            "at_risk_checkpoints": int(counts[o]), "above_null_checkpoints": int(above[o]),
            "fraction_above_null": float(above[o]/counts[o]),
            "metrics": {key: {"estimate": float(mean[o, j]), "cluster_se": float(se[o, j])}
                        for j, key in enumerate(METRICS)},
            "ratio_observed_over_null": float(ratios[o]),
            "ratio_cluster_se": float(np.sqrt(ratio_cov[o, o])),
            "excess_z": float(excess[o]/se[o, 2])}
    return {"N": size, "new_samples": 0,
            "source": {"commit": SOURCE, "path": path, "sha256": digest,
                       "metadata_sha256": hashlib.sha256(metadata_bytes).hexdigest(),
                       "runner": meta["git_commit"], "seed": meta["seed"],
                       "counter_first": meta["replica_counter_first"],
                       "counter_last_exclusive": meta["replica_counter_last_exclusive"],
                       "k0": meta["geometry_pilot_k0"]},
            "dependency_group": f"p334-cooperative-N{size}-20260831",
            "base_permutation_clusters": n, "orientations": rows,
            "vector_order": [f"{o}:{key}" for o in ("first", "second") for key in METRICS],
            "mean_vector": mean.reshape(-1).tolist(), "cluster_covariance": cov.tolist(),
            "ratio_covariance": ratio_cov.tolist(),
            "paired_excess_score": {"chi2": stat, "df": 2,
                                    "asymptotic_log_tail_probability": -stat/2}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=Path("results/p334-trigger-overlap-baseline/latest.json"))
    args = parser.parse_args()
    sizes = [score_size(n) for n in (325, 425)]
    # Descriptive contrast only: not a fitted scaling law or held-out prediction.
    difference = np.array([sizes[1]["orientations"][o]["ratio_observed_over_null"]
                           -sizes[0]["orientations"][o]["ratio_observed_over_null"]
                           for o in ("first", "second")])
    difference_cov = sum(np.array(s["ratio_covariance"]) for s in sizes)
    result = {"schema": "matching-one.p334-trigger-overlap-baseline.v1",
              "lifecycle": "retrospective mechanism analysis of fixed production archives",
              "new_samples": 0, "sizes": sizes,
              "two_size_ratio_difference_descriptive": {
                  "N425_minus_N325": difference.tolist(),
                  "covariance": difference_cov.tolist(),
                  "standard_error": np.sqrt(np.diag(difference_cov)).tolist()},
              "null": "Uniform simple trigger graph G(a,m), conditional on each checkpoint a=b1 and m=C(a,2)-b2_safe; analytic expectation, no graph Monte Carlo.",
              "claim": "Trigger-pair shared-endpoint overlap exceeds this fixed-edge exchangeable baseline. Not a continuum field, causal process memory, or an asymptotic scaling law."}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print(json.dumps({"sizes": [{"N": s["N"], "orientations": s["orientations"],
                                "paired_score": s["paired_excess_score"]} for s in sizes],
                      "size_difference": result["two_size_ratio_difference_descriptive"]}, indent=2))


if __name__ == "__main__":
    main()

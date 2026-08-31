#!/usr/bin/env python3
"""Score exact cooperative branching excess with replica/checkpoint clusters."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

METRICS = (
    "delta_coop_rb", "delta_coop_common", "closure_baseline",
    "branch_success_rb", "s1", "s2", "common_mean_identity",
    "common_square_calibration", "clone1_calibration", "clone2_calibration",
    "clone_product_calibration",
)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def score_size(prefix, contract, commit):
    csv_path = Path(str(prefix) + ".geometry_pilot.csv")
    meta_path = Path(str(prefix) + ".metadata.json")
    meta = json.loads(meta_path.read_text())
    expected = {
        "git_commit": commit, "cooperative_closure": True,
        "branching_clones": True, "samples_per_pair": contract["samples"],
        "batches": contract["batches"], "seed": contract["seed"],
        "replica_counter_first": contract["replica_counter_first"],
        "replica_counter_last_exclusive": contract["replica_counter_last_exclusive"],
        "geometry_pilot_k0": contract["k0"],
    }
    for field, value in expected.items():
        if meta[field] != value:
            raise ValueError(f"{prefix}: metadata {field}: {meta[field]} != {value}")
    n = contract["samples"]
    m = len(METRICS)
    values = np.zeros((n, 2, m))
    risk = np.zeros((n, 2))
    positive = np.zeros(2, dtype=int)
    exact_numerator_max = np.zeros(2, dtype=np.int64)
    audit = {"rows": 0, "identity_failures": 0, "negative_exact_numerators": 0}
    with csv_path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            r = int(row["replica"]) - contract["replica_counter_first"]
            o = ("first", "second").index(row["orientation"])
            if not 0 <= r < n or risk[r, o]:
                raise ValueError("duplicate/out-of-range checkpoint cluster")
            d = int(row["n"]) - int(row["k0"])
            b1 = int(row["checkpoint_b1_safe_count"])
            b2 = int(row["checkpoint_b2_safe_pairs"])
            squares = int(row["checkpoint_sum_child_b1_sq"])
            qnum = int(row["branch_q_after_safe_count"])
            qden = int(row["branch_q_after_denominator"])
            common = int(row["branch_common_safe"])
            y1 = int(row["branch_clone1_survives"])
            y2 = int(row["branch_clone2_survives"])
            both = int(row["branch_both_survive"])
            numerator = b1 * squares - (2 * b2) ** 2
            if (b1 != d-int(row["H2"]) or qden != d-1 or
                    both != y1*y2 or not 0 <= qnum <= qden or
                    (not common and (qnum or y1 or y2)) or
                    abs(float(row["q_after"])-qnum/qden) > 2e-15 or
                    abs(float(row["q_after2"])-(qnum/qden)**2) > 2e-15):
                audit["identity_failures"] += 1
            audit["negative_exact_numerators"] += int(numerator < 0)
            positive[o] += int(numerator > 0)
            exact_numerator_max[o] = max(exact_numerator_max[o], numerator)
            q = qnum / qden
            s1 = b1/d
            s2 = 2*b2/(d*(d-1))
            baseline = s2*s2/s1 if b1 else 0.0
            exact_q2 = squares/(d*(d-1)**2)
            exact_delta = numerator/(d*b1*(d-1)**2) if b1 else 0.0
            values[r, o] = (
                exact_delta, q*q-baseline, baseline, exact_q2, s1, s2,
                q-s2, q*q-exact_q2, y1-q, y2-q, both-q*q,
            )
            risk[r, o] = 1
            audit["rows"] += 1
    if audit["identity_failures"] or audit["negative_exact_numerators"]:
        raise ValueError(f"exact row identities failed: {audit}")
    counts = risk.sum(axis=0)
    means = values.sum(axis=0)/counts[:, None]
    # One influence vector per base permutation. Missing (non-risk) orientation
    # contributes zero, rather than being an independent unpaired sample.
    influence = (values-risk[:, :, None]*means[None, :, :])/(counts/n)[None, :, None]
    flat = influence.reshape(n, 2*m)
    cov = flat.T @ flat / (n*(n-1))
    se = np.sqrt(np.maximum(np.diag(cov), 0.0)).reshape(2, m)
    indices = [0, m]
    delta = means[:, 0]
    block = cov[np.ix_(indices, indices)]
    chi2 = float(delta @ np.linalg.pinv(block, rcond=1e-12) @ delta)
    environments = {}
    for o, label in enumerate(("first", "second")):
        environments[label] = {
            "at_risk_checkpoints": int(counts[o]),
            "risk_rate": float(counts[o]/n),
            "exact_positive_checkpoints": int(positive[o]),
            "positive_fraction_given_risk": float(positive[o]/counts[o]),
            "max_exact_variance_numerator": int(exact_numerator_max[o]),
            "metrics": {name: {"estimate": float(means[o, j]),
                                "se": float(se[o, j]),
                                "z": float(means[o, j]/se[o, j]) if se[o, j] else None}
                        for j, name in enumerate(METRICS)},
            "relative_excess_over_closure_baseline": float(means[o, 0]/means[o, 2]),
        }
    return {
        "inputs": {"csv_sha256": sha256(csv_path), "metadata_sha256": sha256(meta_path)},
        "checkpoint_clusters": n, "audit": audit, "environments": environments,
        "vector_order": [f"{o}:{name}" for o in ("first", "second") for name in METRICS],
        "vector": means.reshape(-1).tolist(), "checkpoint_cluster_covariance": cov.tolist(),
        "primary": {"delta_vector": delta.tolist(), "covariance": block.tolist(),
                    "wald_chi2_2": chi2,
                    "nominal_chi2_tail": float(np.exp(-chi2/2)),
                    "both_positive": bool(np.all(delta > 0)),
                    "gate_passed": bool(np.all(delta > 0) and chi2 > 9.210340371976184)},
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--freeze", type=Path, required=True)
    parser.add_argument("--n325-prefix", type=Path, required=True)
    parser.add_argument("--n425-prefix", type=Path, required=True)
    parser.add_argument("--runner-commit", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    freeze = json.loads(args.freeze.read_text())
    runner = Path(__file__).parents[1]/"src/threshold_rank_integer_period_mc.cpp"
    if sha256(runner) != freeze["runner_source_sha256"]:
        raise ValueError("runner source differs from pre-data freeze")
    sizes = {size: score_size(prefix, freeze["runs"][size], args.runner_commit)
             for size, prefix in (("N325", args.n325_prefix), ("N425", args.n425_prefix))}
    payload = {"schema": "matching-one/p334-cooperative-closure-score/v1",
               "runner_commit": args.runner_commit, "freeze_sha256": sha256(args.freeze),
               "sizes": sizes, "stop_decision": "stop_at_20k_no_extension",
               "both_sizes_pass": all(x["primary"]["gate_passed"] for x in sizes.values()),
               "claim_boundary": freeze["claim_boundary"]}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n")
    print(json.dumps({s: x["primary"] for s, x in sizes.items()}, indent=2))


if __name__ == "__main__":
    main()

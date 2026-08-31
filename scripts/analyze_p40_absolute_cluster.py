#!/usr/bin/env python3
"""Use existing P40 production Gram matrices; never generate a configuration."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import subprocess
import time
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "analysis/p40_absolute_cluster_reanalysis.json"
COORDS = ["q", "S", "K/N", "KK/[N(N-1)]", "T/(2N)", "chi/N"]
METRICS = ["A", "S", "J_raw", "J_clock", "J_full", "var_S", "fraction_clock", "fraction_full"]
H4 = ["A", "S", "J_raw", "J_clock", "J_full"]
ALLOC = [f"{mode}.{part}" for mode in ("raw", "clock", "full") for part in ("readout", "source")]
LABELS = [f"H4.{x}" for x in H4] + [f"H4.{x}" for x in ALLOC] + [f"{d}.{x}" for d in ("first", "second") for x in METRICS]


def digest(content):
    return hashlib.sha256(content).hexdigest()


def git_bytes(commit, path):
    return subprocess.check_output(["git", "show", f"{commit}:{path}"], cwd=ROOT)


def transform(names, n):
    idx = {key: i for i, key in enumerate(names)}
    t = np.zeros((6, len(names)))
    t[0, idx["q"]] = 1
    t[1, idx["C_black"]] = t[1, idx["C_white"]] = 1 / n
    t[2, idx["V"]] = 1 / n
    t[3, idx["E"]] = 1 / (2 * n)
    t[3, idx["E_mc"]] = -1 / (2 * n)
    t[4, idx["E"]] = 1 / (2 * n)
    t[5, idx["V"]] = t[5, idx["F0"]] = 1 / n
    t[5, idx["E"]] = -1 / n
    return t


def load_input(contract, run):
    n, commit = run["N"], contract["input_commit"]
    prefix = f'{contract["input_directory"]}/N{n}/mc'
    inputs = []
    raw = {}
    for suffix, key in (("motifs.jsonl", "motifs_sha256"), ("metadata.json", "metadata_sha256")):
        path = f"{prefix}.{suffix}"
        blob = git_bytes(commit, path)
        if digest(blob) != run[key]:
            raise ValueError(f"source hash differs: {path}")
        raw[suffix] = blob
        inputs.append({"commit": commit, "path": path, "sha256": digest(blob), "bytes": len(blob)})
    meta = json.loads(raw["metadata.json"])
    rows = [json.loads(line) for line in raw["motifs.jsonl"].splitlines()]
    if len(rows) != 100 or [r["batch"] for r in rows] != list(range(100)):
        raise ValueError("expected 100 ordered aligned batches")
    if meta["samples_per_pair"] != 1000000 or not meta["cross_geometry_joint_gram"]:
        raise ValueError("unexpected production metadata")
    t = transform(rows[0]["names"], n)
    counts, sums, grams = [], [], []
    for row in rows:
        if row["names"] != rows[0]["names"] or row["n"] != n or row["samples"] != 10000 or abs(row["p_ref"] - contract["p"]) > 1e-14:
            raise ValueError("changed batch semantics")
        for d in ("first", "second"):
            if [row[d]["a"], row[d]["b"]] != run[d]:
                raise ValueError("changed orientation order")
        sf, ss = [t @ np.asarray(row[d]["sum"], dtype=float) for d in ("first", "second")]
        ff, tt = [t @ np.asarray(row[d]["gram"], dtype=float) @ t.T for d in ("first", "second")]
        ft = t @ np.asarray(row["cross_gram"], dtype=float) @ t.T
        counts.append(row["samples"])
        sums.append(np.concatenate([sf, ss]))
        grams.append(np.block([[ff, ft], [ft.T, tt]]))
    return np.array(counts), np.array(sums), np.array(grams), inputs, meta


def single(mean, cov):
    diagnostics = {}
    vectors = {"raw": np.array([0., 1., 0., 0., 0., 0.])}
    for mode, indices in (("clock", [2, 3]), ("full", [2, 3, 4, 5])):
        zcov = cov[np.ix_(indices, indices)]
        sd = np.sqrt(np.diag(zcov))
        corr = zcov / np.outer(sd, sd)
        eig = np.linalg.eigvalsh(corr)
        if eig[0] <= 1e-10:
            raise ValueError("declared control span numerically unresolved")
        beta = np.linalg.solve(corr, cov[indices, 1] / sd) / sd
        vector = vectors["raw"].copy()
        vector[indices] -= beta
        vectors[mode] = vector
        diagnostics[mode] = {"coefficients": beta.tolist(), "controls": [COORDS[i] for i in indices],
                             "condition_correlation": float(eig[-1] / eig[0]),
                             "source_response_to_controls": (cov[np.ix_(indices, range(6))] @ vector).tolist()}
    response = [cov[0] @ vectors[mode] for mode in ("raw", "clock", "full")]
    fractions = [vectors[mode] @ cov @ vectors[mode] / cov[1, 1] for mode in ("clock", "full")]
    return np.array([mean[0], mean[1], *response, cov[1, 1], *fractions]), vectors, diagnostics


def point(count, sums, gram, delta):
    mean = sums / count
    cov = (gram - np.outer(sums, sums) / count) / (count - 1)
    first, fv, fd = single(mean[:6], cov[:6, :6])
    second, sv, sd = single(mean[6:], cov[6:, 6:])
    h4 = (first[:5] - second[:5]) / delta
    qf, qs = np.eye(12)[0], np.eye(12)[6]
    allocations = []
    for mode in ("raw", "clock", "full"):
        sf = np.concatenate([fv[mode], np.zeros(6)])
        ss = np.concatenate([np.zeros(6), sv[mode]])
        # Exact paired bilinear identity, not two independent causal effects.
        allocations += [float((qf - qs) @ cov @ ((sf + ss) / 2) / delta),
                        float(((qf + qs) / 2) @ cov @ (sf - ss) / delta)]
    diagnostics = {"first": fd, "second": sd,
                   "maximum_allocation_roundoff": float(np.max(np.abs(np.asarray(allocations).reshape(3, 2).sum(axis=1) - h4[2:5])))}
    return np.concatenate([h4, allocations, first, second]), diagnostics


def cos4(a, b):
    return Fraction(a**4 - 6*a*a*b*b + b**4, (a*a + b*b)**2)


def analyze():
    start = time.perf_counter()
    contract = json.loads(MANIFEST.read_text())
    results, inputs = {}, []
    for run in contract["runs"]:
        counts, sums, grams, source, meta = load_input(contract, run)
        total_count, total_sum, total_gram = counts.sum(), sums.sum(axis=0), grams.sum(axis=0)
        delta_exact = cos4(*run["first"]) - cos4(*run["second"])
        delta = float(delta_exact)
        center, diagnostic = point(total_count, total_sum, total_gram, delta)
        loo = np.array([point(total_count - count, total_sum - s, total_gram - g, delta)[0]
                        for count, s, g in zip(counts, sums, grams)])
        deviations = loo - loo.mean(axis=0)
        covariance = 99 / 100 * deviations.T @ deviations
        se = np.sqrt(np.maximum(0, np.diag(covariance)))
        estimates = {key: {"value": float(value), "se": float(error), "z": float(value / error) if error > 0 else None}
                     for key, value, error in zip(LABELS, center, se)}
        results[str(run["N"])] = {"estimates": estimates, "covariance": covariance.tolist(),
            "delete_one_vectors": loo.tolist(), "control_diagnostics": diagnostic,
            "delta_cos4_exact": str(delta_exact), "parent_metadata": meta,
            "transformed_sample_sum": total_sum.tolist(), "transformed_sample_gram": total_gram.tolist()}
        inputs += source
    joint = {}
    for field in ("J_raw", "J_clock", "J_full"):
        stat = sum(r["estimates"][f"H4.{field}"]["z"]**2 for r in results.values())
        joint[field] = {"chi_square": stat, "df": 2, "nominal_p": math.exp(-stat / 2)}
    return {"schema": "matching-one.p40-absolute-cluster-result.v1", "labels": LABELS,
        "coordinate_order_per_orientation": COORDS, "by_N": results, "joint_zero": joint,
        "inputs": inputs, "source_contract": contract,
        "pair_allocation": "Delta Cov(q,S)=Cov(Delta q, mean S)+Cov(mean q, Delta S), also for separately projected sources; an exploratory coupling-dependent coordinate allocation, not a causal decomposition or independent evidence",
        "not_scoreable": {"E_top_source_response": "q squared times S is a third moment absent from the Gram", "N130_norm4_child": "P40 has N65 and N85 only"},
        "uncertainty": "100 aligned delete-one batches per N, all nonlinear source projections refit; N-separated PRNG block assumption for joint zero; nominal Gaussian estimates, not exact confidence certificates",
        "new_samples": 0, "configuration_replays": 0, "scientific_test_suites": [],
        "code": [{"path": str(p.relative_to(ROOT)), "sha256": digest(p.read_bytes())} for p in (Path(__file__), MANIFEST)],
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "machine": platform.machine()},
        "elapsed_seconds": time.perf_counter() - start}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results/p40-absolute-cluster")
    args = parser.parse_args()
    destination = args.output_dir / "latest.json"
    if destination.exists():
        raise ValueError("refusing to overwrite the saved result")
    result = analyze()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"joint_zero": result["joint_zero"], "elapsed_seconds": result["elapsed_seconds"],
                      "by_N": {n: r["estimates"] for n, r in result["by_N"].items()}}, indent=2))


if __name__ == "__main__":
    main()

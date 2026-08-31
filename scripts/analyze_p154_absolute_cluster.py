#!/usr/bin/env python3
"""Absolute cluster source on the inherited Phase-E blocks; no new samples."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
import subprocess
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "analysis/p154_absolute_cluster_replay.json"
KEYS = ["q", "e", "k", "kk", "edges", "s", "chi", "cb", "cw"]
METRICS = ["A", "E", "S", "J_A_S", "J_E_S", "J_A_clock", "J_E_clock",
           "J_A_clock_euler", "J_E_clock_euler", "var_S", "var_after_clock",
           "var_after_clock_euler", "fraction_after_clock", "fraction_after_clock_euler",
           "fraction_after_topology"]
H4 = METRICS[:9]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    with path.open(newline="") as stream:
        raw = list(csv.DictReader(stream))
    rows = {(int(r["batch"]), r["orientation"]):
            {k: v if k == "orientation" else int(v) for k, v in r.items()} for r in raw}
    if len(rows) != len(raw):
        raise ValueError(f"duplicate batch/orientation: {path}")
    return rows


def total(rows, direction, omit=None):
    selected = [r for (b, d), r in rows.items() if d == direction and b != omit]
    return {k: sum(r[k] for r in selected) for k in selected[0]
            if k.startswith("sum_") or k == "samples"}


def point(sums, n):
    count = sums["samples"]
    mean = np.array([sums[f"sum_{k}"] / count for k in KEYS])
    second = np.zeros((len(KEYS), len(KEYS)))
    for i, k in enumerate(KEYS):
        for j in range(i, len(KEYS)):
            second[i, j] = second[j, i] = sums[f"sum_{k}_{KEYS[j]}"] / count
    # Fixed scale changes coordinates, never the declared control span.
    scale = np.array([1, 1, n, n * (n - 1), 2 * n, n, n, n, n])
    mean = mean / scale
    cov = (second - np.outer(mean * scale, mean * scale)) / np.outer(scale, scale)
    cov *= count / (count - 1)
    source = 5
    diagnostics = {}

    def residual(indices, label):
        zcov = cov[np.ix_(indices, indices)]
        sd = np.sqrt(np.diag(zcov))
        corr = zcov / np.outer(sd, sd)
        eigenvalues = np.linalg.eigvalsh(corr)
        if eigenvalues[0] <= 1e-10:
            raise ValueError(f"{label} control span numerically unresolved")
        beta = np.linalg.solve(corr, cov[indices, source] / sd) / sd
        response = cov[:2, source] - cov[np.ix_([0, 1], indices)] @ beta
        variance = cov[source, source] - cov[source, indices] @ beta
        diagnostics[label] = {"indices": [KEYS[i] for i in indices],
                              "coefficients_scaled_coordinates": beta.tolist(),
                              "correlation_condition": float(eigenvalues[-1] / eigenvalues[0])}
        return response, float(variance)

    clock, vclock = residual([2, 3], "clock")
    local, vlocal = residual([2, 3, 4, 6], "clock_euler")
    _, vtop = residual([0, 1], "topology_diagnostic_only")
    # q/E are never included in the clock/Euler residual source.
    values = [*mean[:2], mean[source], *cov[:2, source], *clock, *local,
              cov[source, source], vclock, vlocal, vclock / cov[source, source],
              vlocal / cov[source, source], vtop / cov[source, source]]
    # The exact three-state readout projection limits interpretation of any source.
    probs = np.array([(mean[1] - mean[0]) / 2, 1 - mean[1], (mean[1] + mean[0]) / 2])
    raw_qs = sums["sum_q_s"] / count / n
    raw_es = sums["sum_e_s"] / count / n
    sector_source = np.array([(raw_es - raw_qs) / 2,
                              mean[source] - raw_es, (raw_es + raw_qs) / 2]) / probs
    b = (sector_source[2] - sector_source[0]) / 2
    c = (sector_source[2] + sector_source[0]) / 2 - sector_source[1]
    identity = cov[:2, :2] @ np.array([b, c])
    diagnostics["sector_projection"] = {
        "probabilities_minus_zero_plus": probs.tolist(),
        "mean_S_minus_zero_plus": sector_source.tolist(), "b": float(b), "c": float(c),
        "maximum_response_identity_roundoff": float(np.max(np.abs(identity - cov[:2, source])))
    }
    return np.array(values), diagnostics


def cos4(a, b):
    return Fraction(a**4 - 6*a*a*b*b + b**4, (a*a + b*b)**2)


def vector(rows, n, delta, omit=None):
    first, fd = point(total(rows, "first", omit), n)
    second, sd = point(total(rows, "second", omit), n)
    projected = (first[:len(H4)] - second[:len(H4)]) / delta
    return np.concatenate([projected, first, second]), {"first": fd, "second": sd}


def analyze(output):
    contract = json.loads(MANIFEST.read_text())
    labels = [f"H4.{x}" for x in H4] + [f"{d}.{x}" for d in ("first", "second") for x in METRICS]
    results, inputs = {}, []
    for run in contract["runs"]:
        n = run["N"]
        paths = [output / "raw" / f"n{n}.csv",
                 ROOT / f"results/p154-phase-e-mixed-plane-pilot/raw/n{n}_mixed.batches.csv",
                 ROOT / f"results/p154-fixed-k-interaction/raw/n{n}.csv"]
        rows, original, old_edges = map(read, paths)
        if sha(paths[1]) != run["parent_sha256"] or sha(paths[2]) != run["edge_replay_sha256"]:
            raise ValueError("parent artifact differs from declared source")
        expected = {(b, d) for b in range(100) for d in ("first", "second")}
        if any(set(r) != expected for r in (rows, original, old_edges)):
            raise ValueError("missing aligned batch/orientation")
        identity = ["n", "a", "b", "batch", "samples", "sum_k1", "sum_k2", "sum_i0", "sum_i1", "sum_i2"]
        for key in expected:
            if any(rows[key][k] != original[key][k] for k in identity):
                raise ValueError(f"changed original Phase-E counters: {n}/{key}")
            if any(rows[key][k] != v for k, v in old_edges[key].items()):
                raise ValueError(f"changed original edge replay: {n}/{key}")
            if rows[key]["samples"] != 200:
                raise ValueError("unexpected batch weight")
        for batch in range(100):
            if any(rows[(batch, "first")][f"sum_{k}"] != rows[(batch, "second")][f"sum_{k}"] for k in ("k", "kk")):
                raise ValueError("unpaired occupancy count")
        delta_exact = cos4(*run["first"]) - cos4(*run["second"])
        central, diagnostics = vector(rows, n, float(delta_exact))
        loo = np.array([vector(rows, n, float(delta_exact), b)[0] for b in range(100)])
        centered = loo - loo.mean(axis=0)
        covariance = 99 / 100 * centered.T @ centered
        errors = np.sqrt(np.maximum(0, np.diag(covariance)))
        estimates = {label: {"value": float(value), "se": float(se),
                             "z": float(value/se) if se > 0 else None}
                     for label, value, se in zip(labels, central, errors)}
        results[str(n)] = {"estimates": estimates, "delta_cos4_exact": str(delta_exact),
                           "covariance": covariance.tolist(), "delete_one_vectors": loo.tolist(),
                           "control_diagnostics": diagnostics}
        inputs += [{"path": str(p.relative_to(ROOT)), "sha256": sha(p)} for p in paths]
    joint = {}
    for name in ("J_E_S", "J_E_clock", "J_E_clock_euler"):
        statistic = sum(results[str(n)]["estimates"][f"H4.{name}"]["z"]**2 for n in (65, 130))
        joint[name] = {"chi_square": statistic, "df": 2, "nominal_p": float(chi2.sf(statistic, 2))}
    return {"schema": "matching-one.p154-absolute-cluster-result.v1", "labels": labels,
            "by_N": results, "joint_E_zero": joint, "inputs": inputs,
            "scope": contract["scope"], "limitations": contract["interpretation_limits"],
            "dependency_groups": [f"P154-PhaseE-N{n}-original-20k" for n in (65, 130)],
            "uncertainty": "estimated aligned-delete-one covariance; nominal Gaussian summaries, not exact confidence certificates; full matrix is redundant and is not inverted",
            "code": [{"path": str(p.relative_to(ROOT)), "sha256": sha(p)} for p in
                     (Path(__file__), MANIFEST, ROOT / "src/p154_absolute_cluster_replay.cpp", ROOT / "src/threshold_rank_integer_period_mc.cpp")],
            "environment": {"python": platform.python_version(), "machine": platform.machine(),
                            "numpy": np.__version__, "platform": platform.platform()}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results/p154-absolute-cluster")
    args = parser.parse_args()
    path = args.output_dir / "latest.json"
    if path.exists():
        raise ValueError("refusing to overwrite a saved result")
    backend = subprocess.check_output(["git", "hash-object", "src/threshold_rank_integer_period_mc.cpp"], cwd=ROOT, text=True).strip()
    if backend != json.loads(MANIFEST.read_text())["backend_blob"]:
        raise ValueError("original backend blob changed")
    result = analyze(args.output_dir)
    path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"joint_E_zero": result["joint_E_zero"],
                      "by_N": {n: {k: v for k, v in r["estimates"].items() if k.startswith("H4.J_") or "fraction" in k}
                               for n, r in result["by_N"].items()}}, indent=2))


if __name__ == "__main__":
    main()

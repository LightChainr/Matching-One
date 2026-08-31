#!/usr/bin/env python3
"""One-activation-lagged spatial source from new marks on old permutations."""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import time

import numpy as np
import scipy

from analyze_norm4_source_endpoint_1m import load_profile
from analyze_norm4_global_source_projection import source_response

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis/norm4_lagged_source_contract.json"
OUTPUT = ROOT / "results/norm4-lagged-source"
NS = (65, 85, 130, 170, 260, 340)
FIELDS = ("v", "rootdot", "rank1_rootdot")
TYPES = ("01", "02", "12")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_events(path, n, per_batch, run):
    out = np.zeros((100, 2, n + 1, 6), dtype=np.float64)
    seen = set()
    fields = ["event_count" + t for t in TYPES] + ["sum_s_previous" + t for t in TYPES]
    with gzip.open(path, "rt") as handle:
        for row in csv.DictReader(handle):
            b, k = int(row["batch"]), int(row["k"])
            g = ("first", "second").index(row["orientation"])
            key = (b, g, k)
            if (key in seen or not 0 <= b < 100 or not 0 <= k <= n
                    or int(row["n"]) != n or int(row["samples"]) != per_batch
                    or [int(row["a"]), int(row["b"])] != run[("first", "second")[g]]):
                raise ValueError(f"unaligned lagged event row {path}: {key}")
            seen.add(key)
            out[b, g, k] = [int(row[f]) for f in fields]
    if len(seen) != 100*2*(n+1) or np.any(out < 0) or np.any(out[:, :, 0]):
        raise ValueError(f"incomplete or invalid lagged event archive: {path}")
    return out


def point(baseline, events, samples, n, p0):
    # Existing source loader stores s/N. Early conditional means use bulk s.
    bulk = baseline.copy()
    bulk[..., 2:] *= n
    q, e, s, qs, es = np.moveaxis(bulk, -1, 0)
    den0, den1 = e - q, samples - e
    mu0 = np.divide(es - qs, den0, out=np.zeros_like(q), where=den0 > 0)
    mu1 = np.divide(s - es, den1, out=np.zeros_like(q), where=den1 > 0)
    mu0prev = np.pad(mu0[:, :-1], ((0, 0), (1, 0)))
    mu1prev = np.pad(mu1[:, :-1], ((0, 0), (1, 0)))
    for index, denominator in ((0, den0), (1, den0), (2, den1)):
        unsupported = np.pad(denominator[:, :-1] <= 0, ((0, 0), (1, 0)), constant_values=True)
        if np.any(events[..., index][unsupported] != 0):
            raise ValueError("birth event has no supported previous-rank source layer")
    kernels = events[..., 3:].copy()
    kernels[..., 0] -= mu0prev * events[..., 0]
    kernels[..., 1] -= mu0prev * events[..., 1]
    kernels[..., 2] -= mu1prev * events[..., 2]
    hq = kernels[..., 0] + 2*kernels[..., 1] + kernels[..., 2]
    he = -kernels[..., 0] + kernels[..., 2]
    result = {"p0": float(p0)}
    pairs = [("total", hq, he),
             ("event01", kernels[..., 0], -kernels[..., 0]),
             ("event02", 2*kernels[..., 1], np.zeros_like(hq)),
             ("event12", kernels[..., 2], kernels[..., 2])]
    for name, qsource, esource in pairs:
        profile = bulk.copy()
        profile[..., 2] = 0
        profile[..., 3], profile[..., 4] = qsource, esource
        response, u = source_response(profile, samples, n, p0)
        result.update({name + "." + key: response[key] for key in FIELDS})
        result["U"] = u
    result["first_activation.rootdot"] = result["event01.rootdot"] + .5*result["event02.rootdot"]
    result["second_completion.rootdot"] = result["event12.rootdot"] + .5*result["event02.rootdot"]
    return result, kernels / samples


def vectorize(points):
    return {f"N{n}.{key}": value for n in NS for key, value in points[n].items()}


def addback(points):
    return max(abs(row["total." + field] - sum(row["event" + t + "." + field] for t in TYPES))
               for row in points.values() for field in FIELDS)


def main():
    started = time.perf_counter()
    destination = OUTPUT / "latest.json"
    if destination.exists():
        raise ValueError("Preserve this saved analysis; reproduce in a separate output copy")
    contract = json.loads(CONTRACT.read_text())
    source_path = ROOT / contract["source_result"]
    if digest(source_path) != contract["source_result_sha256"]:
        raise ValueError("changed source result / matching-root reference")
    source = json.loads(source_path.read_text())
    source_index = {label: i for i, label in enumerate(source["labels"])}
    geometry_path = ROOT / "analysis/p40_source_thermal_chain_candidates.json"
    runs = {r["N"]: r for r in json.loads(geometry_path.read_text())["runs"]}
    receipt = json.loads((OUTPUT / "run.json").read_text())
    if receipt["status"] != "completed":
        raise ValueError("lagged source marking has no completed receipt")
    profiles, events, samples, totals, event_totals = {}, {}, {}, {}, {}
    inputs, points, kernels = [], {}, {}
    matching_residual = 0.0
    for n in NS:
        original = ROOT / f"results/norm4-source-thermal/raw/n{n}.csv"
        profiles[n] = load_profile(original, n, 1000, runs[n])
        inputs.append({"path": str(original.relative_to(ROOT)), "sha256": digest(original)})
        samples[n] = 100000
        if n in (260, 340):
            increment = ROOT / f"results/norm4-source-endpoint-1m/increment/raw/n{n}.csv"
            profiles[n] += load_profile(increment, n, 9000, runs[n])
            inputs.append({"path": str(increment.relative_to(ROOT)), "sha256": digest(increment)})
            samples[n] = 1000000
        path = OUTPUT / f"raw/n{n}.csv.gz"
        per_batch = samples[n] // 100
        events[n] = load_events(path, n, per_batch, runs[n])
        inputs.append({"path": str(path.relative_to(ROOT)), "sha256": digest(path)})
        counts = events[n][..., :3]
        rebuilt_q = -per_batch + np.cumsum(counts[..., 0] + 2*counts[..., 1] + counts[..., 2], axis=-1)
        rebuilt_e = per_batch + np.cumsum(-counts[..., 0] + counts[..., 2], axis=-1)
        matching_residual = max(matching_residual,
                                float(np.max(np.abs(rebuilt_q - profiles[n][..., 0]))),
                                float(np.max(np.abs(rebuilt_e - profiles[n][..., 1]))))
        if matching_residual != 0:
            raise ValueError("event marks do not reconstruct the original paired batch q/E profiles")
        totals[n], event_totals[n] = profiles[n].sum(axis=0), events[n].sum(axis=0)
        p0 = source["by_N"][str(n)]["points"]["p0"]
        points[n], kernels[n] = point(totals[n], event_totals[n], samples[n], n, p0)
    central_map = vectorize(points)
    labels, central = list(central_map), np.array(list(central_map.values()))
    covariance = np.zeros((len(labels), len(labels)))
    groups = {}
    addback_residual = addback(points)
    for name, previous in source["covariance_contributions"].items():
        if not name.startswith("source:"):
            continue
        if previous["delete_one_batch_ids"] != list(range(100)):
            raise ValueError("changed source batch deletion order")
        saved = np.asarray(previous["delete_one_vectors"])
        vectors = []
        for b in range(100):
            changed = dict(points)
            for n in previous["Ns"]:
                p0 = saved[b, source_index[f"N{n}.p0"]]
                changed[n], _ = point(totals[n]-profiles[n][b], event_totals[n]-events[n][b],
                                      samples[n]*.99, n, p0)
            addback_residual = max(addback_residual, addback(changed))
            vectors.append(list(vectorize(changed).values()))
        vectors = np.asarray(vectors)
        factor = np.sqrt(.99)*(vectors-vectors.mean(axis=0))
        contribution = factor.T @ factor
        covariance += contribution
        groups[name] = {"Ns": previous["Ns"], "delete_one_batch_ids": list(range(100)),
                        "delete_one_vectors": vectors.tolist(), "factor": factor.tolist(),
                        "operation": "same paired event and early-source omission; reestimate all early-rank means and source responses at the saved retained-sample root"}
    se = np.sqrt(np.maximum(0, covariance.diagonal()))
    estimates = {label: {"value": float(value), "se": float(error),
                         "z": float(value/error) if error else None}
                 for label, value, error in zip(labels, central, se)}
    if not np.isfinite(central).all() or not np.isfinite(covariance).all() or addback_residual > 1e-8:
        raise ValueError("nonfinite lagged response or failed additive reconstruction")
    result = {
        "schema": "matching-one.norm4-lagged-source.v1", "status": "computed_old_permutation_new_joint_marks",
        "execution_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "contract": contract, "inputs": inputs, "samples_per_N": samples,
        "code_sha256": digest(Path(__file__)), "contract_sha256": digest(CONTRACT),
        "source_result_sha256": digest(source_path), "run_receipt_sha256": digest(OUTPUT / "run.json"),
        "labels": labels, "estimates": estimates, "covariance": covariance.tolist(),
        "covariance_contributions": groups,
        "by_N": {str(n): {"points": points[n], "event_kernels_by_geometry_K_type": kernels[n].tolist(),
                          "event_counts_by_geometry_K_type": event_totals[n][..., :3].astype(int).tolist()}
                 for n in NS},
        "identities": {"batch_event_to_original_q_E_max_residual": matching_residual,
                       "central_and_LOO_total_minus_exclusive_events_max_residual": addback_residual},
        "environment": {"python": platform.python_version(), "platform": platform.platform(),
                        "numpy": np.__version__, "scipy": scipy.__version__},
        "elapsed_seconds": time.perf_counter()-started,
        "new_random_samples": 0, "old_permutations_with_new_lagged_marks": sum(samples.values()),
        "root_finders": 0, "server_operations": 0, "test_suites": 0,
    }
    destination.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    print(json.dumps({"seconds": result["elapsed_seconds"], "dimension": len(labels),
                      "identities": result["identities"],
                      "results": {str(n): {key: estimates[f"N{n}.{key}"] for key in
                          ("total.rank1_rootdot", "total.rootdot", "total.v", "first_activation.rootdot", "second_completion.rootdot")}
                                  for n in NS}}, indent=2))


if __name__ == "__main__":
    main()

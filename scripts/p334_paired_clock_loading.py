#!/usr/bin/env python3
"""Whole-archive paired rank-one loading, with exact conditional suffix means."""
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import gzip
import hashlib
import json
from math import comb
from pathlib import Path
import signal
import time

import numpy as np
from scipy.stats import binom

from p334_checkpoint_scalar_collision import archived_permutation
from p334_contracted_birth_network import build
from p334_full_birth_reliability import safety_polynomial
from p334_pair_only_survival import multiply

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis/p334_paired_clock_loading.json"
OUT = ROOT / "results/p334-paired-clock-loading"


def timed_out(_signal, _frame):
    raise TimeoutError("fixed whole-pair prefix-only computation budget")


def read_sources(contract):
    sources = {}
    for n in contract["sizes"]:
        base = ROOT / contract["source_directory"] / f"N{n}"
        meta = json.loads(base.with_suffix(".metadata.json").read_text())
        csv_path = base.with_suffix(".geometry_pilot.csv")
        rows = {}
        with csv_path.open() as stream:
            for raw in csv.DictReader(stream):
                row = {k: (v if k == "orientation" else float(v) if k in ("q_after", "q_after2") else int(v)) for k, v in raw.items()}
                rows[(row["replica"], row["orientation"])] = row
        sources[n] = {"metadata": meta, "rows": rows, "source_sha256": hashlib.sha256(csv_path.read_bytes()).hexdigest()}
    return sources


def clock_for(n, counter, metadata, row, permutation, contract):
    k0 = row["k0"]
    matrix = metadata["designs"][0][row["orientation"] + "_period_matrix"]
    mapping = build({"N": n, "replica_counter": counter, "seed": metadata["seed"], "k0": k0,
                     "occupied_prefix_labels": permutation[:k0], "period_matrix": matrix,
                     "ell": [row["ell_u"], row["ell_v"]]})
    coefficients, core, width, peak = [1], set(), 0, 0
    networks = [c["two_terminal_network"] for c in mapping["port_components"] if "two_terminal_network" in c]
    if not networks or any(len(c["addresses"]) > 2 for c in mapping["port_components"]):
        raise ValueError("no supported whole-event two-port factorization")
    for network in networks:
        sites = set(network["vacant_sites"])
        if sites & core:
            raise ValueError("overlapping random variables across factors")
        core |= sites
        fc, stats = safety_polynomial(network, sites, max_states=contract["maximum_dp_states"])
        coefficients = multiply(coefficients, fc)
        width = max(width, stats["treewidth_upper_bound"])
        peak = max(peak, stats["maximum_states"])
    d = n - k0
    free = d - len(core)
    coefficients = multiply(coefficients, [comb(free, k) for k in range(free + 1)])
    if coefficients[:3] != [1, d - row["H2"], row["checkpoint_b2_safe_pairs"]]:
        raise ValueError("generalized geometry differs from original singleton/pair fields")
    survival = np.array([c / comb(d, k) for k, c in enumerate(coefficients)])
    birth = survival[:-1] - survival[1:]
    thresholds = k0 + np.arange(1, d + 1)
    return {"conditional": [float(birth @ binom.sf(thresholds - 1, n, contract["p_ref"])),
                            float(birth @ ((n - thresholds + 1) / (n + 1)))],
            "safe_coefficients": coefficients, "mean_wait": float(survival[:-1].sum()),
            "factors": len(networks), "core_sites": len(core), "treewidth": width,
            "peak_states": peak, "physical_line": mapping["physical_line"]}


def process_batch(task):
    n, batch, metadata, selected_rows, contract = task
    target = OUT / "batches" / f"N{n}.batch{batch:02d}.json.gz"
    if target.exists():
        return n, batch, "reused"
    signal.signal(signal.SIGALRM, timed_out)
    first = metadata["replica_counter_first"] + batch * 1000
    rows = {(r["replica"], r["orientation"]): r for r in selected_rows}
    records = []
    started = time.monotonic()
    for counter in range(first, first + 1000):
        originals = [rows.get((counter, o)) for o in ("first", "second")]
        baseline = np.zeros(4)
        for o, row in enumerate(originals):
            if row:
                baseline[o] = binom.sf(row["k2"] - 1, n, contract["p_ref"])
                baseline[2 + o] = (n - row["k2"] + 1) / (n + 1)
        record = {"counter": counter, "risk": [int(r is not None) for r in originals],
                  "source_rows": originals, "X": baseline.tolist(), "Y": baseline.tolist(),
                  "status": "outside_rank_one", "clocks": [None, None]}
        if any(originals):
            tick = time.monotonic()
            signal.setitimer(signal.ITIMER_REAL, contract["pair_wall_seconds"])
            try:
                permutation = archived_permutation(n, metadata["seed"], counter)
                conditional = np.zeros(4)
                for o, row in enumerate(originals):
                    if row:
                        clock = clock_for(n, counter, metadata, row, permutation, contract)
                        record["clocks"][o] = clock
                        conditional[o], conditional[2 + o] = clock["conditional"]
                record["Y"] = conditional.tolist()
                record["status"] = "exact_pair"
            except (TimeoutError, RuntimeError, ValueError) as error:
                record["status"] = "whole_pair_fallback"
                record["reason"] = str(error)
            finally:
                signal.setitimer(signal.ITIMER_REAL, 0)
            record["seconds"] = time.monotonic() - tick
        records.append(record)
    payload = {"N": n, "batch": batch, "source_commit": contract["source_commit"],
               "contract_sha256": hashlib.sha256(CONTRACT.read_bytes()).hexdigest(),
               "wall_seconds": time.monotonic() - started, "records": records}
    with target.open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as stream:
            stream.write((json.dumps(payload, separators=(",", ":")) + "\n").encode())
    counts = {name: sum(r["status"] == name for r in records) for name in ("exact_pair", "whole_pair_fallback", "outside_rank_one")}
    print(f"N{n} batch{batch:02d}: {counts}, seconds={payload['wall_seconds']:.2f}", flush=True)
    return n, batch, counts


def h4(rep):
    x, y = rep
    return (x**4 - 6*x*x*y*y + y**4) / (x*x + y*y)**2


def score(contract, sources):
    output = {"contract": contract, "sizes": {}}
    report = ["# Full-archive paired rank-one-stratum loading", "", contract["interpretation"], ""]
    for n, source in sources.items():
        xs, ys, risks, statuses = [], [], [], []
        for batch in range(20):
            with gzip.open(OUT / "batches" / f"N{n}.batch{batch:02d}.json.gz", "rt") as stream:
                rows = json.load(stream)["records"]
            xs.append([r["X"] for r in rows]); ys.append([r["Y"] for r in rows])
            risks.extend(r["risk"] for r in rows); statuses.extend(r["status"] for r in rows)
        x, y = np.array(xs), np.array(ys)
        design = source["metadata"]["designs"][0]
        delta = h4(design["first"]) - h4(design["second"])
        projection = np.array([[1,0,0,0],[0,1,0,0],[1/delta,-1/delta,0,0],
                               [0,0,1,0],[0,0,0,1],[0,0,1/delta,-1/delta]])
        xp, yp = x @ projection.T, y @ projection.T
        joint = np.concatenate([xp.mean(axis=1), yp.mean(axis=1)], axis=1)
        cov = np.cov(joint, rowvar=False, ddof=1) / len(joint)
        residual = (xp - yp).reshape(-1, 6)
        noise = residual.T @ residual / len(residual)
        baseline_var = np.cov(xp.reshape(-1, 6), rowvar=False, ddof=1)
        counts = {s: statuses.count(s) for s in sorted(set(statuses))}
        fields = ["first_canonical", "second_canonical", "H4_normalized_canonical_difference",
                  "first_integrated", "second_integrated", "H4_normalized_integrated_difference"]
        result = {"source_sha256": source["source_sha256"], "delta_cos4": delta, "fields": fields,
                  "risk_frequencies": np.mean(risks, axis=0).tolist(), "pair_status_counts": counts,
                  "baseline_mean": xp.mean(axis=(0,1)).tolist(), "conditional_hybrid_mean": yp.mean(axis=(0,1)).tolist(),
                  "joint_20_batch_means_X_then_Y": joint.tolist(), "joint_covariance_of_means": cov.tolist(),
                  "paired_removed_suffix_noise_second_moment": noise.tolist(), "baseline_individual_covariance": baseline_var.tolist(),
                  "empirical_removed_noise_fraction": (np.diag(noise) / np.diag(baseline_var)).tolist(),
                  "claim_boundary": "All 20000 counters in denominator; only rank-one contribution targeted, no rank-zero/rank-two completion."}
        output["sizes"][str(n)] = result
        report += [f"## N{n}", "", f"Risk frequencies: {result['risk_frequencies']}; paired policy: {counts}.", "",
                   "| Readout | Original mean | Hybrid mean | Hybrid batch SE | Removed suffix-noise fraction |", "|---|---:|---:|---:|---:|"]
        for i, label in enumerate(fields):
            report.append(f"| {label} | {result['baseline_mean'][i]:.10g} | {result['conditional_hybrid_mean'][i]:.10g} | {np.sqrt(cov[6+i,6+i]):.6g} | {result['empirical_removed_noise_fraction'][i]:.6g} |")
        report.append("")
    (OUT / "score.json").write_text(json.dumps(output, indent=2) + "\n")
    (OUT / "REPORT.md").write_text("\n".join(report))
    print("\n".join(report), flush=True)


def main():
    contract = json.loads(CONTRACT.read_text())
    sources = read_sources(contract)
    (OUT / "batches").mkdir(parents=True, exist_ok=True)
    tasks = []
    for batch in range(20):
        for n, source in sources.items():
            first = source["metadata"]["replica_counter_first"] + batch * 1000
            rows = [r for (counter, _), r in source["rows"].items() if first <= counter < first+1000]
            tasks.append((n, batch, source["metadata"], rows, contract))
    with ProcessPoolExecutor(max_workers=contract["workers"]) as pool:
        for done in as_completed([pool.submit(process_batch, task) for task in tasks]):
            done.result()
    score(contract, sources)


if __name__ == "__main__":
    main()

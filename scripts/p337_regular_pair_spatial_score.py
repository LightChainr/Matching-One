#!/usr/bin/env python3
"""Fixed spatial Q-activation readout; no new sampling, fit, or resampling."""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np
from scipy.stats import t

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis/p337_regular_pair_spatial_contract.json"
FREEZE = "3210aeb3"
LABELS = ["total"] + [f"shared{k}" for k in range(5)]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def commit(ref):
    return subprocess.check_output(["git", "rev-parse", ref+"^{commit}"], cwd=ROOT, text=True).strip()


def read_run(directory, spec, contract):
    L = spec["L"]
    path = directory / f"L{L}.csv"
    metadata_path = Path(str(path)+".metadata.json")
    metadata = json.loads(metadata_path.read_text())
    expected = {"L": L, "N": L*L, "r": L//4, "seed": spec["seed"],
                "batches": 200, "samples_per_batch": 1000, "samples": 200000,
                "pairs_per_configuration": 32, "bernoulli_threshold_2pow64": 10934234699625173385}
    if metadata.get("status") != "completed" or any(metadata.get(k) != v for k, v in expected.items()):
        raise ValueError("run metadata differs from the fixed sampling contract")
    with path.open(newline="") as handle:
        rows = [{k: int(v) for k, v in row.items()} for row in csv.DictReader(handle)]
    if len(rows) != 200:
        raise ValueError("only a complete fixed 200-batch block may be scored")
    means = []
    for b, row in enumerate(rows):
        if (row["L"], row["batch"], row["samples"], row["pairs"]) != (L, b, 1000, 32000):
            raise ValueError("sample/pair normalization differs")
        by_shared = [row[f"sum_g16_shared{k}"] for k in range(5)]
        if sum(by_shared) != row["sum_g16"] or by_shared[:2] != [0, 0]:
            raise ValueError("exact kernel channel decomposition differs")
        if sum(row[f"pairs_shared{k}"] for k in range(5)) != row["eligible_pairs"]:
            raise ValueError("eligible-pair decomposition differs")
        means.append(np.asarray([row["sum_g16"], *by_shared], dtype=float)/(16*row["pairs"]))
    block = np.asarray(means)
    covariance = np.cov(block, rowvar=False, ddof=1)/len(block)
    mean = block.mean(axis=0)
    se = np.sqrt(np.maximum(np.diag(covariance), 0))
    critical = float(t.ppf(.995, 199))
    ci = np.column_stack((mean-critical*se, mean+critical*se))
    nonzero_pairs = sum(row["nonzero_pairs"] for row in rows)
    signed_sum = sum(row["sum_g16"] for row in rows)
    # An accounting bound on observed signed integers, not a confidence bound.
    cancellation_floor = max(0., 1-abs(signed_sum)/nonzero_pairs) if nonzero_pairs else None
    return {"L": L, "r": L//4, "mean": mean.tolist(), "mcse": se.tolist(),
            "labels": LABELS, "covariance_of_mean": covariance.tolist(), "ci99": ci.tolist(),
            "observed_zero_variance": (se == 0).tolist(),
            "total_eligible_pairs": sum(row["eligible_pairs"] for row in rows),
            "total_nonzero_pairs": nonzero_pairs, "total_signed_sum_g16": signed_sum,
            "observed_absolute_contribution_cancellation_lower_bound": cancellation_floor,
            "eligible_pairs_by_shared": [sum(row[f"pairs_shared{k}"] for row in rows) for k in range(5)],
            "samples": 200000, "pairs_per_configuration": 32,
            "seed": spec["seed"], "source_file": path.name, "source_sha256": sha(path),
            "metadata_file": metadata_path.name, "metadata_sha256": sha(metadata_path),
            "elapsed_production_seconds": metadata["elapsed_seconds"]}


def fieller(near, far, critical):
    x, y = near["mean"][0], far["mean"][0]
    vx, vy = near["mcse"][0]**2, far["mcse"][0]**2
    a, b, c = x*x-critical**2*vx, -2*x*y, y*y-critical**2*vy
    discriminant = b*b-4*a*c
    result = {"point_ratio": y/x if x != 0 else None,
              "status": "unbounded_or_unresolved", "confidence": .99,
              "covariance_between_sizes": 0, "bounded_interval": None}
    if a > 0 and discriminant >= 0 and vx > 0 and vy > 0:
        result["status"] = "bounded"
        result["bounded_interval"] = [(-b-np.sqrt(discriminant))/(2*a), (-b+np.sqrt(discriminant))/(2*a)]
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--kernel-commit", required=True)
    parser.add_argument("--producer-commit", required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    contract = json.loads(CONTRACT.read_text())
    if [run["L"] for run in contract["runs"]] != [32, 64]:
        raise ValueError("the fixed spatial geometry pair differs")
    outputs = [read_run(args.data_dir.resolve(), run, contract) for run in contract["runs"]]
    near, far = outputs
    mu, se = far["mean"][0], far["mcse"][0]
    pvalue = float(2*t.sf(abs(mu/se), 199)) if se > 0 else None
    decision = "contact_only_zero_rejected_at_L64_r16" if pvalue is not None and pvalue < .01 else "unresolved_stop_at_fixed_budget"
    ratio = fieller(near, far, float(t.ppf(.995, 199)))
    result = {"schema": "p337.regular-pair-spatial.score.v1", "status": "completed_fixed_two_block_readout",
              "decision": decision, "primary_two_sided_p": pvalue, "primary_alpha": .01,
              "contract": contract, "contract_sha256": sha(CONTRACT), "freeze_commit": commit(FREEZE),
              "kernel_commit": commit(args.kernel_commit), "producer_commit": commit(args.producer_commit),
              "code_commit": commit("HEAD"), "created_utc": datetime.now(timezone.utc).isoformat(),
              "runs": outputs, "ratio_C64_over_C32": ratio,
              "uncertainty": "Monte Carlo batch Student-t/Fieller approximations, not exact arithmetic or continuum error bounds",
              "inference_unit": "iid occupation configuration; 32 correlated translations/directions averaged within configuration",
              "new_samples_total": 400000, "independent_size_blocks": True, "resampling_performed": False,
              "field_assignment": None, "exponent_fit_performed": False,
              "elapsed_score_seconds": time.perf_counter()-started}
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    (out/"score.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    lines = ["# Canonical regular-pair spatial Q activation", "", f"Decision: **{decision}**.", "",
             "| L | r | C | Monte Carlo SE | 99% interval |", "|---:|---:|---:|---:|---:|"]
    for row in outputs:
        lo, hi = row["ci99"][0]
        lines.append(f"| {row['L']} | {row['r']} | {row['mean'][0]:.10g} | {row['mcse'][0]:.5g} | [{lo:.10g}, {hi:.10g}] |")
    lines += ["", f"Primary two-sided p={pvalue}. Fixed ratio C64/C32: {ratio}.", "",
              "The kernel is the Q derivative of the actual connected two-insertion colour contraction. "
              "It is not a covariance of separately closed one-site marks. The full shared-component mean/covariance "
              "and input provenance are in score.json. No exponent is fitted and no continuum field is identified. "
              "This fixed-budget result receives no top-up or new completion coefficient."]
    (out/"REPORT.md").write_text("\n".join(lines)+"\n")
    print(json.dumps({"decision": decision, "primary_p": pvalue,
                      "C32": near["mean"][0], "SE32": near["mcse"][0],
                      "C64": mu, "SE64": se, "ratio": ratio,
                      "elapsed_seconds": result["elapsed_score_seconds"]}, allow_nan=False))


if __name__ == "__main__":
    main()

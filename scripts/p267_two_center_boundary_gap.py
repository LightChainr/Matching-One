#!/usr/bin/env python3
"""One affine-invariant rank-two residual from stored three-center LOO vectors."""
from hashlib import sha256
import json
from pathlib import Path
import subprocess

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT/"experiments/p267_two_center_boundary_gap_20260831.json"
OUT = ROOT/"results/p267-two-center-boundary-gap"
LABELS = ["rank_two_gap_Delta2", "early_weight", "early_leverage_Delta2_over_weight"]


def readout(vector):
    centers = np.asarray(vector[1:4], dtype=float)
    weights = np.asarray(vector[4:7], dtype=float)
    weights = weights/weights.sum()
    if np.min(weights) <= 0 or np.min(np.diff(centers)) <= 0:
        raise ValueError("This readout requires the stored positive ordered three-center branch")
    centered = centers-weights@centers
    v = weights@centered**2
    m3, m4 = weights@centered**3, weights@centered**4
    gap_product = np.prod([centers[j]-centers[i] for i in range(3) for j in range(i+1, 3)])
    delta = np.prod(weights)*gap_product**2/v**3
    polynomial = centered**2-(m3/v)*centered-v
    certificate = {
        "center_variance": float(v),
        "minimizing_centered_monic_quadratic": [-float(v), -float(m3/v), 1.],
        "Delta2_from_quadratic_energy": float(weights@polynomial**2/v**2),
        "Delta2_from_Hankel": float(np.linalg.det([[1, 0, v], [0, v, m3], [v, m3, m4]])/v**3),
        "relative_middle_gap": float((centers[1]-centers[0])/(centers[2]-centers[0])),
        "meaning": "Equivalent representations of the one fixed residual, not a distance scan"}
    return np.array([delta, weights[0], delta/weights[0]]), certificate


def covariance(loo):
    centered = loo-loo.mean(axis=0)
    return (len(loo)-1)/len(loo)*centered.T@centered


def main():
    protocol = json.loads(PROTOCOL.read_text())
    rows, hashes = {}, {}
    for path in protocol["source_paths"]:
        blob = subprocess.check_output(["git", "show", protocol["source_commit"]+":"+path], cwd=ROOT)
        hashes[path] = sha256(blob).hexdigest()
        rows.update(json.loads(blob)["sources"])
    sources = {}
    for n, row in rows.items():
        point, certificate = readout(row["estimate"])
        loo = np.array([readout(v)[0] for v in row["leave_one_common_batch_out_vectors"]])
        cov = covariance(loo)
        sources[n] = {"raw_data_commit": row["source_commit"], "common_batches": len(loo),
                      "estimate": point.tolist(), "se": np.sqrt(np.diag(cov)).tolist(),
                      "full_covariance": cov.tolist(), "stored_parameter_LOO_readouts": loo.tolist(),
                      "jackknife_bias_estimate": ((len(loo)-1)*(loo.mean(axis=0)-point)).tolist(),
                      "certificate": certificate,
                      "boundary_calibration": "Regular-branch LOO uncertainty only; zero gap is a singular rank boundary, so no normal-ratio boundary p-value is assigned"}
    changes = {}
    for a, b in (("100", "400"), ("400", "900"), ("100", "900")):
        cov = np.array(sources[a]["full_covariance"])+sources[b]["full_covariance"]
        changes[f"N{b}_minus_N{a}"] = {
            "estimate": (np.array(sources[b]["estimate"])-sources[a]["estimate"]).tolist(),
            "se": np.sqrt(np.diag(cov)).tolist(), "full_covariance": cov.tolist()}
    result = {"schema": protocol["schema"], "definition_freeze": "fa562f47", "protocol": protocol,
              "source_sha256": hashes, "labels": LABELS, "sources": sources, "descriptive_changes": changes,
              "consecutive_increment_cross_covariance": (-np.array(sources["400"]["full_covariance"])).tolist(),
              "scope": "Same three source blocks reused; no new MC, old Gaussian-boundary/center refits, or p-value combination. One gap plus its exact factorization into early weight and leverage. Affine invariance does not mean invariance under arbitrary thermal warps."}
    OUT.mkdir(exist_ok=True)
    (OUT/"score.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    lines = ["# A small early weight does not by itself approach the two-center boundary", "",
             "| N | Fixed rank-two gap Delta2 | LOO SE | early weight | leverage Delta2/w0 | leverage SE |",
             "|---|---:|---:|---:|---:|---:|"]
    for n, row in sources.items():
        point, se = row["estimate"], row["se"]
        lines.append(f"| {n} | {point[0]:.9g} | {se[0]:.7g} | {point[1]:.9g} | {point[2]:.9g} | {se[2]:.7g} |")
    lines += ["", "Changes, not independent of each other:"]
    for name, row in changes.items():
        lines.append(f"- {name}: Delta2 change {row['estimate'][0]:.9g} +/- {row['se'][0]:.7g}.")
    lines += ["", result["scope"], "", "No Gaussian boundary p-value is assigned at the singular rank-two null.", ""]
    (OUT/"REPORT.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()

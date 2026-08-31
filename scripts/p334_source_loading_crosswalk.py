#!/usr/bin/env python3
"""Thin same-batch C/L x direct/collective/fallback crosswalk; no raw reprocessing."""
from hashlib import sha256
import json
from pathlib import Path
import subprocess

import numpy as np

from p334_r1_prevalence_clock_loading import decompose, LABELS

ROOT = Path(__file__).resolve().parents[1]
BASE_COMMIT, BASE_PATH = "3d760b86", "results/p334-r1-prevalence-clock-loading/score.json"
SOURCE_COMMIT = "32270fa2f8c5dfb19bf534b364fde26e2ac117f6"
SOURCE_PATH = "results/p334-direct-collective-population-loading/score.json"
OUT = ROOT/"results/p334-source-loading-crosswalk"
SELECT = [6, 7, 8, 13, 14, 15]


def covariance(loo):
    centered = loo-loo.mean(axis=0)
    return (len(loo)-1)/len(loo)*centered.T@centered


def main():
    hashes = {}
    def read(commit, path):
        blob = subprocess.check_output(["git", "show", commit+":"+path], cwd=ROOT)
        hashes[commit+":"+path] = sha256(blob).hexdigest()
        return json.loads(blob)
    base, source = read(BASE_COMMIT, BASE_PATH), read(SOURCE_COMMIT, SOURCE_PATH)
    if base["source_commit"] != source["source_commit"]:
        raise ValueError("C/L and microscopic-source blocks have different original archives")
    components = source["components"]
    if source["orientations"] != ["first", "second"] or source["readouts"] != ["canonical", "integrated"] or len(components) != 3:
        raise ValueError("Unexpected fixed source layout")
    selected_labels = [f"{component}.{LABELS[k]}" for component in components for k in SELECT]
    sizes = {}
    for n, row in base["sizes"].items():
        src = source["sizes"][n]
        batch = np.array(row["joint_batch_means_risk_then_Y"])
        b = len(batch)
        if b != 20 or src["batch_ids"] != list(range(20)) or src["delta_cos4"] != row["delta_cos4"]:
            raise ValueError("Original size/batch/normalization alignment failed")
        raw = np.array(src["joint_20_batch_means_orientation_readout_source"]).reshape(20, 2, 2, 3)
        source_y = np.stack([raw[:, 0, 0], raw[:, 1, 0], raw[:, 0, 1], raw[:, 1, 1]], axis=1)
        if not np.allclose(source_y.sum(axis=2), batch[:, 2:], atol=1e-12, rtol=0):
            raise ValueError("Microscopic sources do not recover the same batch hybrid Y")
        risk_mean = batch[:, :2].mean(axis=0)
        risk_loo = (b*risk_mean-batch[:, :2])/(b-1)
        estimates, replicas, details = [], [], {}
        for k, component in enumerate(components):
            y_mean = source_y[:, :, k].mean(axis=0)
            y_loo = (b*y_mean-source_y[:, :, k])/(b-1)
            point = decompose(np.r_[risk_mean, y_mean], row["delta_cos4"])
            loo = np.array([decompose(np.r_[r, y], row["delta_cos4"]) for r, y in zip(risk_loo, y_loo)])
            details[component] = {"labels": LABELS, "estimate": point.tolist(),
                                  "se": np.sqrt(np.diag(covariance(loo))).tolist(),
                                  "leave_one_common_batch_out_vectors": loo.tolist(),
                                  "denominator": "The original orientation-specific R1 risk r_i, never a source-event rate"}
            estimates.append(point[SELECT]); replicas.append(loo[:, SELECT])
        source_estimate = np.concatenate(estimates)
        source_loo = np.concatenate(replicas, axis=1)
        old_joint = row["joint_mean_and_variance_readout_covariance"]
        joint_loo = np.column_stack((old_joint["leave_one_common_batch_out_vectors"], source_loo))
        joint_cov = covariance(joint_loo)
        target_loo = np.array(row["leave_one_common_batch_out_vectors"])[:, SELECT]
        closure = sum(replicas)-target_loo
        sizes[n] = {"batch_ids": list(range(b)), "source_specific_rows": details,
                    "source_CLD_labels": selected_labels, "source_CLD_estimate": source_estimate.tolist(),
                    "source_CLD_covariance": covariance(source_loo).tolist(),
                    "source_CLD_se": np.sqrt(np.diag(covariance(source_loo))).tolist(),
                    "joint_labels": old_joint["labels"]+selected_labels,
                    "joint_covariance": joint_cov.tolist(), "joint_LOO_vectors": joint_loo.tolist(),
                    "joint_rank_at_most": b-1,
                    "max_batch_source_Y_sum_error": float(np.max(np.abs(source_y.sum(axis=2)-batch[:, 2:]))),
                    "max_LOO_source_CLD_sum_error": float(np.max(np.abs(closure))),
                    "point_source_CLD_sum_error": (sum(estimates)-np.array(row["estimate"])[SELECT]).tolist()}
    result = {"schema": "matching-one/p334-r1-loading-microscopic-source-crosswalk/v1",
              "original_source_commit": base["source_commit"], "CL_commit": BASE_COMMIT,
              "microscopic_source_commit": SOURCE_COMMIT, "input_sha256": hashes,
              "components": components, "sizes": sizes,
              "identities": ["m_i,s=mean(Y_i,s)/r_i with the same original R1 r_i for every source",
                             "sum_s C_s=C; sum_s L_s=L; sum_s D_s=D"],
              "boundary": "Same 20 original batches, not independent source tests. Full joint covariance is saved, but never inverted and no high-dimensional omnibus is computed. Direct means a final site in the original checkpoint H2 set; later-created triggers belong to collective. All whole-pair fallbacks remain unclassified. Point signs or cancellations are not asserted as statistically established. This is gated R1 only, not full A_top.",
              "new_MC": 0, "new_DP": 0, "raw_batch_reprocessing": 0}
    OUT.mkdir(exist_ok=True)
    (OUT/"score.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    lines = ["# Shared-covariance microscopic sources of prevalence and clock loading", "",
             "| N / readout / source | C | SE | L | SE | D=C+L | SE |", "|---|---:|---:|---:|---:|---:|---:|"]
    for n, row in sizes.items():
        for component, detail in row["source_specific_rows"].items():
            v, se = detail["estimate"], detail["se"]
            for endpoint, start in (("canonical", 6), ("integrated", 13)):
                lines.append(f"| {n} / {endpoint} / {component} | {v[start]:.9g} | {se[start]:.6g} | {v[start+1]:.9g} | {se[start+1]:.6g} | {v[start+2]:.9g} | {se[start+2]:.6g} |")
    lines += ["", result["boundary"], ""]
    (OUT/"REPORT.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()

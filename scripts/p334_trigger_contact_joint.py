#!/usr/bin/env python3
"""Thin same-batch join of exact trigger partition and fixed-degree contact response."""
import gzip
from hashlib import sha256
import json
from pathlib import Path
import subprocess

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PARTITION_COMMIT = "22952d75"
PARTITION_DIR = "results/p334-exact-birth-partition-response"
CONTACT_COMMIT = "b9f79bfb"
CONTACT_PATH = "results/p334-safe-contact-response/score.json"
OUT = ROOT/"results/p334-trigger-contact-joint"


def main():
    hashes = {}
    def read(commit, path):
        blob = subprocess.check_output(["git", "show", commit+":"+path], cwd=ROOT)
        hashes[commit+":"+path] = sha256(blob).hexdigest()
        return blob
    partition = json.loads(read(PARTITION_COMMIT, PARTITION_DIR+"/score.json"))
    contact = json.loads(read(CONTACT_COMMIT, CONTACT_PATH))
    if partition["source_commit"] != contact["source_commit"]:
        raise ValueError("Different underlying random fork blocks")
    OUT.mkdir(parents=True, exist_ok=True)
    sizes = {}
    for n in ("325", "425"):
        p, c = partition["sizes"][n], contact["sizes"][n]
        if p["batch_ids"] != c["batch_ids"] or p["batch_ids"] != list(range(20)):
            raise ValueError("Original batch rows are not aligned")
        blob = read(PARTITION_COMMIT, PARTITION_DIR+"/"+p["complete_covariance_factor_file"])
        saved = json.loads(gzip.decompress(blob))
        factor = np.column_stack((np.array(saved["factor"]), np.array(c["raw_covariance_factor"]), np.array(c["factor"])))
        labels = saved["labels"]+["contact.raw."+x for x in c["raw_labels"]]+["contact."+x for x in c["labels"]]
        point = {"birth_partition."+k: v for k, v in zip(p["gamma_labels"], p["gamma_estimate"]) if v is not None}
        point.update({"four_change_type."+k: v for k, v in zip(p["four_type_labels"], p["four_type_estimate"])})
        point.update({"contact.raw."+k: v for k, v in zip(c["raw_labels"], c["raw_estimate"])})
        point.update({"contact."+k: v for k, v in zip(c["labels"], c["estimate"])})
        # Derived shares use the same original leave-one-batch-out ratios.
        four_batches = np.array(p["four_type_20_batch_means"])
        four_mean = four_batches.mean(axis=0)
        four_loo = (20*four_mean-four_batches)/19
        old_gamma_labels = ["birth_partition."+k for k, v in zip(p["gamma_labels"], p["gamma_estimate"]) if v is not None]
        old_gamma_loo = np.array(saved["gamma_LOO_identified"])
        extras, extra_point, extra_loo = [], [], []
        for group in ("all", "01+10"):
            for ep in ("p_ref", "p_integral"):
                ti = p["four_type_labels"].index(f"{group}.total.{ep}.Gamma")
                bi = p["four_type_labels"].index(f"{group}.between.{ep}.Gamma")
                binary = "birth_partition."+f"{group}.between.{ep}.Gamma_pair"
                bidx = old_gamma_labels.index(binary)
                for mode, value, reps in (("binary", point[binary], old_gamma_loo[:, bidx]),
                                           ("four_type", four_mean[bi], four_loo[:, bi])):
                    name = f"{mode}.{group}.{ep}.between_Gamma_share"
                    extras.append(name); extra_point.append(value/four_mean[ti]); extra_loo.append(reps/four_loo[:, ti])
                name = f"refinement.{group}.{ep}.between_Gamma_increment"
                extras.append(name); extra_point.append(four_mean[bi]-point[binary]); extra_loo.append(four_loo[:, bi]-old_gamma_loo[:, bidx])
        extra_loo = np.column_stack(extra_loo)
        extra_factor = np.sqrt(19/20)*(extra_loo-extra_loo.mean(axis=0))
        factor = np.column_stack((factor, extra_factor)); labels += extras
        point.update(zip(extras, extra_point))
        focused = []
        for group in ("all", "01+10"):
            for ep in ("p_ref", "p_integral"):
                for part in ("within_total", "between"):
                    focused += [f"birth_partition.{group}.{part}.{ep}.Gamma_pair", f"four_change_type.{group}.{part}.{ep}.Gamma"]
        focused += extras
        contact_group = "R0_safe_equal_contact_degree"
        focused += [f"contact.{contact_group}.pooled_slope[contractible_cycles,{k}]" for k in ("K1", "K2", "C", "W")]
        focused += [f"contact.raw.{contact_group}.GX[contractible_cycles,{ep}.{k}]" for ep in ("p_ref", "p_integral") for k in ("F1", "F2")]
        ix = [labels.index(k) for k in focused]
        cov = factor[:, ix].T@factor[:, ix]
        filename = f"N{n}.complete_common_factor.json.gz"
        packed = {"labels": labels, "factor": factor.tolist(), "batch_ids": list(range(20)),
            "convention": "factor.T@factor; original common twenty LOO rows; no inverse", "rank_at_most": 19,
            "contact_raw_labels": c["raw_labels"], "contact_raw_joint_20_batch_means": c["raw_joint_20_batch_means"],
            "contact_derived_labels": c["labels"], "contact_derived_LOO": c["LOO"],
            "additional_labels": extras, "additional_LOO": extra_loo.tolist()}
        blob = gzip.compress((json.dumps(packed, separators=(",", ":"), allow_nan=False)+"\n").encode(), mtime=0)
        (OUT/filename).write_bytes(blob)
        sizes[n] = {"batch_ids": list(range(20)), "labels": focused, "estimate": [point[k] for k in focused],
            "se": np.sqrt(np.diag(cov)).tolist(), "focused_covariance": cov.tolist(),
            "complete_covariance_factor_file": filename, "complete_covariance_factor_sha256": sha256(blob).hexdigest(),
            "complete_coordinate_count": len(labels)}
    result = {"schema": "matching-one/p334-trigger-contact-joint/v1", "source_commit": partition["source_commit"],
        "partition_commit": PARTITION_COMMIT, "contact_commit": CONTACT_COMMIT, "source_sha256": hashes, "sizes": sizes,
        "new_MC": 0, "new_DP": 0, "new_fork_or_contact_raw_reads": 0,
        "boundary": "Trigger partitions use the paired observable and both-orientation labels. Contact response pools own-orientation R0-safe labels at equal contact degree, not the joint-safe class. They share one random block and one common covariance, not independent evidence. No inverse, fitting, clipping or new scan."}
    (OUT/"score.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    lines = ["# Trigger-type mode and fixed-degree contact response: one shared block", ""]
    for n, r in sizes.items():
        lines += [f"## N{n}", "", "| Coordinate | Estimate | Shared-batch SE |", "|---|---:|---:|"]
        for label, mean, se in zip(r["labels"], r["estimate"], r["se"]):
            if label.startswith(("binary.", "four_type.", "refinement.")) or ".pooled_slope" in label:
                lines.append(f"| {label} | {mean:.10g} | {se:.7g} |")
        lines.append("")
    lines += [result["boundary"], ""]
    (OUT/"REPORT.md").write_text("\n".join(lines)); print("\n".join(lines))


if __name__ == "__main__": main()

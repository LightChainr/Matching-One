#!/usr/bin/env python3
"""Reuse saved quartet matrices: one eta magnitude plus the fixed gate handoff."""
import gzip
from hashlib import sha256
import json
from pathlib import Path
import subprocess

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CORE_COMMIT = "24872eef"
CORE_DIR = "results/p334-next-label-doob-quartets"
GATE_COMMIT = "c6ee37a8"
GATE_PATH = "results/p334-next-label-gate-coupling/batch_vectors.json"
OUT = ROOT/"results/p334-next-label-mechanism-joint"


def main():
    hashes = {}
    def read(commit, path):
        blob = subprocess.check_output(["git", "show", commit+":"+path], cwd=ROOT)
        hashes[commit+":"+path] = sha256(blob).hexdigest()
        return blob
    core = json.loads(read(CORE_COMMIT, CORE_DIR+"/score.json"))
    gate = json.loads(read(GATE_COMMIT, GATE_PATH))
    if core["source_commit"] != gate["source_commit"]:
        raise ValueError("The gate and Doob readouts use different new-tail blocks")
    OUT.mkdir(parents=True, exist_ok=True)
    sizes = {}
    for n in ("325", "425"):
        c, g = core["sizes"][n], gate["sizes"][n]
        if c["batch_ids"] != g["batch_ids"] or g["batch_denominators"] != [1000]*20:
            raise ValueError("Original batch alignment differs")
        if c["delta_cos4"] != g["delta_cos4"]:
            raise ValueError("Orientation normalization differs")
        for path, digest in gate["source_sha256"].items():
            if f"/N{n}/" in path and core["source_sha256"][core["source_commit"]+":"+path] != digest:
                raise ValueError("Raw source hash mismatch")
        factor_blob = read(CORE_COMMIT, CORE_DIR+"/"+c["covariance_factor_file"])
        if sha256(factor_blob).hexdigest() != c["covariance_factor_sha256"]:
            raise ValueError("Saved covariance factor differs")
        saved = json.loads(gzip.decompress(factor_blob))
        raw = np.array(saved["raw_joint_20_batch_means"])
        raw_labels = saved["raw_batch_labels"]
        eta, eta_loo = [], []
        for ep in ("p_ref", "p_integral"):
            cols = []
            for first, second in (("F1", "F1"), ("F2", "F2"), ("F1", "F2")):
                ix = [raw_labels.index(f"cell.{a}{b}.Dnext[{ep}.{first},{ep}.{second}]") for a in range(3) for b in range(3)]
                cols.append(raw[:, ix].sum(axis=1))
            batches = np.column_stack(cols)
            mean = batches.mean(axis=0)
            loo = (20*mean-batches)/19
            eta.append(2*mean[2]/(mean[0]+mean[1]))
            eta_loo.append(2*loo[:, 2]/(loo[:, 0]+loo[:, 1]))
        eta_loo = np.column_stack(eta_loo)
        gate_loo = np.array(g["joint_LOO_vectors"])
        added_loo = np.column_stack((eta_loo, gate_loo))
        added_factor = np.sqrt(19/20)*(added_loo-added_loo.mean(axis=0))
        factor = np.column_stack((np.array(saved["factor"]), added_factor))
        eta_labels = ["eta_next.p_ref", "eta_next.p_integral"]
        labels = saved["labels"]+eta_labels+["gate."+name for name in g["labels"]]
        gamma_labels = ["all.next_first_completion_Gamma.p_ref", "all.next_first_completion_Gamma.p_integral"]
        focused = gamma_labels+eta_labels+["gate."+name for name in g["labels"]]
        ix = [labels.index(label) for label in focused]
        values = [c["primary_estimate"][c["primary_labels"].index(name)] for name in gamma_labels]+eta+g["estimate"]
        covariance = factor[:, ix].T@factor[:, ix]
        filename = f"N{n}.complete_joint_factor.json.gz"
        packed = {"labels": labels, "factor": factor.tolist(), "batch_ids": list(range(20)),
                  "convention": "Full same-block covariance = factor.T @ factor; append-only aligned original LOO; no inverse",
                  "rank_at_most": 19}
        blob = gzip.compress((json.dumps(packed, separators=(",", ":"), allow_nan=False)+"\n").encode(), mtime=0)
        (OUT/filename).write_bytes(blob)
        sizes[n] = {"batch_ids": list(range(20)), "labels": focused, "estimate": values,
                    "se": np.sqrt(np.diag(covariance)).tolist(), "focused_covariance": covariance.tolist(),
                    "eta_LOO": eta_loo.tolist(), "gate_joint_20_batch_means": g["joint_20_batch_means"],
                    "complete_covariance_factor_file": filename, "complete_covariance_factor_sha256": sha256(blob).hexdigest()}
    result = {"schema": "matching-one/p334-next-label-eta-gate-joint/v1", "core_commit": CORE_COMMIT,
              "gate_commit": GATE_COMMIT, "source_commit": core["source_commit"], "source_sha256": hashes,
              "definition": "eta_next=2 B12/(B11+B22); B_AA=(B11+B22)(1+eta), B_EE=(B11+B22)(1-eta)",
              "sizes": sizes, "new_MC": 0, "new_DP": 0, "new_raw_csv_reads": 0,
              "boundary": "Eta is the one named magnitude of the already observed positive Gamma, not another independently tested direction. Gate statistics cover only 01/10 and are distinct from the all-prefix continuous first/completion response. All same-source covariance is appended without clipping, refitting or inversion."}
    (OUT/"score.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    lines = ["# Named next-label response magnitude and shared gate covariance", ""]
    for n, r in sizes.items():
        lines += [f"## N{n}", "", "| Named coordinate | Mean | Shared-batch SE |", "|---|---:|---:|"]
        for label, point, se in zip(r["labels"][:4], r["estimate"][:4], r["se"][:4]):
            lines.append(f"| {label} | {point:.10g} | {se:.6g} |")
        lines.append("")
    lines += [result["boundary"], ""]
    (OUT/"REPORT.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""One shared-batch join: full A, marked full-A sources and raw first thermal moment."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
INPUTS = {
    "global": ("3edc785a", "results/p334-full-global-conditional-clock/score.json"),
    "marked": ("2dd865f0b26a4d5d43f52b300293016e6ffd19b8", "results/p334-marked-global-topology-loading/score.json"),
    "thermal": ("e64febe4ff10ca9cfb2f094c1b8ee8f733177fe1", "results/p334-global-first-thermal-moment/batch_vectors.json"),
}
OUT = ROOT/"results/p334-global-source-thermal-joint"


def shared_summary(batches):
    mean = batches.mean(axis=0)
    loo = (20*mean-batches)/19
    centered = loo-loo.mean(axis=0)
    covariance = 19/20*centered.T@centered
    return mean, loo, covariance


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUT)
    args = parser.parse_args()
    data, provenance = {}, {}
    for name, (commit, path) in INPUTS.items():
        commit = subprocess.check_output(["git", "rev-parse", commit+"^{commit}"], cwd=ROOT, text=True).strip()
        blob = subprocess.check_output(["git", "show", commit+":"+path], cwd=ROOT)
        data[name] = json.loads(blob)
        provenance[name] = {"commit": commit, "path": path, "sha256": sha256(blob).hexdigest()}
    full_birth = data["global"]["full_birth_commit"]
    if data["marked"]["birth_commit"] != full_birth or data["thermal"]["source_commit"] != full_birth:
        raise ValueError("Different full-birth dependency blocks")
    if data["marked"]["clock_commit"] != data["global"]["conditional_commit"]:
        raise ValueError("Different conditional dependency blocks")
    sizes = {}
    for n in ("325", "425"):
        g, s, t = (data[k]["sizes"][n] for k in ("global", "marked", "thermal"))
        if any(r["batch_ids"] != list(range(20)) for r in (g, s, t)):
            raise ValueError("Require matching original twenty batches")
        if g["samples_per_batch"] != 1000 or s["counter_count_per_batch"] != 1000 or t["batch_denominators"] != [1000]*20:
            raise ValueError("Different original counter denominators")
        delta = g["delta_cos4"]
        if s["delta_cos4"] != delta:
            raise ValueError("Inconsistent H4 normalization")
        if t["source_sha256"] != s["birth_sha256"]:
            raise ValueError("Different full-birth CSV hashes")
        gb = np.array(g["batch_means"])
        a_orientation = np.array(g["A_orientation_batch_means"])
        marked_orientation = np.array(s["joint_20_batch_means_orientation_readout_source"]).reshape(20, 2, 2, 3)
        positive = np.array(s["positive_F2_20_batch_means"]).reshape(20, 2, 2, 2)
        debt = np.array(s["first_birth_debt_20_batch_means"]).reshape(20, 2, 2, 2)
        tb = np.array(t["joint_20_batch_means"])
        safe_marked = np.array(s["safe_global_hybrid_A_20_batch_means"])
        joins = {
            "source_sum_to_marked_safe": float(np.max(np.abs(marked_orientation.sum(axis=3).reshape(20, 4)-safe_marked))),
            "marked_safe_to_global_safe": float(np.max(np.abs(safe_marked[:, [0, 2, 1, 3]]-a_orientation[:, 4:]))),
            "marked_raw_to_global_baseline": float(np.max(np.abs(np.array(s["raw_A_20_batch_means"])[:, [0, 2, 1, 3]]-a_orientation[:, :4]))),
            "positive_minus_debt_to_marked_DG": float(np.max(np.abs(positive-debt-marked_orientation[:, :, :, :2]))),
            "thermal_J0_to_global_baseline_integral": float(np.max(np.abs(tb[:, [0, 3]]-a_orientation[:, 2:4]))),
        }
        if max(joins.values()) > 1e-12:
            raise ValueError("A source join or exact add-back failed: "+str(joins))
        labels, columns = [], []
        def add(label, vector):
            labels.append(label)
            columns.append(vector)
        for variant in ("baseline", "safe"):
            for ep in ("p_ref", "p_integral"):
                key = f"{variant}.{ep}.A_H4"
                add(key, gb[:, g["labels"].index(key)])
        for k, name in enumerate(data["marked"]["labels"]):
            add("raw_marked."+name, marked_orientation.reshape(20, 12)[:, k])
        for name, array in (("positive_F2", positive), ("first_birth_debt", debt)):
            for k, label in enumerate(data["marked"]["positive_F2_and_first_birth_debt_labels"]):
                add("raw_"+name+"."+label, array.reshape(20, 8)[:, k])
        for j, ep in enumerate(("p_ref", "p_integral")):
            for k, source in enumerate(("original_H2_direct_A", "collective_A", "remainder_A")):
                add(f"marked.{ep}.{source}.H4", (marked_orientation[:, 0, j, k]-marked_orientation[:, 1, j, k])/delta)
        for name, array in (("positive_F2", positive), ("first_birth_debt", debt)):
            for j, ep in enumerate(("p_ref", "p_integral")):
                for k, source in enumerate(("original_H2_direct_A", "collective_A")):
                    add(f"{name}.{ep}.{source}.H4", (array[:, 0, j, k]-array[:, 1, j, k])/delta)
        for k, name in enumerate(data["thermal"]["labels"]):
            add("raw_thermal."+name, tb[:, k])
        for name, value in (("J0", tb[:, 0]-tb[:, 3]),
                            ("J1_center", tb[:, 1]-tb[:, 4]),
                            ("J1_width", tb[:, 2]-tb[:, 5]),
                            ("J1_total", tb[:, 1]+tb[:, 2]-tb[:, 4]-tb[:, 5])):
            add("raw_thermal."+name+".H4", value/delta)
        batches = np.column_stack(columns)
        mean, loo, covariance = shared_summary(batches)
        sizes[n] = {"delta_cos4": delta, "batch_ids": list(range(20)), "counter_denominator_per_batch": 1000,
                    "labels": labels, "mean": mean.tolist(), "se": np.sqrt(np.diag(covariance)).tolist(),
                    "joint_batch_vectors": batches.tolist(), "joint_LOO_vectors": loo.tolist(),
                    "joint_covariance": covariance.tolist(), "joint_covariance_rank_at_most": 19,
                    "join_max_abs_errors": joins,
                    "marked_source_orientation_labels": data["marked"]["labels"],
                    "marked_source_orientation_batch_vectors": marked_orientation.reshape(20, 12).tolist()}
    result = {"schema": "matching-one/p334-global-source-thermal-shared-batches/v1",
              "sources": provenance, "full_birth_dependency_commit": full_birth,
              "conditional_dependency_commit": data["global"]["conditional_commit"],
              "sizes": sizes, "new_MC": 0, "new_DP": 0, "new_raw_replays": 0,
              "boundary": "Same original counter block per N. Marked full-A direct/collective terms exist only on accepted R1 orientations; remainder retains every other full-A term. Thermal J0/J1 use baseline paths, not conditional replacement. No source sign claim, independent confirmation, covariance inversion or omnibus test."}
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output/"score.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    lines = ["# Full topology: microscopic sources and first thermal moment", ""]
    for n, r in sizes.items():
        lines += [f"## N{n}", "", "| Same-batch readout | H4 mean | SE |", "|---|---:|---:|"]
        for label, m, se in zip(r["labels"], r["mean"], r["se"]):
            if label.endswith("H4"):
                lines.append(f"| {label} | {m:.10g} | {se:.6g} |")
        lines += ["", "Maximum batch add-back/join error: "+str(max(r["join_max_abs_errors"].values())), ""]
    lines += [result["boundary"], ""]
    (args.output/"REPORT.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()

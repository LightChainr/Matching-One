#!/usr/bin/env python3
"""Join the final common-label tangent and recover two named cross responses."""
import gzip
from hashlib import sha256
import json
from pathlib import Path
import subprocess

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OLD_COMMIT = "56b383327c834236114be13f2a34f52688803e8c"
OLD_DIR = "results/p334-trigger-contact-joint"
SOURCE_COMMIT = "4db356e1b026853468f94d59d938895a2367ceb7"
SOURCE_PATH = "results/p334-common-label-euler-tangent/score.json"
RANK_COMMIT = "73608ba9d3eef34c6980cb5a049f726cfebdd72d"
RANK_PATH = "results/p334-common-label-response-rank/score.json"
OUT = ROOT/"results/p334-common-label-tangent-joint"


def cross_chart(labels, delta):
    """Columns are raw L_first -> Y_second and L_second -> Y_first."""
    names, columns = [], []
    for group, cells in (("00", ["00"]), ("01+10", ["01", "10"])):
        for ep in ("p_ref", "p_integral"):
            for birth in ("F1", "F2"):
                for direction, coefficients in (
                    ("first_loop->second", (1, 1, -delta/2, -delta/2)),
                    ("second_loop->first", (1, -1, delta/2, -delta/2)),
                ):
                    column = np.zeros(len(labels))
                    for cell in cells:
                        for channel, weight in zip(("plus->S", "minus->S", "plus->D", "minus->D"), coefficients):
                            column[labels.index(f"{cell}.{channel}.{ep}.{birth}")] += weight
                    names.append(f"{group}.{direction}.{ep}.{birth}"); columns.append(column)
    return names, np.column_stack(columns)


def main():
    hashes = {}
    def read(commit, path):
        blob = subprocess.check_output(["git", "show", commit+":"+path], cwd=ROOT)
        hashes[commit+":"+path] = sha256(blob).hexdigest()
        return blob
    old = json.loads(read(OLD_COMMIT, OLD_DIR+"/score.json"))
    source = json.loads(read(SOURCE_COMMIT, SOURCE_PATH))
    rank = json.loads(read(RANK_COMMIT, RANK_PATH))
    if old["source_commit"] != source["source_commit"]:
        raise ValueError("Different underlying random fork sources")
    if rank["source_commit"] != SOURCE_COMMIT:
        raise ValueError("Rank readout does not use this common-label response")
    OUT.mkdir(parents=True, exist_ok=True)
    sizes = {}
    for n in ("325", "425"):
        o, s, rk = old["sizes"][n], source["sizes"][n], rank["sizes"][n]
        if s["batch_ids"] != o["batch_ids"] or s["batch_ids"] != rk["batch_ids"] or s["batch_ids"] != list(range(20)):
            raise ValueError("The original twenty batches must align")
        saved = json.loads(gzip.decompress(read(OLD_COMMIT, OLD_DIR+"/"+o["complete_covariance_factor_file"])))
        raw = np.array(s["joint_20_batch_means"])
        cross_names, transform = cross_chart(s["labels"], s["delta_cos4"])
        cross_batches = raw@transform
        extra_batches = np.column_stack((raw, cross_batches))
        extra_mean = extra_batches.mean(axis=0)
        extra_loo = (20*extra_mean-extra_batches)/19
        extra_factor = np.sqrt(19/20)*(extra_loo-extra_loo.mean(axis=0))
        new_labels = ["common_label."+x for x in s["labels"]]+["cross_response."+x for x in cross_names]
        labels = saved["labels"]+new_labels+["response_rank."+x for x in rk["labels"]]
        factor = np.column_stack((np.array(saved["factor"]), extra_factor, np.array(rk["factor"])))
        point = dict(zip(new_labels, extra_mean))
        point.update(zip(o["labels"], o["estimate"]))
        point.update(zip(["response_rank."+x for x in rk["labels"]], rk["estimate"]))
        selected = ["common_label."+k for k in s["labels"] if k.startswith("all.")]
        selected += ["cross_response."+x for x in cross_names]
        selected += [f"contact.R0_safe_equal_contact_degree.pooled_slope[contractible_cycles,{k}]" for k in ("K1", "K2", "C", "W")]
        selected += [f"response_rank.all.{ep}.A.{k}" for ep in ("p_ref", "p_integral") for k in ("det_J", "oriented_unit_column_area")]
        selected += ["response_rank.all.input_G.correlation"]
        ix = [labels.index(k) for k in selected]
        cov = factor[:, ix].T@factor[:, ix]
        filename = f"N{n}.complete_common_factor.json.gz"
        packed = {"labels": labels, "factor": factor.tolist(), "batch_ids": list(range(20)),
            "convention": "factor.T@factor; original common twenty LOO rows; rank at most 19; no inverse",
            "new_labels": new_labels, "new_joint_20_batch_means": extra_batches.tolist(),
            "new_LOO": extra_loo.tolist(), "cross_response_source_labels": s["labels"],
            "cross_response_labels": cross_names, "cross_response_linear_map": transform.tolist(),
            "rank_derived_labels": rk["labels"], "rank_derived_LOO": rk["LOO"]}
        blob = gzip.compress((json.dumps(packed, separators=(",", ":"), allow_nan=False)+"\n").encode(), mtime=0)
        (OUT/filename).write_bytes(blob)
        ccov = extra_factor[:, len(s["labels"]):].T@extra_factor[:, len(s["labels"]):]
        sizes[n] = {"batch_ids": list(range(20)), "delta_cos4": s["delta_cos4"], "labels": selected,
            "estimate": [point[k] for k in selected], "se": np.sqrt(np.diag(cov)).tolist(), "focused_covariance": cov.tolist(),
            "cross_response_labels": cross_names, "cross_response_estimate": cross_batches.mean(axis=0).tolist(),
            "cross_response_se": np.sqrt(np.diag(ccov)).tolist(), "cross_response_covariance": ccov.tolist(),
            "cross_response_20_batch_means": cross_batches.tolist(),
            "structural_zero_cross_responses": [k for k in cross_names if k.startswith("01+10.") and k.endswith(".F1")],
            "complete_covariance_factor_file": filename, "complete_covariance_factor_sha256": sha256(blob).hexdigest(),
            "complete_coordinate_count": len(labels)}
    result = {"schema": "matching-one/p334-common-label-tangent-joint/v1", "source_commit": source["source_commit"],
        "previous_shared_commit": OLD_COMMIT, "common_label_readout_commit": SOURCE_COMMIT,
        "rank_readout_commit": RANK_COMMIT, "source_sha256": hashes,
        "sizes": sizes, "new_MC": 0, "new_DP": 0, "new_raw_fork_reads": 0,
        "cross_response_convention": "Cfs=RSplus+RSminus-delta/2*(RDplus+RDminus); Csf=RSplus-RSminus+delta/2*(RDplus-RDminus), raw single-orientation scale",
        "boundary": "Common-label S/D responses share one intervention and joint-degree policy. Old own-policy slopes have different domains and normalizations. All share one fork block and the original twenty-batch covariance. No inverse, fitting, clipping or omnibus. The 01+10 F1 cross response is an algebraic zero; stored tiny floating reconstruction residuals are not statistical evidence."}
    (OUT/"score.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    lines = ["# Common-label Euler-invisible 2x2 response and cross-orientation channels", ""]
    for n, r in sizes.items():
        lines += [f"## N{n}", "", "| Observer | plus -> S | plus -> D | minus -> S | minus -> D |", "|---|---:|---:|---:|---:|"]
        for observer in [f"{ep}.{k}" for ep in ("p_ref", "p_integral") for k in ("F1", "F2", "A", "E")]+["K1", "K2", "C", "W"]:
            values = []
            for channel in ("plus->S", "plus->D", "minus->S", "minus->D"):
                k = r["labels"].index(f"common_label.all.{channel}.{observer}")
                values.append(f"{r['estimate'][k]:.9g} +/- {r['se'][k]:.5g}")
            lines.append("| "+observer+" | "+" | ".join(values)+" |")
        lines += ["", "| Cross response | Estimate | Shared-batch SE |", "|---|---:|---:|"]
        for name, value, se in zip(r["cross_response_labels"], r["cross_response_estimate"], r["cross_response_se"]):
            if name in r["structural_zero_cross_responses"]:
                lines.append(f"| {name} | 0 (structural) | not a statistical null |")
            else:
                lines.append(f"| {name} | {value:.10g} | {se:.6g} |")
        lines += ["", "| Previously computed response-rank coordinate | Estimate | Shared-batch SE |", "|---|---:|---:|"]
        for name, value, se in zip(r["labels"], r["estimate"], r["se"]):
            if name.startswith("response_rank."):
                lines.append(f"| {name} | {value:.10g} | {se:.6g} |")
        lines.append("")
    lines += [result["boundary"], ""]
    (OUT/"REPORT.md").write_text("\n".join(lines)); print("\n".join(lines))


if __name__ == "__main__": main()

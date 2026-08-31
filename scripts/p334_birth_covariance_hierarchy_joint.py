#!/usr/bin/env python3
"""Append supplied hierarchy LOO; report signed shares and paired estimator change."""
import gzip
from hashlib import sha256
import json
from pathlib import Path
import subprocess

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OLD_COMMIT = "e2ef9983f426890a299f5a6e1a2eba8b6d072855"
OLD_DIR = "results/p334-euler-dipole-connected-clock"
HIERARCHY_COMMIT = "44dc9e3396e39105cae85a29d04b39d0afc82d84"
HIERARCHY_PATH = "results/p334-birth-covariance-hierarchy/score.json"
RANKCELL_COMMIT = "2bc3529468fbcba589182acaf98fa4855eb0a85e"
RANKCELL_PATH = "results/p334-rankcell-covariance-transport/score.json"
OUT = ROOT/"results/p334-birth-covariance-hierarchy-joint"
DIRECTIONS = ("plus->S", "minus->D")


def main():
    hashes = {}
    def read(commit, path):
        blob = subprocess.check_output(["git", "show", commit+":"+path], cwd=ROOT)
        hashes[commit+":"+path] = sha256(blob).hexdigest()
        return blob
    old = json.loads(read(OLD_COMMIT, OLD_DIR+"/score.json"))
    hierarchy = json.loads(read(HIERARCHY_COMMIT, HIERARCHY_PATH))
    rankcell = json.loads(read(RANKCELL_COMMIT, RANKCELL_PATH))
    if rankcell["source_commit"] != HIERARCHY_COMMIT:
        raise ValueError("Rankcell transport does not use the supplied hierarchy")
    OUT.mkdir(parents=True, exist_ok=True)
    sizes = {}
    for n in ("325", "425"):
        o, h, r = old["sizes"][n], hierarchy["sizes"][n], rankcell["sizes"][n]
        if h["batch_ids"] != o["batch_ids"] or h["batch_ids"] != r["batch_ids"] or h["batch_ids"] != list(range(20)):
            raise ValueError("Original deleted-batch rows are not aligned")
        saved = json.loads(gzip.decompress(read(OLD_COMMIT, OLD_DIR+"/"+o["complete_covariance_factor_file"])))
        raw = np.array(h["raw_batch_means"])
        raw_mean = raw.mean(axis=0)
        raw_loo = (20*raw_mean-raw)/19
        raw_factor = np.sqrt(19/20)*(raw_loo-raw_loo.mean(axis=0))
        h_loo = np.array(h["LOO"])
        rankcell_loo = np.array(r["LOO"])
        old_loo = np.array(o["derived_LOO"])
        extra_labels, extra_values, extra_loo = [], [], []
        for direction in DIRECTIONS:
            total_idx = h["labels"].index(f"all.{direction}.cov_xy.total")
            total = h["estimate"][total_idx]
            total_loo = h_loo[:, total_idx]
            if total == 0 or np.any(total_loo == 0):
                raise ValueError("The named response-share denominator is undefined")
            for part in ("between_prefixes", "within_prefix"):
                i = h["labels"].index(f"all.{direction}.cov_xy.{part}")
                extra_labels.append(f"{direction}.{part}.signed_response_share")
                extra_values.append(h["estimate"][i]/total)
                extra_loo.append(h_loo[:, i]/total_loo)
            rankcell_idx = r["labels"].index(f"all.{direction}.cov_xy.between_rankcells")
            extra_labels.append(f"{direction}.between_rankcells.signed_response_share")
            extra_values.append(r["estimate"][rankcell_idx]/total)
            extra_loo.append(rankcell_loo[:, rankcell_idx]/total_loo)
            old_idx = o["derived_labels"].index(f"{direction}.delta_intrinsic_rank_cov_tau12")
            extra_labels.append(f"{direction}.exactscore_total_minus_matchedmask_intrinsic")
            extra_values.append(total-o["derived_estimate"][old_idx])
            extra_loo.append(total_loo-old_loo[:, old_idx])
        extra_loo = np.column_stack(extra_loo)
        extra_factor = np.sqrt(19/20)*(extra_loo-extra_loo.mean(axis=0))
        factor = np.column_stack((np.array(saved["factor"]), raw_factor, np.array(h["factor"]), np.array(r["factor"]), extra_factor))
        labels = (saved["labels"]+["covariance_hierarchy.raw."+x for x in h["raw_labels"]]
                  +["covariance_hierarchy."+x for x in h["labels"]]+["rankcell_transport."+x for x in r["labels"]]
                  +["hierarchy_joint."+x for x in extra_labels])
        point = {"covariance_hierarchy."+k: v for k, v in zip(h["labels"], h["estimate"])}
        point.update({"rankcell_transport."+k: v for k, v in zip(r["labels"], r["estimate"])})
        point.update({"hierarchy_joint."+k: v for k, v in zip(extra_labels, extra_values)})
        point.update({"connected_clock."+k: v for k, v in zip(o["derived_labels"], o["derived_estimate"])})
        focused = []
        for direction in DIRECTIONS:
            focused += [f"covariance_hierarchy.all.{direction}.cov_xy.{part}" for part in ("total", "between_prefixes", "within_prefix")]
            focused += [f"hierarchy_joint.{direction}.{part}.signed_response_share" for part in ("between_prefixes", "within_prefix")]
            focused += [f"rankcell_transport.all.{direction}.cov_xy.{part}" for part in ("within_rankcell_prefixes", "between_rankcells", "within_rankcell_total")]
            focused += [f"hierarchy_joint.{direction}.between_rankcells.signed_response_share"]
            focused += [f"connected_clock.{direction}.delta_intrinsic_rank_cov_tau12",
                        f"hierarchy_joint.{direction}.exactscore_total_minus_matchedmask_intrinsic"]
        ix = [labels.index(k) for k in focused]
        covariance = factor[:, ix].T@factor[:, ix]
        filename = f"N{n}.complete_common_factor.json.gz"
        packed = {"labels": labels, "factor": factor.tolist(), "batch_ids": list(range(20)),
            "convention": "factor.T@factor; same original twenty deleted-batch rows; rank at most 19; no inverse",
            "new_raw_labels": h["raw_labels"], "new_raw_20_batch_means": raw.tolist(),
            "hierarchy_derived_labels": h["labels"], "hierarchy_derived_LOO": h["LOO"],
            "rankcell_derived_labels": r["labels"], "rankcell_derived_LOO": r["LOO"],
            "joint_derived_labels": extra_labels, "joint_derived_LOO": extra_loo.tolist()}
        blob = gzip.compress((json.dumps(packed, separators=(",", ":"), allow_nan=False)+"\n").encode(), mtime=0)
        (OUT/filename).write_bytes(blob)
        sizes[n] = {"batch_ids": list(range(20)), "prefix_counts": h["prefix_counts"], "labels": focused,
            "estimate": [point[k] for k in focused], "se": np.sqrt(np.diag(covariance)).tolist(),
            "focused_covariance": covariance.tolist(), "new_labels": extra_labels, "new_estimate": extra_values,
            "new_LOO": extra_loo.tolist(), "complete_covariance_factor_file": filename,
            "complete_covariance_factor_sha256": sha256(blob).hexdigest(), "complete_coordinate_count": len(labels)}
    result = {"schema": "matching-one/p334-birth-covariance-hierarchy-joint/v1",
        "fork_source_commit": old["source_commit"], "previous_shared_commit": OLD_COMMIT,
        "hierarchy_readout_commit": HIERARCHY_COMMIT, "exactscore_moment_commit": hierarchy["source_commit"],
        "rankcell_transport_commit": RANKCELL_COMMIT,
        "source_sha256": hashes, "sizes": sizes, "new_MC": 0, "new_DP": 0, "new_raw_prefix_or_fork_reads": 0,
        "share_definition": "Supplied between-prefix or within-prefix covariance RESPONSE divided by total covariance RESPONSE; ratios recomputed in the supplied original-batch LOO",
        "difference_definition": "Exactscore hierarchy total minus old matched-mask intrinsic normalized-rank covariance; paired old and new LOO, not independent standard errors",
        "boundary": "Signed response shares are not baseline variance fractions or probabilities. Prefix law is unchanged; between-prefix transport is the response of relationships among conditional means. Exact-score and matched-mask values reuse the same e32a8593/959a7fa2 block and estimate the same physical target with different estimators/products. No fit, inverse, cell reanalysis, new sampling or independent evidence."}
    (OUT/"score.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    lines = ["# Covariance response hierarchy: shared-batch shares and estimator difference", ""]
    for n, row in sizes.items():
        lines += [f"## N{n}", "", "| Coordinate | Estimate | Shared-batch SE |", "|---|---:|---:|"]
        for label, mean, se in zip(row["labels"], row["estimate"], row["se"]):
            lines.append(f"| {label} | {mean:.10g} | {se:.6g} |")
        lines.append("")
    lines += [result["boundary"], ""]
    (OUT/"REPORT.md").write_text("\n".join(lines)); print("\n".join(lines))


if __name__ == "__main__": main()

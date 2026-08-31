#!/usr/bin/env python3
"""Append the fixed original-eight-quartet projection to the shared batch factor.

No raw replay, feature selection, projection solve, or simulation is performed.
The supplied producer LOO contains every nonlinear refit and ratio operation.
"""
import gzip
from hashlib import sha256
import json
from pathlib import Path
import subprocess

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OLD_COMMIT = "ce20158a5928e55b67324cba7ed3a18a5c163b39"
OLD_DIR = "results/p334-birth-covariance-hierarchy-joint"
PROJECTION_COMMIT = "9022659843ff0e9c2919c37e9468b0e7b5307268"
PROJECTION_PATH = "results/p334-prefix-response-projection/score.json"
OUT = ROOT / "results/p334-prefix-response-projection-joint"
RECEIVERS = ("first", "second")
RESPONSES = tuple(f"source_{s}.{y}" for s in RECEIVERS for y in ("C", "W"))
CONTACTS = ("joint_safe_mass", "own_score_energy", "own_safe_degree", "own_safe_loop")


def jk_factor(loo):
    return np.sqrt(19 / 20) * (loo - loo.mean(axis=0))


def main():
    hashes = {}

    def read(commit, path):
        data = subprocess.check_output(["git", "show", f"{commit}:{path}"], cwd=ROOT)
        hashes[f"{commit}:{path}"] = sha256(data).hexdigest()
        return data

    old = json.loads(read(OLD_COMMIT, OLD_DIR + "/score.json"))
    projection = json.loads(read(PROJECTION_COMMIT, PROJECTION_PATH))
    if projection["moment_source_commit"] != old["exactscore_moment_commit"]:
        raise ValueError("Projection and previous joint archive use different moment sources")
    OUT.mkdir(parents=True, exist_ok=True)
    sizes = {}
    report = ["# Prefix-response projection: original-block covariance join", "",
              "All entries are estimates ± one original-20-batch jackknife SE. The source",
              "is the original eight-quartet block, not the separately collected new64 block.", ""]
    for n in ("325", "425"):
        o, p = old["sizes"][n], projection["sizes"][n]
        if p["batch_ids"] != o["batch_ids"] or p["batch_ids"] != list(range(20)):
            raise ValueError("Original deleted-batch rows are not aligned")
        if p["prefix_counts"] != o["prefix_counts"] or p["prefix_counts"] != [1000] * 20:
            raise ValueError("Full-population batch denominators differ")
        saved = json.loads(gzip.decompress(read(OLD_COMMIT, OLD_DIR + "/" + o["complete_covariance_factor_file"])))
        raw = np.asarray(p["raw_batch_means"])
        raw_loo = (20 * raw.mean(axis=0) - raw) / 19
        supplied_loo = np.asarray(p["LOO"])
        supplied_factor = np.asarray(p["factor"])
        label_index = {label: i for i, label in enumerate(p["labels"])}

        # Only one new linear contrast: what the four contacts add over score energy.
        gains, gain_labels, gain_loo = [], [], []
        for receiver in RECEIVERS:
            for suffix in ("", "_share"):
                contact = label_index[f"{receiver}.contact.own_intrinsic_cov_loading{suffix}"]
                strength = label_index[f"{receiver}.strength.own_intrinsic_cov_loading{suffix}"]
                gain_labels.append(f"{receiver}.contact_minus_strength.loading{suffix}")
                gains.append(p["estimate"][contact] - p["estimate"][strength])
                gain_loo.append(supplied_loo[:, contact] - supplied_loo[:, strength])
        gain_loo = np.column_stack(gain_loo)
        factor = np.column_stack((np.asarray(saved["factor"]), jk_factor(raw_loo),
                                  supplied_factor, jk_factor(gain_loo)))
        labels = (saved["labels"] + ["prefix_projection.raw." + x for x in p["raw_labels"]]
                  + ["prefix_projection." + x for x in p["labels"]]
                  + ["projection_joint." + x for x in gain_labels])
        # Preserve inherited factor order, including any historical duplicate aliases.
        previous_columns = len(saved["labels"])
        projection_offset = previous_columns + len(p["raw_labels"])
        gain_offset = projection_offset + len(p["labels"])
        old_ix = [saved["labels"].index(k) for k in o["labels"]]
        focused_projection = []
        for receiver in RECEIVERS:
            focused_projection += [f"{receiver}.receiver_R0_mass", f"{receiver}.own_source.intrinsic_cov_loading"]
            for model in ("strength", "contact"):
                focused_projection += [f"{receiver}.{model}.own_intrinsic_cov_loading{s}"
                                       for s in ("", "_share", "_residual")]
            for response in RESPONSES:
                focused_projection += [f"{receiver}.{response}.contact_increment_after_clock"]
                focused_projection += [f"{receiver}.{model}.{response}.projected_response_variance"
                                       for model in ("clock", "contact_clock")]
            focused_projection += [f"{receiver}.clock_partial_cov.{feature}|{response}"
                                   for feature in CONTACTS for response in RESPONSES]
        ix = old_ix + [projection_offset + label_index[k] for k in focused_projection] + list(range(gain_offset, len(labels)))
        focused_factor = factor[:, ix]
        focused_covariance = focused_factor.T @ focused_factor
        focused_labels = o["labels"] + ["prefix_projection." + x for x in focused_projection] + ["projection_joint." + x for x in gain_labels]
        focused_values = o["estimate"] + [p["estimate"][label_index[k]] for k in focused_projection] + gains
        packed = {
            "batch_ids": p["batch_ids"], "prefix_counts": p["prefix_counts"],
            "labels": labels, "factor": factor.tolist(),
            "convention": "factor.T @ factor; original 20 deleted-batch rows, rank <=19; no inverse",
            "blocks": {"previous": [0, previous_columns], "projection_raw": [previous_columns, projection_offset],
                       "projection_derived": [projection_offset, gain_offset], "linear_gains": [gain_offset, len(labels)]},
            "projection_raw_labels": p["raw_labels"], "projection_raw_20_batch_means": p["raw_batch_means"],
            "projection_labels": p["labels"], "projection_estimate": p["estimate"],
            "projection_LOO": p["LOO"], "linear_gain_labels": gain_labels,
            "linear_gain_estimate": gains, "linear_gain_LOO": gain_loo.tolist(),
        }
        blob = gzip.compress((json.dumps(packed, separators=(",", ":"), allow_nan=False) + "\n").encode(), mtime=0)
        filename = f"N{n}.complete_common_factor.json.gz"
        (OUT / filename).write_bytes(blob)
        sizes[n] = {
            "batch_ids": p["batch_ids"], "prefix_counts": p["prefix_counts"],
            "labels": focused_labels, "estimate": focused_values,
            "se": np.sqrt(np.diag(focused_covariance)).tolist(),
            "focused_covariance": focused_covariance.tolist(),
            "projection_labels": p["labels"], "projection_estimate": p["estimate"], "projection_se": p["se"],
            "linear_gain_labels": gain_labels, "linear_gain_estimate": gains,
            "linear_gain_se": np.linalg.norm(jk_factor(gain_loo), axis=0).tolist(),
            "complete_covariance_factor_file": filename, "complete_covariance_factor_sha256": sha256(blob).hexdigest(),
            "complete_coordinate_count": len(labels),
        }

        def show(label, scale=1):
            k = label_index[label]
            return f"{p['estimate'][k] * scale:.6g} ± {np.linalg.norm(supplied_factor[:, k]) * scale:.4g}"

        report += [f"## N{n}", "", "### Own-source signed covariance loading", "",
                   "| Receiver | Strength share (%) | Contact share (%) | Strength residual | Contact residual |",
                   "|---|---:|---:|---:|---:|"]
        for receiver in RECEIVERS:
            report.append(f"| {receiver} | " + " | ".join([
                show(f"{receiver}.strength.own_intrinsic_cov_loading_share", 100),
                show(f"{receiver}.contact.own_intrinsic_cov_loading_share", 100),
                show(f"{receiver}.strength.own_intrinsic_cov_loading_residual"),
                show(f"{receiver}.contact.own_intrinsic_cov_loading_residual")]) + " |")
        report += ["", "### Additional projected response variance after both baseline clocks", "",
                   "| Receiver | Physical source and response | Clock projection | Contact + clock projection | Increment |",
                   "|---|---|---:|---:|---:|"]
        for receiver in RECEIVERS:
            for response in RESPONSES:
                report.append(f"| {receiver} | {response} | " + " | ".join([
                    show(f"{receiver}.clock.{response}.projected_response_variance"),
                    show(f"{receiver}.contact_clock.{response}.projected_response_variance"),
                    show(f"{receiver}.{response}.contact_increment_after_clock")]) + " |")
        report += ["", "### Contact-response cross-moments after projecting out both baseline clocks", "",
                   "| Receiver / contact | source first C | source first W | source second C | source second W |",
                   "|---|---:|---:|---:|---:|"]
        for receiver in RECEIVERS:
            for feature in CONTACTS:
                report.append(f"| {receiver} / {feature} | " + " | ".join(
                    show(f"{receiver}.clock_partial_cov.{feature}|{response}") for response in RESPONSES) + " |")
        report += ["", "### Paired contact-minus-strength gain", "",
                   "| Coordinate | Estimate | Shared-batch SE |", "|---|---:|---:|"]
        for k, mean, se in zip(gain_labels, gains, sizes[n]["linear_gain_se"]):
            report.append(f"| {k} | {mean:.9g} | {se:.6g} |")
        report.append("")
    boundary = (
        "Exploratory fixed six-predictor projections of receiver-rank0, within-rank-cell covariance, "
        "zero padded to each full population. Loading shares are signed response shares, not variance fractions, "
        "probabilities, causal attribution, or an out-of-sample closure result. Contact-after-clock projection "
        "increments are plug-in estimates with supplied refit LOO errors, not confidence lower bounds or R-squared. "
        "The clock-only loading identity is tautological; partial contact cross-moments retain their sampling uncertainty. "
        "All physical receivers, physical sources and C/W responses stay jointly correlated with ce20158a. "
        "The inherited covariance factor has rank at most19 per size; no high-dimensional inverse or omnibus test is used. "
        "The separate new64 block is neither read nor pooled."
    )
    result = {"schema": "matching-one/p334-prefix-response-projection-joint/v1",
              "allocation_commit": "93ee4e9841d15fa57cdd7732f67e2dd232d18f47", "previous_shared_commit": OLD_COMMIT,
              "projection_commit": PROJECTION_COMMIT, "projection_path": PROJECTION_PATH,
              "fork_source_commit": old["fork_source_commit"],
              "exactscore_moment_commit": projection["moment_source_commit"],
              "descriptor_commit": projection["descriptor_commit"], "estimand": projection["estimand"],
              "predictors": projection["predictors"], "responses": projection["responses"], "models": projection["models"],
              "source_sha256": hashes, "script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
              "sizes": sizes, "new_MC": 0, "new_DP": 0, "new_raw_prefix_or_fork_reads": 0,
              "projection_refits_by_this_reader": 0, "boundary": boundary}
    (OUT / "score.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    report += ["## Interpretation boundary", "", boundary, ""]
    (OUT / "REPORT.md").write_text("\n".join(report))
    print(f"Wrote {OUT}; complete common-factor coordinates: " + ", ".join(f"N{n}={r['complete_coordinate_count']}" for n, r in sizes.items()))


if __name__ == "__main__":
    main()

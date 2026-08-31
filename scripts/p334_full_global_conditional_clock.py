#!/usr/bin/env python3
"""Full A_top, same counters and paired prefix-safe conditional replacement."""
import argparse
import csv
import gzip
from hashlib import sha256
import io
import json
from pathlib import Path
import subprocess

import numpy as np
from scipy.stats import binom

from p334_r1_prevalence_clock_loading import four_state_variance, risk_pair_sums, VARIANCE_LABELS

ROOT = Path(__file__).resolve().parents[1]
CONDITIONAL_COMMIT = "0d1e586dafbade5e7d1f9bfc598170d0c881e337"
COND_DIR = "results/p334-paired-clock-loading"
FULL_DIR = "results/p334-full-birth-archive"
OUT = ROOT/"results/p334-full-global-conditional-clock"
OBSERVERS = ("old_gated_R1_F2", "baseline_full_A", "safe_full_A")


def observation_fields(ranks, f1, f2, hybrid, old_y, delta, alpha):
    fields = {}
    for orientation, o in (("first", 0), ("second", 1)):
        for r in range(3):
            fields[f"rank{r}.{orientation}"] = (ranks[:, o] == r).astype(float)
    for name, f2v in (("baseline", f2), ("safe", hybrid)):
        a = f1+f2v-1
        for endpoint, start in (("p_ref", 0), ("p_integral", 2)):
            fields[f"{name}.{endpoint}.F1_H4"] = (f1[:, start]-f1[:, start+1])/delta
            fields[f"{name}.{endpoint}.F2_H4"] = (f2v[:, start]-f2v[:, start+1])/delta
            fields[f"{name}.{endpoint}.A_H4"] = (a[:, start]-a[:, start+1])/delta
            for r in range(3):
                contributions = a[:, start:start+2]*(ranks == r)
                for orientation, o in (("first", 0), ("second", 1)):
                    fields[f"{name}.{endpoint}.A_R{r}_{orientation}"] = contributions[:, o]
                fields[f"{name}.{endpoint}.A_R{r}_H4"] = (contributions[:, 0]-contributions[:, 1])/delta
                if start == 2:
                    fields[f"{name}.{endpoint}.centered_A_R{r}_H4"] = fields[f"{name}.{endpoint}.A_R{r}_H4"]-alpha*((ranks[:, 0] == r).astype(float)-(ranks[:, 1] == r))/delta
    for endpoint, start in (("p_ref", 0), ("p_integral", 2)):
        fields[f"old_gated.{endpoint}.F2_R1_H4"] = (old_y[:, start]-old_y[:, start+1])/delta
    return fields


def aux_readout(sums, residual_second_mean):
    vectors, details = [], {}
    for name, block in zip(OBSERVERS, sums):
        vector, info = four_state_variance(block)
        vectors.append(vector)
        details[name] = info
    baseline_var = np.diag(details["baseline_full_A"]["total_individual_covariance"])
    safe_var = np.diag(details["safe_full_A"]["total_individual_covariance"])
    ratios = np.r_[safe_var/baseline_var, residual_second_mean/baseline_var]
    return np.r_[np.concatenate(vectors), ratios], details


def jackknife_covariance(loo):
    centered = loo-loo.mean(axis=0)
    return (len(loo)-1)/len(loo)*centered.T@centered


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full-commit", required=True)
    parser.add_argument("--output", type=Path, default=OUT)
    args = parser.parse_args()
    full_commit = subprocess.check_output(["git", "rev-parse", args.full_commit+"^{commit}"], cwd=ROOT, text=True).strip()
    hashes = {}
    def read(commit, path):
        blob = subprocess.check_output(["git", "show", commit+":"+path], cwd=ROOT)
        hashes[commit+":"+path] = sha256(blob).hexdigest()
        return blob
    conditional = json.loads(read(CONDITIONAL_COMMIT, COND_DIR+"/score.json"))
    p = conditional["contract"]["p_ref"]
    alpha = 1-2*p
    sizes = {}
    for n in (325, 425):
        metadata = json.loads(read(full_commit, f"{FULL_DIR}/N{n}.csv.metadata.json"))
        full_rows = [{k: int(v) for k, v in row.items()} for row in csv.DictReader(io.StringIO(read(full_commit, f"{FULL_DIR}/N{n}.csv").decode()))]
        if len(full_rows) != 20000 or len({r["counter"] for r in full_rows}) != 20000:
            raise ValueError("Need the complete original 20000-counter archive")
        delta = conditional["sizes"][str(n)]["delta_cos4"]
        by_batch = {b: sorted([r for r in full_rows if r["batch"] == b], key=lambda r:r["counter"]) for b in range(20)}
        batch_means, batch_grams, pair_batches, residual_batches, a_batches = [], [], [], [], []
        statuses = {"eligible_with_R1_replacement": 0, "both_R2_no_change": 0,
                    "blocked_any_R0": 0, "blocked_solver_failure": 0}
        for b in range(20):
            rows = by_batch[b]
            payload = json.loads(gzip.decompress(read(CONDITIONAL_COMMIT, f"{COND_DIR}/batches/N{n}.batch{b:02d}.json.gz")))
            records = {r["counter"]: r for r in payload["records"]}
            if len(rows) != 1000 or set(records) != {r["counter"] for r in rows}:
                raise ValueError("Full and conditional original counter batches do not align")
            k1 = np.array([[r["first_k1"], r["second_k1"]] for r in rows])
            k2 = np.array([[r["first_k2"], r["second_k2"]] for r in rows])
            ranks = np.array([[r["first_rank"], r["second_rank"]] for r in rows])
            k0 = np.array([r["k0"] for r in rows])[:, None]
            if not np.array_equal(ranks, np.where(k0 < k1, 0, np.where(k0 < k2, 1, 2))):
                raise ValueError("Replayed rank is inconsistent with the full ordered birth times")
            ordered = [records[r["counter"]] for r in rows]
            risks = (ranks == 1).astype(float)
            if not np.array_equal(risks, np.array([r["risk"] for r in ordered])):
                raise ValueError("Replayed R1 membership differs from the original conditional archive")
            f1 = np.column_stack((binom.sf(k1-1, n, p), (n+1-k1)/(n+1)))
            f2 = np.column_stack((binom.sf(k2-1, n, p), (n+1-k2)/(n+1)))
            hybrid = f2.copy()
            old_y = np.array([r["Y"] for r in ordered])
            for i, record in enumerate(ordered):
                if np.any(ranks[i] == 0):
                    statuses["blocked_any_R0"] += 1
                elif np.all(ranks[i] == 2):
                    statuses["both_R2_no_change"] += 1
                elif record["status"] == "exact_pair" and all(record["clocks"][o] is not None for o in range(2) if ranks[i, o] == 1):
                    statuses["eligible_with_R1_replacement"] += 1
                    for o in range(2):
                        if ranks[i, o] == 1:
                            hybrid[i, o], hybrid[i, 2+o] = record["clocks"][o]["conditional"]
                else:
                    statuses["blocked_solver_failure"] += 1
            fields = observation_fields(ranks, f1, f2, hybrid, old_y, delta, alpha)
            labels = list(fields)
            matrix = np.column_stack(list(fields.values()))
            batch_means.append(matrix.mean(axis=0))
            batch_grams.append(matrix.T@matrix/len(rows))
            baseline_a, safe_a = f1+f2-1, f1+hybrid-1
            pair_batches.append([risk_pair_sums(risks, y, delta) for y in (old_y, baseline_a, safe_a)])
            residual = np.column_stack(((baseline_a[:, 0]-baseline_a[:, 1])-(safe_a[:, 0]-safe_a[:, 1]),
                                        (baseline_a[:, 2]-baseline_a[:, 3])-(safe_a[:, 2]-safe_a[:, 3])))/delta
            residual_batches.append(np.mean(residual**2, axis=0))
            a_batches.append(np.r_[baseline_a.mean(axis=0), safe_a.mean(axis=0)])
        batch_means, pair_batches = np.array(batch_means), np.array(pair_batches)
        residual_batches = np.array(residual_batches)
        mean = batch_means.mean(axis=0)
        loo = (20*mean-batch_means)/19
        empirical_cov = np.mean(batch_grams, axis=0)-np.outer(mean, mean)
        aux, info = aux_readout(pair_batches.sum(axis=0), residual_batches.mean(axis=0))
        aux_loo = np.array([aux_readout(pair_batches.sum(axis=0)-pair_batches[b],
                         (residual_batches.sum(axis=0)-residual_batches[b])/19)[0] for b in range(20)])
        aux_labels = [name+"."+field for name in OBSERVERS for field in VARIANCE_LABELS]
        aux_labels += ["safe_over_baseline_variance.p_ref", "safe_over_baseline_variance.p_integral",
                       "paired_removed_noise_fraction.p_ref", "paired_removed_noise_fraction.p_integral"]
        joint = np.column_stack((loo, aux_loo))
        joint_cov = jackknife_covariance(joint)
        centered = {}
        for variant in ("baseline", "safe"):
            keys = ([f"{variant}.p_integral.A_R{r}_H4" for r in range(3)]
                    +[f"{variant}.p_integral.centered_A_R{r}_H4" for r in range(3)]
                    +[f"{variant}.p_integral.A_H4"])
            ix = [labels.index(k) for k in keys]
            centered[variant] = {"labels": keys, "mean": mean[ix].tolist(),
                "individual_covariance": empirical_cov[np.ix_(ix, ix)].tolist(),
                "mean_covariance": joint_cov[np.ix_(ix, ix)].tolist(),
                "sum_uncentered_minus_global": float(mean[ix[:3]].sum()-mean[ix[-1]]),
                "sum_centered_minus_global": float(mean[ix[3:6]].sum()-mean[ix[-1]])}
        sizes[str(n)] = {"delta_cos4": delta, "metadata": metadata, "policy_counts": statuses,
            "batch_ids": list(range(20)), "samples_per_batch": 1000, "labels": labels,
            "mean": mean.tolist(), "mean_se": np.sqrt(np.diag(joint_cov)[:len(labels)]).tolist(),
            "batch_means": batch_means.tolist(), "individual_variance": np.diag(empirical_cov).tolist(),
            "A_orientation_fields": [v+"."+ep+"."+o for v in ("baseline", "safe") for ep in ("p_ref", "p_integral") for o in ("first", "second")],
            "A_orientation_batch_means": np.array(a_batches).tolist(),
            "four_R1_state_variance": info, "auxiliary_labels": aux_labels, "auxiliary_estimate": aux.tolist(),
            "auxiliary_se": np.sqrt(np.diag(joint_cov)[len(labels):]).tolist(),
            "joint_labels": labels+aux_labels, "joint_LOO_vectors": joint.tolist(),
            "joint_covariance": joint_cov.tolist(), "joint_rank_at_most": 19,
            "fixed_origin_centered_integral_strata": centered}
    result = {"schema": "matching-one/p334-full-global-prefix-safe-clock/v1",
              "full_birth_commit": full_commit, "conditional_commit": CONDITIONAL_COMMIT,
              "source_sha256": hashes, "p_ref": p, "fixed_integral_origin_alpha": alpha,
              "sizes": sizes, "new_MC": 0, "new_DP": 0, "new_birth_replays": 0,
              "policy": "Replace R1 F2 only if both orientations have checkpoint rank >=1 and every R1 clock is exact_pair; any R0 or failure keeps the entire baseline pair. Both R2 require no change.",
              "identities": ["A=F1+F2-1", "integral A=1-(K1+K2)/(N+1)",
                             "sum_r gated A_r=A", "sum_r centered A_r=sum_r gated A_r; alpha=1-2*p_ref"],
              "boundary": "Full A_top on the original paired counter blocks. Baseline, safe hybrid, gated strata and centering are correlated readouts of the same archive, never independent evidence or separately summed errors. Four-state fractions condition on the old R1 flags only. No high-dimensional covariance inverse or new omnibus test."}
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output/"score.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    lines = ["# Complete A_top under the paired prefix-safe conditional policy", ""]
    for n, row in sizes.items():
        lines += [f"## N{n}", "", str(row["policy_counts"]), ""]
        for variant in ("baseline", "safe"):
            for ep in ("p_ref", "p_integral"):
                k = row["labels"].index(f"{variant}.{ep}.A_H4")
                lines.append(f"- {variant} {ep} A_H4 = {row['mean'][k]:.10g} +/- {row['mean_se'][k]:.6g}")
        lines += ["", "| Observer | canonical between-R1-state fraction | integral between-R1-state fraction |", "|---|---:|---:|"]
        for k, name in enumerate(OBSERVERS):
            a, se = row["auxiliary_estimate"], row["auxiliary_se"]
            lines.append(f"| {name} | {100*a[k*11+3]:.7g}% +/- {100*se[k*11+3]:.5g}pp | {100*a[k*11+7]:.7g}% +/- {100*se[k*11+7]:.5g}pp |")
        lines += ["", f"Safe/baseline variance ratios (canonical, integral): {row['auxiliary_estimate'][-4:-2]}",
                  f"Paired residual-noise fractions: {row['auxiliary_estimate'][-2:]}", ""]
    lines += [result["boundary"], ""]
    (args.output/"REPORT.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()

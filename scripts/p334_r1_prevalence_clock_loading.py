#!/usr/bin/env python3
"""Symmetric prevalence/conditional-clock decomposition of a final R1 archive.

Reads only immutable git blobs from an explicitly supplied final commit.
It never accesses the producer's in-progress working directory.
"""
import argparse
import gzip
from hashlib import sha256
import json
from pathlib import Path
import subprocess

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT/"results/p334-r1-prevalence-clock-loading"
ENDPOINTS = ("p_ref", "p_integral")
LABELS = ["risk_first", "risk_second"] + [f"{endpoint}.{field}"
    for endpoint in ENDPOINTS for field in ("Y_first", "Y_second", "m_first", "m_second",
                                          "prevalence_C", "conditional_clock_L", "H4_difference_D")]
VARIANCE_LABELS = [f"{endpoint}.{field}" for endpoint in ENDPOINTS
                  for field in ("between_Rpair_variance", "within_Rpair_variance",
                                "total_individual_variance", "between_Rpair_fraction")]
VARIANCE_LABELS += ["cross_endpoint.between_covariance", "cross_endpoint.within_covariance",
                   "cross_endpoint.total_covariance"]


def decompose(mean, delta_cos4):
    """Input: rf,rs,Yf_pref,Ys_pref,Yf_integral,Ys_integral; all unconditioned means."""
    mean = np.asarray(mean, dtype=float)
    rf, rs = mean[:2]
    if mean.shape != (6,) or min(rf, rs) <= 0 or not np.isfinite(delta_cos4) or delta_cos4 == 0:
        raise ValueError("Both R1 risks and a nonzero source delta_cos4 are required")
    values = [rf, rs]
    for start in (2, 4):
        yf, ys = mean[start:start+2]
        mf, ms = yf/rf, ys/rs
        prevalence = .5*(rf-rs)*(mf+ms)/delta_cos4
        loading = .5*(rf+rs)*(mf-ms)/delta_cos4
        difference = (yf-ys)/delta_cos4
        values.extend((yf, ys, mf, ms, prevalence, loading, difference))
    return np.array(values)


def score_batch_means(batch_means, delta_cos4):
    b = len(batch_means)
    mean = np.mean(batch_means, axis=0)
    point = decompose(mean, delta_cos4)
    loo = np.array([decompose((b*mean-row)/(b-1), delta_cos4) for row in batch_means])
    centered = loo-loo.mean(axis=0)
    cov = (b-1)/b*centered.T@centered
    residuals = np.stack((loo[:, 6]+loo[:, 7]-loo[:, 8], loo[:, 13]+loo[:, 14]-loo[:, 15]))
    return {"labels": LABELS, "estimate": point.tolist(), "se": np.sqrt(np.diag(cov)).tolist(),
            "full_covariance": cov.tolist(), "leave_one_common_batch_out_vectors": loo.tolist(),
            "joint_batch_means_risk_then_Y": np.asarray(batch_means).tolist(),
            "max_LOO_additive_identity_residual": float(np.max(np.abs(residuals))),
            "jackknife_bias_estimate": ((b-1)*(loo.mean(axis=0)-point)).tolist()}


def risk_pair_sums(risks, y, delta):
    codes = (2*risks[:, 0]+risks[:, 1]).astype(int)
    contrast = np.column_stack((y[:, 0]-y[:, 1], y[:, 2]-y[:, 3]))/delta
    columns = [np.ones(len(codes)), contrast[:, 0], contrast[:, 1],
               contrast[:, 0]**2, contrast[:, 0]*contrast[:, 1], contrast[:, 1]**2]
    return np.array([np.bincount(codes, weights=x, minlength=4) for x in columns]).T


def four_state_variance(sums):
    """Finite empirical law: each original counter receives weight 1/n."""
    count = sums[:, 0]
    total_count = count.sum()
    probabilities = count/total_count
    conditional_mean = np.divide(sums[:, 1:3], count[:, None],
                                 out=np.zeros((4, 2)), where=count[:, None] > 0)
    mean = sums[:, 1:3].sum(axis=0)/total_count
    within, between, states = np.zeros((2, 2)), np.zeros((2, 2)), []
    for k in range(4):
        second = np.array([[sums[k, 3], sums[k, 4]], [sums[k, 4], sums[k, 5]]])
        covariance = second/count[k]-np.outer(conditional_mean[k], conditional_mean[k]) if count[k] else np.zeros((2, 2))
        deviation = conditional_mean[k]-mean
        between_contribution = probabilities[k]*np.outer(deviation, deviation)
        within += probabilities[k]*covariance
        between += between_contribution
        states.append({"Rpair": [k//2, k%2], "counter_count": int(count[k]),
                       "probability": float(probabilities[k]),
                       "conditional_mean": conditional_mean[k].tolist() if count[k] else None,
                       "conditional_covariance": covariance.tolist() if count[k] else None,
                       "between_covariance_contribution": between_contribution.tolist()})
    total = within+between
    if np.any(np.diag(total) <= 0):
        raise ValueError("A zero-variance endpoint has no prevalence variance fraction")
    fraction = np.diag(between)/np.diag(total)
    vector = np.r_[[x for k in range(2) for x in
                    (between[k, k], within[k, k], total[k, k], fraction[k])],
                   between[0, 1], within[0, 1], total[0, 1]]
    return vector, {"mean_hybrid_contrast": mean.tolist(), "states": states,
                    "between_covariance": between.tolist(), "within_covariance": within.tolist(),
                    "total_individual_covariance": total.tolist(),
                    "variance_denominator": int(total_count),
                    "scope": "Exact ddof=0 decomposition of this finite empirical counter distribution; not an unbiased population variance estimator or a causal mediation fraction"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-commit", required=True, help="Completed producer commit supplied by its owner")
    parser.add_argument("--source-directory", default="results/p334-paired-clock-loading")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    commit = subprocess.check_output(["git", "rev-parse", "--verify", args.source_commit+"^{commit}"], cwd=ROOT, text=True).strip()
    hashes = {}
    def read(path):
        blob = subprocess.check_output(["git", "show", commit+":"+path], cwd=ROOT)
        hashes[path] = sha256(blob).hexdigest()
        return blob
    # Final score is mandatory; an archive with only some completed batches is
    # intentionally not a valid input, even if those batches have been committed.
    source = json.loads(read(args.source_directory+"/score.json"))
    contract = source["contract"]
    count, samples = contract["batches_per_size"], contract["samples_per_size"]
    if count != 20 or samples != 20000 or samples % count:
        raise ValueError("This fixed analysis expects the declared 20 x 1000 paired archive")
    sizes = {}
    for n in contract["sizes"]:
        reference = source["sizes"][str(n)]
        delta = reference["delta_cos4"]
        batch_means, pair_sums, seen, contract_hashes, k0s = [], [], set(), set(), set()
        for batch in range(count):
            path = f"{args.source_directory}/batches/N{n}.batch{batch:02d}.json.gz"
            payload = json.loads(gzip.decompress(read(path)))
            if (payload["N"], payload["batch"], payload["source_commit"]) != (n, batch, contract["source_commit"]):
                raise ValueError("Batch/source identity does not match final producer contract")
            contract_hashes.add(payload["contract_sha256"])
            records = payload["records"]
            counters = {row["counter"] for row in records}
            if len(records) != samples//count or len(counters) != len(records) or seen & counters:
                raise ValueError("Original equal-size counter batches are not complete and disjoint")
            seen |= counters
            risks = np.array([row["risk"] for row in records], dtype=float)
            y = np.array([row["Y"] for row in records], dtype=float)
            if risks.shape != (samples//count, 2) or y.shape != (samples//count, 4):
                raise ValueError("Unexpected paired risk/Y layout")
            if not np.all(np.isin(risks, [0., 1.])) or not np.all(np.isfinite(y)) or np.any(y[np.tile(risks == 0, (1, 2))] != 0):
                raise ValueError("Y must be a finite R1-weighted contribution, zero outside its orientation's R1 risk")
            batch_means.append(np.r_[risks.mean(axis=0), y.mean(axis=0)])
            pair_sums.append(risk_pair_sums(risks, y, delta))
            k0s.update(int(row["k0"]) for record in records for row in record["source_rows"] if row is not None)
        if len(contract_hashes) != 1:
            raise ValueError("Mixed producer contracts inside one size")
        sizes[str(n)] = {"delta_cos4": delta, "original_batch_count": count,
                         "original_counter_count": len(seen), "producer_contract_sha256": contract_hashes.pop(),
                         "raw_source_sha256": reference["source_sha256"],
                         "paired_hybrid_status_counts": reference["pair_status_counts"],
                         **score_batch_means(np.array(batch_means), delta)}
        row = sizes[str(n)]
        pair_sums = np.array(pair_sums)
        summed = pair_sums.sum(axis=0)
        variance_point, variance_details = four_state_variance(summed)
        variance_loo = np.array([four_state_variance(summed-batch)[0] for batch in pair_sums])
        variance_centered = variance_loo-variance_loo.mean(axis=0)
        variance_cov = (count-1)/count*variance_centered.T@variance_centered
        joint_loo = np.column_stack((row["leave_one_common_batch_out_vectors"], variance_loo))
        joint_centered = joint_loo-joint_loo.mean(axis=0)
        row["four_state_hybrid_contrast_variance"] = {
            "labels": VARIANCE_LABELS, "estimate": variance_point.tolist(),
            "se": np.sqrt(np.diag(variance_cov)).tolist(), "full_covariance": variance_cov.tolist(),
            "leave_one_common_batch_out_vectors": variance_loo.tolist(),
            "batch_sufficient_sums_count_pref_int_pref2_prefint_int2": pair_sums.tolist(),
            **variance_details}
        row["joint_mean_and_variance_readout_covariance"] = {
            "labels": LABELS+VARIANCE_LABELS,
            "full_covariance": ((count-1)/count*joint_centered.T@joint_centered).tolist(),
            "leave_one_common_batch_out_vectors": joint_loo.tolist(),
            "rank_at_most": count-1,
            "scope": "C, L, and the four-state variance decomposition share these same 20 batches; they are not independent tests"}
        if len(k0s) == 1:
            k0 = next(iter(k0s))
            m = np.array(row["estimate"])[[11, 12]]
            m_cov = np.array(row["full_covariance"])[np.ix_([11, 12], [11, 12])]
            row["integrated_waiting_interpretation"] = {
                "k0": k0, "d": n-k0,
                "hybrid_mean_wait_given_R1_first_second": (n-k0+1-(n+1)*m).tolist(),
                "mean_wait_covariance": ((n+1)**2*m_cov).tolist(),
                "identity": "m_i_int=(d+1-E_hybrid[T_i|R1])/(N+1)",
                "scope": "Conditional means include the declared whole-pair fallbacks"}
    result = {"schema": "matching-one/p334-r1-prevalence-clock-loading/v1", "source_commit": commit,
              "source_directory": args.source_directory, "source_sha256": hashes, "source_contract": contract,
              "producer_schema_commit": "5b81e2f9", "sizes": sizes,
              "identities": ["D=(Yf-Ys)/delta_cos4", "C=(rf-rs)*(mf+ms)/(2*delta_cos4)",
                             "L=(rf+rs)*(mf-ms)/(2*delta_cos4)", "mi=Yi/ri; D=C+L"],
              "boundary": "R1-layer weighted hybrid Y only: conditional suffix means on solved pairs and original suffix readouts on whole-pair fallbacks. Not full F2/A_top, not a causal decomposition, not an independent random block, and no attribution percentages when C and L cancel. Sizes are reported separately; no cross-size independence is assumed here.",
              "new_MC": 0, "network_reruns": 0}
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output/"score.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    lines = ["# R1 prevalence versus conditional-clock loading", "",
             "| N / readout | C: prevalence | SE | L: conditional clock | SE | C+L: H4 difference | SE |",
             "|---|---:|---:|---:|---:|---:|---:|"]
    for n, row in sizes.items():
        value, se = row["estimate"], row["se"]
        for endpoint, start in zip(ENDPOINTS, (6, 13)):
            lines.append(f"| {n} / {endpoint} | {value[start]:.10g} | {se[start]:.6g} | {value[start+1]:.10g} | {se[start+1]:.6g} | {value[start+2]:.10g} | {se[start+2]:.6g} |")
    lines += ["", "## Variance explained by the four R1 pair states", "",
              "| N / readout | Between-state variance | Within-state variance | Between fraction | Batch SE |",
              "|---|---:|---:|---:|---:|"]
    for n, row in sizes.items():
        variance = row["four_state_hybrid_contrast_variance"]
        value, se = variance["estimate"], variance["se"]
        for endpoint, start in zip(ENDPOINTS, (0, 4)):
            lines.append(f"| {n} / {endpoint} | {value[start]:.10g} | {value[start+1]:.10g} | {100*value[start+3]:.7g}% | {100*se[start+3]:.5g} percentage points |")
    lines += ["", "These variance terms describe individual hybrid contrasts, not their mean-estimator covariance. The 20-batch joint uncertainty for both decompositions is saved separately.", "", result["boundary"], ""]
    (args.output/"REPORT.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()

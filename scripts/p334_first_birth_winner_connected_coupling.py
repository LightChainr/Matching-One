#!/usr/bin/env python3
"""Connected first-birth/winner coupling from fixed same-gate sufficient means."""
from collections import Counter
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import subprocess

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results/p334-marked-global-topology-loading/score.json"
OUT = ROOT / "results/p334-first-birth-winner-connected-coupling"
MARKED_COMMIT = "2dd865f0b26a4d5d43f52b300293016e6ffd19b8"
BIRTH_COMMIT = "9c495ab13e65f2bc93dc0849ee3b73f88724c4b1"
ORIENTATIONS = ("first", "second")
SUFFICIENT = ("r", "E_I_K1", "E_I_piD", "E_I_mu", "E_I_hmu", "E_I_K1_piD")
ORIENTATION_FIELDS = ("r", "conditional_K1", "conditional_piD", "conditional_mu",
                      "conditional_hmu", "conditional_K1_piD", "conditional_cov_K1_piD",
                      "direct_positive_F2", "direct_marginal_debt", "direct_connected_debt",
                      "direct_A_without_connected", "direct_A_full",
                      "collective_positive_F2", "collective_marginal_debt", "collective_connected_debt",
                      "collective_A_without_connected", "collective_A_full")
CONTRAST_FIELDS = ("positive_F2", "marginal_debt_prevalence", "marginal_debt_K1",
                   "marginal_debt_winner", "marginal_debt_total", "connected_debt",
                   "A_without_connected", "A_full")


def derive(sufficient, n, d, delta):
    """Pool first, then form all products/ratios; repeated for every deleted batch."""
    raw = np.asarray(sufficient).reshape(2, 6)
    r, k, p, mu, hmu, kp = raw.T
    kb, pb, mub = k/r, p/r, mu/r
    cov = kp/r-kb*pb
    positive = np.column_stack([hmu, (d+1)*r-mu-hmu])/(n+1)
    marginal = np.column_stack([k*p/r, k-k*p/r])/(n+1)
    connected = np.column_stack([r*cov, -r*cov])/(n+1)
    without = positive-marginal
    full = without-connected
    orientations = np.column_stack([r, kb, pb, mub, hmu/r, kp/r, cov,
        positive[:,0], marginal[:,0], connected[:,0], without[:,0], full[:,0],
        positive[:,1], marginal[:,1], connected[:,1], without[:,1], full[:,1]])
    contrasts = []
    for source in (0,1):
        winner = pb if source == 0 else 1-pb
        prevalence = (r[0]-r[1]) * np.mean(kb*winner) / ((n+1)*delta)
        k1_term = np.mean(r)*np.mean(winner)*(kb[0]-kb[1])/((n+1)*delta)
        winner_term = np.mean(r)*np.mean(kb)*(winner[0]-winner[1])/((n+1)*delta)
        diff = lambda values: float((values[0,source]-values[1,source])/delta)
        contrasts.append([diff(positive), prevalence, k1_term, winner_term,
                          diff(marginal), diff(connected), diff(without), diff(full)])
    vector = np.r_[orientations.reshape(-1), np.array(contrasts).reshape(-1)]
    return vector


def main():
    source_blob = SOURCE.read_bytes()
    old = json.loads(source_blob)
    labels = ([f"{o}_{f}" for o in ORIENTATIONS for f in ORIENTATION_FIELDS] +
              [f"H4_{s}_{f}" for s in ("direct", "collective") for f in CONTRAST_FIELDS])
    output = {"marked_source_commit": MARKED_COMMIT, "marked_source_path": str(SOURCE.relative_to(ROOT)),
              "marked_source_sha256": hashlib.sha256(source_blob).hexdigest(),
              "birth_commit": BIRTH_COMMIT, "sufficient_statistic_labels": [f"{o}_{s}" for o in ORIENTATIONS for s in SUFFICIENT],
              "derived_labels": labels,
              "selection": "I=that orientation is R1 and the original global whole-pair gate accepts; same population and policy as 2dd865f0.",
              "denominator": "All original20000 counters per size. Conditional moments are pooled weighted means divided by the same I frequency.",
              "LOO_rule": "Delete one original1000-counter batch, pool the remaining19 batch sufficient means, and then recompute every conditional ratio/product/covariance. Never average batchwise products.",
              "new_DP": 0, "new_MC": 0, "new_replay": 0, "new_clock_evaluation": 0,
              "sizes": {}}
    report = ["# First-birth/winner connected coupling on the fixed global source gate", "",
              "Point decompositions; the common covariance coordinator receives the original20 sufficient-statistic batches and all correctly re-pooled LOO vectors.", ""]
    for n in (325,425):
        row = old["sizes"][str(n)]
        blob = subprocess.check_output(["git", "show", BIRTH_COMMIT+f":results/p334-full-birth-archive/N{n}.csv"], cwd=ROOT)
        births = {int(r["counter"]): {k:int(v) for k,v in r.items()}
                  for r in csv.DictReader(io.StringIO(blob.decode()))}
        q = np.zeros((20,2,6))
        clock_means = np.array(row["source_clock_20_batch_means"]).reshape(20,2,4)
        positive = np.array(row["positive_F2_20_batch_means"]).reshape(20,2,2,2)
        debt = np.array(row["first_birth_debt_20_batch_means"]).reshape(20,2,2,2)
        q[:,:,2] = clock_means[:,:,0]
        q[:,:,3] = clock_means[:,:,2]+clock_means[:,:,3]
        q[:,:,4] = (n+1)*positive[:,:,1,0]
        q[:,:,5] = (n+1)*debt[:,:,1,0]
        selected = Counter()
        # Only gate and K1 are read here. Saved coefficients are never evaluated.
        for batch in range(20):
            path = ROOT / f"results/p334-paired-clock-loading/batches/N{n}.batch{batch:02d}.json.gz"
            with gzip.open(path, "rt") as stream:
                records = json.load(stream)["records"]
            for record in records:
                b = births[record["counter"]]
                if record["status"] != "exact_pair" or min(b["first_rank"],b["second_rank"]) < 1:
                    continue
                for o, orientation in enumerate(ORIENTATIONS):
                    if b[orientation+"_rank"] == 1:
                        q[batch,o,0] += 1/1000
                        q[batch,o,1] += b[orientation+"_k1"]/1000
                        selected[orientation] += 1
        point = derive(q.mean(axis=0), n, row["d"], row["delta_cos4"])
        loo = np.array([derive((q.sum(axis=0)-q[b])/19,n,row["d"],row["delta_cos4"])
                        for b in range(20)])
        orientation_point = point[:34].reshape(2,17)
        contrast_point = point[34:].reshape(2,8)
        saved_marked_h4 = np.array(row["H4_means"])[1,:2]
        result = {"batch_ids":list(range(20)), "batch_counter_count":1000, "counter_count":20000,
                  "selected_orientation_counts":dict(selected),"d":row["d"],"delta_cos4":row["delta_cos4"],
                  "birth_sha256":hashlib.sha256(blob).hexdigest(),
                  "joint_20_batch_sufficient_means":q.reshape(20,12).tolist(),
                  "pooled_sufficient_means":q.mean(axis=0).reshape(-1).tolist(),
                  "derived_point":point.tolist(),"leave_one_original_batch_out":loo.tolist(),
                  "orientation_fields":ORIENTATION_FIELDS,"orientation_point":orientation_point.tolist(),
                  "H4_contrast_fields":CONTRAST_FIELDS,"H4_contrast_point":contrast_point.tolist(),
                  "marginal_only_direction_vs_completion":
                    {s: bool(contrast_point[i,0]*contrast_point[i,6]<0) for i,s in enumerate(("direct","collective"))},
                  "connected_changes_marginal_only_direction":
                    {s: bool(contrast_point[i,6]*contrast_point[i,7]<0) for i,s in enumerate(("direct","collective"))},
                  "max_saved_marked_contrast_difference":float(np.max(np.abs(contrast_point[:,7]-saved_marked_h4))),
                  "point_connected_source_sum":float(contrast_point[:,5].sum()),
                  "max_LOO_connected_source_sum_abs":float(np.max(np.abs(loo[:,39]+loo[:,47])))}
        output["sizes"][str(n)] = result
        report += [f"## N{n}","",f"Selected counts: {dict(selected)}; original denominator20000.","",
                   "| Orientation | r | E[K1|I] | E[piD|I] | Cov(K1,piD|I) |",
                   "|---|---:|---:|---:|---:|"]
        for o,orientation in enumerate(ORIENTATIONS):
            v=orientation_point[o]
            report.append(f"| {orientation} | {v[0]:.8g} | {v[1]:.10g} | {v[2]:.10g} | {v[6]:.10g} |")
        report += ["","| H4 integral source | Completion | Marginal debt | Connected debt | A without connected | Full A |",
                   "|---|---:|---:|---:|---:|---:|"]
        for i,s in enumerate(("direct","collective")):
            v=contrast_point[i]
            report.append(f"| {s} | {v[0]:.10g} | {v[4]:.10g} | {v[5]:.10g} | {v[6]:.10g} | {v[7]:.10g} |")
        report += ["","| H4 marginal first-birth debt | Prevalence | Conditional K1 | Winner share | Sum |",
                   "|---|---:|---:|---:|---:|"]
        for i,s in enumerate(("direct","collective")):
            v=contrast_point[i]
            report.append(f"| {s} | {v[1]:.10g} | {v[2]:.10g} | {v[3]:.10g} | {v[4]:.10g} |")
        report.append("")
        print(f"N{n} conditional={orientation_point[:,:7].tolist()} H4={contrast_point.tolist()}",flush=True)
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/"score.json").write_text(json.dumps(output,indent=2,allow_nan=False)+"\n")
    (OUT/"REPORT.md").write_text("\n".join(report))


if __name__ == "__main__":
    main()

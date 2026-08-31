#!/usr/bin/env python3
"""Complete A/E in nine joint checkpoint cells on saved original paired paths."""
import csv
import gzip
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np
from scipy.stats import binom


ROOT = Path(__file__).resolve().parents[1]
COND = ROOT / "results/p334-paired-clock-loading"
OUT = ROOT / "results/p334-nine-layer-complete-ae"
FIELDS = [f"{variant}.{endpoint}.{observer}" for variant in ("baseline", "safe")
          for endpoint in ("p_ref", "p_integral") for observer in ("A", "E")]


def main():
    contract = json.loads((COND / "score.json").read_text())
    p = contract["contract"]["p_ref"]
    result = {
        "full_birth_commit": "9c495ab13e65f2bc93dc0849ee3b73f88724c4b1",
        "conditional_commit": "0d1e586dafbade5e7d1f9bfc598170d0c881e337",
        "global_policy_commit": "3edc785a0312e4dce688bc6966593780907abc51",
        "analysis_code_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "p_ref": p,
        "new_MC": 0, "new_DP": 0, "new_path_replays": 0,
        "effect_fields": FIELDS,
        "cell_fields": ["mass"] + FIELDS,
        "labels": [f"R{r}{s}.{field}" for r in range(3) for s in range(3)
                   for field in ["mass"] + FIELDS],
        "source_sha256": {}, "sizes": {},
        "covariance_handoff": "Original20 batch vectors to the global coordinator. Cells are disjoint: all cross-cell individual raw second moments are exactly zero; subtract outer means for covariance, not for independent errors.",
        "policy": "A shared ordered-prefix gate: both orientations rank>=1 and every R1 clock exact_pair; otherwise retain whole original pair. BothR2 identity. Exactly the old complete-global policy, no new clock solve.",
        "identity": "A=F1+F2-1; E=1-F1+F2. Shared joint-cell indicator cancels constants cellwise. Integral A=1-2C/(N+1); integral E=1-W/(N+1).",
    }
    for n in (325, 425):
        path = ROOT / "results/p334-full-birth-archive" / f"N{n}.csv"
        result["source_sha256"][str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
        with path.open() as stream:
            rows = [{k: int(v) for k, v in row.items()} for row in csv.DictReader(stream)]
        delta = contract["sizes"][str(n)]["delta_cos4"]
        batch_means = []
        second = np.zeros((9, 9, 9))
        counts = np.zeros(9, dtype=int)
        statuses = {"accepted_R1_pair": 0, "both_R2_identity": 0, "any_R0_original": 0,
                    "clock_fallback_original": 0}
        for batch in range(20):
            chosen = [r for r in rows if r["batch"] == batch]
            if len(chosen) != 1000:
                raise ValueError("incomplete original batch")
            blob_path = COND / "batches" / f"N{n}.batch{batch:02d}.json.gz"
            raw = blob_path.read_bytes()
            result["source_sha256"][str(blob_path.relative_to(ROOT))] = hashlib.sha256(raw).hexdigest()
            records = {r["counter"]: r for r in json.loads(gzip.decompress(raw))["records"]}
            k1 = np.array([[r["first_k1"], r["second_k1"]] for r in chosen])
            k2 = np.array([[r["first_k2"], r["second_k2"]] for r in chosen])
            ranks = np.array([[r["first_rank"], r["second_rank"]] for r in chosen])
            f1 = np.column_stack((binom.sf(k1-1, n, p), (n+1-k1)/(n+1)))
            f2 = np.column_stack((binom.sf(k2-1, n, p), (n+1-k2)/(n+1)))
            hybrid = f2.copy()
            for i, row in enumerate(chosen):
                rec = records[row["counter"]]
                if np.any(ranks[i] == 0):
                    statuses["any_R0_original"] += 1
                elif np.all(ranks[i] == 2):
                    statuses["both_R2_identity"] += 1
                elif rec["status"] == "exact_pair" and all(
                    rec["clocks"][o] is not None for o in range(2) if ranks[i, o] == 1
                ):
                    statuses["accepted_R1_pair"] += 1
                    for o in range(2):
                        if ranks[i, o] == 1:
                            hybrid[i, o], hybrid[i, 2+o] = rec["clocks"][o]["conditional"]
                else:
                    statuses["clock_fallback_original"] += 1
            effects = []
            for f2v in (f2, hybrid):
                for start in (0, 2):
                    a, e = f1+f2v-1, 1-f1+f2v
                    effects.extend([(a[:, start]-a[:, start+1])/delta,
                                    (e[:, start]-e[:, start+1])/delta])
            vector = np.column_stack((np.ones(len(chosen)), np.column_stack(effects)))
            cell = 3*ranks[:, 0] + ranks[:, 1]
            batch_vector = np.zeros((9, 9))
            for c in range(9):
                selected = vector[cell == c]
                counts[c] += len(selected)
                batch_vector[c] = selected.sum(axis=0)/1000
                second[c] += selected.T @ selected / 20000
            batch_means.append(batch_vector.ravel())
        batch_means = np.array(batch_means)
        totals = batch_means.reshape(20, 9, 9)[:, :, 1:].sum(axis=1)
        result["sizes"][str(n)] = {
            "delta_cos4": delta, "batch_ids": list(range(20)), "samples_per_batch": 1000,
            "cell_counts": counts.tolist(), "policy_counts": statuses,
            "joint_20_batch_means": batch_means.tolist(),
            "cell_individual_raw_second_moments": second.tolist(),
            "mean": batch_means.mean(axis=0).tolist(),
            "full_AE_total_batch_means": totals.tolist(),
            "full_AE_total_mean": totals.mean(axis=0).tolist(),
        }
        print(n, "counts", counts.tolist(), "total A/E", totals.mean(axis=0).tolist())
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "batch_vectors.json").write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()

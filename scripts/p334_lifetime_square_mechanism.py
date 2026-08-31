#!/usr/bin/env python3
"""Population-centered mean-lifetime / lifetime-spread decomposition, old paths only."""
import argparse
import csv
from hashlib import sha256
import io
import json
from pathlib import Path
import platform
import subprocess

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
SOURCE = "9c495ab13e65f2bc93dc0849ee3b73f88724c4b1"
NORMALIZATION_SOURCE = "3edc785a"
RAW_NAMES = ("W", "W2", "C", "C2", "CW")
NAMES = ("W_mean", "W_mean_squared", "W_variance", "W_second_moment",
         "C_mean", "C_variance", "CW_covariance", "mixture_Y_variance", "K1K2_covariance")
OUT = ROOT/"results/p334-lifetime-square-mechanism"


def descriptors(raw, n, delta):
    """Center at full population means (or full retained population in each LOO)."""
    raw = raw.reshape(2, 5)
    rows = []
    for w, w2, c, c2, cw in raw:
        var_w, var_c = w2-w*w, c2-c*c
        rows.append([w, w*w, var_w, w2, c, var_c, cw-c*w,
                     (var_c+w2/4)/(n+1)**2, var_c-var_w/4])
    rows = np.array(rows)
    difference = rows[0]-rows[1]
    h4 = difference/delta
    width_factor = -1/(4*(n+1)*(n+2))
    width = width_factor*h4[[1, 2, 3]]
    return np.r_[rows.ravel(), difference, h4, width]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUT)
    args = parser.parse_args()
    provenance = {}
    def read(commit, path):
        blob = subprocess.check_output(["git", "show", commit+":"+path], cwd=ROOT)
        provenance[commit+":"+path] = sha256(blob).hexdigest()
        return blob
    norm = json.loads(read(NORMALIZATION_SOURCE, "results/p334-full-global-conditional-clock/score.json"))
    labels = [prefix+"."+name for prefix in ("first", "second", "difference", "H4") for name in NAMES]
    labels += ["J1_width.H4.mean_lifetime_squared", "J1_width.H4.lifetime_variance", "J1_width.H4.total"]
    sizes = {}
    for n in (325, 425):
        path = f"results/p334-full-birth-archive/N{n}.csv"
        rows = list(csv.DictReader(io.StringIO(read(SOURCE, path).decode())))
        metadata = json.loads(read(SOURCE, path+".metadata.json"))
        delta = norm["sizes"][str(n)]["delta_cos4"]
        raw_batches = []
        for b in range(20):
            part = [row for row in rows if int(row["batch"]) == b]
            if len(part) != 1000:
                raise ValueError("Need each complete original 1000-counter batch")
            both = []
            for orientation in ("first", "second"):
                k1 = np.array([int(row[orientation+"_k1"]) for row in part])
                k2 = np.array([int(row[orientation+"_k2"]) for row in part])
                w, c = k2-k1, (k1+k2)/2
                both.extend([w.mean(), np.mean(w*w), c.mean(), np.mean(c*c), np.mean(c*w)])
            raw_batches.append(both)
        raw_batches = np.array(raw_batches)
        raw_mean = raw_batches.mean(axis=0)
        estimate = descriptors(raw_mean, n, delta)
        raw_loo = (20*raw_mean-raw_batches)/19
        loo = np.array([descriptors(r, n, delta) for r in raw_loo])
        centered = loo-loo.mean(axis=0)
        covariance = 19/20*centered.T@centered
        closure = {}
        for prefix in ("first", "second", "difference", "H4"):
            ix = [labels.index(prefix+"."+name) for name in ("W_mean_squared", "W_variance", "W_second_moment")]
            closure[prefix] = {"estimate_residual": float(estimate[ix[0]]+estimate[ix[1]]-estimate[ix[2]]),
                               "LOO_max_abs_residual": float(np.max(np.abs(loo[:, ix[0]]+loo[:, ix[1]]-loo[:, ix[2]])))}
        sizes[str(n)] = {"metadata": metadata, "delta_cos4": delta, "batch_ids": list(range(20)),
            "batch_denominators": [1000]*20, "raw_labels": [o+"."+v for o in ("first", "second") for v in RAW_NAMES],
            "raw_batch_means": raw_batches.tolist(), "full_population_raw_mean": raw_mean.tolist(),
            "raw_LOO_mean": raw_loo.tolist(), "labels": labels, "estimate": estimate.tolist(),
            "se": np.sqrt(np.diag(covariance)).tolist(), "joint_LOO_vectors": loo.tolist(),
            "joint_covariance": covariance.tolist(), "joint_rank_at_most": 19, "exact_add_back": closure}
    result = {"schema": "matching-one/p334-lifetime-square-population-decomposition/v1",
        "source_commit": SOURCE, "source_sha256": provenance, "normalization_source": NORMALIZATION_SOURCE,
        "python": platform.python_version(), "numpy": np.__version__, "sizes": sizes,
        "definitions": {"W": "K2-K1", "C": "(K1+K2)/2", "W2_split": "E W^2 = (E W)^2 + Var W",
            "J1_width": "-E W^2/[4(N+1)(N+2)]",
            "mixture_Y": "Y selects K1/(N+1) or K2/(N+1) with equal probability, then averages over original paths",
            "mixture_Y_variance": "[Var C + E W^2/4]/(N+1)^2",
            "K1K2_covariance": "Var C-Var W/4",
            "full_even_integral": "integral E_top=1-E W/(N+1); H4 contrast=-H4(E W)/(N+1)"},
        "centering": "Full 20000-path population raw means per orientation; each delete-one-batch replicate recenters the retained 19000 paths. No within-batch centering or average of batch variances.",
        "boundary": "Post-reveal decomposition of the previously observed W^2 direction, not an independent test or new mechanism selection. Empirical plug-in moments with ddof=0 keep the exact mean-square/variance identity; twenty original paired batches determine the common uncertainty. No covariance inverse, new MC, DP, replay or resampled paths.",
        "new_MC": 0, "new_DP": 0, "new_birth_replay": 0}
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output/"score.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    lines = ["# P334 lifetime-square mechanism", ""]
    for n, r in sizes.items():
        lines += [f"## N{n}", "", "| Full-population readout | Estimate | Paired-batch SE |", "|---|---:|---:|"]
        for name, value, se in zip(r["labels"], r["estimate"], r["se"]):
            lines.append(f"| {name} | {value:.11g} | {se:.7g} |")
        lines.append("")
    lines += [result["centering"], "", result["boundary"], ""]
    (args.output/"REPORT.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Final-commit scorer for fixed eight-quartet nested next-label production."""
import argparse
import gzip
from hashlib import sha256
import io
import json
from pathlib import Path
import re
import subprocess

import numpy as np
from scipy.stats import binom

ROOT = Path(__file__).resolve().parents[1]
OLD_COMMIT = "bb79fd47"
OLD_PATH = "results/p334-nine-layer-complete-ae/batch_vectors.json"
OUT = ROOT/"results/p334-next-label-doob-quartets"
HEADER = ["N", "batch", "counter", "k0", "first_rank", "second_rank", "quartet", "group", "replica", "next_label",
          "first_next_rank", "second_next_rank", "first_k1", "first_k2", "second_k1", "second_k2"]
LABELS = [ep+"."+obs for ep in ("p_ref", "p_integral") for obs in ("F1", "F2", "A", "E")]+["K1", "K2", "W"]
MATRIX_NAMES = ("Vtot", "Dnext", "Vafter")
CELLS = [str(i)+str(j) for i in range(3) for j in range(3)]
GROUPS = {"00": [0], "11": [4], "22": [8], "01+10": [1, 3], "02+20": [2, 6], "12+21": [5, 7]}
TRI = np.triu_indices(len(LABELS))


def paired_observer(k1f, k2f, k1s, k2s, n, delta, tail):
    k1, k2 = (k1f-k1s)/delta, (k2f-k2s)/delta
    f1, f2 = (tail[k1f]-tail[k1s])/delta, (tail[k2f]-tail[k2s])/delta
    g1, g2 = -k1/(n+1), -k2/(n+1)
    return np.stack((f1, f2, f1+f2, -f1+f2, g1, g2, g1+g2, -g1+g2, k1, k2, k2-k1), axis=-1)


def old_observer(ae, n, offset):
    a, e, ai, ei = (ae[:, offset+j] for j in range(4))
    f1, f2, g1, g2 = (a-e)/2, (a+e)/2, (ai-ei)/2, (ai+ei)/2
    k1, k2 = -(n+1)*g1, -(n+1)*g2
    return np.column_stack((f1, f2, a, e, g1, g2, ai, ei, k1, k2, k2-k1))


def score_batch(blob, n, batch, delta, tail):
    text = gzip.decompress(blob).decode()
    if text.splitlines()[0].split(",") != HEADER:
        raise ValueError("Quartet CSV schema differs from the frozen 16 fields")
    rows = np.loadtxt(io.StringIO(text), delimiter=",", skiprows=1, dtype=np.int64)
    if rows.shape != (32000, len(HEADER)):
        raise ValueError("Need 1000 complete prefixes times eight quartets times four tails")
    rows = rows[np.lexsort((rows[:, 8], rows[:, 7], rows[:, 6], rows[:, 2]))]
    x = rows.reshape(1000, 8, 2, 2, len(HEADER))
    counters = x[:, 0, 0, 0, 2]
    if len(np.unique(counters)) != 1000 or not np.all(x[..., 2] == counters[:, None, None, None]):
        raise ValueError("Prefix counter blocks are not complete")
    if not np.all(x[..., 0] == n) or not np.all(x[..., 1] == batch):
        raise ValueError("Wrong N or original batch")
    if (not np.all(x[..., 6] == np.arange(8)[None, :, None, None])
            or not np.all(x[..., 7] == np.arange(2)[None, None, :, None])
            or not np.all(x[..., 8] == np.arange(2)[None, None, None, :])):
        raise ValueError("Missing or duplicated quartet/group/replica rows")
    for column in (3, 4, 5):
        if not np.all(x[..., column] == x[:, 0, 0, 0, column][:, None, None, None]):
            raise ValueError("Checkpoint state changes within a fixed prefix")
    for column in (9, 10, 11):
        if not np.array_equal(x[:, :, :, 0, column], x[:, :, :, 1, column]):
            raise ValueError("The two suffix replicas do not share their next-label checkpoint")
    rank = x[:, 0, 0, 0, 4:6]
    if np.any((rank < 0) | (rank > 2)):
        raise ValueError("Invalid joint-rank cell")
    cell = 3*rank[:, 0]+rank[:, 1]
    obs = paired_observer(x[..., 12], x[..., 13], x[..., 14], x[..., 15], n, delta, tail)
    a, b = obs[:, :, 0, 0]-obs[:, :, 1, 0], obs[:, :, 0, 1]-obs[:, :, 1, 1]
    means, matrices, mass = [], [], []
    for k in range(9):
        mask = cell == k
        av, bv = a[mask].reshape(-1, 11), b[mask].reshape(-1, 11)
        total = (av.T@av+bv.T@bv)/(4*8000)
        next_cov = (av.T@bv+bv.T@av)/(4*8000)
        after = ((av-bv).T@(av-bv))/(4*8000)
        matrices.append([total, next_cov, after])
        means.append(obs[mask].sum(axis=(0, 1, 2, 3))/32000)
        mass.append(mask.mean())
    return np.array(mass), np.array(means), np.array(matrices), counters


def block_summary(mass, means, matrices):
    def summary(b):
        return {"mean": b.mean(axis=0).tolist(), "se": (b.std(axis=0, ddof=1)/np.sqrt(20)).tolist()}
    return {"mass": summary(mass), "observer_contribution": summary(means),
            "conditional_covariance_contribution": {name: summary(matrices[:, i]) for i, name in enumerate(MATRIX_NAMES)},
            "new_32_tail_mean_conditional_covariance": summary(matrices[:, 1]/16+matrices[:, 2]/32)}


def flatten_raw(mass, means, matrices, old_base, old_safe):
    labels, arrays = [], []
    for c, cell in enumerate(CELLS):
        labels.append("cell."+cell+".mass")
        arrays.append(mass[:, c:c+1])
        labels += ["cell."+cell+".mean."+name for name in LABELS]
        arrays.append(means[:, c])
        for k, name in enumerate(MATRIX_NAMES):
            labels += [f"cell.{cell}.{name}[{LABELS[i]},{LABELS[j]}]" for i, j in zip(*TRI)]
            arrays.append(matrices[:, c, k][:, TRI[0], TRI[1]])
    for name, values in (("old_baseline", old_base), ("old_safe", old_safe)):
        labels += [name+"."+label for label in LABELS]
        arrays.append(values)
    return labels, np.column_stack(arrays)


def primary_readout(matrices, new_mean, old_base, old_safe):
    """Only named fraction targets; retain negative next-label estimates verbatim."""
    global_m = matrices.sum(axis=0)
    diag = np.diagonal(global_m, axis1=-2, axis2=-1)
    labels, values = [], []
    for name, x in (("new_32_tail", new_mean), ("old_baseline", old_base), ("old_safe", old_safe),
                    ("new_minus_old_baseline", new_mean-old_base), ("new_minus_old_safe", new_mean-old_safe)):
        labels += [name+"."+label for label in LABELS]
        values.extend(x)
    # All-global covariance entries and their common uncertainty, not a
    # purported covariance matrix of next-label eigenstates.
    for k, name in enumerate(MATRIX_NAMES):
        labels += [f"all.{name}[{LABELS[i]},{LABELS[j]}]" for i, j in zip(*TRI)]
        values.extend(global_m[k][TRI])
    for name, numerator in (("next_fraction", diag[1]), ("after_fraction", diag[2]),
                            ("32tail_remaining_fraction", diag[1]/16+diag[2]/32)):
        labels += ["all."+name+"."+label for label in LABELS]
        values.extend(np.divide(numerator, diag[0], out=np.full(11, np.nan), where=diag[0] != 0))
    group_m = matrices[[1, 3]].sum(axis=0)
    gd = np.diagonal(group_m, axis1=-2, axis2=-1)
    for coordinate in ("K1", "p_ref.E", "p_integral.E"):
        j = LABELS.index(coordinate)
        for name, numerator, denominator in (("total_share_of_all_suffix", gd[0, j], diag[0, j]),
                ("next_share_of_all_suffix", gd[1, j], diag[0, j]),
                ("after_share_of_all_suffix", gd[2, j], diag[0, j]),
                ("next_fraction_within_group", gd[1, j], gd[0, j])):
            labels.append("01+10."+name+"."+coordinate)
            values.append(numerator/denominator if denominator != 0 else np.nan)
    return labels, np.array(values)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-commit", required=True, help="Final completed production commit; no working directory input")
    parser.add_argument("--source-directory", required=True)
    parser.add_argument("--output", type=Path, default=OUT)
    args = parser.parse_args()
    source = subprocess.check_output(["git", "rev-parse", args.source_commit+"^{commit}"], cwd=ROOT, text=True).strip()
    paths = subprocess.check_output(["git", "ls-tree", "-r", "--name-only", source, "--", args.source_directory], cwd=ROOT, text=True).splitlines()
    files = {}
    for path in paths:
        match = re.fullmatch(r"N(325|425)\.batch(\d{2})\.csv\.gz", Path(path).name)
        if match:
            key = tuple(map(int, match.groups()))
            if key in files:
                raise ValueError("Duplicate production batch file")
            files[key] = path
    if set(files) != {(n, b) for n in (325, 425) for b in range(20)}:
        raise ValueError("Final source must contain all forty complete batch gzip files")
    hashes = {}
    def read(commit, path):
        blob = subprocess.check_output(["git", "show", commit+":"+path], cwd=ROOT)
        hashes[commit+":"+path] = sha256(blob).hexdigest()
        return blob
    old = json.loads(read(OLD_COMMIT, OLD_PATH))
    p = old["p_ref"]
    args.output.mkdir(parents=True, exist_ok=True)
    sizes = {}
    for n in (325, 425):
        old_n = old["sizes"][str(n)]
        delta = old_n["delta_cos4"]
        tail = binom.sf(np.arange(n+1)-1, n, p)
        masses, means, matrices, counters = [], [], [], []
        for b in range(20):
            mass, mean, matrix, ids = score_batch(read(source, files[n, b]), n, b, delta, tail)
            masses.append(mass); means.append(mean); matrices.append(matrix); counters.extend(ids.tolist())
        if len(set(counters)) != 20000:
            raise ValueError("Production does not contain 20000 unique original prefixes")
        masses, means, matrices = np.array(masses), np.array(means), np.array(matrices)
        old_ae = np.array(old_n["full_AE_total_batch_means"])
        old_base, old_safe = old_observer(old_ae, n, 0), old_observer(old_ae, n, 4)
        raw_labels, raw_batches = flatten_raw(masses, means, matrices, old_base, old_safe)
        raw_mean = raw_batches.mean(axis=0)
        raw_loo = (20*raw_mean-raw_batches)/19
        labels, point = primary_readout(matrices.mean(axis=0), means.mean(axis=0).sum(axis=0), old_base.mean(axis=0), old_safe.mean(axis=0))
        loo = np.array([primary_readout((matrices.sum(axis=0)-matrices[b])/19,
                    (means.sum(axis=(0, 1))-means[b].sum(axis=0))/19,
                    (old_base.sum(axis=0)-old_base[b])/19, (old_safe.sum(axis=0)-old_safe[b])/19)[1] for b in range(20)])
        identified = np.isfinite(point) & np.isfinite(loo).all(axis=0)
        joint = np.column_stack((raw_loo, loo[:, identified]))
        factor = np.sqrt(19/20)*(joint-joint.mean(axis=0))
        se = np.full(len(point), np.nan)
        se[identified] = np.sqrt(np.sum(factor[:, len(raw_labels):]**2, axis=0))
        covariance_path = f"N{n}.joint_covariance_factor.json.gz"
        packed = {"labels": raw_labels+[label for label, keep in zip(labels, identified) if keep],
                  "covariance_convention": "Full common covariance equals factor.T @ factor; rows are aligned original delete-one-batch replicates; never inverted",
                  "factor": factor.tolist(), "raw_batch_labels": raw_labels,
                  "raw_joint_20_batch_means": raw_batches.tolist(), "primary_LOO_identified": loo[:, identified].tolist(),
                  "batch_ids": list(range(20)), "rank_at_most": 19}
        compressed = gzip.compress((json.dumps(packed, separators=(",", ":"), allow_nan=False)+"\n").encode(), mtime=0)
        (args.output/covariance_path).write_bytes(compressed)
        blocks = {"all": block_summary(masses.sum(axis=1), means.sum(axis=1), matrices.sum(axis=1))}
        for cell, k in zip(CELLS, range(9)):
            blocks["cell."+cell] = block_summary(masses[:, k], means[:, k], matrices[:, k])
        for name, ix in GROUPS.items():
            blocks["group."+name] = block_summary(masses[:, ix].sum(axis=1), means[:, ix].sum(axis=1), matrices[:, ix].sum(axis=1))
        first55 = factor[:, len(raw_labels):len(raw_labels)+55]
        sizes[str(n)] = {"delta_cos4": delta, "samples": 20000, "new_tails": 640000, "quartets_per_prefix": 8,
            "batch_ids": list(range(20)), "prefixes_per_batch": 1000, "blocks": blocks,
            "primary_labels": labels, "primary_estimate": [float(v) if ok else None for v, ok in zip(point, identified)],
            "primary_se": [float(v) if ok else None for v, ok in zip(se, identified)],
            "unidentified_primary": [label for label, keep in zip(labels, identified) if not keep],
            "new_old_mean_joint_covariance": (first55.T@first55).tolist(),
            "covariance_factor_file": covariance_path, "covariance_factor_sha256": sha256(compressed).hexdigest(),
            "total_next_after_max_batch_closure": float(np.max(np.abs(matrices[:, :, 0]-matrices[:, :, 1]-matrices[:, :, 2]))),
            "Dnext_clipping": "none"}
    result = {"schema": "matching-one/p334-next-label-nested-quartet/v1", "source_commit": source,
        "source_directory": args.source_directory, "source_sha256": hashes, "old_complete_AE_commit": OLD_COMMIT,
        "p_ref": p, "observer_labels": LABELS, "sizes": sizes,
        "matrix_targets": {"Vtot": "E_Z Cov(X|Z)", "Dnext": "E_Z Cov_U(E[X|Z,U]|Z)",
            "Vafter": "E_Z E_U Cov(X|Z,U)", "new32_mean_noise": "Dnext/16+Vafter/32; 16 independent next labels with two suffix replicas each"},
        "boundary": "Fresh suffix production on the original prefix population, not old-tail replay or independent new prefixes. All errors use twenty original paired batches. Cell/group covariance terms use whole-population denominators, not conditional cell denominators. Dnext remains signed and may be non-PSD at finite samples; no clipping, covariance inversion, new DP or simulation by the scorer."}
    (args.output/"score.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    lines = ["# Next-label versus after-next suffix covariance", ""]
    for n, r in sizes.items():
        lines += [f"## N{n}", "", "| Named readout | Estimate | Shared-batch SE |", "|---|---:|---:|"]
        for name, value, error in zip(r["primary_labels"], r["primary_estimate"], r["primary_se"]):
            if (name.startswith("01+10.") or name.startswith("all.next_fraction.")
                    or (name.startswith(("new_32_tail.", "new_minus_old_baseline.", "new_minus_old_safe.")) and name.split(".", 1)[1] in ("K1", "K2", "W", "p_integral.A", "p_integral.E"))):
                cells = (f"{value:.9g}", f"{error:.5g}") if value is not None else ("not_scoreable", "not_scoreable")
                lines.append(f"| {name} | {cells[0]} | {cells[1]} |")
        lines.append("")
    lines += [result["boundary"], ""]
    (args.output/"REPORT.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()

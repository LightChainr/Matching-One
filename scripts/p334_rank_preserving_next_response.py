#!/usr/bin/env python3
"""Split the saved nested-fork next-label response by rank preservation."""
import gzip
from hashlib import sha256
import io
import json
from pathlib import Path
import subprocess

import numpy as np
from scipy.stats import binom

from p334_next_label_doob_quartets import HEADER, LABELS as PAIRED_LABELS, CELLS, paired_observer

ROOT = Path(__file__).resolve().parents[1]
SOURCE = "e32a85939279b8574278024d647b56d2d1485247"
SOURCE_DIR = "results/p334-nested-next-label-forks"
CORE = "24872eef"
CORE_DIR = "results/p334-next-label-doob-quartets"
OLD_JOINT = "e0494fdf"
OLD_DIR = "results/p334-next-label-mechanism-joint"
OUT = ROOT/"results/p334-rank-preserving-next-response"
MASKS = ("neither", "mixed", "both")
LABELS = PAIRED_LABELS+[ep+"."+name for ep in ("p_ref", "p_integral") for name in ("F1f_raw", "F2f_raw", "F1s_raw", "F2s_raw")]
TRI = np.triu_indices(19)


def batch_response(blob, n, delta, tail):
    text = gzip.decompress(blob).decode()
    if text.splitlines()[0].split(",") != HEADER:
        raise ValueError("Wrong fixed quartet CSV schema")
    rows = np.loadtxt(io.StringIO(text), delimiter=",", skiprows=1, dtype=np.int64)
    rows = rows[np.lexsort((rows[:, 8], rows[:, 7], rows[:, 6], rows[:, 2]))]
    x = rows.reshape(1000, 8, 2, 2, 16)
    ranks = x[:, 0, 0, 0, 4:6]
    cell = 3*ranks[:, 0]+ranks[:, 1]
    changed = np.any(x[:, :, :, 0, 10:12] != ranks[:, None, None, :], axis=-1)
    mask = changed.sum(axis=2)  # 0 neither, 1 mixed, 2 both, once per quartet.
    k1f, k2f, k1s, k2s = (x[..., j] for j in (12, 13, 14, 15))
    paired = paired_observer(k1f, k2f, k1s, k2s, n, delta, tail)
    raw = np.stack([tail[k] for k in (k1f, k2f, k1s, k2s)]+[(n+1-k)/(n+1) for k in (k1f, k2f, k1s, k2s)], axis=-1)
    obs = np.concatenate((paired, raw), axis=-1)
    a, b = obs[:, :, 0, 0]-obs[:, :, 1, 0], obs[:, :, 0, 1]-obs[:, :, 1, 1]
    matrices = np.zeros((9, 3, 19, 19))
    mass = np.zeros((9, 3))
    for c in range(9):
        for m in range(3):
            selected = (cell[:, None] == c) & (mask == m)
            av, bv = a[selected], b[selected]
            matrices[c, m] = (av.T@bv+bv.T@av)/(4*8000)
            mass[c, m] = selected.sum()/8000
    return mass, matrices


def gamma_coordinates(matrices, delta):
    labels, values = [], []
    for name, ix in (("all", list(range(9))), ("01+10", [1, 3])):
        masks = matrices[ix].sum(axis=0)
        for mask, matrix in [("all", masks.sum(axis=0))]+list(zip(MASKS, masks)):
            for ep, start in (("p_ref", 11), ("p_integral", 15)):
                first, second = matrix[start, start+1], matrix[start+2, start+3]
                cross1, cross2 = matrix[start, start+3], matrix[start+1, start+2]
                within, cross = (first+second)/delta**2, -(cross1+cross2)/delta**2
                names = ("Gamma_pair", "within_first_raw", "within_second_raw", "crn_F1f_F2s_raw", "crn_F2f_F1s_raw", "within_sum_H4", "crn_sum_H4")
                labels += [f"{name}.{mask}.{ep}.{field}" for field in names]
                values.extend([within+cross, first, second, cross1, cross2, within, cross])
        for ep, j in (("p_ref", 0), ("p_integral", 4)):
            total = masks[:, j, j+1].sum()
            for mask, matrix in zip(MASKS, masks):
                labels.append(f"{name}.{mask}.{ep}.signed_Gamma_share")
                values.append(matrix[j, j+1]/total if total != 0 else np.nan)
    return labels, np.array(values)


def block_summary(mass, matrix):
    blocks = {}
    for name, ix in [("all", list(range(9))), ("01+10", [1, 3])]+[("cell."+c, [k]) for k, c in enumerate(CELLS)]:
        m, b = mass[:, ix].sum(axis=1), matrix[:, ix].sum(axis=1)
        blocks[name] = {}
        for label, masses, matrices in [("all", m.sum(axis=1), b.sum(axis=1))]+[(label, m[:, k], b[:, k]) for k, label in enumerate(MASKS)]:
            blocks[name][label] = {"quartet_mass": float(masses.mean()), "quartet_mass_se": float(masses.std(ddof=1)/np.sqrt(20)),
                                  "B": matrices.mean(axis=0).tolist(), "B_se": (matrices.std(axis=0, ddof=1)/np.sqrt(20)).tolist()}
    return blocks


def main():
    hashes = {}
    def read(commit, path):
        blob = subprocess.check_output(["git", "show", commit+":"+path], cwd=ROOT)
        hashes[commit+":"+path] = sha256(blob).hexdigest()
        return blob
    core = json.loads(read(CORE, CORE_DIR+"/score.json"))
    old = json.loads(read(OLD_JOINT, OLD_DIR+"/score.json"))
    OUT.mkdir(parents=True, exist_ok=True)
    sizes = {}
    for n in (325, 425):
        c, delta = core["sizes"][str(n)], core["sizes"][str(n)]["delta_cos4"]
        tail = binom.sf(np.arange(n+1)-1, n, core["p_ref"])
        masses, matrices = [], []
        for batch in range(20):
            blob = read(SOURCE, f"{SOURCE_DIR}/N{n}/N{n}.batch{batch:02d}.csv.gz")
            mass, matrix = batch_response(blob, n, delta, tail)
            masses.append(mass); matrices.append(matrix)
        masses, matrices = np.array(masses), np.array(matrices)
        mean = matrices.mean(axis=0)
        loo_b = (20*mean-matrices)/19
        labels, point = gamma_coordinates(mean, delta)
        loo = np.array([gamma_coordinates(b, delta)[1] for b in loo_b])
        identified = np.isfinite(point) & np.isfinite(loo).all(axis=0)
        raw_labels, raw_cols = [], []
        for k, cell in enumerate(CELLS):
            for m, mask in enumerate(MASKS):
                prefix = f"rank_preserving.cell.{cell}.{mask}"
                raw_labels.append(prefix+".quartet_mass")
                raw_cols.append(masses[:, k, m:m+1])
                raw_labels += [f"{prefix}.B[{LABELS[i]},{LABELS[j]}]" for i, j in zip(*TRI)]
                raw_cols.append(matrices[:, k, m][:, TRI[0], TRI[1]])
        raw_batches = np.column_stack(raw_cols)
        raw_loo = (20*raw_batches.mean(axis=0)-raw_batches)/19
        added_loo = np.column_stack((raw_loo, loo[:, identified]))
        added_factor = np.sqrt(19/20)*(added_loo-added_loo.mean(axis=0))
        old_info = old["sizes"][str(n)]
        old_blob = read(OLD_JOINT, OLD_DIR+"/"+old_info["complete_covariance_factor_file"])
        saved = json.loads(gzip.decompress(old_blob))
        factor = np.column_stack((np.array(saved["factor"]), added_factor))
        factor_labels = saved["labels"]+raw_labels+["rank_preserving."+label for label, keep in zip(labels, identified) if keep]
        filename = f"N{n}.complete_common_factor.json.gz"
        packed = {"labels": factor_labels, "factor": factor.tolist(), "covariance": "factor.T @ factor; same twenty original LOO rows; no inverse",
                  "new_raw_labels": raw_labels, "new_raw_20_batch_means": raw_batches.tolist(),
                  "gamma_labels_identified": [label for label, keep in zip(labels, identified) if keep],
                  "gamma_LOO_identified": loo[:, identified].tolist(), "batch_ids": list(range(20)), "rank_at_most": 19}
        blob = gzip.compress((json.dumps(packed, separators=(",", ":"), allow_nan=False)+"\n").encode(), mtime=0)
        (OUT/filename).write_bytes(blob)
        gamma_factor = added_factor[:, len(raw_labels):]
        gamma_cov = gamma_factor.T@gamma_factor
        se = np.full(len(point), np.nan); se[identified] = np.sqrt(np.diag(gamma_cov))
        closure = 0.
        for ep, j, start in (("p_ref", 0, 11), ("p_integral", 4, 15)):
            implied = (matrices[..., start, start+1]+matrices[..., start+2, start+3]
                       -matrices[..., start, start+3]-matrices[..., start+1, start+2])/delta**2
            closure = max(closure, float(np.max(np.abs(implied-matrices[..., j, j+1]))))
        sizes[str(n)] = {"delta_cos4": delta, "batch_ids": list(range(20)), "prefixes_per_batch": 1000,
            "blocks": block_summary(masses, matrices), "gamma_labels": labels,
            "gamma_estimate": [float(v) if ok else None for v, ok in zip(point, identified)],
            "gamma_se": [float(v) if ok else None for v, ok in zip(se, identified)],
            "gamma_covariance_identified": gamma_cov.tolist(), "unidentified_labels": [label for label, keep in zip(labels, identified) if not keep],
            "complete_covariance_factor_file": filename, "complete_covariance_factor_sha256": sha256(blob).hexdigest(),
            "raw_vs_paired_Gamma_max_batch_mask_cell_residual": closure,
            "all_masks_original_Gamma_difference": {ep: float(mean[:, :, j, j+1].sum()-c["primary_estimate"][c["primary_labels"].index("all.next_first_completion_Gamma."+ep)]) for ep, j in (("p_ref", 0), ("p_integral", 4))}}
    result = {"schema": "matching-one/p334-rank-preserving-next-response/v1", "source_commit": SOURCE,
        "source_sha256": hashes, "previous_complete_joint": OLD_JOINT, "observer_labels": LABELS, "masks": list(MASKS),
        "sizes": sizes, "new_MC": 0, "new_DP": 0, "new_path_replay": 0,
        "estimand": "For each prefix Z, neither/safe-safe contributes pi_safe(Z)^2 Cov_U(m(U)|Z,safe), averaged over prefixes; no division by a pooled safety rate.",
        "Gamma_identity": "(B_F1f,F2f+B_F1s,F2s-B_F1f,F2s-B_F2f,F1s)/delta_cos4^2",
        "boundary": "Three masks partition original quartets by next-rank changes in either orientation; raw orientation blocks use this same paired mask. Signed Gamma allocations are not mean fractions or complete within-safe variance fractions. All original cells/prefix denominators and cross-provider covariance remain common; no PSD clipping, test suite, covariance inversion or additional production."}
    (OUT/"score.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    lines = ["# Rank-preserving next-label first/completion response", ""]
    for n, r in sizes.items():
        lines += [f"## N{n}", "", "| Sector / mask | canonical Gamma +/- SE | integrated Gamma +/- SE |", "|---|---:|---:|"]
        for group in ("all", "01+10"):
            for mask in ("all",)+MASKS:
                columns = []
                for ep in ("p_ref", "p_integral"):
                    ix = r["gamma_labels"].index(f"{group}.{mask}.{ep}.Gamma_pair")
                    columns.append(f"{r['gamma_estimate'][ix]:.10g} +/- {r['gamma_se'][ix]:.6g}")
                lines.append("| "+group+" / "+mask+" | "+" | ".join(columns)+" |")
        lines.append("")
    lines += [result["estimand"], "", result["boundary"], ""]
    (OUT/"REPORT.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()

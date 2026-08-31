#!/usr/bin/env python3
"""Use an exact prefix census to separate within-type and between-type response."""
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

from p334_rank_preserving_next_response import SOURCE, SOURCE_DIR, CORE, CORE_DIR, LABELS, TRI
from p334_next_label_doob_quartets import HEADER, CELLS, paired_observer

ROOT = Path(__file__).resolve().parents[1]
OLD = "30c7ddb0"
OLD_DIR = "results/p334-rank-preserving-next-response"
OUT = ROOT/"results/p334-exact-birth-partition-response"
PARTS = ("total", "within_safe", "within_birth", "within_total", "between")
FOUR_PARTS = ("total", "within_00", "within_01", "within_10", "within_11", "within_total", "between")
READOUTS = ("Gamma", "A_variance", "E_variance")


def weighted_batch(blob, census, n, delta, tail):
    text = gzip.decompress(blob).decode()
    if text.splitlines()[0].split(",") != HEADER:
        raise ValueError("Unexpected nested-fork schema")
    a = np.loadtxt(io.StringIO(text), delimiter=",", skiprows=1, dtype=np.int64)
    a = a[np.lexsort((a[:, 8], a[:, 7], a[:, 6], a[:, 2]))].reshape(1000, 8, 2, 2, 16)
    first = a[:, 0, 0, 0]
    cr = [census[int(row[2])] for row in first]
    safe = np.array([row["joint_safe_count"]/row["d"] for row in cr])
    for row, ref in zip(first, cr):
        if (ref["N"] != n or ref["batch"] != row[1] or ref["k0"] != row[3]
                or ref["d"] != n-row[3] or ref["first_oldrank"] != row[4] or ref["second_oldrank"] != row[5]):
            raise ValueError("Census and original fork prefix disagree")
    if np.any((safe < 0) | (safe > 1)):
        raise ValueError("Invalid exact census probability")
    rank = first[:, 4:6]
    cell = 3*rank[:, 0]+rank[:, 1]
    changed_bits = a[:, :, :, 0, 10:12] != rank[:, None, None, :]
    change_type = 2*changed_bits[..., 0]+changed_bits[..., 1]
    changed = np.any(changed_bits, axis=-1)
    mask = changed.sum(axis=2)
    if np.any((mask == 0) & (safe[:, None] == 0)) or np.any((mask == 2) & (safe[:, None] == 1)):
        raise ValueError("A sampled label contradicts an exact empty census class")
    ws = np.divide((mask == 0).astype(float), safe[:, None], out=np.zeros_like(mask, dtype=float), where=safe[:, None] > 0)
    wb = np.divide((mask == 2).astype(float), 1-safe[:, None], out=np.zeros_like(mask, dtype=float), where=safe[:, None] < 1)
    k1f, k2f, k1s, k2s = (a[..., j] for j in (12, 13, 14, 15))
    paired = paired_observer(k1f, k2f, k1s, k2s, n, delta, tail)
    raw = np.stack([tail[k] for k in (k1f, k2f, k1s, k2s)]+[(n+1-k)/(n+1) for k in (k1f, k2f, k1s, k2s)], axis=-1)
    x = np.concatenate((paired, raw), axis=-1)
    da, db = x[:, :, 0, 0]-x[:, :, 1, 0], x[:, :, 0, 1]-x[:, :, 1, 1]
    # Only the named Gamma/A/E readouts are retained for the four observed types.
    probabilities = np.array([[r["joint_safe_count"], r["first_safe_count"]-r["joint_safe_count"],
        r["second_safe_count"]-r["joint_safe_count"], r["d"]-r["first_safe_count"]-r["second_safe_count"]+r["joint_safe_count"]]
        for r in cr], dtype=float)/np.array([r["d"] for r in cr])[:, None]
    if np.any(probabilities < 0): raise ValueError("Invalid four-type census")
    sample = np.stack([np.stack(((da[..., j]*db[..., j+1]+da[..., j+1]*db[..., j])/4,
        da[..., j+2]*db[..., j+2]/2, da[..., j+3]*db[..., j+3]/2), axis=-1) for j in (0, 4)], axis=-2)
    four = np.zeros((2, len(FOUR_PARTS), 2, 3))
    four_weights = []
    for c in range(4):
        same = (change_type[:, :, 0] == c) & (change_type[:, :, 1] == c)
        if np.any(same & (probabilities[:, c, None] == 0)):
            raise ValueError("Sampled four-type class has zero exact probability")
        four_weights.append(np.divide(same, probabilities[:, c, None], out=np.zeros_like(same, dtype=float), where=probabilities[:, c, None] > 0))
    for g, select in enumerate((np.ones(1000, dtype=bool), np.isin(cell, [1, 3]))):
        four[g, 0] = sample[select].sum(axis=(0, 1))/8000
        for c, w in enumerate(four_weights):
            four[g, c+1] = (sample[select]*w[select, :, None, None]).sum(axis=(0, 1))/8000
        four[g, 5] = four[g, 1:5].sum(axis=0)
        four[g, 6] = four[g, 0]-four[g, 5]
    mats = np.zeros((9, 5, 19, 19))
    for c in range(9):
        ix = cell == c
        av, bv = da[ix].reshape(-1, 19), db[ix].reshape(-1, 19)
        for k, w in enumerate((np.ones_like(ws), ws, wb)):
            bv_w = bv*w[ix].reshape(-1, 1)
            mats[c, k] = (av.T@bv_w+bv_w.T@av)/(4*8000)
        mats[c, 3] = mats[c, 1]+mats[c, 2]
        mats[c, 4] = mats[c, 0]-mats[c, 3]
    return mats, four, {"safe_zero": int(np.sum(safe == 0)), "safe_one": int(np.sum(safe == 1)),
                  "nontrivial": int(np.sum((safe > 0) & (safe < 1))),
                  "largest_inverse_probability": float(max(ws.max(), wb.max())),
                  "four_type_empty_classes": (probabilities == 0).sum(axis=0).tolist(),
                  "largest_four_type_inverse_probability": float(max(w.max() for w in four_weights))}


def gamma_readout(mats, delta):
    labels, values = [], []
    for group, cells in (("all", list(range(9))), ("01+10", [1, 3])):
        total = mats[cells].sum(axis=0)
        for k, part in enumerate(PARTS):
            for ep, start, j in (("p_ref", 11, 0), ("p_integral", 15, 4)):
                b = total[k]
                within = (b[start, start+1]+b[start+2, start+3])/delta**2
                cross = -(b[start, start+3]+b[start+1, start+2])/delta**2
                for name, v in (("Gamma_pair", within+cross), ("within_orientation", within), ("cross_orientation", cross)):
                    labels.append(f"{group}.{part}.{ep}.{name}"); values.append(v)
                if part != "total":
                    den = total[0, j, j+1]
                    labels.append(f"{group}.{part}.{ep}.signed_Gamma_share")
                    values.append(b[j, j+1]/den if den != 0 else np.nan)
    return labels, np.array(values)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--census-commit", required=True)
    parser.add_argument("--census-template", required=True, help="Committed CSV path containing {N}; optionally gzip")
    args = parser.parse_args()
    census_commit = subprocess.check_output(["git", "rev-parse", args.census_commit+"^{commit}"], cwd=ROOT, text=True).strip()
    hashes = {}
    def read(commit, path):
        blob = subprocess.check_output(["git", "show", commit+":"+path], cwd=ROOT)
        hashes[commit+":"+path] = sha256(blob).hexdigest(); return blob
    core = json.loads(read(CORE, CORE_DIR+"/score.json"))
    old = json.loads(read(OLD, OLD_DIR+"/score.json"))
    OUT.mkdir(parents=True, exist_ok=True)
    sizes = {}
    for n in (325, 425):
        path = args.census_template.format(N=n)
        blob = read(census_commit, path)
        if path.endswith(".gz"): blob = gzip.decompress(blob)
        census_rows = [{k: int(v) for k, v in row.items()} for row in csv.DictReader(io.StringIO(blob.decode()))]
        census = {row["counter"]: row for row in census_rows}
        if len(census) != 20000: raise ValueError("Need all 20000 original prefix census rows")
        delta = core["sizes"][str(n)]["delta_cos4"]
        tail = binom.sf(np.arange(n+1)-1, n, core["p_ref"])
        batches, four_batches, counts = [], [], []
        for batch in range(20):
            matrix, four, count = weighted_batch(read(SOURCE, f"{SOURCE_DIR}/N{n}/N{n}.batch{batch:02d}.csv.gz"), census, n, delta, tail)
            batches.append(matrix); four_batches.append(four); counts.append(count)
        batches = np.array(batches); mean = batches.mean(axis=0)
        b_loo = (20*mean-batches)/19
        labels, point = gamma_readout(mean, delta)
        loo = np.array([gamma_readout(b, delta)[1] for b in b_loo])
        identified = np.isfinite(point) & np.isfinite(loo).all(axis=0)
        raw_labels = [f"birth_partition.cell.{cell}.{part}.B[{LABELS[i]},{LABELS[j]}]" for cell in CELLS for part in PARTS for i, j in zip(*TRI)]
        raw_batches = batches[..., TRI[0], TRI[1]].reshape(20, -1)
        raw_loo = (20*raw_batches.mean(axis=0)-raw_batches)/19
        four_labels = [f"{g}.{part}.{ep}.{readout}" for g in ("all", "01+10") for part in FOUR_PARTS
                       for ep in ("p_ref", "p_integral") for readout in READOUTS]
        four_batches = np.array(four_batches).reshape(20, -1)
        four_point = four_batches.mean(axis=0)
        four_loo = (20*four_point-four_batches)/19
        new_loo = np.column_stack((raw_loo, loo[:, identified], four_loo))
        new_factor = np.sqrt(19/20)*(new_loo-new_loo.mean(axis=0))
        old_row = old["sizes"][str(n)]
        old_factor = json.loads(gzip.decompress(read(OLD, OLD_DIR+"/"+old_row["complete_covariance_factor_file"])))
        factor = np.column_stack((np.array(old_factor["factor"]), new_factor))
        fn = f"N{n}.complete_common_factor.json.gz"
        packed = {"labels": old_factor["labels"]+raw_labels+["birth_partition."+label for label, keep in zip(labels, identified) if keep]+["four_change_type."+label for label in four_labels],
                  "factor": factor.tolist(), "covariance": "factor.T@factor; original twenty LOO rows, no inverse", "rank_at_most": 19,
                  "new_raw_labels": raw_labels, "new_raw_20_batch_means": raw_batches.tolist(), "gamma_LOO_identified": loo[:, identified].tolist(),
                  "four_type_labels": four_labels, "four_type_20_batch_means": four_batches.tolist(),
                  "batch_ids": list(range(20))}
        compressed = gzip.compress((json.dumps(packed, separators=(",", ":"), allow_nan=False)+"\n").encode(), mtime=0)
        (OUT/fn).write_bytes(compressed)
        gf = new_factor[:, len(raw_labels):len(raw_labels)+int(identified.sum())]; cov = gf.T@gf
        four_factor = np.sqrt(19/20)*(four_loo-four_loo.mean(axis=0)); four_cov = four_factor.T@four_factor
        se = np.full(len(point), np.nan); se[identified] = np.sqrt(np.diag(cov))
        blocks = {}
        for group, ix in [("all", list(range(9))), ("01+10", [1, 3])]+[("cell."+c, [j]) for j, c in enumerate(CELLS)]:
            b = batches[:, ix].sum(axis=1)
            blocks[group] = {part: {"B": b[:, k].mean(axis=0).tolist(), "B_se": (b[:, k].std(axis=0, ddof=1)/np.sqrt(20)).tolist()} for k, part in enumerate(PARTS)}
        sizes[str(n)] = {"delta_cos4": delta, "batch_ids": list(range(20)), "census_probability_counts_by_batch": counts,
            "blocks": blocks, "gamma_labels": labels, "gamma_estimate": [float(x) if ok else None for x, ok in zip(point, identified)],
            "gamma_se": [float(x) if ok else None for x, ok in zip(se, identified)], "gamma_covariance_identified": cov.tolist(),
            "four_type_labels": four_labels, "four_type_estimate": four_point.tolist(), "four_type_se": np.sqrt(np.diag(four_cov)).tolist(),
            "four_type_covariance": four_cov.tolist(), "four_type_20_batch_means": four_batches.tolist(),
            "complete_covariance_factor_file": fn, "complete_covariance_factor_sha256": sha256(compressed).hexdigest()}
    result = {"schema": "matching-one/p334-exact-birth-partition-response/v2", "source_commit": SOURCE,
        "census_commit": census_commit, "census_template": args.census_template, "source_sha256": hashes,
        "observer_labels": LABELS, "parts": list(PARTS), "sizes": sizes,
        "target": "Within=E_Z[ps Bs+pb Bb], Between=E_Z[ps pb (mu_s-mu_b)(mu_s-mu_b)^T]; exact ps is applied separately before prefix averaging",
        "four_type_target": "c=(first_change,second_change); Within=sum_c E_Z[p_c Cov(m|c,Z)]; Between=E_Z Cov_c(E[m|c,Z]); same-c quartet B divided by exact per-prefix p_c, empty classes zero",
        "boundary": "Prefix census is deterministic information, not new suffix data. Empty classes at ps=0/1 contribute zero. Finite-sample matrices and signed Gamma shares are not PSD-clipped. Same original twenty batches and common covariance; no new MC/DP/replay/test suite."}
    (OUT/"score.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    lines = ["# Exact-probability within/between immediate-birth response", ""]
    for n, r in sizes.items():
        lines += [f"## N{n}", "", "| Group / part | canonical Gamma +/- SE | integrated Gamma +/- SE |", "|---|---:|---:|"]
        for group in ("all", "01+10"):
            for part in PARTS:
                columns = []
                for ep in ("p_ref", "p_integral"):
                    k = r["gamma_labels"].index(f"{group}.{part}.{ep}.Gamma_pair")
                    columns.append(f"{r['gamma_estimate'][k]:.10g} +/- {r['gamma_se'][k]:.6g}")
                lines.append("| "+group+" / "+part+" | "+" | ".join(columns)+" |")
        lines.append("")
        lines += ["| Four-change group / part | canonical Gamma +/- SE | integrated Gamma +/- SE |", "|---|---:|---:|"]
        for group in ("all", "01+10"):
            for part in FOUR_PARTS:
                columns = []
                for ep in ("p_ref", "p_integral"):
                    k = r["four_type_labels"].index(f"{group}.{part}.{ep}.Gamma")
                    columns.append(f"{r['four_type_estimate'][k]:.10g} +/- {r['four_type_se'][k]:.6g}")
                lines.append("| "+group+" / "+part+" | "+" | ".join(columns)+" |")
        lines.append("")
    lines += [result["target"], "", result["boundary"], ""]
    (OUT/"REPORT.md").write_text("\n".join(lines)); print("\n".join(lines))


if __name__ == "__main__": main()

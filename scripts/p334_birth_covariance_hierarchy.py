#!/usr/bin/env python3
"""Separate a birth covariance response into suffix, label and prefix parts.

No paths are sampled. Cross-quartet U-products remove conditional mean-product
bias; cross-prefix U-products remove the global mean-product bias. One original
batch is deleted at a time, retaining the parent block's covariance convention.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
CELLS = ("all", "00", "01", "02", "10", "20")
ORIS = ("first", "second")
MARKS = ("plus", "minus")
MOMENTS = ("x", "y", "xx", "xy", "yy", "cross_xx", "cross_xy", "cross_yy")
PAIRS = ((0, 0), (0, 1), (1, 1))
PAIR_NAMES = ("xx", "xy", "yy")
DELTA = {325: -0.7634556213017751, 425: -0.8928996539792388}


def prefix_products(b, h):
    """b: (prefix,quartet,ori,moment), h: (prefix,quartet,mark,ori,moment)."""
    q = b.shape[1]
    bs, hs = b.sum(axis=1), h.sum(axis=1)
    bm, hm = bs/q, hs/q
    ub, diagb, uh, diagh = [], [], [], []
    for i, j in PAIRS:
        ub.append((bs[..., i]*bs[..., j]-(b[..., i]*b[..., j]).sum(axis=1))/(q*(q-1)))
        diagb.append(bm[..., i]*bm[..., j])
        # Both product derivatives are needed, including factor two for i=j.
        cross = (bs[:, None, :, i]*hs[..., j]+bs[:, None, :, j]*hs[..., i])
        same = (b[:, :, None, :, i]*h[..., j]+b[:, :, None, :, j]*h[..., i]).sum(axis=1)
        uh.append((cross-same)/(q*(q-1)))
        diagh.append(bm[:, None, :, i]*hm[..., j]+bm[:, None, :, j]*hm[..., i])
    return bm, hm, np.stack(ub, -1), np.stack(diagb, -1), np.stack(uh, -1), np.stack(diagh, -1)


def batch_coordinates(b, h, batches, cells):
    bm, hm, ub, db, uh, dh = prefix_products(b, h)
    fields = {}
    for cell in CELLS:
        mask = np.ones(len(b), dtype=bool) if cell == "all" else cells == 3*int(cell[0])+int(cell[1])
        fields[f"{cell}.mass"] = mask.astype(float)
        for oi, ori in enumerate(ORIS):
            for j, name in enumerate(MOMENTS):
                fields[f"{cell}.B.{ori}.{name}"] = bm[:, oi, j]*mask
            for j, name in enumerate(PAIR_NAMES):
                fields[f"{cell}.B.{ori}.offquartet.{name}"] = ub[:, oi, j]*mask
                fields[f"{cell}.B.{ori}.prefixdiag.{name}"] = db[:, oi, j]*mask
            for mi, mark in enumerate(MARKS):
                for j, name in enumerate(MOMENTS):
                    fields[f"{cell}.H.{mark}.{ori}.{name}"] = hm[:, mi, oi, j]*mask
                for j, name in enumerate(PAIR_NAMES):
                    fields[f"{cell}.H.{mark}.{ori}.offquartet.{name}"] = uh[:, mi, oi, j]*mask
                    fields[f"{cell}.H.{mark}.{ori}.prefixdiag.{name}"] = dh[:, mi, oi, j]*mask
    labels = list(fields)
    x = np.column_stack(list(fields.values()))
    batch_ids = np.unique(batches)
    counts = np.array([(batches == k).sum() for k in batch_ids])
    if not np.all(counts == counts[0]):
        raise ValueError("This original-batch scorer requires equal prefix counts")
    return labels, np.array([x[batches == k].mean(axis=0) for k in batch_ids]), batch_ids, counts


def shape_coordinates(xx, xy, yy):
    return {"var_x": xx, "cov_xy": xy, "var_y": yy,
            "var_C": (xx+2*xy+yy)/4, "var_W": xx-2*xy+yy,
            "cov_CW": (yy-xx)/2}


def derive(row, labels, prefix_count, delta):
    d, out = dict(zip(labels, row)), {}
    p = prefix_count
    for cell in CELLS:
        out[f"{cell}.mass"] = d[f"{cell}.mass"]
        # Tangent components are zero-padded to the complete population; use
        # the global baseline means when assigning the centering to each cell.
        orows = {}
        for ori in ORIS:
            mu = [d[f"all.B.{ori}.{v}"] for v in ("x", "y")]
            if cell == "all":
                base_parts = {part: [] for part in ("total", "within_prefix", "within_suffix", "between_labels", "between_prefixes")}
                for (i, j), pair in zip(PAIRS, PAIR_NAMES):
                    raw = d[f"all.B.{ori}.{pair}"]
                    cross = d[f"all.B.{ori}.cross_{pair}"]
                    off = d[f"all.B.{ori}.offquartet.{pair}"]
                    diag = d[f"all.B.{ori}.prefixdiag.{pair}"]
                    global_product = (p*mu[i]*mu[j]-diag)/(p-1)
                    for part, v in {"total": raw-global_product,
                                    "within_prefix": raw-off,
                                    "within_suffix": raw-cross,
                                    "between_labels": cross-off,
                                    "between_prefixes": off-global_product}.items():
                        base_parts[part].append(v)
                for part, v in base_parts.items():
                    for k, value in shape_coordinates(*v).items():
                        out[f"baseline.{ori}.{k}.{part}"] = value
            for mark in MARKS:
                hm = [d[f"{cell}.H.{mark}.{ori}.{v}"] for v in ("x", "y")]
                parts = {part: [] for part in ("total", "within_prefix", "within_suffix", "between_labels", "between_prefixes")}
                for (i, j), pair in zip(PAIRS, PAIR_NAMES):
                    raw = d[f"{cell}.H.{mark}.{ori}.{pair}"]
                    cross = d[f"{cell}.H.{mark}.{ori}.cross_{pair}"]
                    off = d[f"{cell}.H.{mark}.{ori}.offquartet.{pair}"]
                    diag = d[f"{cell}.H.{mark}.{ori}.prefixdiag.{pair}"]
                    global_product = (p*(mu[i]*hm[j]+mu[j]*hm[i])-diag)/(p-1)
                    for part, v in {"total": raw-global_product,
                                    "within_prefix": raw-off,
                                    "within_suffix": raw-cross,
                                    "between_labels": cross-off,
                                    "between_prefixes": off-global_product}.items():
                        parts[part].append(v)
                one = {}
                for part, v in parts.items():
                    for k, value in shape_coordinates(*v).items():
                        one[f"{k}.{part}"] = value
                        out[f"{cell}.{mark}->{ori}.{k}.{part}"] = value
                orows[(mark, ori)] = one
        for mark in MARKS:
            a, b = orows[(mark, "first")], orows[(mark, "second")]
            for k in a:
                out[f"{cell}.{mark}->S.{k}"] = (a[k]+b[k])/2
                out[f"{cell}.{mark}->D.{k}"] = (a[k]-b[k])/delta
    return out


def score_arrays(n, b, h, batches, cells):
    labels, raw, batch_ids, counts = batch_coordinates(b, h, batches, cells)
    k = len(batch_ids)
    raw_mean = raw.mean(axis=0)
    raw_loo = (k*raw_mean-raw)/(k-1)
    point = derive(raw_mean, labels, int(counts.sum()), DELTA[n])
    loo = np.array([list(derive(row, labels, int(counts.sum()-counts[j]), DELTA[n]).values())
                    for j, row in enumerate(raw_loo)])
    factor = np.sqrt((k-1)/k)*(loo-loo.mean(axis=0))
    return {"N": n, "batch_ids": batch_ids.tolist(), "prefix_counts": counts.tolist(),
            "delta_cos4": DELTA[n], "raw_labels": labels, "raw_batch_means": raw.tolist(),
            "labels": list(point), "estimate": list(point.values()),
            "se": np.linalg.norm(factor, axis=0).tolist(), "LOO": loo.tolist(),
            "factor": factor.tolist(), "factor_convention": "factor.T@factor; original deleted-batch order"}


def load_arrays(path):
    # The extractor's schema is explicit; command flags avoid positional or
    # shape-based guesses if a future archive uses different field names.
    with np.load(path, allow_pickle=False) as z:
        return {k: z[k] for k in z.files}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--output", type=Path, default=ROOT/"results/p334-birth-covariance-hierarchy")
    parser.add_argument("--baseline-key", default="baseline")
    parser.add_argument("--tangent-key", default="tangent")
    parser.add_argument("--batch-key", default="batch")
    parser.add_argument("--cell-key", default="cell")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    result = {"schema": "matching-one/p334-birth-covariance-hierarchy/v1",
              "source_commit": args.source_commit,
              "reader_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
              "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "moment_order": MOMENTS, "cell_order": CELLS,
              "estimand": "Covariance of X=K1/(N+1),Y=K2/(N+1), with fixed original prefix distribution. Per-geometry components before S/D.",
              "identity": "total=within_suffix+between_labels+between_prefixes; within_prefix=within_suffix+between_labels",
              "products": "Different quartet means for within-prefix products; different prefixes for global mean products. Each original-batch deletion recomputes the latter with its retained prefix count.",
              "dependency": "Original e32a8593/959a7fa2 prefix/fork/contact block; exact census from PR509. No new prefixes, quartets, suffixes, DP or model fitting.",
              "input_sha256": {}, "sizes": {}}
    for n in (325, 425):
        path = args.input_dir/f"N{n}.npz"
        z = load_arrays(path)
        r = score_arrays(n, z[args.baseline_key], z[args.tangent_key], z[args.batch_key], z[args.cell_key])
        result["sizes"][str(n)] = r
        result["input_sha256"][path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
        for mark, axis in (("plus", "S"), ("minus", "D")):
            for part in ("total", "within_prefix", "within_suffix", "between_labels", "between_prefixes"):
                name = f"all.{mark}->{axis}.cov_xy.{part}"
                j = r["labels"].index(name)
                print(n, name, f"{r['estimate'][j]:+.10g} +/- {r['se'][j]:.6g}", flush=True)
    (args.output/"score.json").write_text(json.dumps(result, separators=(",", ":"), allow_nan=False)+"\n")


if __name__ == "__main__":
    main()

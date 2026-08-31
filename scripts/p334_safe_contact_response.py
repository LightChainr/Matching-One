#!/usr/bin/env python3
"""Read how topology-safe R0 insertions load future birth responses."""
import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import subprocess

import numpy as np
from scipy.stats import binom

ROOT = Path(__file__).resolve().parents[1]
SOURCE = "e32a85939279b8574278024d647b56d2d1485247"
FORK_PATH = "results/p334-nested-next-label-forks"
P_REF = 0.59274605079
FEATURES = ["contractible_cycles", "component_mergers", "isolated_site"]
RESPONSES = ["p_ref.F1", "p_ref.F2", "p_integral.F1", "p_integral.F2"]
GROUPS = ["all_R0_orientations", "01+10_R0_orientations", "R0_safe_equal_contact_degree"]
GTRI = np.triu_indices(3)
XTRI = np.triu_indices(4)


def blob(commit, path):
    return subprocess.check_output(["git", "show", f"{commit}:{path}"], cwd=ROOT)


def array_csv(data):
    with gzip.GzipFile(fileobj=io.BytesIO(data)) as stream:
        header = stream.readline().decode().strip().split(",")
        rows = np.loadtxt(stream, delimiter=",", dtype=np.int64)
    return header, rows


def labels():
    out = []
    for group in GROUPS:
        out.append(group + ".safe_pair_mass")
        out.extend(f"{group}.GG[{FEATURES[i]},{FEATURES[j]}]" for i, j in zip(*GTRI))
        out.extend(f"{group}.GX[{g},{x}]" for g in FEATURES for x in RESPONSES)
        out.extend(f"{group}.XX[{RESPONSES[i]},{RESPONSES[j]}]" for i, j in zip(*XTRI))
    return out


def score_batch(n, batch, contact_commit, contact_path, tail):
    raw_path = f"{FORK_PATH}/N{n}/N{n}.batch{batch:02}.csv.gz"
    mark_path = f"{contact_path}/N{n}/N{n}.batch{batch:02}.csv.gz"
    raw_blob, mark_blob = blob(SOURCE, raw_path), blob(contact_commit, mark_path)
    h, raw = array_csv(raw_blob)
    ch, marks = array_csv(mark_blob)
    hi, ci = {x: i for i, x in enumerate(h)}, {x: i for i, x in enumerate(ch)}
    raw = raw[np.lexsort(tuple(raw[:, hi[k]] for k in ("replica", "group", "quartet", "counter")))]
    marks = marks[np.lexsort(tuple(marks[:, ci[k]] for k in ("group", "quartet", "counter")))]
    raw = raw.reshape(1000, 8, 2, 2, -1)
    marks = marks.reshape(1000, 8, 2, -1)
    for key in ("counter", "quartet", "group", "next_label"):
        if not np.array_equal(raw[..., 0, hi[key]], marks[..., ci[key]]):
            raise ValueError(f"Contact/source key disagreement: N{n} batch{batch} {key}")
    old = raw[:, 0, 0, 0, [hi["first_rank"], hi["second_rank"]]]
    cell = 3 * old[:, 0] + old[:, 1]
    accum = {name: [0.0, np.zeros((3, 3)), np.zeros((3, 4)), np.zeros((4, 4))]
             for name in GROUPS}
    for orientation, rank_index in (("first", 0), ("second", 1)):
        e, c = (marks[..., ci[f"{orientation}_{key}"]] for key in ("e", "c"))
        next_rank = raw[..., 0, hi[f"{orientation}_next_rank"]]
        safe_pair = (next_rank[:, :, 0] == 0) & (next_rank[:, :, 1] == 0)
        active = (old[:, rank_index] == 0)[:, None] & safe_pair
        geom = np.stack((e-c, np.maximum(c-1, 0), (c == 0).astype(float)), axis=-1)
        dg = geom[:, :, 0]-geom[:, :, 1]
        k1, k2 = (raw[..., hi[f"{orientation}_{key}"]] for key in ("k1", "k2"))
        obs = np.stack((tail[k1], tail[k2], 1-k1/(n+1), 1-k2/(n+1)), axis=-1)
        a, b = obs[:, :, 0, 0]-obs[:, :, 1, 0], obs[:, :, 0, 1]-obs[:, :, 1, 1]
        dm = (a+b)/2
        for name in GROUPS:
            if name == GROUPS[0]:
                chosen = active
            elif name == GROUPS[1]:
                chosen = active & np.isin(cell, [1, 3])[:, None]
            else:
                chosen = active & (e[:, :, 0] == e[:, :, 1])
            gg, xx, aa, bb = dg[chosen], dm[chosen], a[chosen], b[chosen]
            # 1/2 for the orientation mixture, 1/2 for the iid-label identity.
            accum[name][0] += chosen.sum()/(2*8000)
            accum[name][1] += gg.T @ gg/(4*8000)
            accum[name][2] += gg.T @ xx/(4*8000)
            accum[name][3] += (aa.T @ bb + bb.T @ aa)/(8*8000)
    values = []
    for name in GROUPS:
        mass, gg, gx, xx = accum[name]
        values.extend([mass, *gg[GTRI], *gx.ravel(), *xx[XTRI]])
    hashes = {raw_path: hashlib.sha256(raw_blob).hexdigest(), mark_path: hashlib.sha256(mark_blob).hexdigest()}
    return np.array(values), hashes


def named_readout(raw, names, n):
    d = dict(zip(names, raw))
    out = {}
    for group in GROUPS:
        for feature in FEATURES:
            v = d[f"{group}.GG[{feature},{feature}]"]
            for response in RESPONSES:
                cv = d[f"{group}.GX[{feature},{response}]"]
                out[f"{group}.cov[{feature},{response}]"] = cv
                out[f"{group}.pooled_slope[{feature},{response}]"] = cv/v if v > 0 else np.nan
            for i in (1, 2):
                cv = -(n+1)*d[f"{group}.GX[{feature},p_integral.F{i}]"]
                out[f"{group}.cov[{feature},K{i}]"] = cv
                out[f"{group}.pooled_slope[{feature},K{i}]"] = cv/v if v > 0 else np.nan
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contact-commit", required=True)
    parser.add_argument("--contact-path", default="results/p334-next-label-contact-coordinates")
    parser.add_argument("--output", type=Path, default=ROOT/"results/p334-safe-contact-response")
    args = parser.parse_args()
    contact_commit = subprocess.check_output(["git", "rev-parse", args.contact_commit+"^{commit}"], cwd=ROOT, text=True).strip()
    if args.output.exists():
        raise FileExistsError(f"Will not replace an existing scientific readout: {args.output}")
    args.output.mkdir(parents=True)
    out = {"source_commit": SOURCE, "contact_commit": contact_commit,
           "p_ref": P_REF, "features": FEATURES, "responses": RESPONSES,
           "estimand": "Equal orientation mixture, full original-prefix denominator: E[1_R0 pi_safe^2 Cov_U(g,m|Z,own-rank-safe)]. Equal-degree mask instead sums pi_safe,e^2 Cov_U(g,m|Z,own-rank-safe,e). This is not paired H4 or a causal feature intervention.",
           "source_sha256": {}, "sizes": {}}
    names = labels()
    for n in (325, 425):
        tail = binom.sf(np.arange(n+2)-1, n, P_REF)
        rows = []
        for batch in range(20):
            row, hashes = score_batch(n, batch, contact_commit, args.contact_path, tail)
            rows.append(row)
            out["source_sha256"].update(hashes)
        raw = np.stack(rows)
        estimate = raw.mean(axis=0)
        raw_factor = (raw-estimate)/np.sqrt(20*19)
        primary = named_readout(estimate, names, n)
        pnames = list(primary)
        values = np.array(list(primary.values()))
        loo = np.array([list(named_readout((20*estimate-row)/19, names, n).values()) for row in raw])
        factor = np.sqrt(19/20)*(loo-loo.mean(axis=0))
        out["sizes"][str(n)] = {"batch_ids": list(range(20)), "raw_labels": names,
            "raw_joint_20_batch_means": raw.tolist(), "raw_estimate": estimate.tolist(),
            "raw_covariance_factor": raw_factor.tolist(), "labels": pnames,
            "estimate": values.tolist(), "se": np.linalg.norm(factor, axis=0).tolist(),
            "LOO": loo.tolist(), "factor": factor.tolist()}
        print(f"N{n}: complete own-safe contact-response readout", flush=True)
        for name, value, se in zip(pnames, values, np.linalg.norm(factor, axis=0)):
            if name.startswith((GROUPS[0], GROUPS[2])) and ".pooled_slope" in name and name.endswith((",K1]", ",K2]")):
                print(f"  {name}: {value:.9g} +/- {se:.5g}", flush=True)
    (args.output/"score.json").write_text(json.dumps(out, indent=2)+"\n")
    print(args.output, flush=True)


if __name__ == "__main__":
    main()

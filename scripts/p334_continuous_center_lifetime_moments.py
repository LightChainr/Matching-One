#!/usr/bin/env python3
"""One pass over frozen common-label forks: continuous C/W raw moments and tangents."""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np

from p334_safe_contact_response import ROOT, SOURCE, FORK_PATH, blob, array_csv
from p334_common_label_euler_tangent import CONTACT, DELTA, DENOMINATOR

CELLS = ["00", "01", "02", "10", "20", "all"]
MARKS = ["plus", "minus"]
ORIENTATIONS = ["first", "second"]
MOMENTS = ["C", "W", "CW", "C2", "W2", "tau1_cubed", "tau2_cubed", "R1_plateau_M2"]
POLICY_SOURCE = "4db356e1b026853468f94d59d938895a2367ceb7"
BASELINE_TAILS = 32000  # 1000 prefixes x 8 quartets x 2 labels x 2 tails.


def moment_numerators(k1, k2):
    """Joint uniform order statistics, not products of conditional means."""
    q1, q2, q12 = k1*(k1+1), k2*(k2+1), k1*(k2+1)
    t1, t2 = q1*(k1+2), q2*(k2+2)
    return np.stack((k1+k2, k2-k1, q2-q1, q1+2*q12+q2,
                     q1-2*q12+q2, t1, t2, t2-t1), axis=-1)


def moment_denominators(n):
    d2, d3 = (n+1)*(n+2), (n+1)*(n+2)*(n+3)
    return np.array([2*(n+1), n+1, 2*d2, 4*d2, d2, d3, d3, 3*d3], dtype=np.int64)


def read_batch(n, batch):
    paths = [f"{FORK_PATH}/N{n}/N{n}.batch{batch:02}.csv.gz",
             f"results/p334-next-label-contact-coordinates/N{n}/N{n}.batch{batch:02}.csv.gz"]
    bb = [blob(SOURCE, paths[0]), blob(CONTACT, paths[1])]
    h, raw = array_csv(bb[0])
    ch, marks = array_csv(bb[1])
    hi, ci = {x:i for i,x in enumerate(h)}, {x:i for i,x in enumerate(ch)}
    raw = raw[np.lexsort(tuple(raw[:,hi[k]] for k in ("replica", "group", "quartet", "counter")))].reshape(1000,8,2,2,-1)
    marks = marks[np.lexsort(tuple(marks[:,ci[k]] for k in ("group", "quartet", "counter")))].reshape(1000,8,2,-1)
    for key in ("counter", "quartet", "group", "next_label"):
        if not np.array_equal(raw[...,0,hi[key]], marks[...,ci[key]]):
            raise ValueError(f"Source/contact key disagreement: N{n}/{batch}/{key}")
    ranks = [raw[:,0,0,0,hi[f"{o}_rank"]] for o in ORIENTATIONS]
    cell = ranks[0]*3+ranks[1]
    selected = np.ones((1000,8), dtype=bool)
    loops = []
    for o, old in zip(ORIENTATIONS, ranks):
        nr = marks[...,ci[f"{o}_rank_after"]]
        e, c = (marks[...,ci[f"{o}_{x}"]] for x in ("e", "c"))
        selected &= np.all(nr == old[:,None,None], axis=2)
        selected &= e[:,:,0] == e[:,:,1]
        loops.append((old == 0)[:,None,None]*(e-c))
    doubled = [loops[0]+loops[1], loops[0]-loops[1]]
    dh = np.stack([g[:,:,0]-g[:,:,1] for g in doubled], axis=-1)
    baseline = np.zeros((6,2,len(MOMENTS)), dtype=np.int64)
    tangent = np.zeros((6,2,2,len(MOMENTS)), dtype=np.int64)
    prevalence, selected_counts = np.zeros(6, dtype=np.int64), np.zeros(6, dtype=np.int64)
    for ic, name in enumerate(CELLS):
        in_cell = np.ones(1000, dtype=bool) if name == "all" else cell == 3*int(name[0])+int(name[1])
        chosen = selected & in_cell[:,None]
        prevalence[ic], selected_counts[ic] = in_cell.sum(), chosen.sum()
        for io, o in enumerate(ORIENTATIONS):
            k1, k2 = (raw[...,hi[f"{o}_{x}"]] for x in ("k1", "k2"))
            moment = moment_numerators(k1, k2)
            baseline[ic,io] = moment[in_cell].sum(axis=(0,1,2,3))
            # One label's two tails estimate its conditional mean. The existing
            # 64000 denominator includes both half factors and the suffix mean.
            difference = moment[:,:,0].sum(axis=2)-moment[:,:,1].sum(axis=2)
            tangent[ic,:,io] = dh[chosen].T @ difference[chosen]
    # Five cells exhaust tangent support, but not the unperturbed population.
    if not np.array_equal(tangent[-1], tangent[:-1].sum(axis=0)):
        raise ValueError("An inactive rank cell unexpectedly contributes to the tangent")
    hashes = {p:hashlib.sha256(b).hexdigest() for p,b in zip(paths,bb)}
    return baseline, tangent, prevalence, selected_counts, hashes


def pack_size(n, baseline, tangent, prevalence, selected_counts):
    denominators = moment_denominators(n)
    b = baseline/(BASELINE_TAILS*denominators)
    h = tangent/(DENOMINATOR*denominators)
    fields = {}
    for ic, cell in enumerate(CELLS):
        for io, orientation in enumerate(ORIENTATIONS):
            for im, moment in enumerate(MOMENTS):
                fields[f"{cell}.baseline.{orientation}.{moment}"] = b[:,ic,io,im]
            for ig, mark in enumerate(MARKS):
                for im, moment in enumerate(MOMENTS):
                    fields[f"{cell}.H.{mark}.{orientation}.{moment}"] = h[:,ic,ig,io,im]
    x = np.column_stack(list(fields.values()))
    mean, se = x.mean(axis=0), x.std(axis=0, ddof=1)/np.sqrt(20)
    return {"batch_ids":list(range(20)), "delta_cos4":DELTA[n],
            "moment_denominators":denominators.tolist(),
            "batch_baseline_integer_numerators":baseline.tolist(),
            "batch_tangent_integer_numerators":tangent.tolist(),
            "batch_prefix_counts":prevalence.tolist(),
            "batch_selected_pair_counts":selected_counts.tolist(),
            "labels":list(fields), "joint_20_batch_means":x.tolist(),
            "estimate":mean.tolist(), "se":se.tolist()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT/"results/p334-continuous-center-lifetime-moments")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    result = {"source_commit":SOURCE, "contact_commit":CONTACT, "policy_commit":POLICY_SOURCE,
        "reader_commit":subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
        "reader_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "new_samples":0, "new_DP":0, "source_read_passes":1,
        "baseline_batch_denominator":BASELINE_TAILS, "tangent_batch_denominator":DENOMINATOR,
        "cell_order":CELLS, "mark_order":MARKS, "orientation_order":ORIENTATIONS,
        "moment_order":MOMENTS, "baseline_axes":["batch","cell","orientation","moment"],
        "tangent_axes":["batch","cell","mark","orientation","moment"],
        "baseline_semantics":"Unperturbed 32-tail conditional average per original prefix; every cell zero-padded to full1000/batch. all includes all nine rank cells, so five displayed cells do not exhaust baseline.",
        "tangent_semantics":"Shared labels; both orientations rank-safe; U/V equal joint contact degrees. g_plus/minus=half sum/difference of R0-only (e-c). Existing pi_a-preserving exp(t*pi_a*g) policy. Five cells exhaust tangent.",
        "dependence":"Same original20 paired prefix batches per N, all existing tails and contacts; no independent evidence. Retain batch IDs across old/new observables; downstream connected products use pooled means and delete-one recomputation.",
        "continuous_semantics":"Joint uniform order statistics conditional on integer K1,K2; includes their common Beta/order-statistic variance, not squares of K/(N+1).",
        "source_sha256":{}, "sizes":{}}
    for n in (325,425):
        baseline, tangent, prevalence, selected_counts = [],[],[],[]
        for batch in range(20):
            b,h,p,c,hashes = read_batch(n,batch)
            baseline.append(b); tangent.append(h); prevalence.append(p); selected_counts.append(c)
            result["source_sha256"].update(hashes)
        s = pack_size(n, np.array(baseline), np.array(tangent), np.array(prevalence), np.array(selected_counts))
        result["sizes"][str(n)] = s
        for name, value, se in zip(s["labels"],s["estimate"],s["se"]):
            if name.startswith("all.") and name.endswith((".C",".W",".CW",".C2",".W2")):
                print(n,name,f"{value:.11g} +/- {se:.5g}",flush=True)
        print(f"N{n}: one complete read, 20 original batches",flush=True)
    (args.output/"batch_moments.json").write_text(json.dumps(result,separators=(",",":"),allow_nan=False)+"\n")
    print(args.output,flush=True)


if __name__ == "__main__":
    main()

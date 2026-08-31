#!/usr/bin/env python3
"""Read a shared-label tangent preserving both immediate Euler/rank states."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.stats import binom

from p334_safe_contact_response import ROOT, SOURCE, FORK_PATH, P_REF, blob, array_csv

CONTACT = "959a7fa26677c416b874d272f1ba66523fb38f73"
CELLS = ["00", "01", "02", "10", "20"]
MARKS = ["plus", "minus"]
ORIENTATIONS = ["first", "second"]
DENOMINATOR = 64000  # 8000 quartets x 2 half-difference x 2 mark x 2 suffix mean.
DELTA = {325: -0.7634556213017751, 425: -0.8928996539792388}


def read_batch(n, batch):
    paths = [f"{FORK_PATH}/N{n}/N{n}.batch{batch:02}.csv.gz",
             f"results/p334-next-label-contact-coordinates/N{n}/N{n}.batch{batch:02}.csv.gz"]
    bb = [blob(SOURCE, paths[0]), blob(CONTACT, paths[1])]
    h, raw = array_csv(bb[0])
    ch, marks = array_csv(bb[1])
    hi, ci = {x:i for i,x in enumerate(h)}, {x:i for i,x in enumerate(ch)}
    raw = raw[np.lexsort(tuple(raw[:, hi[k]] for k in ("replica", "group", "quartet", "counter")))].reshape(1000,8,2,2,-1)
    marks = marks[np.lexsort(tuple(marks[:,ci[k]] for k in ("group", "quartet", "counter")))].reshape(1000,8,2,-1)
    for key in ("counter", "quartet", "group", "next_label"):
        if not np.array_equal(raw[...,0,hi[key]], marks[...,ci[key]]):
            raise ValueError(f"Mismatched saved source keys N{n}/{batch}/{key}")
    ranks = [raw[:,0,0,0,hi[f"{o}_rank"]] for o in ORIENTATIONS]
    cell = ranks[0]*3+ranks[1]
    selected = np.ones((1000,8), dtype=bool)
    loops = []
    for o, old in zip(ORIENTATIONS, ranks):
        nr = marks[...,ci[f"{o}_rank_after"]]
        e = marks[...,ci[f"{o}_e"]]
        c = marks[...,ci[f"{o}_c"]]
        selected &= np.all(nr==old[:,None,None], axis=2)
        selected &= e[:,:,0]==e[:,:,1]
        loops.append((old==0)[:,None,None]*(e-c))
    # Doubled marks keep the complete signed histograms exactly integral.
    doubled = [loops[0]+loops[1], loops[0]-loops[1]]
    dh = np.stack([x[:,:,0]-x[:,:,1] for x in doubled], axis=-1)
    hist = np.zeros((5,2,2,2,n+1), dtype=np.int64)
    gg = np.zeros((5,2,2))
    counts = np.zeros((5,3), dtype=np.int64)
    for ic, cname in enumerate(CELLS):
        chosen = selected & (cell == 3*int(cname[0])+int(cname[1]))[:,None]
        dmark = dh[chosen]
        counts[ic] = [chosen.sum(), np.count_nonzero(dmark[:,0]), np.count_nonzero(dmark[:,1])]
        gg[ic] = dmark.T @ dmark / DENOMINATOR
        for im in range(2):
            for io, o in enumerate(ORIENTATIONS):
                for ib, clock in enumerate(("k1", "k2")):
                    ks = raw[...,hi[f"{o}_{clock}"]][chosen]
                    for g, sign in ((0,1),(1,-1)):
                        for replica in (0,1):
                            np.add.at(hist[ic,im,io,ib], ks[:,g,replica], sign*dmark[:,im])
    return hist, gg, counts, {p:hashlib.sha256(b).hexdigest() for p,b in zip(paths,bb)}


def summarize(hist, gg, counts, n):
    # Integer cumulative coefficients avoid a spurious floating endpoint tail.
    coeff = np.cumsum(hist, axis=-1)
    pref = (coeff @ binom.pmf(np.arange(n+1), n, P_REF))/DENOMINATOR
    integrated = coeff.sum(axis=-1)/(DENOMINATOR*(n+1))
    clock = (hist @ np.arange(n+1))/DENOMINATOR
    fields = {}
    for group, ix in [("all", list(range(5))), *[(name,[i]) for i,name in enumerate(CELLS)]]:
        hp, he, hk = (x[:,ix].sum(axis=1) for x in (pref, integrated, clock))
        gm = gg[:,ix].sum(axis=1)
        cc = counts[:,ix].sum(axis=1)
        for j, mark in enumerate(MARKS):
            for output in ("S", "D"):
                w = np.array([.5,.5]) if output=="S" else np.array([1.,-1.])/DELTA[n]
                p, e, k = (np.einsum('boj,o->bj', x[:,j], w) for x in (hp,he,hk))
                label = f"{group}.{mark}->{output}."
                for ep, val in (("p_ref",p),("p_integral",e)):
                    fields[label+ep+".F1"] = val[:,0]
                    fields[label+ep+".F2"] = val[:,1]
                    fields[label+ep+".A"] = val[:,0]+val[:,1]
                    fields[label+ep+".E"] = val[:,1]-val[:,0]
                for name, val in (("K1",k[:,0]),("K2",k[:,1]),("C",k.sum(1)/2),("W",k[:,1]-k[:,0])):
                    fields[label+name] = val
            fields[f"{group}.GG[{mark},{mark}]"] = gm[:,j,j]
            fields[f"{group}.nonzero_{mark}_pair_mass"] = cc[:,j+1]/8000
        fields[f"{group}.GG[plus,minus]"] = gm[:,0,1]
        fields[f"{group}.selected_pair_mass"] = cc[:,0]/8000
    x = np.column_stack(list(fields.values()))
    mean = x.mean(axis=0)
    factor = (x-mean)/np.sqrt(20*19)
    return {"batch_ids":list(range(20)),"delta_cos4":DELTA[n],"labels":list(fields),
            "joint_20_batch_means":x.tolist(),"estimate":mean.tolist(),
            "se":np.linalg.norm(factor,axis=0).tolist(),"factor":factor.tolist()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=ROOT/"results/p334-common-label-euler-tangent")
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    result={"source_commit":SOURCE,"contact_commit":CONTACT,"new_samples":0,
            "p_ref":P_REF,"batch_denominator":DENOMINATOR,"cell_order":CELLS,
            "mark_order":MARKS,"orientation_order":ORIENTATIONS,"birth_order":["F1","F2"],
            "histogram_axes":["batch","cell","mark","orientation","birth","k"],
            "policy":"Shared label; both-rank-safe joint degree classes; pi_a fixed; tilt exp(t*pi_a*g_plus/minus).",
            "dependence":"Same20 original paired prefix batches per N; all existing suffixes and contact marks; no new independent population.",
            "source_sha256":{},"sizes":{}}
    score={k:v for k,v in result.items() if k not in ("sizes","source_sha256")}
    score["sizes"]={}
    for n in (325,425):
        hist,gg,counts=[],[],[]
        for batch in range(20):
            h,g,c,hashes=read_batch(n,batch)
            hist.append(h);gg.append(g);counts.append(c)
            result["source_sha256"].update(hashes)
        hist,gg,counts=np.array(hist),np.array(gg),np.array(counts)
        result["sizes"][str(n)]={"batch_ids":list(range(20)),"delta_cos4":DELTA[n],
            "batch_integer_histograms":hist.tolist(),"batch_mark_covariances":gg.tolist(),
            "batch_selected_nonzero_counts":counts.tolist()}
        scored=summarize(hist,gg,counts,n)
        score["sizes"][str(n)]=scored
        for name,value,se in zip(scored["labels"],scored["estimate"],scored["se"]):
            if name.startswith("all.") and name.endswith(("p_ref.A","p_ref.E","p_integral.A","p_integral.E",".W")):
                print(n,name,f"{value:.10g} +/- {se:.6g}",flush=True)
    (args.output/"signed_birth_histograms.json").write_text(json.dumps(result,separators=(",",":"),allow_nan=False)+"\n")
    (args.output/"score.json").write_text(json.dumps(score,indent=2,allow_nan=False)+"\n")
    print(args.output,flush=True)


if __name__=="__main__":
    main()

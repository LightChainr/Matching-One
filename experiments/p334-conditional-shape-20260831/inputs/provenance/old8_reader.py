#!/usr/bin/env python3
"""Prefix/quartet moments under exact census scores; no new tails or DP."""
import argparse
import hashlib
import json
import platform
import subprocess
import time
from pathlib import Path

import numpy as np

from p334_safe_contact_response import ROOT, SOURCE, FORK_PATH, blob, array_csv
from p334_common_label_euler_tangent import CONTACT, DELTA

CENSUS_COMMIT = "ac5761ce504c3cd170fa42c86c17d6fb87f0375b"
CENSUS_ROOT = "experiments/p334-finite-source-20260831"
FEATURES = ["x", "y", "xx", "xy", "yy", "cross_xx", "cross_xy", "cross_yy"]
ORIENTATIONS = ["first", "second"]
MARKS = ["plus", "minus"]


def census(n):
    path = f"{CENSUS_ROOT}/census/N{n}/census.csv.gz"
    data = blob(CENSUS_COMMIT, path)
    names, rows = array_csv(data)
    col = {name:i for i,name in enumerate(names)}
    counters, ids = np.unique(rows[:,col["counter"]], return_inverse=True)
    if len(counters) != 20000:
        raise ValueError("Census must represent all original 20000 prefixes")
    ix = (ids, rows[:,col["first_e"]], rows[:,col["second_e"]])
    count = np.zeros((20000,5,5), dtype=np.int32)
    loop_sum = np.zeros((20000,5,5,2), dtype=np.int32)
    np.add.at(count, ix, rows[:,col["count"]])
    for io, o in enumerate(ORIENTATIONS):
        np.add.at(loop_sum[...,io], ix, rows[:,col["count"]]*rows[:,col[f"L_{o}"]])
    k0 = np.zeros(20000, dtype=np.int16)
    k0[ids] = rows[:,col["k0"]]
    return counters, count, loop_sum, k0, {path:hashlib.sha256(data).hexdigest()}


def read_batch(n, batch, counters, class_count, class_loop_sum, census_k0):
    raw_path = f"{FORK_PATH}/N{n}/N{n}.batch{batch:02}.csv.gz"
    contact_path = f"results/p334-next-label-contact-coordinates/N{n}/N{n}.batch{batch:02}.csv.gz"
    raw_bytes, contact_bytes = blob(SOURCE,raw_path), blob(CONTACT,contact_path)
    names, raw = array_csv(raw_bytes)
    contact_names, contact = array_csv(contact_bytes)
    hi, ci = {s:i for i,s in enumerate(names)}, {s:i for i,s in enumerate(contact_names)}
    raw = raw[np.lexsort(tuple(raw[:,hi[k]] for k in ("replica","group","quartet","counter")))].reshape(1000,8,2,2,-1)
    contact = contact[np.lexsort(tuple(contact[:,ci[k]] for k in ("group","quartet","counter")))].reshape(1000,8,2,-1)
    for key in ("counter","quartet","group","next_label"):
        if not np.array_equal(raw[...,0,hi[key]], contact[...,ci[key]]):
            raise ValueError(f"Frozen key disagreement N{n}/{batch}/{key}")
    counter = raw[:,0,0,0,hi["counter"]]
    ids = np.searchsorted(counters,counter)
    if not np.array_equal(counters[ids],counter):
        raise ValueError("Prefix absent from census")
    k0 = raw[:,0,0,0,hi["k0"]]
    if not np.array_equal(k0,census_k0[ids]):
        raise ValueError("Census and fork checkpoint disagree")
    old_rank = np.stack([raw[:,0,0,0,hi[f"{o}_rank"]] for o in ORIENTATIONS],axis=-1)
    rank_after = np.stack([contact[...,ci[f"{o}_rank_after"]] for o in ORIENTATIONS],axis=-1)
    es = np.stack([contact[...,ci[f"{o}_e"]] for o in ORIENTATIONS],axis=-1)
    cs = np.stack([contact[...,ci[f"{o}_c"]] for o in ORIENTATIONS],axis=-1)
    loops = (es-cs)*(old_rank[:,None,None,:] == 0)
    safe = np.all(rank_after == old_rank[:,None,None,:],axis=-1)
    ix = (np.broadcast_to(ids[:,None,None],safe.shape), es[...,0], es[...,1])
    mass_count, total_loop = class_count[ix], class_loop_sum[ix]
    # Exact centered score numerator. H = score_numerator/(2*d), zero off safe.
    # This avoids constructing finite-q_t weights and carries no fitted means.
    centered2 = mass_count[...,None]*loops-total_loop
    score_numerator = np.stack((centered2.sum(axis=-1),
                               centered2[...,0]-centered2[...,1]),axis=-1)
    score_numerator *= safe[...,None]
    d = n-k0
    score = score_numerator/(2*d[:,None,None,None])
    birth_k = np.stack([np.stack([raw[...,hi[f"{o}_{k}"]] for k in ("k1","k2")],axis=-1)
                        for o in ORIENTATIONS],axis=-2)
    x, y = birth_k[...,0]/(n+1), birth_k[...,1]/(n+1)
    # x,y axes: prefix,quartet,label,suffix,orientation.
    label_moments = np.stack((x,y,x*x,x*y,y*y),axis=-1).mean(axis=3)
    x0,x1,y0,y1 = x[:,:,:,0,:],x[:,:,:,1,:],y[:,:,:,0,:],y[:,:,:,1,:]
    label_cross = np.stack((x0*x1,(x0*y1+x1*y0)/2,y0*y1),axis=-1)
    label_features = np.concatenate((label_moments,label_cross),axis=-1)
    b = label_features.mean(axis=2)
    h = 0.5*(score[:,:,0]-score[:,:,1])[:,:,:,None,None]*(label_features[:,:,0]-label_features[:,:,1])[:,:,None,:,:]
    inactive = np.all(old_rank > 0,axis=-1)
    if np.count_nonzero(h[inactive]) != 0:
        raise ValueError("Exact R0-loop score outside its five-cell support")
    arrays = {"b":b,"h":h,"score":score,"score_numerator":score_numerator.astype(np.int32),
              "birth_k":birth_k.astype(np.uint16),"next_label":contact[...,ci["next_label"]].astype(np.uint16),
              "joint_safe":safe,"contact_e":es.astype(np.int8),"contact_c":cs.astype(np.int8),
              "r0_loop":loops.astype(np.int8),"rank_after":rank_after.astype(np.int8),
              "old_rank":old_rank.astype(np.int8),"rankcell":(3*old_rank[:,0]+old_rank[:,1]).astype(np.int8),
              "batch":np.full(1000,batch,dtype=np.int16),"counter":counter.astype(np.int64),
              "k0":k0.astype(np.int16),"d":d.astype(np.int16),
              "census_class_count":class_count[ids],"census_class_loop_sum":class_loop_sum[ids]}
    return arrays, {raw_path:hashlib.sha256(raw_bytes).hexdigest(),contact_path:hashlib.sha256(contact_bytes).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=ROOT/"results/p334-exact-score-quartet-moments")
    args = parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    started = time.perf_counter()
    result = {"schema":"p334.exact-score-quartet-moments.v1","source_commit":SOURCE,
        "contact_commit":CONTACT,"census_commit":CENSUS_COMMIT,
        "policy_reference":f"{CENSUS_ROOT}/analyze_finite_source.py@{CENSUS_COMMIT}",
        "reader_commit":subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
        "reader_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "new_samples":0,"new_DP":0,"raw_passes":1,"feature_order":FEATURES,
        "orientation_order":ORIENTATIONS,"mark_order":MARKS,"label_order":["U","V"],
        "suffix_order":[0,1],"birth_order":["K1","K2"],
        "axes":{"b":["prefix","quartet","orientation","feature"],
                "h":["prefix","quartet","mark","orientation","feature"],
                "score":["prefix","quartet","label","mark"],
                "score_numerator":["prefix","quartet","label","mark"],
                "birth_k":["prefix","quartet","label","suffix","orientation","birth"],
                "next_label":["prefix","quartet","label"],
                "joint_safe":["prefix","quartet","label"],
                "contact_e":["prefix","quartet","label","orientation"],
                "contact_c":["prefix","quartet","label","orientation"],
                "r0_loop":["prefix","quartet","label","orientation"],
                "rank_after":["prefix","quartet","label","orientation"],
                "old_rank":["prefix","orientation"],
                "census_class_count":["prefix","first_e","second_e"],
                "census_class_loop_sum":["prefix","first_e","second_e","orientation"]},
        "definitions":{"x":"K1/(N+1)","y":"K2/(N+1)",
            "xx_xy_yy":"x*x, x*y, y*y in a single saved suffix; not unconditional Beta second moments",
            "label_first_five":"Mean of two independent saved suffixes at that label",
            "cross_xx":"x0*x1 at one label",
            "cross_xy":"(x0*y1+x1*y0)/2 at one label",
            "cross_yy":"y0*y1 at one label",
            "b":"(f_U+f_V)/2: first five average all four tails; final three average the two label cross-products, no additional suffix mean",
            "h":"0.5*(H(U)-H(V))*(f_U-f_V)",
            "H":"joint-safe: pi_a*(g_plus/minus-mean_a g_plus/minus); otherwise zero",
            "g_plus_minus":"half sum/difference of R0-only (e-c) loop marks",
            "a":"joint contact-degree class (e_first,e_second) among joint rank-safe vacant labels",
            "pi_a":"census_class_count[a]/d; exact prefix-specific class mass",
            "score_numerator":"2*d*H; exact int32, zero outside joint-safe",
            "rankcell":"3*first_oldrank+second_oldrank, all nine cells retained; tangent support 0,1,2,3,6",
            "beta_recovery":"For either b or h: continuous xx=((N+1)*xx+x)/(N+2), xy=((N+1)*xy+x)/(N+2), yy=((N+1)*yy+y)/(N+2). Cross independent-suffix products unchanged.",
            "baseline_population":"Each original prefix has equal weight; each of eight quartets retained, no cell conditioning or prevalence renormalization"},
        "dependence":"Same frozen paired prefixes/labels/tails as e32/959; exact-score and matched-mask estimates share an estimand and are not independent evidence. Use distinct quartets for products estimating conditional means squared, retain original20 batch IDs.",
        "environment":{"python":platform.python_version(),"numpy":np.__version__,"machine":platform.machine()},
        "input_sha256":{},"sizes":{}}
    for n in (325,425):
        counters,count,loop_sum,k0,hashes = census(n)
        result["input_sha256"].update(hashes)
        chunks = []
        for batch in range(20):
            arrays,hashes = read_batch(n,batch,counters,count,loop_sum,k0)
            chunks.append(arrays)
            result["input_sha256"].update(hashes)
        arrays = {key:np.concatenate([a[key] for a in chunks],axis=0) for key in chunks[0]}
        path = args.output/f"N{n}.npz"
        np.savez_compressed(path,**arrays)
        result["sizes"][str(n)] = {"file":path.name,"sha256":hashlib.sha256(path.read_bytes()).hexdigest(),
            "file_bytes":path.stat().st_size,"prefixes":20000,"quartets_per_prefix":8,
            "delta_cos4":DELTA[n],"rankcell_counts":np.bincount(arrays["rankcell"],minlength=9).tolist(),
            "arrays":{key:{"shape":list(value.shape),"dtype":str(value.dtype)} for key,value in arrays.items()},
            "nonzero_score_label_counts":np.count_nonzero(arrays["score"],axis=(0,1,2)).tolist()}
        print(f"N{n}: {path.stat().st_size/1e6:.3f} MB; 20000 prefixes x 8 quartets; source pass complete",flush=True)
        del arrays,chunks,count,loop_sum
    result["elapsed_seconds"] = time.perf_counter()-started
    (args.output/"metadata.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(f"Complete in {result['elapsed_seconds']:.3f}s",flush=True)


if __name__ == "__main__":
    main()

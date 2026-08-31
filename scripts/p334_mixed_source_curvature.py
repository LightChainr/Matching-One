#!/usr/bin/env python3
"""Exact ff/fs/ss density-score responses, from old8 NPZ and one New64 gzip pass."""
import argparse
import hashlib
import io
import json
import subprocess
import time
from pathlib import Path

import numpy as np
from scipy.stats import binom

from p334_safe_contact_response import ROOT, blob, array_csv

OLD = "375cd3a12b2b7a87d79148a59f62b95898f9e471"
CLASS = "1cfa4ae892a2f7f4168e9a71690efd7a5560d4cd"
NEW = "8ad30617b0a3076a5c01a208eb213096d8879b32"
EXT = "experiments/p334-mechanism-response-20260831/extension"
P_REF = 0.59274605079
ORIS = ("first","second")
OBS = ("A_ref","E_ref","C","W")
COMP = ("ff","fs","ss")
PAIRS = ((0,0),(0,1),(1,1))


def load_old(n):
    paths = (f"results/p334-exact-score-quartet-moments/N{n}.npz",
             f"results/p334-exact-prefix-structure/N{n}.npz")
    bb = (blob(OLD,paths[0]),blob(CLASS,paths[1]))
    with np.load(io.BytesIO(bb[0]),allow_pickle=False) as z, np.load(io.BytesIO(bb[1]),allow_pickle=False) as c:
        if not np.array_equal(z["counter"],c["counter"]) or not np.array_equal(z["batch"],c["batch"]):
            raise ValueError("Old8 and exact class products disagree on prefix identity")
        old = {k:z[k].copy() for k in ("counter","batch","rankcell","old_rank","k0","d","birth_k","contact_e","contact_c","rank_after")}
        classes = {k:c[k].copy() for k in ("class_count","class_loop_sum","class_loop_product_sum")}
    return old,classes,{f"{commit}:{path}":hashlib.sha256(data).hexdigest() for commit,path,data in zip((OLD,CLASS),paths,bb)}


def second_score(ids,old,classes,es,cs,after):
    ranks = old["old_rank"][ids]
    safe = np.all(after == ranks[:,None,None,:],axis=-1)
    ls = (es.astype(np.int64)-cs)*(ranks[:,None,None,:] == 0)
    ix = (np.broadcast_to(ids[:,None,None],safe.shape),es[...,0],es[...,1])
    count = classes["class_count"][ix]
    sums = classes["class_loop_sum"][ix]
    products = classes["class_loop_product_sum"][ix]
    centered_num = count[...,None]*ls-sums
    numerator = np.stack([centered_num[...,i]*centered_num[...,j]-(count*products[...,k]-sums[...,i]*sums[...,j])
                          for k,(i,j) in enumerate(PAIRS)],axis=-1)
    numerator *= safe[...,None]
    if np.count_nonzero(numerator[old["rankcell"][ids] != 0,...,1]) != 0:
        raise ValueError("Mixed exact density score must vanish outside original00")
    d = old["d"][ids].astype(np.int64)
    return numerator/(d*d)[:,None,None,None]


def means(n,birth_k,score):
    k1,k2 = birth_k[...,0],birth_k[...,1]
    tail = binom.sf(np.arange(n+1)-1,n,P_REF)
    f = np.stack((tail[k1]+tail[k2]-1,1-tail[k1]+tail[k2],
                  (k1.astype(np.float64)+k2)/(2*(n+1)),(k2.astype(np.float64)-k1)/(n+1)),axis=-1).mean(axis=3)
    # f: prefix,quartet,label,receiver,observable. No extra mixed factor.
    dy = f[:,:,0]-f[:,:,1]
    dt = score[:,:,0]-score[:,:,1]
    response = dy[...,None]*dt[:,:,None,None,:]/2
    return f.mean(axis=(1,2)),response.mean(axis=1)


def read_new_batch(n,batch,old,classes):
    ids = np.flatnonzero((old["batch"] == batch)&(old["rankcell"] == 0))
    path = f"{EXT}/N{n}.batch{batch:02}.csv.gz"
    data = blob(NEW,path)
    header,raw = array_csv(data)
    col = {s:i for i,s in enumerate(header)}
    raw = raw[np.lexsort(tuple(raw[:,col[k]] for k in ("replica","group","quartet","counter")))].reshape(len(ids),64,2,2,-1)
    if not np.array_equal(raw[:,0,0,0,col["counter"]],old["counter"][ids]):
        raise ValueError("New64 does not match the original00 prefix order")
    if not np.array_equal(raw[0,:,0,0,col["quartet"]],np.arange(8,72)):
        raise ValueError("Unexpected frozen New64 quartet domain")
    es = np.stack([raw[...,0,col[f"{o}_e"]] for o in ORIS],axis=-1)
    cs = np.stack([raw[...,0,col[f"{o}_c"]] for o in ORIS],axis=-1)
    after = np.stack([raw[...,0,col[f"{o}_next_rank"]] for o in ORIS],axis=-1)
    birth = np.stack([np.stack([raw[...,col[f"{o}_{k}"]] for k in ("k1","k2")],axis=-1) for o in ORIS],axis=-2)
    score = second_score(ids,old,classes,es,cs,after)
    baseline,response = means(n,birth,score)
    return ids,baseline,response,{f"{NEW}:{path}":hashlib.sha256(data).hexdigest()}


def batch_result(streams):
    fields = {}
    for stream,z in streams.items():
        groups = ["all",*[f"{a}{b}" for a in range(3) for b in range(3)]] if stream == "old8" else ["00"]
        for group in groups:
            keep = np.ones(len(z["counter"]),dtype=bool) if group == "all" else z["rankcell"] == 3*int(group[0])+int(group[1])
            for oi,ori in enumerate(ORIS):
                for fi,observable in enumerate(OBS):
                    for ci,component in enumerate(COMP):
                        value = z["mean_response2"][:,oi,fi,ci]
                        fields[f"{stream}.{group}.{ori}.{observable}.{component}"] = np.array([value[keep&(z["batch"] == b)].sum()/1000 for b in range(20)])
    raw = np.column_stack(list(fields.values()))
    point = raw.mean(axis=0)
    factor = (raw-point)/np.sqrt(20*19)
    return {"batch_ids":list(range(20)),"labels":list(fields),"raw_batch_means":raw.tolist(),
            "estimate":point.tolist(),"se":np.linalg.norm(factor,axis=0).tolist(),
            "factor":factor.tolist(),"mean_covariance":(factor.T@factor).tolist(),
            "LOO":((20*point-raw)/19).tolist(),"LOO_factor":(-factor).tolist()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=ROOT/"results/p334-mixed-source-curvature")
    args = parser.parse_args(); args.output.mkdir(parents=True,exist_ok=False)
    started = time.perf_counter()
    result = {"schema":"p334.mixed-source-curvature.v1","allocation_commit":"6bace935",
        "old8_commit":OLD,"class_products_commit":CLASS,"new64_commit":NEW,"p_ref":P_REF,
        "reader_commit":subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
        "reader_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "new_samples":0,"new_DP":0,"old_fork_gzip_reads":0,"new64_gzip_passes":1,
        "finite_weight_evaluations":0,"determinant_evaluations":0,
        "orientation_order":ORIS,"observable_order":OBS,"component_order":COMP,
        "prefix_axes":{"mean_observable":["prefix","receiver","observable"],
                       "mean_response2":["prefix","receiver","observable","component"]},
        "definitions":{"A_ref":"Pr(Binomial(N,p_ref)>=K1)+Pr(Binomial(N,p_ref)>=K2)-1",
            "E_ref":"1-Pr(Binomial(N,p_ref)>=K1)+Pr(Binomial(N,p_ref)>=K2)",
            "C":"(K1+K2)/(2*(N+1))","W":"(K2-K1)/(N+1)",
            "second_density_score":"[(n*L_i-S_i)*(n*L_j-S_j)-(n*Q_ij-S_i*S_j)]/d^2, zero outside joint-safe",
            "class_symbols":"n=joint-safe degree-class count; S_i=sum_class L_i; Q_ij=sum_class L_i L_j; d=N-k0",
            "L":"R0-only e-c; source coordinates are physical first and second, without plus/minus half factors",
            "paired_response":"0.5*(t_ij(U)-t_ij(V))*(F_U-F_V), F_label=two-suffix average",
            "mixed_normalization":"fs is the actual d_f d_s response, no extra factor2; ff and ss are actual pure second derivatives",
            "prefix_means":"Average8 old quartets or64 new quartets separately; identities and original batches retained",
            "population":"old8 all9 cells; new64 only original00. Every batch contribution divided by original1000; new64.00 is not renormalized by00 prevalence",
            "factor_convention":"factor=(batch-mean)/sqrt(20*19); LOO_factor=-factor. mean_covariance=factor.T@factor; align signs when joining other delete-one results"},
        "interpretation":"Nonzero mixed response is nonadditivity in the specified commuting exponential source coordinates, not path memory, noncommuting source actions or field identification.",
        "dependence":"Old8 and New64 share original prefixes and20batch blocks; streams kept separate, no independent population replication.",
        "input_sha256":{},"sizes":{}}
    for n in (325,425):
        old,classes,hashes = load_old(n)
        result["input_sha256"].update(hashes)
        ids = np.arange(20000)
        score2 = second_score(ids,old,classes,old["contact_e"],old["contact_c"],old["rank_after"])
        baseline,response = means(n,old["birth_k"],score2)
        identity = {k:old[k] for k in ("counter","batch","rankcell","old_rank","k0","d")}
        streams = {"old8":{**identity,"mean_observable":baseline,"mean_response2":response}}
        ni,nb,nr = [],[],[]
        for batch in range(20):
            ids,b,r,hashes = read_new_batch(n,batch,old,classes)
            ni.append(ids);nb.append(b);nr.append(r);result["input_sha256"].update(hashes)
        ids = np.concatenate(ni)
        streams["new64"] = {**{k:v[ids] for k,v in identity.items()},
                            "mean_observable":np.concatenate(nb),"mean_response2":np.concatenate(nr)}
        files = {}
        for stream,z in streams.items():
            path = args.output/f"{stream}_N{n}.npz"
            np.savez_compressed(path,**z)
            files[stream] = {"file":path.name,"prefixes":len(z["counter"]),"quartets_per_prefix":8 if stream == "old8" else 64,
                "sha256":hashlib.sha256(path.read_bytes()).hexdigest(),"file_bytes":path.stat().st_size}
        summary = batch_result(streams)
        summary["prefix_files"] = files
        result["sizes"][str(n)] = summary
        for name,value,se in zip(summary["labels"],summary["estimate"],summary["se"]):
            if name.startswith("new64."):
                print(n,name,f"{value:.11g} +/- {se:.6g}",flush=True)
    result["elapsed_seconds"] = time.perf_counter()-started
    (args.output/"score.json").write_text(json.dumps(result,separators=(",",":"),allow_nan=False)+"\n")
    print(f"Complete in {result['elapsed_seconds']:.3f}s",flush=True)


if __name__ == "__main__":
    main()

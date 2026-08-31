#!/usr/bin/env python3
"""Exact census descriptors aligned to the saved prefix-quartet moment rows."""
import argparse
import hashlib
import io
import json
import platform
import subprocess
import time
from pathlib import Path

import numpy as np

from p334_safe_contact_response import ROOT, blob, array_csv

PREFIX_COMMIT = "375cd3a12b2b7a87d79148a59f62b95898f9e471"
CENSUS_COMMIT = "ac5761ce504c3cd170fa42c86c17d6fb87f0375b"
FEATURES = ["joint_safe_mass", "class_collision", "safe_degree_first", "safe_degree_second",
            "safe_loop_first", "safe_loop_second", "score_gram_ff", "score_gram_fs", "score_gram_ss",
            "score_energy_plus", "score_energy_minus", "score_gram_plus_minus"]
MOMENT_PAIRS = [(0,0),(0,1),(1,1)]


def read_size(n):
    prefix_path = f"results/p334-exact-score-quartet-moments/N{n}.npz"
    census_path = f"experiments/p334-finite-source-20260831/census/N{n}/census.csv.gz"
    prefix_bytes, census_bytes = blob(PREFIX_COMMIT,prefix_path), blob(CENSUS_COMMIT,census_path)
    with np.load(io.BytesIO(prefix_bytes),allow_pickle=False) as source:
        identity = {k:source[k].copy() for k in ("counter","batch","rankcell","old_rank","k0","d")}
    names, rows = array_csv(census_bytes)
    col = {name:i for i,name in enumerate(names)}
    counters = identity["counter"]
    order = np.argsort(counters)
    ids = order[np.searchsorted(counters[order],rows[:,col["counter"]])]
    for name in ("counter","batch","k0"):
        if not np.array_equal(identity[name][ids],rows[:,col[name]]):
            raise ValueError(f"Prefix/census key disagreement N{n}/{name}")
    if not np.array_equal(identity["rankcell"][ids],3*rows[:,col["first_rank"]]+rows[:,col["second_rank"]]):
        raise ValueError("Census rankcell differs from fixed prefix identity")
    # Preserve complete class sufficient counts, with every sum integral.
    count = np.zeros((20000,5,5),dtype=np.int64)
    loop_sum = np.zeros((20000,5,5,2),dtype=np.int64)
    loop_product_sum = np.zeros((20000,5,5,3),dtype=np.int64)
    ix = (ids,rows[:,col["first_e"]],rows[:,col["second_e"]])
    weight = rows[:,col["count"]]
    loops = rows[:,[col["L_first"],col["L_second"]]]
    np.add.at(count,ix,weight)
    for i in range(2):
        np.add.at(loop_sum[...,i],ix,weight*loops[:,i])
    for j,(a,b) in enumerate(MOMENT_PAIRS):
        np.add.at(loop_product_sum[...,j],ix,weight*loops[:,a]*loops[:,b])
    # E_uniform[s_i s_j] = sum_a [n_a² Q_ij - n_a S_i S_j]/d³.
    class_gram_num = np.stack([count*count*loop_product_sum[...,j]-count*loop_sum[...,a]*loop_sum[...,b]
                               for j,(a,b) in enumerate(MOMENT_PAIRS)],axis=-1)
    gram_num = class_gram_num.sum(axis=(1,2))
    d = identity["d"].astype(np.int64)
    d3 = d*d*d
    gram = gram_num/d3[:,None]
    pm_num = np.stack((gram_num[:,0]+2*gram_num[:,1]+gram_num[:,2],
                       gram_num[:,0]-gram_num[:,2],
                       gram_num[:,0]-2*gram_num[:,1]+gram_num[:,2]),axis=-1)
    pm = pm_num/(4*d3[:,None])
    safe_count = count.sum(axis=(1,2))
    collision_num = (count*count).sum(axis=(1,2))
    degree_num = np.stack(((count*np.arange(5)[None,:,None]).sum(axis=(1,2)),
                           (count*np.arange(5)[None,None,:]).sum(axis=(1,2))),axis=-1)
    loop_num = loop_sum.sum(axis=(1,2))
    features = np.column_stack((safe_count/d,collision_num/(d*d),degree_num/d[:,None],
                                loop_num/d[:,None],gram,pm[:,0],pm[:,2],pm[:,1]))
    physical_matrix = np.stack((gram[:,0],gram[:,1],gram[:,1],gram[:,2]),axis=-1).reshape(20000,2,2)
    pm_matrix = np.stack((pm[:,0],pm[:,1],pm[:,1],pm[:,2]),axis=-1).reshape(20000,2,2)
    loop_mean = np.divide(loop_sum,count[...,None],out=np.zeros_like(loop_sum,dtype=np.float64),where=count[...,None]>0)
    arrays = {**identity,"features":features,"score_gram_physical":physical_matrix,
        "score_gram_plus_minus":pm_matrix,"class_count":count,"class_loop_sum":loop_sum,
        "class_loop_product_sum":loop_product_sum,"class_mass":count/d[:,None,None],
        "class_loop_mean":loop_mean,"class_score_gram_numerator":class_gram_num,
        "score_gram_numerator":gram_num,"score_gram_denominator":d3,
        "plus_minus_gram_numerator":pm_num,"plus_minus_gram_denominator":4*d3,
        "safe_count":safe_count,"class_collision_numerator":collision_num,
        "safe_degree_numerator":degree_num,"safe_loop_numerator":loop_num}
    return arrays,{prefix_path:hashlib.sha256(prefix_bytes).hexdigest(),census_path:hashlib.sha256(census_bytes).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=ROOT/"results/p334-exact-prefix-structure")
    args = parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=False)
    started = time.perf_counter()
    meta = {"schema":"p334.exact-prefix-structure.v1","prefix_commit":PREFIX_COMMIT,"census_commit":CENSUS_COMMIT,
        "reader_commit":subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
        "reader_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "new_samples":0,"new_DP":0,"fork_gzip_reads":0,"recomputed_quartet_responses":0,
        "feature_order":FEATURES,"orientation_order":["first","second"],"mark_order":["plus","minus"],
        "physical_gram_component_order":["ff","fs","ss"],"plus_minus_gram_component_order":["++","+-","--"],
        "feature_blocks":{"contact_structure":FEATURES[:6],"physical_source_gram":FEATURES[6:9],
                          "derived_plus_minus_gram":FEATURES[9:]},
        "definitions":{"class":"Joint-safe contact degrees (e_first,e_second), each degree 0..4",
            "d":"N-k0, all vacant labels, including labels outside joint safety",
            "L":"R0-only e-c mark; non-R0 orientation has L=0. Census contains only joint-safe labels.",
            "joint_safe_mass":"sum_a n_a/d", "class_collision":"sum_a (n_a/d)^2; excludes outside-safe category",
            "safe_degree_first_second":"sum_a n_a*e_first_second/d; not divided by joint_safe_mass",
            "safe_loop_first_second":"sum_a S_first_second[a]/d; not divided by joint_safe_mass",
            "S_i_Q_ij":"S_i=sum_in_class L_i; Q_ij=sum_in_class L_i L_j, from full label census counts",
            "s_Li":"pi_a*(L_i-S_i/n_a), outside joint-safe zero, pi_a=n_a/d",
            "physical_gram":"E_uniform(s_Li*s_Lj)=sum_a(n_a^2 Q_ij-n_a S_i S_j)/d^3",
            "plus_minus_scores":"s_plus=(s_first+s_second)/2; s_minus=(s_first-s_second)/2",
            "class_score_gram_numerator":"Unconditional per-class contribution numerator (ff,fs,ss), denominator d^3; not the within-class conditional Gram",
            "class_loop_mean":"S_i/n_a, set to zero for empty classes; class_mass=0 identifies them",
            "redundancy":"Last three feature columns are exact transforms of physical Gram columns 6:9, not additional independent predictors",
            "population":"All20000 original prefixes in the exact NPZ row order, nine rank cells preserved; no finite-quartet score or descriptor estimation"},
        "environment":{"python":platform.python_version(),"numpy":np.__version__,"machine":platform.machine()},
        "input_sha256":{},"sizes":{}}
    for n in (325,425):
        arrays,hashes = read_size(n)
        meta["input_sha256"].update(hashes)
        path = args.output/f"N{n}.npz"
        np.savez_compressed(path,**arrays)
        meta["sizes"][str(n)] = {"file":path.name,"file_bytes":path.stat().st_size,
            "sha256":hashlib.sha256(path.read_bytes()).hexdigest(),
            "arrays":{key:{"shape":list(a.shape),"dtype":str(a.dtype)} for key,a in arrays.items()},
            "rankcell_counts":np.bincount(arrays["rankcell"],minlength=9).tolist()}
        print(f"N{n}: exact descriptors {arrays['features'].shape}, {path.stat().st_size/1e6:.3f} MB",flush=True)
    meta["elapsed_seconds"] = time.perf_counter()-started
    (args.output/"metadata.json").write_text(json.dumps(meta,indent=2,allow_nan=False)+"\n")
    print(f"Complete in {meta['elapsed_seconds']:.3f}s",flush=True)


if __name__ == "__main__":
    main()

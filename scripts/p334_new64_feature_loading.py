#!/usr/bin/env python3
"""Original00 feature loading with old8 clocks and independent New64 tangents."""
import argparse
import hashlib
import io
import json
import subprocess
import time
from pathlib import Path

import numpy as np

from p334_safe_contact_response import ROOT, blob

OLD = "375cd3a12b2b7a87d79148a59f62b95898f9e471"
DESC = "1cfa4ae892a2f7f4168e9a71690efd7a5560d4cd"
NEW = "8ad30617b0a3076a5c01a208eb213096d8879b32"
MODEL = "011f50e3"
ORIS = ("first", "second")
STREAMS = ("old8", "new64")
PRED = ("joint_safe_mass", "own_score_energy", "own_safe_degree", "own_safe_loop", "mu_C", "mu_W")
RESP = ("source_first.C", "source_first.W", "source_second.C", "source_second.W")
MODELS = {"strength":(1,), "contact":(0,1,2,3)}
TRI = np.triu_indices(6)


def load(n):
    specs = [(OLD,f"results/p334-exact-score-quartet-moments/N{n}.npz"),
             (DESC,f"results/p334-exact-prefix-structure/N{n}.npz"),
             (NEW,f"experiments/p334-mechanism-response-20260831/results-extension/prefix_statistics_N{n}.npz")]
    bytearrays = [blob(commit,path) for commit,path in specs]
    with np.load(io.BytesIO(bytearrays[0]),allow_pickle=False) as z, np.load(io.BytesIO(bytearrays[1]),allow_pickle=False) as f, np.load(io.BytesIO(bytearrays[2]),allow_pickle=False) as new:
        if not np.array_equal(z["counter"],f["counter"]) or not np.array_equal(z["batch"],f["batch"]):
            raise ValueError("Exact features and original moments have different prefix identities")
        take = z["rankcell"] == 0
        counter,batch = z["counter"][take],z["batch"][take]
        if not np.array_equal(counter,new["counter"]) or not np.array_equal(batch,new["batch"]):
            raise ValueError("New64 must be exactly the same original00 counters/batches")
        b,h = z["b"][take],z["h"][take]
        feature = np.stack((f["features"][take][:,[0,6,2,4]],f["features"][take][:,[0,8,3,5]]),axis=1)
        clock = np.stack(((b[...,0]+b[...,1])/2,b[...,1]-b[...,0]),axis=-1)
        x = np.concatenate((np.broadcast_to(feature[:,None],(*clock.shape[:-1],4)),clock),axis=-1)
        # Physical sources, not the plus/minus half-sum coordinates.
        physical = np.stack((h[:,:,0]+h[:,:,1],h[:,:,0]-h[:,:,1]),axis=2)
        response = np.stack(((physical[...,0]+physical[...,1])/2,physical[...,1]-physical[...,0]),axis=-1)
        old_y = response.transpose(0,1,3,2,4).reshape(len(b),8,2,4)
        names,values = list(new["labels"]),new["new64"]
        new_y = np.empty((len(b),2,4))
        for oi,ori in enumerate(ORIS):
            for si,source in enumerate(ORIS):
                for ti,(observable,scale) in enumerate((("A",-.5),("E",-1.))):
                    key = f"p_integral.{observable}.mean_J[{ori},{source}]"
                    new_y[:,oi,2*si+ti] = values[:,names.index(key)]*scale
    hashes = {f"{commit}:{path}":hashlib.sha256(data).hexdigest() for (commit,path),data in zip(specs,bytearrays)}
    return x,old_y,new_y,counter,batch,hashes


def sufficient_batches(x,old_y,new_y,batch):
    q = x.shape[1]
    sx,sy = x.sum(axis=1),old_y.sum(axis=1)
    mx,my = sx/q,sy/q
    ku = (np.einsum('noi,noj->noij',sx,sx)-np.einsum('nqoi,nqoj->noij',x,x))/(q*(q-1))
    kd = np.einsum('noi,noj->noij',mx,mx)
    vu_old = (np.einsum('noi,noj->noij',sx,sy)-np.einsum('nqoi,nqoj->noij',x,old_y))/(q*(q-1))
    vd_old = np.einsum('noi,noj->noij',mx,my)
    cross = np.einsum('noi,noj->noij',mx,new_y)
    fields = {"cell00.mass":np.ones(len(x))}
    for oi,ori in enumerate(ORIS):
        for j,name in enumerate(PRED): fields[f"{ori}.meanX.{name}"] = mx[:,oi,j]
        for i,j in zip(*TRI):
            for kind,matrix in (("U",ku),("diag",kd)):
                fields[f"{ori}.K.{kind}.{i}.{j}"] = matrix[:,oi,i,j]
        for stream,y,vu,vd in (("old8",my,vu_old,vd_old),("new64",new_y,cross,cross)):
            for j,name in enumerate(RESP): fields[f"{ori}.{stream}.meanY.{name}"] = y[:,oi,j]
            for i in range(6):
                for j in range(4):
                    fields[f"{ori}.{stream}.V.U.{i}.{j}"] = vu[:,oi,i,j]
                    fields[f"{ori}.{stream}.V.diag.{i}.{j}"] = vd[:,oi,i,j]
    a = np.column_stack(list(fields.values()))
    ids = np.arange(20)
    # Explicit zero padding to all1000 original prefixes of each batch.
    raw = np.array([a[batch==j].sum(axis=0)/1000 for j in ids])
    counts = np.array([(batch==j).sum() for j in ids])
    return list(fields),raw,counts


def solve(k,v):
    scale = np.sqrt(np.diag(k))
    if not np.all(np.isfinite(scale)) or np.any(scale <= 0):
        raise ValueError("Named predictor variance not identified; no repair")
    cor = k/scale[:,None]/scale[None,:]
    eigen = np.linalg.eigvalsh(cor)
    if eigen[0] <= 0: raise ValueError("Named latent Gram not positive definite; no ridge")
    return np.linalg.solve(cor,v/scale[:,None])/scale[:,None],eigen[0],eigen[-1]/eigen[0]


def derive(row,labels,population):
    d,out = dict(zip(labels,row)),{}
    pi = d["cell00.mass"]
    m = int(round(population*pi))
    out["cell00.mass"] = pi
    for oi,ori in enumerate(ORIS):
        mx = np.array([d[f"{ori}.meanX.{name}"]/pi for name in PRED])
        matrices = {}
        for kind in ("U","diag"):
            k = np.zeros((6,6))
            for i,j in zip(*TRI): k[i,j] = k[j,i] = d[f"{ori}.K.{kind}.{i}.{j}"]
            matrices[kind] = k
        k = matrices["U"]-(m*pi*np.outer(mx,mx)-matrices["diag"])/(m-1)
        for i,j in zip(*TRI): out[f"{ori}.K.{PRED[i]}|{PRED[j]}"] = k[i,j]
        for stream in STREAMS:
            stem = f"{ori}.{stream}."
            my = np.array([d[f"{ori}.{stream}.meanY.{name}"]/pi for name in RESP])
            vu,vd = [np.array([[d[f"{ori}.{stream}.V.{kind}.{i}.{j}"] for j in range(4)] for i in range(6)]) for kind in ("U","diag")]
            v = vu-(m*pi*np.outer(mx,my)-vd)/(m-1)
            for i,pred in enumerate(PRED):
                for j,response in enumerate(RESP): out[stem+f"cov.{pred}|{response}"] = v[i,j]
            own = (2*oi,2*oi+1)
            observed = 2*v[4,own[0]]-.5*v[5,own[1]]
            out[stem+"own_signed_loading"] = observed
            for model,cols in MODELS.items():
                ix = np.array(cols)
                beta,eig,condition = solve(k[np.ix_(ix,ix)],v[ix])
                loading = k[4:6,ix]@beta
                captured = 2*loading[0,own[0]]-.5*loading[1,own[1]]
                out[stem+model+".captured_signed_loading"] = captured
                out[stem+model+".residual_signed_loading"] = observed-captured
                out[stem+model+".signed_loading_share"] = captured/observed
                out[stem+model+".scaled_Gram_min_eigenvalue"] = eig
                out[stem+model+".scaled_Gram_condition"] = condition
                for ii,col in enumerate(cols):
                    for j,response in enumerate(RESP): out[stem+model+f".beta.{PRED[col]}|{response}"] = beta[ii,j]
            beta_clock,eig,condition = solve(k[4:6,4:6],v[4:6])
            partial = v[:4]-k[:4,4:6]@beta_clock
            out[stem+"clock.scaled_Gram_min_eigenvalue"] = eig
            out[stem+"clock.scaled_Gram_condition"] = condition
            for i,pred in enumerate(PRED[:4]):
                for j,response in enumerate(RESP): out[stem+f"clock_partial_cov.{pred}|{response}"] = partial[i,j]
        # Paired stream contrast; shared prefixes, descriptors and old clocks.
        suffixes = ("own_signed_loading",*[f"{model}.{item}" for model in MODELS for item in
                    ("captured_signed_loading","residual_signed_loading","signed_loading_share")],
                    *[f"clock_partial_cov.{pred}|{response}" for pred in PRED[:4] for response in RESP])
        for suffix in suffixes:
            out[f"{ori}.new64_minus_old8.{suffix}"] = out[f"{ori}.new64.{suffix}"]-out[f"{ori}.old8.{suffix}"]
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=ROOT/"results/p334-new64-feature-loading")
    args = parser.parse_args(); args.output.mkdir(parents=True,exist_ok=False)
    started = time.perf_counter()
    result = {"schema":"p334.new64-feature-loading.v1","old8_commit":OLD,"descriptor_commit":DESC,
        "new64_commit":NEW,"frozen_feature_model_commit":MODEL,"predictors":PRED,"responses":RESP,"models":MODELS,
        "reader_commit":subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
        "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "new_samples":0,"fork_gzip_reads":0,"determinant_evaluations":0,
        "estimand":"Original00 within-cell latent prefix covariance, weighted by original00 prevalence in all20000 original prefixes. Separate old8 and new64 tangents, same four-feature and source-energy-only model definitions; coefficients are estimated separately, not treated as frozen out-of-sample predictions.",
        "cross_stream":"Old8 clock products use ordered distinct quartets. old8 clock means x new64 response means are conditionally independent-stream products; no same-quartet subtraction is applied to this cross-stream moment.",
        "units":"New64 integralA response /(-2)=H_C; integralE response /(-1)=H_W. Both already normalized by N+1. Physical sources first/second, physical receivers first/second.",
        "signed_loading":"2 Cov(mu_C,H_C_own)-0.5 Cov(mu_W,H_W_own). Captured loading substitutes the named linear predictor; its share is signed, not explained variance or a closure probability.",
        "uncertainty":"Fulloriginal20batch delete-one refits with original population20000/19000. All stream/source/receiver/contact directions share one factor; no response variance claim or independent population replication.",
        "input_sha256":{},"sizes":{}}
    for n in (325,425):
        x,old_y,new_y,counter,batch,hashes = load(n)
        labels,raw,counts = sufficient_batches(x,old_y,new_y,batch)
        mean = raw.mean(axis=0)
        point = derive(mean,labels,20000)
        loo = np.array([list(derive((20*mean-r)/19,labels,19000).values()) for r in raw])
        factor = np.sqrt(19/20)*(loo-loo.mean(axis=0))
        se = np.linalg.norm(factor,axis=0)
        result["input_sha256"].update(hashes)
        result["sizes"][str(n)] = {"original00_prefixes":len(counter),"original00_batch_counts":counts.tolist(),
            "batch_ids":list(range(20)),"raw_labels":labels,"raw_batch_means":raw.tolist(),
            "labels":list(point),"estimate":list(point.values()),"se":se.tolist(),"LOO":loo.tolist(),"factor":factor.tolist()}
        for name,value,error in zip(point,point.values(),se):
            if name.endswith(("own_signed_loading","signed_loading_share","residual_signed_loading")):
                print(n,name,f"{value:.11g} +/- {error:.6g}",flush=True)
    result["elapsed_seconds"] = time.perf_counter()-started
    (args.output/"score.json").write_text(json.dumps(result,separators=(",",":"),allow_nan=False)+"\n")
    print(f"Complete in {result['elapsed_seconds']:.3f}s",flush=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Exact thermal dipole/quadrupole of the existing common-label tangent."""
import json
from pathlib import Path
import subprocess

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
SOURCE="4db356e1b026853468f94d59d938895a2367ceb7"
HIST_PATH="results/p334-common-label-euler-tangent/signed_birth_histograms.json"


def rising(x,q):
    out=np.ones_like(x,dtype=np.float64)
    for j in range(q):
        out*=x+j
    return out


source=json.loads(subprocess.check_output(["git","show",SOURCE+":"+HIST_PATH],cwd=ROOT))
result={"source_commit":SOURCE,"source_path":HIST_PATH,"new_samples":0,
    "p_ref":source["p_ref"],
    "identity":"For q>=1, E[tau_K^q|K]=rising(K,q)/rising(N+1,q). Tangent integral p^m E=(H_tau1^(m+1)-H_tau2^(m+1))/(m+1); A uses negative sum.",
    "scope":"Linear exact functional of existing integer histograms; same20 original paired batches. Centering uses fixed p_ref, not a data-selected root or divided small integral.",
    "sizes":{}}
for ns,s in source["sizes"].items():
    n=int(ns)
    hist=np.asarray(s["batch_integer_histograms"],dtype=np.int64)
    denom=source["batch_denominator"]
    k=np.arange(n+1,dtype=np.float64)
    # Axes batch,cell,mark,orientation,birth,power(q=1,2,3).
    tau=np.stack([hist @ (rising(k,q)/rising(np.array(float(n+1)),q))/denom for q in (1,2,3)],axis=-1)
    cols={}
    groups=[("all",list(range(5))), *[(c,[i]) for i,c in enumerate(source["cell_order"])]]
    for group,ix in groups:
        h=tau[:,ix].sum(1)
        for im,mark in enumerate(source["mark_order"]):
            for io,orientation in enumerate(source["orientation_order"]):
                for ib,birth in enumerate(("tau1","tau2")):
                    for iq in range(3):
                        cols[f"{group}.{mark}->{orientation}.{birth}_power{iq+1}"]=h[:,im,io,ib,iq]
                h1,h2=h[:,im,io,0],h[:,im,io,1]
                for name,val in (("C",(h1[:,0]+h2[:,0])/2),("W",h2[:,0]-h1[:,0]),("CW",(h2[:,1]-h1[:,1])/2)):
                    cols[f"{group}.{mark}->{orientation}.{name}"]=val
            for channel in ("S","D"):
                w=np.array([.5,.5]) if channel=="S" else np.array([1.,-1.])/s["delta_cos4"]
                hj=np.einsum("bojq,o->bjq",h[:,im],w)
                for obs in ("A","E"):
                    moments=-(hj[:,0]+hj[:,1])/np.arange(1,4) if obs=="A" else (hj[:,0]-hj[:,1])/np.arange(1,4)
                    base=f"{group}.{mark}->{channel}.{obs}."
                    for m in range(3):
                        cols[base+f"I{m}"]=moments[:,m]
                    p=source["p_ref"]
                    cols[base+"dipole_at_p_ref"]=moments[:,1]-p*moments[:,0]
                    cols[base+"quadrupole_at_p_ref"]=moments[:,2]-2*p*moments[:,1]+p*p*moments[:,0]
    x=np.column_stack(list(cols.values()))
    mean=x.mean(0)
    factor=(x-mean)/np.sqrt(20*19)
    se=np.linalg.norm(factor,axis=0)
    result["sizes"][ns]={"batch_ids":s["batch_ids"],"delta_cos4":s["delta_cos4"],
        "labels":list(cols),"joint_20_batch_means":x.tolist(),"estimate":mean.tolist(),
        "se":se.tolist(),"factor":factor.tolist()}
    for name,value,error in zip(cols,mean,se):
        if name.startswith("all.") and ".E." in name:
            print(ns,name,f"{value:.11g} +/- {error:.6g}")
target=ROOT/"results/p334-euler-thermal-moments"
target.mkdir(exist_ok=False)
(target/"score.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")

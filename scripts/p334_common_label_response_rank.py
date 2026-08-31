#!/usr/bin/env python3
"""Two-control response rank, using only the saved common-label batch vectors."""
import json
from pathlib import Path
import subprocess

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
SOURCE="4db356e1b026853468f94d59d938895a2367ceb7"
PATH="results/p334-common-label-euler-tangent/score.json"


def readout(v,names,delta):
    d=dict(zip(names,v))
    result={}
    for group,cells in (("all",["all"]),("00",["00"]),("01+10",["01","10"])):
        get=lambda suffix: sum(d[c+"."+suffix] for c in cells)
        gp,gm,gpm=(get(f"GG[{a},{b}]") for a,b in (("plus","plus"),("minus","minus"),("plus","minus")))
        gff,gss,gfs=gp+gm+2*gpm,gp+gm-2*gpm,gp-gm
        result[group+".input_G.det"]=gff*gss-gfs*gfs
        result[group+".input_G.correlation"]=gfs/np.sqrt(gff*gss)
        for endpoint in ("p_ref","p_integral"):
            for obs in ("F1","F2","A","E"):
                response=lambda mark,out: get(f"{mark}->{out}.{endpoint}.{obs}")
                sp,sm,dp,dm=(response(m,o) for m,o in (("plus","S"),("minus","S"),("plus","D"),("minus","D")))
                # Rows=future first/second geometry, columns=input L_first/L_second.
                j=np.array([[sp+sm+delta*(dp+dm)/2,sp-sm+delta*(dp-dm)/2],
                            [sp+sm-delta*(dp+dm)/2,sp-sm-delta*(dp-dm)/2]])
                det=j[0,0]*j[1,1]-j[0,1]*j[1,0]
                base=f"{group}.{endpoint}.{obs}."
                for r,o in enumerate(("first","second")):
                    for c,inp in enumerate(("first","second")):
                        result[base+f"J[output_{o},input_{inp}]"]=j[r,c]
                result[base+"det_J"]=det
                result[base+"oriented_unit_column_area"]=det/(np.linalg.norm(j[:,0])*np.linalg.norm(j[:,1]))
    return result


source=json.loads(subprocess.check_output(["git","show",SOURCE+":"+PATH],cwd=ROOT))
out={"source_commit":SOURCE,"source_path":PATH,"new_samples":0,
     "interpretation":"Exploratory 2x2 mean-response rank. Nonzero determinant excludes a single fixed mean-response direction for these two controls; not a proof against arbitrary prefix-dependent scalar latent-state models or a CFT field count.",
     "uncertainty":"Delete one of the same20 original batches. No inverse, fitted model, independent-source claim, or multiplicity-adjusted test.","sizes":{}}
for ns,s in source["sizes"].items():
    x=np.asarray(s["joint_20_batch_means"])
    mean=x.mean(0)
    primary=readout(mean,s["labels"],s["delta_cos4"])
    loo=np.array([list(readout((20*mean-row)/19,s["labels"],s["delta_cos4"]).values()) for row in x])
    factor=np.sqrt(19/20)*(loo-loo.mean(0))
    names=list(primary)
    estimate=np.array(list(primary.values()))
    se=np.linalg.norm(factor,axis=0)
    out["sizes"][ns]={"batch_ids":s["batch_ids"],"labels":names,"estimate":estimate.tolist(),
        "se":se.tolist(),"LOO":loo.tolist(),"factor":factor.tolist()}
    for name,val,err in zip(names,estimate,se):
        if name.startswith("all.") and (".A." in name or "input_G" in name):
            print(ns,name,f"{val:.10g} +/- {err:.6g}")
target=ROOT/"results/p334-common-label-response-rank"
target.mkdir(exist_ok=False)
(target/"score.json").write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")

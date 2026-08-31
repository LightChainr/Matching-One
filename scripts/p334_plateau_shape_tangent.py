#!/usr/bin/env python3
"""Lifetime-normalized plateau location/width from existing joint moments."""
import argparse
import json
from pathlib import Path
import subprocess

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
SOURCE_PATH="results/p334-continuous-center-lifetime-moments/batch_moments.json"
DELTA={325:-0.7634556213017751,425:-0.8928996539792388}


def readout(row,labels,n):
    d=dict(zip(labels,row))
    out={}
    base={}
    for ori in ("first","second"):
        b=lambda name:d[f"all.baseline.{ori}.{name}"]
        mass=b("W")
        eta=b("CW")/mass
        var=b("R1_plateau_M2")/mass-eta*eta
        if mass<=0 or var<=0:
            raise ValueError("Plateau baseline mass/variance is not positive")
        base[ori]=(mass,eta,var)
        for name,value in (("mass",mass),("centroid",eta),("variance",var),("width",np.sqrt(var))):
            out[f"baseline.{ori}.{name}"]=value
    for mark in ("plus","minus"):
        vals={}
        for ori in ("first","second"):
            h=lambda name:d[f"all.H.{mark}.{ori}.{name}"]
            mass,eta,var=base[ori]
            h_eta=(h("CW")-eta*h("W"))/mass
            h_var=(h("R1_plateau_M2")-(var+eta*eta)*h("W"))/mass-2*eta*h_eta
            vals[ori]={"mass":h("W"),"centroid":h_eta,"variance":h_var,
                       "width":h_var/(2*np.sqrt(var)),
                       "centroid_minus_unweighted_C":h_eta-h("C")}
            out.update({f"{mark}->{ori}.{name}":value for name,value in vals[ori].items()})
        for channel in ("S","D"):
            weights=(.5,.5) if channel=="S" else (1/DELTA[n],-1/DELTA[n])
            for name in vals["first"]:
                out[f"{mark}->{channel}.{name}"]=sum(w*vals[ori][name] for w,ori in zip(weights,("first","second")))
    return out


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--moment-commit",required=True)
    args=parser.parse_args()
    sha=subprocess.check_output(["git","rev-parse",args.moment_commit+"^{commit}"],cwd=ROOT,text=True).strip()
    source=json.loads(subprocess.check_output(["git","show",sha+":"+SOURCE_PATH],cwd=ROOT))
    out={"source_commit":sha,"source_path":SOURCE_PATH,"new_samples":0,
         "definition":"Normalize the R1 plateau measure by its unperturbed positive mean lifetime, separately in each geometry: eta=E(CW)/E W, V=M2/E W-eta^2. Differentiate both ratios, then take S/D.",
         "scope":"Not a response centroid divided by a near-zero signed integral. V mixes between-window center spread and within-window width. Same20 original batches, plugin derivatives with all expressions reformed in LOO.",
         "sizes":{}}
    for ns,s in source["sizes"].items():
        n=int(ns)
        x=np.array(s["joint_20_batch_means"])
        mean=x.mean(0)
        value=readout(mean,s["labels"],n)
        loo=np.array([list(readout((20*mean-row)/19,s["labels"],n).values()) for row in x])
        factor=np.sqrt(19/20)*(loo-loo.mean(0))
        se=np.linalg.norm(factor,axis=0)
        out["sizes"][ns]={"batch_ids":s["batch_ids"],"labels":list(value),
            "estimate":list(value.values()),"se":se.tolist(),"LOO":loo.tolist(),"factor":factor.tolist()}
        for name,v,e in zip(value,value.values(),se):
            if name.startswith(("baseline.","plus->S.","minus->D.")):
                print(ns,name,f"{v:.11g} +/- {e:.6g}")
    target=ROOT/"results/p334-plateau-shape-tangent"
    target.mkdir(exist_ok=False)
    (target/"score.json").write_text(json.dumps(out,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    main()

#!/usr/bin/env python3
"""Pole versus weight effects of the saved T4 Schur bridge: small blocks only."""
from __future__ import annotations

import os
for _key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_key] = "1"

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import linalg
from threadpoolctl import threadpool_limits

ROOT=Path(__file__).resolve().parents[1]
PROTOCOL=ROOT/"analysis/p398_width8_T4_pole_weight.json"
BLOCKS=ROOT/"results/p398-width8-T4-schur-bridge/latest.json"
FULL=ROOT/"results/p398-width8-reversible-current-control/latest.json"
OUT=ROOT/"results/p398-width8-T4-pole-weight/latest.json"


def encode(value):
    x=np.asarray(value)
    return np.stack((x.real,x.imag),axis=-1).tolist()


def decode(value):
    x=np.asarray(value)
    return x[...,0]+1j*x[...,1]


def endpoint(block,parameter,t):
    h=block.copy()
    h[7,:7]*=parameter
    perturbation=np.zeros_like(h)
    perturbation[7,:7]=block[7,:7]
    values,vectors=linalg.eig(h)
    inverse=linalg.inv(vectors)
    residues=vectors[0,:]*inverse[:,0]
    visible=np.flatnonzero(abs(residues)>1e-10*max(abs(residues)))
    selected=int(visible[np.argmax(values[visible].real)])
    z,r=values[selected],residues[selected]
    projector=np.outer(vectors[:,selected],inverse[selected,:])
    reduced_resolvent=sum((np.outer(vectors[:,j],inverse[j,:])/(z-values[j])
                          for j in range(len(values)) if j!=selected),np.zeros_like(h))
    dz=inverse[selected,:]@perturbation@vectors[:,selected]
    dp=reduced_resolvent@perturbation@projector+projector@perturbation@reduced_resolvent
    dr=dp[0,0]
    total=complex(np.sum(residues*np.exp(values*t)))
    slow=r*np.exp(z*t)
    return {"lambda":parameter,"slow_generator_pole_re_im":encode(z),"slow_mass":float(-z.real),
            "slow_residue_re_im":encode(r),"slow_contribution_t4":float(slow.real),
            "total_correlation_t4":float(total.real),"other_modes_t4":float((total-slow).real),
            "generator_pole_derivative_re_im":encode(dz),"mass_derivative":float(-dz.real),
            "residue_derivative_re_im":encode(dr),
            "log_slow_derivative_pole_part":float((t*dz).real),
            "log_slow_derivative_weight_part":float((dr/r).real),
            "nearest_other_eigenvalue_distance":float(min(abs(z-values[j]) for j in range(len(values)) if j!=selected)),
            "all_generator_poles_re_im":encode(values),"all_source_residues_re_im":encode(residues)}


def analyze(block,t,reference):
    old,new=endpoint(block,0,t),endpoint(block,1,t)
    z0,z1=decode(old["slow_generator_pole_re_im"]),decode(new["slow_generator_pole_re_im"])
    r0,r1=decode(old["slow_residue_re_im"]),decode(new["slow_residue_re_im"])
    exp0,exp1=np.exp(z0*t),np.exp(z1*t)
    pole_part=complex((r0+r1)*(exp1-exp0)/2)
    weight_part=complex((exp0+exp1)*(r1-r0)/2)
    delta=complex(new["total_correlation_t4"]-old["total_correlation_t4"])
    other=delta-pole_part-weight_part
    a,b,c,d=block[:7,:7],block[:7,7],block[7,:7],block[7,7]
    resolvent=linalg.inv(z1*np.eye(7)-a)
    left=(resolvent@b)[0]
    right=(c@resolvent)[0]
    loop=c@resolvent@b
    derivative=1+c@resolvent@resolvent@b
    schur={"z_minus_d_re_im":encode(z1-d),"c_R7_b_re_im":encode(loop),
           "denominator_at_new_pole_re_im":encode(z1-d-loop),
           "denominator_z_derivative_re_im":encode(derivative),
           "left_source_resolvent_factor_re_im":encode(left),
           "right_source_resolvent_factor_re_im":encode(right),
           "residue_from_Schur_re_im":encode(left*right/derivative),
           "pole_derivative_from_Schur_re_im":encode(loop/derivative)}
    return {"old7":old,"new8":new,"schur_at_lambda1":schur,
            "slow_log_change":{"pole":float((t*(z1-z0)).real),
                               "weight":float(np.log(r1/r0).real),
                               "total":float(np.log((r1*exp1)/(r0*exp0)).real)},
            "t4_symmetric_additive_decomposition":{
                "full_G_reference":reference,
                "total_change":float(delta.real),"slow_pole_part":float(pole_part.real),
                "slow_weight_part":float(weight_part.real),"other_modes_part":float(other.real),
                "parts_in_full_reference_percentage_points":{
                    "slow_pole":float(100*pole_part.real/reference),
                    "slow_weight":float(100*weight_part.real/reference),
                    "other_modes":float(100*other.real/reference),
                    "total":float(100*delta.real/reference)},
                "pole_signed_share_of_total_change":float((pole_part/delta).real),
                "weight_signed_share_of_total_change":float((weight_part/delta).real),
                "other_modes_signed_share_of_total_change":float((other/delta).real)}}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json",type=Path,default=OUT)
    args=parser.parse_args()
    protocol=json.loads(PROTOCOL.read_text())
    blocks=json.loads(BLOCKS.read_text())
    full=json.loads(FULL.read_text())
    rows=[]
    with threadpool_limits(limits=1):
        for row,baseline in zip(blocks["ray_rows"],full["ray_rows"]):
            reference=next(s["u"] for s in baseline["original"]["samples"] if s["s"]==protocol["fixed_distance"])
            result=analyze(decode(row["original"]["generator_8_re_im"]),protocol["fixed_distance"],reference)
            result.update({"ray":row["ray"],"sign":row["sign"]})
            rows.append(result)
    result={"schema":protocol["schema"],"protocol":str(PROTOCOL.relative_to(ROOT)),"ray_rows":rows,
            "computation":"only saved 8x8 matrices at lambda=0,1; no full-state model imports or reconstruction",
            "boundary":"a descriptive decomposition of the same named bridge, not a coupling fit or an independent physical process"}
    result["input_sha256"]={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (PROTOCOL,Path(__file__),BLOCKS,FULL)}
    args.json.parent.mkdir(parents=True,exist_ok=True)
    args.json.write_text(json.dumps(result,indent=2)+"\n")
    for row in rows:
        print(row["ray"],row["t4_symmetric_additive_decomposition"])
        print("log",row["slow_log_change"])
        for name in ("old7","new8"):
            print(name,{k:row[name][k] for k in ("slow_mass","slow_residue_re_im","mass_derivative","residue_derivative_re_im","log_slow_derivative_pole_part","log_slow_derivative_weight_part")})
        print("Schur",row["schur_at_lambda1"])


if __name__=="__main__":
    main()

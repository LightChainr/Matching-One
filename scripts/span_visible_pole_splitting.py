#!/usr/bin/env python3
"""Projected two-pole diagnostic for complete-component span spectra.

Reads the already generated span-spectrum JSON and does not rebuild a transfer
matrix. Fits d_{h+2}=S d_{h+1}-P d_h on a declared late window after dividing
each row by d_{h+1}, then fits the two exponential residues.
"""
from __future__ import annotations
import argparse, json, math
from pathlib import Path


def solve_2x2(a11,a12,a22,b1,b2):
    det=a11*a22-a12*a12
    if abs(det)<1e-30: raise ArithmeticError("singular normal equations")
    return ((b1*a22-b2*a12)/det,(a11*b2-a12*b1)/det)


def fit_recurrence(d,start_h=13):
    rows=[]
    for h in range(start_h,len(d)-1):
        a,b,c=d[h-1],d[h],d[h+1]
        if min(a,b,c)<=0: continue
        rows.append((1.0,-a/b,c/b,h))
    if len(rows)<4: raise ValueError("too few late-window triples")
    a11=sum(x*x for x,z,y,h in rows); a12=sum(x*z for x,z,y,h in rows); a22=sum(z*z for x,z,y,h in rows)
    b1=sum(x*y for x,z,y,h in rows); b2=sum(z*y for x,z,y,h in rows)
    S,P=solve_2x2(a11,a12,a22,b1,b2)
    disc=S*S-4*P
    if disc<=0: raise ArithmeticError("nonpositive recurrence discriminant")
    rr=math.sqrt(disc); rho1=(S+rr)/2; rho2=(S-rr)/2
    if not (rho1>rho2>0): raise ArithmeticError("unexpected recurrence roots")
    rec=[]
    for x,z,y,h in rows:
        pred=S+z*P; rec.append(abs(pred-y)/abs(y))

    q=rho2/rho1; xs=[]; ys=[]; hs=[]
    for h in range(start_h,len(d)+1):
        val=d[h-1]
        if val<=0: continue
        xs.append(q**(h-1)); ys.append(val/(rho1**(h-1))); hs.append(h)
    n=len(xs); sx=sum(xs); sxx=sum(x*x for x in xs); sy=sum(ys); sxy=sum(x*y for x,y in zip(xs,ys))
    c1,c2=solve_2x2(float(n),sx,sxx,sy,sxy)
    fit=[]
    for h in hs:
        val=d[h-1]; pred=c1*rho1**(h-1)+c2*rho2**(h-1); fit.append(abs(pred-val)/val)
    dg=math.log(rho1/rho2)
    return {"start_h":start_h,"rho1":rho1,"rho2":rho2,"gamma1":-math.log(rho1),"gamma2":-math.log(rho2),
            "delta_gamma":dg,"relaxation_length":1/dg,"c1":c1,"c2":c2,"c2_over_c1":c2/c1,
            "max_recurrence_relative_residual":max(rec),"max_two_exp_relative_residual":max(fit)}


def tail_hazard(d,tail,h):
    if h<1 or h>=len(d)+1: return None
    t0=sum(d[h-1:])+tail; t1=sum(d[h:])+tail
    return -math.log(t1/t0) if t0>0 and t1>0 else None


def run_record(rec,start_h):
    d=[float(x) for x in rec["d_h_float"]]; fit=fit_recurrence(d,start_h); w=int(rec["width"]); nu=float(rec.get("nu_total_float",sum(d)))
    return {"label":rec.get("label"),"width":w,"matching":bool(rec.get("matching")),"p":rec.get("p"),"d_max":len(d),"nu":nu,**fit,
            "relaxation_length_over_width":fit["relaxation_length"]/w,
            "delta_gamma_over_w_nu":fit["delta_gamma"]/(w*nu) if nu>0 else None,
            "gamma_eff_at_w":tail_hazard(d,float(rec.get("tail_bin_float",0.0)),w),
            "gamma_eff_at_2w":tail_hazard(d,float(rec.get("tail_bin_float",0.0)),2*w)}


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--input",required=True); ap.add_argument("--start-h",type=int,default=13); ap.add_argument("--output")
    a=ap.parse_args(); payload=json.loads(Path(a.input).read_text()); rows=[]; errors=[]
    for rec in payload["runs"]:
        try: rows.append(run_record(rec,a.start_h))
        except Exception as exc: errors.append({"label":rec.get("label"),"width":rec.get("width"),"matching":rec.get("matching"),"p":rec.get("p"),"error":str(exc)})
    out={"schema":"matching-one/span-visible-pole-splitting/v1","source":a.input,"start_h":a.start_h,
         "method":"normalized second-order recurrence + two-exponential residue fit","rows":rows,"errors":errors,
         "claim_boundary":"Projected/Prony diagnostic only; not a certified full transfer eigendecomposition."}
    txt=json.dumps(out,indent=2)
    if a.output: Path(a.output).write_text(txt+"\n")
    print(txt)
if __name__=="__main__": main()

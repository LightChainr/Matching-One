#!/usr/bin/env python3
"""Fast threshold-rank analyzer for matching-parity spin-4 derivatives.

This is the production companion to `analyze_matching_parity_derivatives.py`.
It uses O(N) binomial recurrences rather than one incomplete-beta call per rank,
which matters for delete-one-batch intrinsic-center jackknifes.

Input columns: n,a,b,orientation,batch,samples,kind,k,count, with
orientation in {first,second}, kind in {minus,plus}.
"""
from __future__ import annotations

import argparse, csv, json
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple
import mpmath as mp

METRICS=("P4_S","P4_D","P4_S_prime","P4_D_prime")
Key=Tuple[int,str,int]

@dataclass
class H:
    n:int; a:int; b:int; orientation:str; batch:int; samples:int
    minus:List[int]; plus:List[int]

def read(path:Path)->Dict[Key,H]:
    out={}
    with path.open(newline="",encoding="utf-8") as f:
        r=csv.DictReader(f)
        req={"n","a","b","orientation","batch","samples","kind","k","count"}
        miss=req-set(r.fieldnames or [])
        if miss: raise ValueError("missing columns: "+",".join(sorted(miss)))
        for x in r:
            n=int(x["n"]); o=x["orientation"]; bch=int(x["batch"]); k=int(x["k"])
            if o not in ("first","second") or x["kind"] not in ("minus","plus"):
                raise ValueError("bad orientation/kind")
            key=(n,o,bch)
            if key not in out:
                out[key]=H(n,int(x["a"]),int(x["b"]),o,bch,int(x["samples"]),[0]*(n+1),[0]*(n+1))
            row=out[key]
            if (row.a,row.b,row.samples)!=(int(x["a"]),int(x["b"]),int(x["samples"])):
                raise ValueError("inconsistent batch metadata")
            getattr(row,x["kind"])[k]+=int(x["count"])
    for row in out.values():
        if sum(row.minus)!=row.samples or sum(row.plus)!=row.samples:
            raise ValueError("histogram total mismatch")
    if not out: raise ValueError("empty input")
    return out

def combine(rows:List[H])->H:
    r0=rows[0]; z=H(r0.n,r0.a,r0.b,r0.orientation,-1,0,[0]*(r0.n+1),[0]*(r0.n+1))
    for r in rows:
        if (r.n,r.a,r.b,r.orientation)!=(r0.n,r0.a,r0.b,r0.orientation): raise ValueError("incompatible rows")
        z.samples+=r.samples
        for k in range(1,r.n+1): z.minus[k]+=r.minus[k]; z.plus[k]+=r.plus[k]
    return z

def remove(total:H,row:H)->H:
    return H(total.n,total.a,total.b,total.orientation,-1,total.samples-row.samples,
             [a-b for a,b in zip(total.minus,row.minus)],
             [a-b for a,b in zip(total.plus,row.plus)])

def cos4(a:int,b:int)->mp.mpf:
    n=a*a+b*b
    return mp.mpf(a**4-6*a*a*b*b+b**4)/(n*n)

def tail_expectation(n:int,counts:List[int],samples:int,p:mp.mpf)->mp.mpf:
    """E[ P(Bin(N,p)>=K) ] averaged over a threshold-rank histogram."""
    if p<=0: return mp.mpf(0)
    if p>=1: return mp.mpf(1)
    q=1-p; prob=q**n; cumulative=0; total=mp.mpf(0)
    for occupied in range(n+1):
        if occupied: cumulative+=counts[occupied]
        total+=cumulative*prob
        if occupied<n: prob*=mp.mpf(n-occupied)*p/((occupied+1)*q)
    return total/samples

def tail_derivative(n:int,counts:List[int],samples:int,p:mp.mpf)->mp.mpf:
    """Derivative of E[P(Bin(N,p)>=K)] by the beta-density recurrence."""
    if not 0<p<1: return mp.mpf(0)
    q=1-p; density=n*q**(n-1); total=mp.mpf(0)
    for k in range(1,n+1):
        total+=counts[k]*density
        if k<n: density*=mp.mpf(n-k)*p/(k*q)
    return total/samples

def obs(h:H,p:mp.mpf):
    rg=tail_expectation(h.n,h.plus,h.samples,p)
    minus_cdf=tail_expectation(h.n,h.minus,h.samples,p)
    rh=1-minus_cdf
    rgp=tail_derivative(h.n,h.plus,h.samples,p)
    rhp=-tail_derivative(h.n,h.minus,h.samples,p)
    return {"S":(rg+rh)/2,"D":(rg-rh)/2,"Sp":(rgp+rhp)/2,"Dp":(rgp-rhp)/2,"M":rg-rh}

def center(a:H,b:H,it:int=120)->mp.mpf:
    lo=mp.mpf(0); hi=mp.mpf(1)
    for _ in range(it):
        p=(lo+hi)/2
        m=(obs(a,p)["M"]+obs(b,p)["M"])/2
        if m<0: lo=p
        else: hi=p
    return (lo+hi)/2

def project(a:H,b:H):
    p=center(a,b); x=obs(a,p); y=obs(b,p); dc=cos4(a.a,a.b)-cos4(b.a,b.b)
    if dc==0: raise ValueError("Delta cos4 is zero")
    v={"P4_S":(x["S"]-y["S"])/dc,"P4_D":(x["D"]-y["D"])/dc,
       "P4_S_prime":(x["Sp"]-y["Sp"])/dc,"P4_D_prime":(x["Dp"]-y["Dp"])/dc}
    return p,dc,v

def jkcov(rows):
    m=len(rows); means=[mp.fsum(r[k] for r in rows)/m for k in METRICS]; fac=mp.mpf(m-1)/m
    return [[fac*mp.fsum((r[a]-means[i])*(r[b]-means[j]) for r in rows)
             for j,b in enumerate(METRICS)] for i,a in enumerate(METRICS)]

def analyze(data):
    ans={"format_version":1,"metrics":list(METRICS),"by_N":{}}
    for n in sorted({k[0] for k in data}):
        g={o:[data[k] for k in sorted(data) if k[0]==n and k[1]==o] for o in ("first","second")}
        if len(g["first"])!=len(g["second"]) or len(g["first"])<2: raise ValueError(f"N={n}: bad batch pairing")
        A=combine(g["first"]); B=combine(g["second"]); p,dc,v=project(A,B)
        jk=[]; jp=[]
        for ra,rb in zip(g["first"],g["second"]):
            pi,_,vi=project(remove(A,ra),remove(B,rb)); jp.append(pi); jk.append(vi)
        cov=jkcov(jk); se={k:mp.sqrt(cov[i][i]) for i,k in enumerate(METRICS)}
        pm=mp.fsum(jp)/len(jp); pse=mp.sqrt(mp.mpf(len(jp)-1)/len(jp)*mp.fsum((x-pm)**2 for x in jp))
        scaled={"A_S_N1":n*v["P4_S"],"A_D_N13_8":n**(mp.mpf(13)/8)*v["P4_D"],
                "A_Dprime_N5_8":n**(mp.mpf(5)/8)*v["P4_D_prime"],"A_Sprime_N5_4":n**(mp.mpf(5)/4)*v["P4_S_prime"]}
        ans["by_N"][str(n)]={"first_rep":[A.a,A.b],"second_rep":[B.a,B.b],"samples":A.samples,"batches":len(g["first"]),
          "p0":mp.nstr(p,30),"p0_se":mp.nstr(pse,15),"delta_cos4":mp.nstr(dc,30),
          "point":{k:mp.nstr(v[k],30) for k in METRICS},"se":{k:mp.nstr(se[k],15) for k in METRICS},
          "scaled":{k:mp.nstr(x,30) for k,x in scaled.items()},"covariance":[[mp.nstr(x,15) for x in r] for r in cov]}
    return ans

def write_csv(path:Path,payload):
    rows=[]
    for n,r in payload["by_N"].items():
        row={"N":n,"first_rep":",".join(map(str,r["first_rep"])),"second_rep":",".join(map(str,r["second_rep"])),
             "samples":r["samples"],"p0":r["p0"],"p0_se":r["p0_se"],"delta_cos4":r["delta_cos4"]}
        for k in METRICS: row[k]=r["point"][k]; row[k+"_se"]=r["se"][k]
        row.update(r["scaled"]); rows.append(row)
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator="\n"); w.writeheader(); w.writerows(rows)

def main():
    ap=argparse.ArgumentParser(description=__doc__); ap.add_argument("--histograms",type=Path,required=True); ap.add_argument("--dps",type=int,default=40)
    ap.add_argument("--json",type=Path,required=True); ap.add_argument("--csv",type=Path,required=True); a=ap.parse_args()
    mp.mp.dps=a.dps
    try: result=analyze(read(a.histograms))
    except (ValueError,ArithmeticError) as e: raise SystemExit(str(e))
    a.json.parent.mkdir(parents=True,exist_ok=True); a.csv.parent.mkdir(parents=True,exist_ok=True)
    a.json.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8"); write_csv(a.csv,result)
    print(a.json); print(a.csv); return 0
if __name__=="__main__": raise SystemExit(main())

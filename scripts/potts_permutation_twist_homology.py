#!/usr/bin/env python3
"""Finite-Q commuting-permutation twist controls for torus FK homology."""
from __future__ import annotations
import argparse, json, math
from pathlib import Path

def parse_perm(text):
    p=tuple(int(x) for x in text.split(","))
    if sorted(p)!=list(range(len(p))): raise ValueError("not a permutation")
    return p

def compose(p,q): return tuple(p[q[i]] for i in range(len(p)))
def identity(n): return tuple(range(n))
def inverse(p):
    out=[0]*len(p)
    for i,j in enumerate(p): out[j]=i
    return tuple(out)
def power(p,k):
    if k<0: return power(inverse(p),-k)
    out=identity(len(p)); base=p
    while k:
        if k&1: out=compose(base,out)
        base=compose(base,base); k>>=1
    return out
def fixed_points(p): return tuple(i for i,j in enumerate(p) if i==j)

def analyze(sigma,tau,a,b):
    if len(sigma)!=len(tau): raise ValueError("size mismatch")
    if compose(sigma,tau)!=compose(tau,sigma): raise ValueError("twists must commute")
    if math.gcd(abs(a),abs(b))!=1: raise ValueError("slope must be primitive")
    Q=len(sigma); hol=compose(power(sigma,a),power(tau,b)); fu=fixed_points(hol)
    common=tuple(sorted(set(fixed_points(sigma)) & set(fixed_points(tau))))
    return {"Q":Q,"slope":[a,b],"sigma":list(sigma),"tau":list(tau),"slope_holonomy":list(hol),
            "fixed_slope_colors":list(fu),"rank1_neutral_count_fugacity":len(fu)/Q,
            "common_fixed_colors_rank2":list(common),"rank2_cross_fugacity":len(common)/Q,
            "rank0_cluster_fugacity":1.0}

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--sigma",required=True); ap.add_argument("--tau",required=True)
    ap.add_argument("--a",type=int,required=True); ap.add_argument("--b",type=int,required=True); ap.add_argument("--output")
    x=ap.parse_args(); out=analyze(parse_perm(x.sigma),parse_perm(x.tau),x.a,x.b); txt=json.dumps(out,indent=2)
    if x.output: Path(x.output).write_text(txt+"\n")
    print(txt)
if __name__=="__main__": main()

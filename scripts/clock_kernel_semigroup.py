#!/usr/bin/env python3
"""Exact unmarked Poisson intensity-clock gap semigroup."""
from __future__ import annotations
import argparse, json, math
from pathlib import Path

def joint_laplace(r: float, s: float, t: float) -> float:
    if r < 1: raise ValueError("r must be >= 1")
    return (r + s) / ((1.0 + s) * (s + r * (1.0 + t)))

def summary(r: float) -> dict:
    if r < 1: raise ValueError("r must be >= 1")
    return {
        "intensity_ratio": r,
        "delta_log_intensity": math.log(r),
        "one_side_mean": 1.0,
        "one_side_variance": 1.0,
        "one_side_covariance": 1.0 / r,
        "fixed_location_gap_mean": 2.0,
        "fixed_location_gap_variance": 2.0,
        "fixed_location_gap_covariance": 2.0 / r,
        "fixed_location_gap_correlation": 1.0 / r,
        "no_record_one_side": 1.0 / r,
        "no_cut_two_sides": 1.0 / (r * r),
        "generator": "L f(u)=u f'(u)+u int_0^1 [f(uv)-f(u)] dv",
        "warning": "Fixed-location sampling, not once-per-component Palm sampling."
    }

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--ratio",type=float,required=True)
    ap.add_argument("--s",type=float,default=.7); ap.add_argument("--t",type=float,default=.4); ap.add_argument("--output")
    a=ap.parse_args(); out=summary(a.ratio); j=joint_laplace(a.ratio,a.s,a.t)
    out["joint_laplace_example"]={"s":a.s,"t":a.t,"one_side":j,"two_side":j*j}
    txt=json.dumps(out,indent=2)
    if a.output: Path(a.output).write_text(txt+"\n")
    print(txt)
if __name__=="__main__": main()

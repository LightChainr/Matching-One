#!/usr/bin/env python3
"""High-precision Q=1 critical torus rank-shape controls."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import mpmath as mp

def eta_tau(tau, terms=140):
    q=mp.e**(2*mp.pi*1j*tau); prod=mp.mpc(1)
    for n in range(1,terms+1): prod*=1-q**n
    return mp.e**(mp.pi*1j*tau/12)*prod

def theta_unit_area(g,tau,cutoff=14):
    y=mp.im(tau); eta=eta_tau(tau); total=mp.mpf(0)
    for m in range(-cutoff,cutoff+1):
        for n in range(-cutoff,cutoff+1):
            total += mp.e**(-mp.pi*g*abs(m*tau-n)**2/y)
    return mp.sqrt(g/y)*total/(abs(eta)**2)

def p0_xy(x,y):
    tau=mp.mpf(x)+1j*mp.mpf(y)
    return (theta_unit_area(mp.mpf(8)/3,tau)-theta_unit_area(mp.mpf(2)/3,tau))/2

def c_xy(x,y):
    p=p0_xy(x,y); return 1/(2*p)-1

def point(x,y):
    c=c_xy(x,y)
    gx=mp.diff(lambda xx:c_xy(xx,y),x); gy=mp.diff(lambda yy:c_xy(x,yy),y)
    hxx=mp.diff(lambda xx:c_xy(xx,y),x,2); hyy=mp.diff(lambda yy:c_xy(x,yy),y,2)
    hxy=mp.diff(lambda xx:mp.diff(lambda yy:c_xy(xx,yy),y),x)
    return {"x":mp.nstr(x,30),"y":mp.nstr(y,30),"P0":mp.nstr(p0_xy(x,y),45),"c":mp.nstr(c,45),
            "gradient":[mp.nstr(gx,28),mp.nstr(gy,28)],
            "hessian":[[mp.nstr(hxx,32),mp.nstr(hxy,32)],[mp.nstr(hxy,32),mp.nstr(hyy,32)]]}

def collision_y(x): return mp.findroot(lambda y:c_xy(x,y)-1,(mp.mpf("1.5"),mp.mpf("2.1")))

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--dps",type=int,default=70); ap.add_argument("--output")
    a=ap.parse_args(); mp.mp.dps=a.dps
    hx,hy=mp.mpf("0.5"),mp.sqrt(3)/2; sx,sy=mp.mpf(0),mp.mpf(1)
    xs=[mp.mpf(k)/20 for k in range(11)]; ys=[collision_y(x) for x in xs]; y0=ys[0]
    basis=[1-mp.cos(2*mp.pi*x) for x in xs]
    amp=mp.fsum(basis[i]*(ys[i]-y0) for i in range(len(xs)))/mp.fsum(v*v for v in basis)
    resid=max(abs(ys[i]-(y0+amp*basis[i])) for i in range(len(xs)))
    out={"hex":point(hx,hy),"square":point(sx,sy),
         "collision_curve":[{"x":mp.nstr(x,12),"y":mp.nstr(y,40)} for x,y in zip(xs,ys)],
         "first_harmonic_collision_fit":{"y0":mp.nstr(y0,40),"amplitude":mp.nstr(amp,28),"max_grid_residual":mp.nstr(resid,24)},
         "claim_boundary":"Numerical Morse/collision controls only; no global hex-minimum theorem."}
    txt=json.dumps(out,indent=2)
    if a.output: Path(a.output).write_text(txt+"\n")
    print(txt)
if __name__=="__main__": main()

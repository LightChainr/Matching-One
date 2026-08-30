#!/usr/bin/env python3
"""Exact square-lattice angular harmonics for Gaussian-integer tori.

For z=a+i b, N=|z|^2 and theta=arg(z),

    exp(i*4*m*theta) = z^(4m) / N^(2m).

Thus cos/sin of every square-lattice harmonic 4m is rational.  This module is
used to freeze angular design matrices and Gaussian-multiplier genealogies
without floating trigonometry.
"""
from __future__ import annotations

import argparse, json, math
from fractions import Fraction
from typing import Tuple

Pair=Tuple[int,int]

def gmul(z:Pair,w:Pair)->Pair:
    a,b=z; c,d=w
    return a*c-b*d, a*d+b*c

def gpow(z:Pair,n:int)->Pair:
    if n<0: raise ValueError("nonnegative power required")
    out=(1,0); base=z
    while n:
        if n&1: out=gmul(out,base)
        base=gmul(base,base); n//=2
    return out

def norm(z:Pair)->int:
    return z[0]*z[0]+z[1]*z[1]

def harmonic(z:Pair,m:int)->Tuple[Fraction,Fraction]:
    if m<=0: raise ValueError("m must be positive; harmonic spin is 4m")
    n=norm(z)
    if n==0: raise ValueError("zero Gaussian integer")
    re,im=gpow(z,4*m); den=n**(2*m)
    return Fraction(re,den),Fraction(im,den)

def canonical_square(z:Pair)->Pair:
    """Canonical representative under sign, conjugation and coordinate swap."""
    a,b=z
    candidates=set()
    for x,y in ((a,b),(a,-b),(-a,b),(-a,-b),(b,a),(b,-a),(-b,a),(-b,-a)):
        candidates.add((abs(x),abs(y)))
    return max(candidates,key=lambda p:(max(p),min(p)))

def check(z:Pair)->None:
    c4,s4=harmonic(z,1); c8,s8=harmonic(z,2); c12,s12=harmonic(z,3)
    assert c4*c4+s4*s4==1
    assert c8*c8+s8*s8==1
    assert c12*c12+s12*s12==1
    assert c8==2*c4*c4-1
    assert c12==4*c4*c4*c4-3*c4
    assert s8==2*s4*c4
    assert s12==s4*(4*c4*c4-1)

def payload(z:Pair,max_m:int=3)->dict:
    check(z)
    return {
        "z":list(z),"N":norm(z),"primitive":math.gcd(abs(z[0]),abs(z[1]))==1,
        "harmonics":{
            str(4*m):{
                "cos_fraction":str(harmonic(z,m)[0]),
                "sin_fraction":str(harmonic(z,m)[1]),
                "cos":float(harmonic(z,m)[0]),
                "sin":float(harmonic(z,m)[1]),
            } for m in range(1,max_m+1)
        }
    }
def main()->int:
    p=argparse.ArgumentParser(description=__doc__); p.add_argument("a",type=int); p.add_argument("b",type=int)
    p.add_argument("--max-m",type=int,default=3); p.add_argument("--multiply",nargs=2,type=int,metavar=("U","V"))
    a=p.parse_args(); z=(a.a,a.b); out={"parent":payload(z,a.max_m)}
    if a.multiply:
        h=tuple(a.multiply); child=gmul(z,h)  # type: ignore[arg-type]
        out["multiplier"]={"z":list(h),"norm":norm(h)}
        out["child_raw"]=payload(child,a.max_m)
        # exact consistency of harmonic multiplication
        for m in range(1,a.max_m+1):
            pc,ps=harmonic(z,m); hc,hs=harmonic(h,m); cc,cs=harmonic(child,m)
            assert cc==pc*hc-ps*hs and cs==ps*hc+pc*hs
    print(json.dumps(out,indent=2)); return 0
if __name__=="__main__": raise SystemExit(main())

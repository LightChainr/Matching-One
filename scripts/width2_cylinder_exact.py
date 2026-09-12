#!/usr/bin/env python3
"""Exact 2-by-m square-site torus probabilities and their cylinder limit.

No new percolation production. Bernoulli-polynomial coefficients are integers.
A width-two occupied row makes a transverse cycle; both periodic bonds must
be retained even though their endpoint pairs coincide. m >= 2 is required.
"""
from __future__ import annotations
import argparse
from fractions import Fraction
import json
from math import comb
from pathlib import Path


def add(a, b):
    out = [0] * max(len(a), len(b))
    for i, v in enumerate(a): out[i] += v
    for i, v in enumerate(b): out[i] += v
    while len(out)>1 and out[-1]==0: out.pop()
    return out


def mul(a, b):
    out = [0] * (len(a)+len(b)-1)
    for i, u in enumerate(a):
        for j, v in enumerate(b): out[i+j] += u*v
    while len(out)>1 and out[-1]==0: out.pop()
    return out


def power(a, m):
    out = [1]
    for _ in range(m): out = mul(out,a)
    return out


def evaluate(a, p):
    out = 0
    for v in reversed(a): out = out*p+v
    return out


def matching_polynomial(m):
    """Integer power coefficients, low degree first, from a 3-state trace."""
    if m < 2: raise ValueError('honest two-by-m torus requires m >= 2')
    # Symmetric 2-state block: trace=p, determinant=-p^3(1-p).
    # tr(block^m)=p*tr(block^(m-1))+p^3(1-p)*tr(block^(m-2)).
    previous, current = [2], [0,1]
    for _ in range(2,m+1):
        previous,current = current,add(mul([0,1],current),mul([0,0,0,1,-1],previous))
    return add(add(current,power([0,1,-1],m)),[-x for x in power([1,0,-1],m)])


def ambient_rank(mask, m):
    """Independent lifted-edge graph traversal; not a row compatibility test."""
    if m < 2: raise ValueError('m >= 2 required')
    n = 2*m
    positions, span = {}, []
    for root in range(n):
        if not (mask>>root)&1 or root in positions: continue
        positions[root] = (0,0)
        stack = [root]
        while stack:
            v=stack.pop(); x,y=v%2,v//2; px,py=positions[v]
            for dx,dy in ((1,0),(-1,0),(0,1),(0,-1)):
                u=((y+dy)%m)*2+(x+dx)%2
                if not (mask>>u)&1: continue
                proposed=(px+dx,py+dy)
                if u not in positions:
                    positions[u]=proposed;stack.append(u)
                else:
                    wx,wy=proposed[0]-positions[u][0],proposed[1]-positions[u][1]
                    if wx%2 or wy%m: raise AssertionError('nonperiodic cycle displacement')
                    wx,wy=wx//2,wy//m
                    if wx or wy:
                        if not span: span.append((wx,wy))
                        elif span[0][0]*wy-span[0][1]*wx: return 2
    return len(span)


def enumerated_polynomial(m):
    """Tiny independent exact census, collapsed by occupation count."""
    n=2*m; bern=[0]*(n+1)
    for mask in range(1<<n): bern[mask.bit_count()] += ambient_rank(mask,m)-1
    poly=[0]*(n+1)
    for k,b in enumerate(bern):
        for j in range(n-k+1): poly[k+j]+=b*comb(n-k,j)*(-1)**j
    while len(poly)>1 and poly[-1]==0: poly.pop()
    return poly,bern


def root_bracket(poly, steps=100):
    lo,hi=Fraction(0),Fraction(1)
    for _ in range(steps):
        mid=(lo+hi)/2
        if evaluate(poly,mid)<0: lo=mid
        else: hi=mid
    return [str(lo),str(hi)]


def report(max_check=6):
    from mpmath import mp
    with mp.workdps(70):
        q=mp.findroot(lambda p:2*p**3+2*p**2-1,('0.55','0.58'))
        x=q*(1-q); lc=1-q*q
        lp_derivative=(lc+3*q*q-4*q**3)/(2*lc-q)
        hprime=lp_derivative/lc+2*q/lc
        checks={}
        for m in range(2,max_check+1):
            got,bern=enumerated_polynomial(m); expected=matching_polynomial(m)
            if got!=expected: raise AssertionError(f'coefficient disagreement m={m}')
            checks[str(m)]={'configurations':1<<(2*m),'coefficient_identity':True,
                            'power_coefficients':expected,'bernstein_counts':bern}
        roots={}
        for m in (2,3,4,6,8,12,20):
            poly=matching_polynomial(m); bracket=root_bracket(poly)
            lo,hi=map(Fraction,bracket)
            root=(mp.mpf(lo.numerator)/lo.denominator+mp.mpf(hi.numerator)/hi.denominator)/2
            predicted=-(x/lc)**m/(m*hprime)
            roots[str(m)]={'exact_rational_bracket':bracket,'root_diagnostic':mp.nstr(root,28),
                           'root_minus_cylinder':mp.nstr(root-q,20),
                           'leading_shift':mp.nstr(predicted,20),
                           'shift_over_leading':mp.nstr((root-q)/predicted,15)}
        return {'schema':'matching-one.width2-cylinder-exact.v1',
                'scope':'axis square-site 2-by-m honest torus, m>=2; row transfer, not an all-width pTL intertwiner',
                'cylinder_minimal_polynomial_low_first':[-1,0,2,2],
                'cylinder_root_diagnostic':mp.nstr(q,60),
                'cylinder_root_isolation':root_bracket([-1,0,2,2]),
                'published_n2_agrees_within_1e_minus_40':bool(abs(q-mp.mpf('0.5651977173836393964375280132470308160984'))<mp.mpf('1e-40')),
                'hprime_diagnostic':mp.nstr(hprime,30),
                'decay_ratio_diagnostic':mp.nstr(x/lc,30),
                'finite_m_defect_at_crossing':'[p(1-p)]^m + lambda_minus^m > 0',
                'equal_leading_coefficients':True,'small_exact_checks':checks,'finite_roots':roots}


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path);args=ap.parse_args()
    text=json.dumps(report(),indent=2,allow_nan=False)+'\n'
    if args.out:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        with args.out.open('x') as f:f.write(text)
    else: print(text,end='')

#!/usr/bin/env python3
"""Geometrically anchored thermodynamic integration of finite rank-sector odds.

A fixed-size exact-polynomial oracle validates the conditional-mean formula.
This does not claim that obtaining independent conditional samples is free.
"""
from __future__ import annotations
import argparse,json
from fractions import Fraction as F
from math import comb
from pathlib import Path
import time
from mpmath import mp

def value(c,t):
    y=0
    for x in reversed(c):y=y*t+x
    return y

def deriv(c):return [i*c[i] for i in range(1,len(c))]

def normalize_sector(c):
    k=next(i for i,a in enumerate(c) if a)
    a=c[k]
    return k,a,[F(x,a) for x in c[k:]]

def mp_coeff(c):return [mp.mpf(x.numerator)/x.denominator if isinstance(x,F) else mp.mpf(x) for x in c]

def occupation_moments(c,t):
    z=value(c,t)
    if z<=0:raise ValueError('positive fugacity and nonempty sector required')
    first=sum(i*a*t**i for i,a in enumerate(c))/z
    second=sum(i*i*a*t**i for i,a in enumerate(c))/z
    return first,second-first*first

def anchor_bound(N,k,a,t0):
    """Coefficient-free finite lower-end correction bound using binomial majorants."""
    lower=(1+t0)**N-1
    upper=sum(F(comb(N,k+j),a)*t0**j for j in range(1,N-k+1))
    return max(lower,upper)

def report(hessianpath,dps=70):
    started=time.perf_counter();h=json.loads(hessianpath.read_text())
    counts=h['rank_bernstein_counts'];N=len(counts[0])-1
    k0,a0,r0=normalize_sector(counts[0]);k2,a2,r2=normalize_sector(counts[2])
    assert k0==0 and a0==1 and k2==7 and a2==16
    with mp.workdps(dps):
        R0,R2=mp_coeff(r0),mp_coeff(r2)
        D0,D2=deriv(R0),deriv(R2)
        def integrand(t):return value(D2,t)/value(R2,t)-value(D0,t)/value(R0,t)
        def integrated_g(t):
            return mp.log(mp.mpf(a2)/a0)+(k2-k0)*mp.log(t)+mp.quad(integrand,[0,t/4,t/2,t])
        controls=[]
        for pp in ('0.2','0.5','0.5906721123310283','0.8'):
            p=mp.mpf(pp);t=p/(1-p)
            g=integrated_g(t);direct=mp.log(value(counts[2],t)/value(counts[0],t))
            mu0,var0=occupation_moments(counts[0],t);mu2,var2=occupation_moments(counts[2],t)
            identity=abs((k2-k0)+t*integrand(t)-(mu2-mu0))
            if abs(g-direct)>mp.mpf('1e-55') or identity>mp.mpf('1e-55'):
                raise ArithmeticError('integration or conditional-mean identity failed')
            controls.append({'p':pp,'log_odds_by_integration':mp.nstr(g,45),
                             'log_odds_direct':mp.nstr(direct,45),
                             'absolute_error':mp.nstr(abs(g-direct),10),
                             'conditional_K_mean_rank0':mp.nstr(mu0,30),'conditional_K_mean_rank2':mp.nstr(mu2,30),
                             'conditional_K_variance_rank0':mp.nstr(var0,30),'conditional_K_variance_rank2':mp.nstr(var2,30)})
        root=mp.findroot(lambda p:integrated_g(p/(1-p)),('.58','.60'))
        direct=mp.findroot(lambda p:value(counts[2],p/(1-p))-value(counts[0],p/(1-p)),('.58','.60'))
        assert abs(root-direct)<mp.mpf('1e-55')
        tail=[]
        for den in (10**4,10**8,10**12):
            t0=F(1,den);bound=anchor_bound(N,k2,a2,t0)
            x=mp.mpf(t0.numerator)/t0.denominator
            actual=mp.log(value(R2,x)/value(R0,x))
            if abs(actual)>mp.mpf(bound.numerator)/bound.denominator:raise AssertionError('anchor bound failed')
            tail.append({'t0':str(t0),'rigorous_absolute_anchor_error_bound':str(bound),
                         'bound_diagnostic':float(bound),'actual_anchor_correction':mp.nstr(actual,35)})
    return {'schema':'matching-one.conditional-odds-integration.v1','sites':N,'rank2_onset':k2,
            'rank2_onset_count':a2,'rank0_onset':k0,'rank0_onset_count':a0,
            'dps':dps,'controls':controls,'root_via_integrated_conditional_means':mp.nstr(root,50),
            'direct_polynomial_root':mp.nstr(direct,50),'low_fugacity_anchor_bounds':tail,
            'scope':'fixed exact 4x4 oracle; path-sampling identity not a new sampler or complexity claim',
            'elapsed_seconds':time.perf_counter()-started}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--hessian',type=Path,required=True)
    p.add_argument('--out',type=Path);p.add_argument('--dps',type=int,default=70);a=p.parse_args()
    text=json.dumps(report(a.hessian,a.dps),indent=2)+'\n'
    if a.out:
        with a.out.open('x') as f:f.write(text)
    else:print(text,end='')
if __name__=='__main__':main()

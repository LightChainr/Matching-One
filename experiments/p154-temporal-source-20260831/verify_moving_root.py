#!/usr/bin/env python3
"""Separate finite-difference check of the lagged source's moving-root derivative.

The auxiliary measure uses (1+t*z), i.e. only a local first-order path tangent.
This is not a finite exp(t*z) production measurement or a second sample.
"""
import json,math
from pathlib import Path
from fractions import Fraction
import numpy as np
from scipy.optimize import brentq
from analyze import ROOT,NS,load,moments

def main():
    result=json.loads((ROOT/'results/latest.json').read_text());checks=[]
    for n in NS:
        total=load(n).sum(axis=0);pool=total.sum(axis=2);m=1000000 if n in (260,340) else 100000
        count=total[...,0];mu=np.divide(total[...,3],count,out=np.zeros_like(count),where=count>0)
        raw=pool[...,[1,2,3,4,5]].copy();centered=raw.copy();centered[...,2]=0
        centered[...,3]=(total[...,4]-mu*total[...,1]).sum(axis=-1)
        centered[...,4]=(total[...,5]-mu*total[...,2]).sum(axis=-1)
        delta=float(Fraction('1152/845' if n in (65,130,260) else '2304/1445'))
        p0=result['estimates'][f'N{n}.p0']['value']
        for name,profile in (('early_raw',raw),('early_centered',centered)):
            def u_at_t(t):
                def get(p):
                    out=[]
                    for g in range(2):
                        v,d,_=moments(profile[g],m,p,n)
                        den=1+t*v[2];denp=t*d[2]
                        numerator=v[:2]+t*v[3:];nump=d[:2]+t*d[3:]
                        value=numerator/den;slope=(nump*den-numerator*denp)/den**2
                        out.append((value,slope))
                    return out
                p=brentq(lambda p:sum(x[0][0] for x in get(p))/2,p0-.01,p0+.01,xtol=5e-15,rtol=5e-15)
                packet=get(p)
                return n**(13/8)/2*(packet[0][1][1]-packet[1][1][1])/delta/np.mean([x[1][0] for x in packet])
            analytical=result['estimates'][f'N{n}.{name}.v']['value']
            values=[]
            for h in (1e-5,5e-6):values.append((u_at_t(h)-u_at_t(-h))/(2*h))
            err=abs(values[-1]-analytical)
            assert err<2e-5,(n,name,values,analytical)
            checks.append({'N':n,'source':name,'analytic_v':analytical,'symmetric_steps':[1e-5,5e-6],
                           'finite_differences':values,'last_absolute_error':err})
    out={'status':'passed','method':'finite-difference U at newly solved roots of the first-order (1+t*z) auxiliary measure; L(K) and its conditional source means stay fixed',
         'checks':checks,'root_solves':48,'production_replays':0,'new_random_samples':0,
         'scope':'algebraic derivative check on same archived moments, not finite-strength experimental evidence'}
    (ROOT/'results/moving_root_verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print('passed',len(checks),'max error',max(x['last_absolute_error'] for x in checks))
if __name__=='__main__':main()

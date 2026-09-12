#!/usr/bin/env python3
"""Bounded rank-source zeros and extensive-source crossover, width-four controls.

The zero-free strip is an elementary probability theorem, not a new Lee-Yang
circle theorem. Fixed-width spectral numbers use the exact PR710 definition.
"""
from __future__ import annotations
import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib,json
from pathlib import Path
import time
from mpmath import mp

DEF_SHA='5edc624377399ad0878c7a608e722ebb65a7be454a6f022a1b6382a0661e91b0'

def peval(c,x):
    v=0
    for a in reversed(c):v=v*x+a
    return v

def read_definition(path):
    raw=path.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=DEF_SHA:raise ValueError('wrong PR710 definition')
    return json.loads(raw)

def traces(co,maximum):
    d=len(co)-1;out=[mp.mpf(d)]
    for m in range(1,maximum+1):
        if m<=d:out.append(-sum(co[i]*out[m-i] for i in range(1,m))-m*co[m])
        else:out.append(-sum(co[i]*out[m-i] for i in range(1,d+1)))
    return out

class Spectra:
    def __init__(self,definition):self.d=definition
    def co(self,name,p):
        t=p/(1-p)
        return [peval(a,t)/(1+t)**(4*i) for i,a in enumerate(self.d['factors'][name]['coefficients_desc_x_ascending_t'])]
    def probabilities(self,p,maximum):
        seq={name:traces(self.co(name,p),maximum) for name in self.d['factors']}
        blocks={key:[sum(mult*seq[name][m] for name,mult in data['factorization']) for m in range(maximum+1)] for key,data in self.d['blocks'].items()}
        out={}
        for m in range(2,maximum+1):
            correction=blocks['16'][m]-2*(p*(1-p))**(2*m)
            p0,p2=blocks['5'][m]-correction,blocks['15'][m]-correction
            if not (p0>0 and p2>0 and p0+p2<1):raise ArithmeticError('working precision or invalid probabilities')
            out[m]=(p0,1-p0-p2,p2)
        return out
    def crossing(self):
        t=mp.findroot(lambda t:peval(self.d['crossing_polynomial_ascending_t'],t),('1.44','1.45'))
        return t/(1+t)
    def leading(self,p):
        c=self.co('closed_2',p)
        a=(-c[1]+mp.sqrt(c[1]*c[1]-4*c[2]))/2
        c=self.co('open_5',p)
        # Verify the selected positive root is maximal among independent companion eigenvalues.
        n=len(c)-1;A=mp.matrix(n)
        for i in range(n):A[0,i]=-c[i+1]
        for i in range(1,n):A[i,i-1]=1
        vals=mp.eig(A,left=False,right=False)
        b=max(vals,key=lambda x:mp.re(x))
        if abs(mp.im(b))>mp.mpf('1e-70'):raise ArithmeticError('nonreal leading root')
        return a,mp.re(b)

def source_zeros(probabilities):
    a,b,c=probabilities
    if min(a,b,c)<=0:raise ValueError('positive rank probabilities required')
    disc=b*b-4*a*c
    if disc>=0:
        temp=b+mp.sqrt(disc)
        roots=[-temp/(2*c),-2*a/temp]
    else:roots=[(-b+1j*mp.sqrt(-disc))/(2*c),(-b-1j*mp.sqrt(-disc))/(2*c)]
    values=[]
    for u in roots:
        scale=abs(c*u*u)+abs(b*u)+a
        residual=abs(c*u*u+b*u+a)/scale
        if residual>mp.mpf('1e-65'):raise ArithmeticError('root residual too large')
        z=mp.log(u)
        values.append({'u_real':mp.nstr(mp.re(u),35),'u_imag':mp.nstr(mp.im(u),35),
                       's_real':mp.nstr(mp.re(z),35),'s_imag':mp.nstr(mp.im(z),35),
                       'relative_polynomial_residual':mp.nstr(residual,8)})
    return {'discriminant':mp.nstr(disc,35),'type':'negative_real_pair' if disc>0 else ('double_negative_root' if disc==0 else 'complex_conjugate_pair'),
            'balance_source':mp.nstr(mp.log(a/c)/2,35),
            'centered_kappa':mp.nstr(b/(2*mp.sqrt(a*c)),35),'zeros':values}

def source_jet_identity_controls():
    # Polynomial-independent exact mgf controls, including a double source zero.
    examples=[(F(1,4),F(1,2),F(1,4)),(F(1,3),F(1,3),F(1,3)),(F(1,100),F(49,50),F(1,100))]
    out=[]
    for a,b,c in examples:
        disc=b*b-4*a*c
        out.append({'probabilities':list(map(str,(a,b,c))),'discriminant':str(disc),
                    'u_product':str(a/c),'u_sum':str(-b/c)})
    assert examples[0][1]**2-4*examples[0][0]*examples[0][2]==0
    return out

def report(path,dps=100):
    start=time.perf_counter();d=read_definition(path)
    rows=[]
    with mp.workdps(dps):
        model=Spectra(d);q=model.crossing()
        for label,p in [('half',mp.mpf('.5')),('q4',q),('seven_tenths',mp.mpf('.7'))]:
            l0,l2=model.leading(p)
            sigma0,sigma2=mp.log(l0),-mp.log(l2)
            probs=model.probabilities(p,256)
            z=[];phase=[]
            for m in (4,8,16,64,256):
                ps=probs[m]
                row=source_zeros(ps);row['length']=m;row['probabilities']=list(map(lambda x:mp.nstr(x,38),ps))
                row['scaled_source_zero_real_parts']=sorted((mp.nstr(mp.mpf(a['s_real'])/m,30) for a in row['zeros']), key=mp.mpf)
                z.append(row)
                for sigma in (sigma0-mp.mpf('.1'),mp.mpf(0),sigma2+mp.mpf('.1')):
                    logs=[mp.log(ps[0])-m*sigma,mp.log(ps[1]),mp.log(ps[2])+m*sigma]
                    mx=max(logs);zsum=sum(mp.exp(x-mx) for x in logs);logZ=mx+mp.log(zsum)
                    weights=[mp.exp(x-mx)/zsum for x in logs]
                    limit=max(mp.mpf(0),mp.log(l0)-sigma,mp.log(l2)+sigma)
                    phase.append({'length':m,'sigma':mp.nstr(sigma,30),'scaled_log_Z':mp.nstr(logZ/m,30),
                                  'limiting_max_formula':mp.nstr(limit,30),'absolute_error':mp.nstr(abs(logZ/m-limit),12),
                                  'tilted_rank_probabilities':list(map(lambda x:mp.nstr(x,30),weights))})
            rows.append({'parameter':label,'p':mp.nstr(p,50),'lambda0':mp.nstr(l0,40),'lambda2':mp.nstr(l2,40),
                         'negative_source_threshold':mp.nstr(sigma0,35),'positive_source_threshold':mp.nstr(sigma2,35),
                         'zero_controls':z,'extensive_source_controls':phase})
    return {'schema':'matching-one.rank-source-zeros.v1','definition_sha256':DEF_SHA,'dps':dps,
            'exact_quadratic_controls':source_jet_identity_controls(),'parameters':rows,
            'theorem':'For every law on {-1,0,1}, Re Z(x+iy)>0 when |y|<pi/2.',
            'nonclaims':['No real-axis zero pinch at bounded unscaled topological source.',
                         'Does not exclude a transition in occupation probability p.',
                         'Extensive source s=m*sigma defines a different tilted ensemble.',
                         'Numerical zeros and small-block Perron values are high-precision diagnostics, not interval certificates.'],
            'elapsed_seconds':time.perf_counter()-start}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--definition',type=Path,required=True)
    p.add_argument('--out',type=Path);p.add_argument('--dps',type=int,default=100);a=p.parse_args()
    text=json.dumps(report(a.definition,a.dps),indent=2)+'\n'
    if a.out:
        with a.out.open('x') as f:f.write(text)
    else:print(text,end='')
if __name__=='__main__':main()

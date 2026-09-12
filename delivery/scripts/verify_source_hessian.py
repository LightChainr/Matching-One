#!/usr/bin/env python3
"""Independent automaton + polynomial checks of full spatial source response.

Does not import the physical rank enumerator or Hessian generator. Compares
all 65,536 ranks, exact score derivatives at rational p, and finite-amplitude
bivariate field roots. Optional mpmath diagnostics are distinctly labelled.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib,json
from math import comb
from pathlib import Path
import time

EXPECTED='50b7297deefe7c50215aea2ed534ca5810461af3'

def read_cert(path):
    data=path.read_bytes()
    actual=hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
    if actual!=EXPECTED:raise ValueError('wrong PR708 certificate blob')
    return json.loads(data)

def all_ranks(cert):
    out=[]
    for mask in range(1<<16):
        s=cert['quotient_initial'][mask&15]
        for y in (1,2,3):s=cert['quotient_transitions'][s][(mask>>(4*y))&15]
        out.append(cert['quotient_rank_output'][s])
    return out

def eval_b(c,p):
    n=len(c)-1
    return sum((v*p**k*(1-p)**(n-k) for k,v in enumerate(c)),F(0))

def pattern(kx,ky):
    return [(1,1,-1,-1)[(kx*(v%4)+ky*(v//4))%4] for v in range(16)]

def group_counts(ranks,field):
    plus=sum(1<<i for i,x in enumerate(field) if x==1)
    minus=((1<<16)-1)^plus
    out=[[0]*9 for _ in range(9)]
    for mask,r in enumerate(ranks):out[(mask&plus).bit_count()][(mask&minus).bit_count()]+=r-1
    return out

def bivariate(counts,p,e):
    u,v=p+e,p-e
    return sum(c*u**a*(1-u)**(8-a)*v**b*(1-v)**(8-b)
               for a,row in enumerate(counts) for b,c in enumerate(row) if c)

def verify(certpath,hessianpath,dps=65):
    from mpmath import mp
    started=time.perf_counter()
    cert=read_cert(certpath); h=json.loads(hessianpath.read_text())
    ranks=all_ranks(cert)
    if hashlib.sha256(bytes(ranks)).hexdigest()!=h['physical_rank_sha256']:
        raise AssertionError('physical and automaton rank arrays differ')
    checks=0
    # Derive all off-diagonal derivatives by an independent likelihood-score sum.
    score_coeff=[]
    for d in range(1,16):
        rows=[[0]*17 for _ in range(4)]
        for mask,r in enumerate(ranks):
            i=2*(mask&1)+((mask>>d)&1)
            rows[i][mask.bit_count()]+=r-1
        score_coeff.append(rows)
        for p in (F(1,3),F(1,2),F(2,3)):
            actual=sum((eval_b(rows[i],p)*(i//2-p)*(i%2-p)/(p*p*(1-p)**2) for i in range(4)),F(0))
            if actual!=eval_b(h['origin_M_hessian_coefficients'][d],p):
                raise AssertionError('independent score derivative mismatch')
            checks+=1
    finite=[]
    with mp.workdps(dps):
        co=h['M_coefficients']
        def mp_b(c,p):return sum(mp.mpf(v)*p**k*(1-p)**(len(c)-1-k) for k,v in enumerate(c))
        root=mp.findroot(lambda p:mp_b(co,p),('.58','.60'))
        for kx,ky in ((1,0),(2,0),(1,1),(2,1),(2,2)):
            f=pattern(kx,ky)
            assert sum(f)==0 and sum(x*x for x in f)==16
            counts=group_counts(ranks,f)
            # Translation flips each selected field, so this equality is exact.
            assert counts==list(map(list,zip(*counts)))
            c=next(r['M_hessian_eigenvalue_coefficients'] for r in h['fourier_modes'] if r['wavevector']==[kx,ky])
            curvature=-16*mp_b(c,root)/mp_b(h['Mprime_coefficients'],root)
            values=[]
            for denominator in (64,128,256):
                eps=mp.mpf(1)/denominator
                shifted=mp.findroot(lambda p:bivariate(counts,p,eps),(root-mp.mpf('.01'),root+mp.mpf('.01')))
                estimate=2*(shifted-root)/eps**2
                values.append({'epsilon':f'1/{denominator}','root':mp.nstr(shifted,45),
                               'central_curvature':mp.nstr(estimate,35),
                               'absolute_error':mp.nstr(abs(estimate-curvature),15)})
            finite.append({'wavevector':[kx,ky],'field':f,'bivariate_M_counts':counts,
                           'exact_Hessian_curvature_diagnostic':mp.nstr(curvature,40),
                           'finite_amplitude_roots':values})
        # Full two-site disorder law. Remaining site parameters are p.
        # E_{xi0,xi1} M(p+sigma xi0,p+sigma xi1,p,..)=M(p) exactly by multiaffinity.
        d=1; rows=score_coeff[d-1]
        def two_site(p,e0,e1):
            a,b=p+e0,p+e1
            out=mp.mpf(0)
            for typ in range(4):
                i,j=typ//2,typ%2
                for k,c in enumerate(rows[typ]):
                    if c:out+=c*a**i*(1-a)**(1-i)*b**j*(1-b)**(1-j)*p**(k-i-j)*(1-p)**(14-k+i+j)
            return out
        # Rational exact validation of annealed equality (not root averaging).
        for p in (F(2,5),F(3,5)):
            for eps in (F(1,20),F(1,50)):
                # Rational independent expression, no floating cancellation.
                def rational_two(e0,e1):
                    a,b=p+e0,p+e1
                    return sum((c*a**(typ//2)*(1-a)**(1-typ//2)*b**(typ%2)*(1-b)**(1-typ%2)*p**(k-typ//2-typ%2)*(1-p)**(14-k+typ//2+typ%2)
                                for typ in range(4) for k,c in enumerate(rows[typ]) if c),F(0))
                assert sum(rational_two(i*eps,j*eps) for i in (-1,1) for j in (-1,1))/4==eval_b(co,p)
        # Two independently fluctuating sites: half trace on those two coordinates = R_00.
        bias=mp_b(h['Msecond_coefficients'],root)/(16**2*mp_b(h['Mprime_coefficients'],root))
        disorder=[]
        for denominator in (64,128,256):
            eps=mp.mpf(1)/denominator
            roots=[mp.findroot(lambda p:two_site(p,i*eps,j*eps),(root-mp.mpf('.01'),root+mp.mpf('.01'))) for i in (-1,1) for j in (-1,1)]
            actual=(sum(roots)/4-root)/eps**2
            disorder.append({'sigma':f'1/{denominator}','roots':list(map(lambda x:mp.nstr(x,40),roots)),
                             'mean_root_bias_over_sigma_squared':mp.nstr(actual,35),
                             'predicted_limit':mp.nstr(bias,35),'absolute_error':mp.nstr(abs(actual-bias),15)})
    return {'schema':'matching-one.full-site-hessian-independent.v1', 'source_certificate_git_blob':EXPECTED,
            'all_configurations_compared':len(ranks),'rational_score_derivative_checks':checks,
            'finite_amplitude_dps':dps,'field_controls':finite,
            'annealed_exact_controls':4,'two_site_disorder':disorder,
            'warning':'finite amplitude roots are high-precision numerical controls, not rational interval certificates',
            'elapsed_seconds':time.perf_counter()-started}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--certificate',type=Path,required=True);p.add_argument('--hessian',type=Path,required=True)
    p.add_argument('--out',type=Path);p.add_argument('--dps',type=int,default=65)
    a=p.parse_args();d=verify(a.certificate,a.hessian,a.dps);text=json.dumps(d,indent=2)+'\n'
    if a.out:
        with a.out.open('x') as f:f.write(text)
    else:print(text,end='')
if __name__=='__main__':main()

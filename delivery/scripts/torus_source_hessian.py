#!/usr/bin/env python3
"""Exact 4x4 site-source Hessians, Fourier signatures, and conditional covariance.

Standalone standard-library reference. Polynomial coefficients use the
unnormalized Bernstein convention sum c[k]*p**k*(1-p)**(n-k).
The exhaustive physical lifted graph is independent of the width-four automaton.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import hashlib
import json
from math import comb
from pathlib import Path
import time

N = 16
L = 4

def physical_rank(mask: int) -> int:
    if not isinstance(mask, int) or not 0 <= mask < (1 << N):
        raise ValueError('expected a 16-bit occupation mask')
    potentials = {}
    first = None
    for root in range(N):
        if root in potentials or not (mask >> root) & 1:
            continue
        potentials[root] = (0, 0)
        todo = [root]
        while todo:
            v = todo.pop()
            x, y = v % L, v // L
            ax, ay = potentials[v]
            for dx, dy in ((1,0),(-1,0),(0,1),(0,-1)):
                u = ((y+dy) % L)*L+(x+dx) % L
                if not (mask >> u) & 1:
                    continue
                candidate = (ax+dx, ay+dy)
                if u not in potentials:
                    potentials[u] = candidate
                    todo.append(u)
                else:
                    wx = candidate[0]-potentials[u][0]
                    wy = candidate[1]-potentials[u][1]
                    if wx % L or wy % L:
                        raise AssertionError('nonperiodic closing displacement')
                    if wx or wy:
                        if first is None:
                            first = (wx, wy)
                        elif first[0]*wy-first[1]*wx:
                            return 2
    return int(first is not None)

def bernstein(c, p):
    n = len(c)-1
    return sum((a*p**i*(1-p)**(n-i) for i,a in enumerate(c)), F(0))

def derivative(c):
    n = len(c)-1
    return [(i+1)*c[i+1]-(n-i)*c[i] for i in range(n)]

def bounds(c, lo, hi):
    n = len(c)-1
    low = high = F(0)
    for i,a in enumerate(c):
        u = lo**i*(1-hi)**(n-i)
        v = hi**i*(1-lo)**(n-i)
        low += a*(u if a>=0 else v)
        high += a*(v if a>=0 else u)
    return low, high

def signed_interval_div(a, b):
    if not b[0]>0:
        raise ValueError('denominator interval must be strictly positive')
    values=[x/y for x in a for y in b]
    return min(values), max(values)

def submasks(mask):
    current=mask
    while True:
        yield current
        if current==0:break
        current=(current-1)&mask

def build_counts(ranks):
    sectors=[[0]*(N+1) for _ in range(3)]
    pair=[[[0]*(N+1) for _ in range(N)] for _ in range(3)]
    for mask,r in enumerate(ranks):
        k=mask.bit_count()
        sectors[r][k]+=1
        if mask&1:
            for d in range(N):
                if (mask>>d)&1:pair[r][d][k]+=1
    hessian=[[0]*(N-1) for _ in range(N)]
    for d in range(1,N):
        rest=((1<<N)-1)^1^(1<<d)
        for mask in submasks(rest):
            k=mask.bit_count()
            delta=ranks[mask|1|(1<<d)]-ranks[mask|1]-ranks[mask|(1<<d)]+ranks[mask]
            hessian[d][k]+=delta
    matching=[b-a for a,b in zip(sectors[0],sectors[2])]
    return sectors,pair,hessian,matching

def cos4(k):return (1,0,-1,0)[k%4]

def modes(hessian):
    out=[]
    for ky in range(L):
        for kx in range(L):
            c=[sum(hessian[d][i]*cos4(kx*(d%L)+ky*(d//L)) for d in range(N)) for i in range(N-1)]
            out.append(((kx,ky),c))
    return out

def rational_root(matching, bits=100):
    lo,hi=F(0),F(1)
    for _ in range(bits):
        mid=(lo+hi)/2
        if bernstein(matching,mid)<0:lo=mid
        else:hi=mid
    assert bernstein(matching,lo)<0<bernstein(matching,hi)
    return lo,hi

def analyze(certificate: Path | None=None):
    start=time.perf_counter()
    ranks=[physical_rank(mask) for mask in range(1<<N)]
    checksum=hashlib.sha256(bytes(ranks)).hexdigest()
    automaton_checked=0
    if certificate is not None:
        cert=json.loads(certificate.read_text())
        trans=cert['quotient_transitions']
        initial=cert['quotient_initial'] if 'quotient_initial' in cert else cert['quotient_initial_classes']
        output=cert['quotient_rank_output']
        for mask,r in enumerate(ranks):
            word=[(mask>>(4*y))&15 for y in range(4)]
            s=initial[word[0]]
            for b in word[1:]:s=trans[s][b]
            if output[s]!=r:raise AssertionError('automaton comparison failed')
        automaton_checked=len(ranks)
    sectors,pairs,hs,matching=build_counts(ranks)
    assert all(sum(sectors[r][k] for r in range(3))==comb(N,k) for k in range(N+1))
    mp=derivative(matching); mpp=derivative(mp)
    assert [N*sum(hs[d][k] for d in range(N)) for k in range(N-1)]==mpp
    lo,hi=rational_root(matching)
    p=(lo+hi)/2
    mpI=bounds(mp,lo,hi)
    assert mpI[0]>0
    sp=[bernstein(c,p) for c in sectors]
    means=[bernstein([k*a for k,a in enumerate(c)],p)/sp[r] for r,c in enumerate(sectors)]
    deltaK=means[2]-means[0]
    covs=[]
    for r in (0,2):
        mu=means[r]/N
        covs.append([bernstein(pairs[r][d],p)/sp[r]-mu*mu for d in range(N)])
    dc=[b-a for a,b in zip(*covs)]
    rows=[]
    for (kx,ky),c in modes(hs):
        fi=bounds(c,lo,hi)
        eig=bernstein(c,p)
        additive=F(0) if (kx,ky)==(0,0) else -N*eig/bernstein(mp,p)
        if (kx,ky)==(0,0):curvI=(F(0),F(0))
        else:curvI=signed_interval_div((-N*fi[1],-N*fi[0]),mpI)
        # Zero-mean modes have z-first-derivative zero. Uniform source is a z gauge,
        # but its root reported in p has logistic curvature pq(1-2p).
        dceig=sum(dc[d]*cos4(kx*(d%L)+ky*(d//L)) for d in range(N))
        logit_z=(F(0) if (kx,ky)==(0,0) else -N*dceig/deltaK)
        logit=(p*(1-p)*(1-2*p) if (kx,ky)==(0,0) else p*(1-p)*logit_z)
        if (kx,ky)!=(0,0):
            # Equality is exact at M=0; rational midpoint has tiny M, so diagnostic tolerance.
            transformed=p*p*(1-p)**2*additive+p*(1-p)*(2*p-1)
            assert abs(float(logit-transformed))<1e-25
        if (kx,ky)!=(0,0):
            # Exact coordinate-conversion enclosure for the logit-field curvature.
            pqI=(lo*(1-hi),hi*(1-lo))
            products=[x*y*y for x in curvI for y in pqI]
            contact=(pqI[0]*(2*lo-1),pqI[1]*(2*hi-1))
            logitI=(min(products)+contact[0],max(products)+contact[1])
            if (logitI[0]>0)!=(curvI[0]>0) or (logitI[1]<0)!=(curvI[1]<0):
                raise AssertionError('coordinate change sign is not preserved in this control')
        else:logitI=(hi*(1-lo)*(1-2*hi),lo*(1-hi)*(1-2*lo))
        rows.append({'wavevector':[kx,ky], 'M_hessian_eigenvalue_coefficients':c,
                     'eigenvalue_at_half':str(bernstein(c,F(1,2))),
                     'eigenvalue_at_root':float(eig),'eigenvalue_root_interval':list(map(str,fi)),
                     'additive_p_root_second_derivative_unit_RMS':float(additive),
                     'additive_root_curvature_interval':list(map(str,curvI)),
                     'logit_source_p_root_second_derivative_unit_RMS':float(logit),
                     'logit_source_logit_root_second_derivative_unit_RMS':float(logit_z),
                     'logit_root_curvature_interval':list(map(str,logitI)),
                     'conditional_covariance_difference_eigenvalue':float(dceig)})
    trace_curv=bernstein(mpp,p)/(N*bernstein(mp,p))
    assert abs(sum(row['additive_p_root_second_derivative_unit_RMS']/N for row in rows)-float(trace_curv))<1e-14
    output={'schema':'matching-one.full-site-hessian.v1','scope':'4x4 NN square-site rank; p is occupation probability',
            'physical_configurations':len(ranks),'physical_rank_sha256':checksum,'automaton_configurations_checked':automaton_checked,
            'rank_bernstein_counts':sectors,'rank_totals':list(map(sum,sectors)),
            'origin_pair_sector_counts':pairs,'origin_M_hessian_coefficients':hs,'M_coefficients':matching,
            'Mprime_coefficients':mp,'Msecond_coefficients':mpp,
            'root_bracket':list(map(str,(lo,hi))),'root_diagnostic':float(p),'Mprime_root_interval':list(map(str,mpI)),
            'fourier_modes':rows,
            'conditional_K_means_at_root':list(map(float,means)),
            'delta_conditional_K_at_root':float(deltaK),
            'independent_additive_disorder_bias_per_variance':float(trace_curv/2),
            'root_Hessian_trace':float(trace_curv),
            'limits':['All Fourier amplitudes refer to unit-RMS real modes.',
                      'Nonzero wavevectors have zero first-order root response by translation invariance.',
                      'Logit source uses p_v=logistic(logit(p)+epsilon*h_v), not p+epsilon*h_v.',
                      'Root diagnostics are rational midpoint evaluations; signs use rational enclosures.',
                      'Finite 4x4 result; no infinite-volume exponent or original-U identification.'],
            'elapsed_seconds':time.perf_counter()-start}
    return output

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--certificate',type=Path)
    ap.add_argument('--out',type=Path)
    a=ap.parse_args(); result=analyze(a.certificate)
    text=json.dumps(result,indent=2,allow_nan=False)+'\n'
    if a.out:
        a.out.parent.mkdir(parents=True,exist_ok=True)
        with a.out.open('x') as f:f.write(text)
    else:print(text,end='')
if __name__=='__main__':main()

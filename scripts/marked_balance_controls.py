#!/usr/bin/env python3
"""Positive source reweightings preserve odds bounds but need not a unique root.

An exact 4x4 percolation example gives at least three balance roots under one
positive, p-independent, NONPRODUCT configuration weight. This is a counterexample
to finite monotonicity under arbitrary marks, not to the asymptotic root theorem.
"""
from __future__ import annotations
from fractions import Fraction as F
from pathlib import Path
import argparse,json,time
from oblique_torus_balance import Torus,NN,rank_lift,probability


def report():
    start=time.perf_counter();t=Torus(4,0,4);edges=t.edges(NN);n=t.n
    cross=sum(1<<i for i,(x,y) in enumerate(t.points) if x==0 or y==0)
    block=sum(1<<i for i,(x,y) in enumerate(t.points) if x!=0 and y!=0)
    assert cross.bit_count()==7 and rank_lift(cross,edges)==2
    assert block.bit_count()==9 and rank_lift(block,edges)==0
    bins=[[0]*(n+1) for _ in range(3)]
    for mask in range(1<<n):bins[rank_lift(mask,edges)][mask.bit_count()]+=1
    coeff=[bins[2][k]-bins[0][k] for k in range(n+1)]
    assert list(map(sum,bins))==[36559,19932,9045]
    A=10**12;checks=[];sgns=[]
    for p in (F(1,1000),F(1,4),F(3,4),F(999,1000)):
        P0=probability(bins[0],p);P2=probability(bins[2],p)
        w0=P0+A*p**9*(1-p)**7
        w2=P2+A*p**7*(1-p)**9
        numerator=w2-w0
        den=1+A*(p**7*(1-p)**9+p**9*(1-p)**7)
        assert w0>0 and w2>0 and den>0
        ratio=w2/w0;original=P2/P0
        assert original/F(A+1)<=ratio<=(A+1)*original
        sign=(numerator>0)-(numerator<0);sgns.append(sign)
        checks.append({'p':str(p),'unnormalized_balance':str(numerator),
                       'normalizer':str(den),'normalized_balance':str(numerator/den),
                       'sign':sign,'multiplicative_odds_bound_checked':True})
    assert sgns==[-1,1,-1,1]
    return {'schema':'matching-one.marked-balance-controls.v1','N':n,
            'scope':'axis 4x4, uniform base Bernoulli law, one positive nonproduct p-independent mark',
            'physical_configurations':1<<n,'cross_mask':cross,'rank0_block_mask':block,
            'cross_k_rank':[7,2],'block_k_rank':[9,0],
            'source_weight':'1+10^12*(1{omega=cross}+1{omega=3x3_block})',
            'weight_min':'1','weight_max':str(A+1),
            'base_matching_Bernstein_counts':coeff,'rank_counts_by_occupation':bins,
            'marked_numerator':'M(p)+10^12*p^7*(1-p)^7*(1-2p)',
            'rational_sign_checks':checks,'at_least_three_distinct_roots':True,
            'root_intervals':[['1/1000','1/4'],['1/4','3/4'],['3/4','999/1000']],
            'not_claimed':'This is not a product local-field source; arbitrary marks do not inherit finite monotonicity. All roots can still concentrate asymptotically under an oscillation bound.',
            'elapsed_seconds':time.perf_counter()-start}

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--out',type=Path)
    a=ap.parse_args();text=json.dumps(report(),indent=2,allow_nan=False)+'\n'
    if a.out:
        a.out.parent.mkdir(parents=True,exist_ok=True)
        with a.out.open('x') as f:f.write(text)
    else:print(text,end='')

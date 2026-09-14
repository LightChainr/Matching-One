#!/usr/bin/env python3
"""Small exact SITE interfaces for the no-merger proof, not a scaling experiment.

The cylinder has periodic horizontal and FREE vertical boundaries. A pair of
nested masks is the exact three-colour common-label ensemble. No plane/cylinder
asymptotic, Markov approximation, or Poisson fit is inferred from this census.
"""
from __future__ import annotations
import argparse
from collections import Counter
from fractions import Fraction
from math import comb
from pathlib import Path
import json


def edges(w: int, h: int):
    if w < 3 or h < 1:
        raise ValueError('Use an honest width >= 3 and positive height')
    out = [[] for _ in range(w*h)]
    for y in range(h):
        for x in range(w):
            v=y*w+x
            for dx,dy in ((1,0),(-1,0),(0,1),(0,-1)):
                if 0 <= y+dy < h:
                    out[v].append(((y+dy)*w+(x+dx)%w, dx))
    return out


def components_bfs(mask: int, adj):
    """Return (vertex mask, nonzero horizontal homology) for each component."""
    unseen=mask; result=[]
    while unseen:
        start=(unseen & -unseen).bit_length()-1
        stack=[start]; lift={start:0}; bits=0; winding=False
        unseen &= ~(1<<start)
        while stack:
            u=stack.pop(); bits |= 1<<u
            for v,dx in adj[u]:
                if not (mask>>v)&1: continue
                proposed=lift[u]+dx
                if v in lift:
                    winding |= lift[v] != proposed
                else:
                    lift[v]=proposed; unseen &= ~(1<<v); stack.append(v)
        result.append((bits,bool(winding)))
    return tuple(sorted(result))


def components_dsu(mask: int, adj):
    """Independent displacement union-find check."""
    n=len(adj); parent=list(range(n)); potential=[0]*n; wind=[False]*n
    def find(x):
        if parent[x] != x:
            r,d=find(parent[x]); potential[x]+=d; parent[x]=r
        return parent[x],potential[x]
    for u in range(n):
        if not (mask>>u)&1: continue
        for v,dx in adj[u]:
            if u>=v or not (mask>>v)&1: continue
            ru,du=find(u); rv,dv=find(v)
            if ru==rv: wind[ru] |= dv-du != dx
            else:
                parent[rv]=ru; potential[rv]=dx+du-dv
                wind[ru] |= wind[rv]
    groups={}
    for u in range(n):
        if (mask>>u)&1:
            r,_=find(u); groups[r]=groups.get(r,0)|(1<<u)
    return tuple(sorted((b,bool(wind[r])) for r,b in groups.items()))


def fmt(x: Fraction):
    return {'fraction':str(x),'decimal':float(x)}


def census(w=3,h=3):
    adj=edges(w,h); n=w*h
    if n>12: raise ValueError('This finite interface census is limited to 12 sites')
    total=1<<n
    cc=[components_bfs(m,adj) for m in range(total)]
    for m in range(total):
        assert cc[m]==components_dsu(m,adj)
        assert all(b.bit_count()>=w for b,wind in cc[m] if wind)
    essential=[tuple(b for b,wind in c if wind) for c in cc]
    # Inclusion-minimal increasing winding witnesses; they need not be entire clusters.
    minimal=[]
    for m in range(1,total):
        if essential[m] and all(not essential[m^(1<<v)] for v in range(n) if (m>>v)&1):
            minimal.append(m)
    double=[]
    for m in range(total):
        witnesses=[s for s in minimal if s&m==s]
        double.append(any(not(a&b) for i,a in enumerate(witnesses) for b in witnesses[i+1:]))
    hist=Counter(); configs=merging=0; witness=None
    # Sum over A subset B: each site is early / added / still closed.
    for late in range(total):
        early=late
        while True:
            nums=[sum(a&c==a for a in essential[early]) for c in essential[late]]
            mcount=sum(comb(k,2) for k in nums)
            loss=len(essential[early])-sum(k>0 for k in nums)
            assert sum(nums)==len(essential[early])
            assert 0<=loss<=mcount
            assert (not mcount) or double[late]
            assert mcount<=comb(h,2)
            key=(early.bit_count(),(late^early).bit_count(),n-late.bit_count())
            hist[key,'normalizer']+=1
            hist[key,'pairs']+=mcount
            hist[key,'loss']+=loss
            hist[key,'event']+=bool(mcount)
            if mcount:
                merging+=1
                if witness is None or late.bit_count()<witness['late'].bit_count():
                    witness={'early':early,'late':late,'ancestor_counts':nums,'pairs':mcount}
            configs+=1
            if early==0: break
            early=(early-1)&late
    probabilities=[]
    for lo,hi in ((Fraction(1,5),Fraction(1,4)),(Fraction(1,4),Fraction(1,3)),(Fraction(2,5),Fraction(3,5))):
        sums={kind:Fraction(0) for kind in ('normalizer','pairs','loss','event')}
        for (key,kind),count in hist.items():
            e,a,c=key; sums[kind]+=count*lo**e*(hi-lo)**a*(1-hi)**c
        assert sums['normalizer']==1
        q=sum(hi**m.bit_count()*(1-hi)**(n-m.bit_count()) for m in range(total) if essential[m])
        d=sum(hi**m.bit_count()*(1-hi)**(n-m.bit_count()) for m in range(total) if double[m])
        assert d<=q*q
        assert sums['event']<=sums['pairs']<=comb(h,2)*d
        probabilities.append({'p_minus':fmt(lo),'p_plus':fmt(hi),
            'merger_probability':fmt(sums['event']),
            'factorial_pair_expectation':fmt(sums['pairs']),
            'lineage_loss_expectation':fmt(sums['loss']),
            'one_witness_probability_at_p_plus':fmt(q),
            'two_disjoint_witnesses_probability_at_p_plus':fmt(d),
            'BK_upper':fmt(q*q)})
    return {'width':w,'height':h,'vertical_boundary':'free',
        'physical_masks_crosschecked':total,'common_label_pairs':configs,
        'pairs_with_a_merger':merging,'minimal_winding_witnesses':len(minimal),
        'explicit_merger':witness,'exact_probabilities':probabilities,
        'scope':'Finite common-label and BK interfaces only; no estimate of an infinite-cylinder density.'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    data=census()
    text=json.dumps(data,ensure_ascii=False,indent=2)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True); args.output.write_text(text,encoding='utf8')
    else: print(text,end='')

if __name__=='__main__': main()

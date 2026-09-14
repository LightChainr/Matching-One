#!/usr/bin/env python3
"""Independent finite controls for shielded query noise and width-eight cavities.

All probabilities/moments used in identities are Fractions. No previous
transfer engine, sampling, or large-width fit is imported. The enumerations
check finite interfaces; the accompanying arguments prove the limits.
"""
from __future__ import annotations
import argparse
from collections import Counter, deque
from fractions import Fraction as F
from itertools import product, combinations
import json
from pathlib import Path
from typing import Iterable

XY = tuple[int, int]
NN = ((1,0),(-1,0),(0,1),(0,-1))
KING = tuple((x,y) for x in (-1,0,1) for y in (-1,0,1) if x or y)
P_COEFF = (9,18,27,36,45,54,63,8,7,6,5,4,3,2,1)

def packed(x: F | int) -> dict[str, str]:
    f = F(x)
    return {'fraction':str(f), 'decimal':format(float(f),'.16g')}

def covariance(rows: list[tuple[F, tuple[F, ...]]]) -> tuple[list[F], list[list[F]]]:
    assert sum((p for p,_ in rows),F()) == 1
    d=len(rows[0][1])
    mu=[sum((p*x[j] for p,x in rows),F()) for j in range(d)]
    cv=[[sum((p*(x[i]-mu[i])*(x[j]-mu[j]) for p,x in rows),F())
         for j in range(d)] for i in range(d)]
    return mu,cv

def wired_statistics(mask: int) -> tuple[int,int,int]:
    """Actual KING BFS: 3x3 random interior, radius-two all-white shell."""
    inside=tuple(product(range(-1,2),repeat=2))
    shell={v for v in product(range(-2,3),repeat=2) if max(map(abs,v))==2}
    white=shell|{v for i,v in enumerate(inside) if mask>>i&1}
    reached=set(shell); queue=deque(shell)
    while queue:
        x,y=queue.popleft()
        for dx,dy in KING:
            z=x+dx,y+dy
            if z in white and z not in reached:
                reached.add(z); queue.append(z)
    k=len(set(inside)&reached)
    boundary={v for v in inside if v not in white
              and any((v[0]+dx,v[1]+dy) in reached for dx,dy in KING)}
    return k,len(boundary),9-k-len(boundary)

def shield_control(q: F) -> dict:
    if not 0<q<1: raise ValueError('0<q<1 required')
    p=1-q; rows=[]
    inside=tuple(product(range(-1,2),repeat=2)); centre=inside.index((0,0))
    for mask in range(512):
        n=mask.bit_count(); z=(mask>>centre)&1; t=n-z
        k,b,h=wired_statistics(mask)
        assert (k,b,h)==(t+z*(t>0),8-t+(1-z)*(t>0),int(t==0))
        score=F(k)/q-F(b)/p
        rows.append((q**n*p**(9-n),(F(k),score,F(k+b),F(h))))
    mu,cv=covariance(rows); r=p**8
    vS=(9-r)/(p*q); vK=p*q*(9-r)+q*q*r*(17-r)
    kS=9-r+8*q*p**7
    residual=vK-kS*kS/vS
    poly=sum((F(a)*p**i for i,a in enumerate(P_COEFF)),F())
    positive=p**8*q**4*poly/(9-p**8)
    assert cv[1][1]==vS and cv[0][0]==vK and cv[0][1]==kS
    assert mu[1]==0 and mu[3]==r and residual==positive>0
    return {'q':packed(q),'p':packed(p),'configurations':512,
            'mean_K':packed(mu[0]),'variance_K':packed(vK),
            'variance_score':packed(vS),'covariance_K_score':packed(kS),
            'conditional_residual':packed(residual),
            'global_lower_bound_per_theta':packed(q**16*residual/25),
            'local_hole_probability':packed(r)}

def canon(shape: Iterable[XY], width: int) -> tuple[XY,...]:
    pts=tuple(shape); lo=min(y for x,y in pts)
    return min(tuple(sorted((((x+s)%width,y-lo) for x,y in pts)))
               for s in range(width))

def grow_animals(width: int, maximum: int) -> list[set[tuple[XY,...]]]:
    out=[set(),{((0,0),)}]
    for n in range(2,maximum+1):
        nxt=set()
        for animal in out[-1]:
            a=set(animal)
            for x,y in animal:
                for dx,dy in NN:
                    v=((x+dx)%width,y+dy)
                    if v not in a: nxt.add(canon(a|{v},width))
        out.append(nxt)
    return out

def winds(shape: Iterable[XY], width: int) -> bool:
    pts={(x%width,y) for x,y in shape}; lifts={}; winding=False
    for root in pts:
        if root in lifts: continue
        lifts[root]=root[0]; todo=[root]
        while todo:
            u=todo.pop(); x,y=u
            for dx,dy in NN:
                v=((x+dx)%width,y+dy)
                if v not in pts: continue
                target=lifts[u]+dx
                if v in lifts: winding |= lifts[v]!=target
                else: lifts[v]=target;todo.append(v)
    return winding

def hidden_sites(shape: Iterable[XY], width: int) -> tuple[XY,...]:
    """Black shape in otherwise white cylinder; exterior wired at both ends."""
    black={(x%width,y) for x,y in shape}; hi=max(y for x,y in black);lo=min(y for x,y in black)
    allowed={(x,y) for x in range(width) for y in range(lo-1,hi+2)}-black
    reached={v for v in allowed if v[1] in (lo-1,hi+1)};todo=list(reached)
    while todo:
        x,y=todo.pop()
        for dx,dy in KING:
            v=((x+dx)%width,y+dy)
            if v in allowed and v not in reached:
                reached.add(v);todo.append(v)
    queried=set(reached)
    for x,y in reached:
        queried.update(((x+dx)%width,y+dy) for dx,dy in KING)
    return tuple(sorted((x,y) for x in range(width) for y in range(lo,hi+1)
                        if (x,y) not in queried))

def animal_control(width: int=8, maximum: int=8) -> dict:
    animals=grow_animals(width,maximum);summary=[];wit=[]
    for n in range(1,maximum+1):
        essential=[];holes=[]
        for a in sorted(animals[n]):
            if winds(a,width): essential.append(a)
            else:
                h=hidden_sites(a,width)
                if h: holes.append((a,h))
        if n<8: assert not holes
        if n<width: assert not essential
        if width==8 and n==8:
            assert len(essential)==1 and set(essential[0])=={(x,0) for x in range(8)}
            ring=canon(KING,width)
            assert len(holes)==1 and holes[0][0]==ring and len(holes[0][1])==1
            wit=[{'black_shape':a,'hidden_sites':h} for a,h in holes]
        summary.append({'black_sites':n,'translation_classes':len(animals[n]),
                        'essential_classes':len(essential),'cavity_classes':len(holes)})
    return {'width':width,'maximum_black_sites':maximum,'counts':summary,
            'hole_witnesses':wit,'checked_classes':sum(len(a) for a in animals)}

def pattern_control() -> dict:
    w=8
    types=[('barrier',frozenset((x,0) for x in range(w)))]
    types += [('hole',frozenset(((x+dx)%w,dy) for dx,dy in KING)) for x in range(w)]
    overlap=Counter();sizes=Counter()
    for i,(ti,si) in enumerate(types):
        for j,(tj,sj0) in enumerate(types):
            for dy in range(-2,3):
                if i==j and dy==0: continue
                sj={(x,y+dy) for x,y in sj0}
                if si & sj:
                    k=len(si|sj);overlap[(ti,tj,k)]+=1;sizes[k]+=1
    assert min(sizes)>=12
    return {'types_per_row':{'barrier':1,'hole':8},
            'ordered_overlap_pairs':[{'first':a,'second':b,'union_black_sites':k,'count':n}
                for (a,b,k),n in sorted(overlap.items())],
            'least_distinct_union':min(sizes),
            'finite_window_pattern_poisson_error_order':'O(p^4) for length O(p^-8)',
            'physical_reduction_bad_cluster_order':'O(p), on a fixed scaled window'}

def limiting_laws() -> dict:
    r=F(8)
    probabilities=[F(1,9)*F(8,9)**k for k in range(9)]
    mixed=[]
    for u in (F(0),F(1,2),F(2)):
        for z in (F(0),F(1,2),F(1)):
            f=1/(1+u+r*(1-z))
            # Competing exponential clocks / factorial-series algebra.
            assert f==1/(9+u-8*z)
            mixed.append({'u':packed(u),'z':packed(z),'transform':packed(f)})
    return {'interpretation':'Exact moments of limiting two-pattern model; not a finite-p moment certificate.',
            'H_probabilities_k_0_to_8':[packed(x) for x in probabilities],
            'mean_H':packed(r),'variance_H':packed(r*(1+r)),
            'covariance_E_H':packed(r),'correlation_E_H_squared':packed(r/(1+r)),
            'mean_limit_query_excess_16_minus_H':packed(16-r),
            'cavity_fugacity_pole_limiting_model':packed(F(9,8)),
            'joint_transforms':mixed}

def boundary_of(sites: set[XY]) -> frozenset[XY]:
    return frozenset({(x+dx,y+dy) for x,y in sites for dx,dy in KING}-sites)

def defect_boolean_terms(v: XY, degree: int=11) -> dict[frozenset[XY],int]:
    """OR of the one singleton cage and four axial two-site cages.

Return its square-free black-indicator polynomial, truncated by number of
required black sites. This is not claimed to contain cages of perimeter >=12.
"""
    patterns=[boundary_of({v})]
    patterns += [boundary_of({v,(v[0]+dx,v[1]+dy)}) for dx,dy in NN]
    out=Counter()
    for n in range(1,6):
        for inds in combinations(range(5),n):
            union=frozenset().union(*(patterns[i] for i in inds))
            if len(union)<=degree: out[union]+=(-1)**(n+1)
    return {k:v for k,v in out.items() if v}

def sparse_defect_series() -> dict:
    degree=11; left=defect_boolean_terms((0,0),degree)
    mean=Counter()
    for sites,coef in left.items(): mean[len(sites)]+=coef
    spectra=Counter(); offsets=[]
    # Outside this range supports are disjoint. Their product has degree >=16.
    for dx,dy in product(range(-5,6),repeat=2):
        right=defect_boolean_terms((dx,dy),degree);cov=Counter()
        for u,a in left.items():
            for v,b in right.items():
                if len(u|v)<=degree:cov[len(u|v)]+=a*b
        cov={k:v for k,v in cov.items() if v}
        if cov:
            offsets.append({'offset':(dx,dy),'coefficients':dict(sorted(cov.items()))})
            spectra.update(cov)
    assert dict(mean)=={8:1,10:4,11:-4}
    assert len(offsets)==5 and dict(spectra)=={8:1,10:8,11:-4}
    residual=Counter()
    for k,v in spectra.items():
        for j,b in enumerate((1,-2,1)):
            if k+j<=degree:residual[k+j]+=v*b
    assert dict(residual)=={8:1,9:-2,10:9,11:-20}
    return {'model':'Exact Boolean motif algebra through black degree 11. '
                    'The proof of the O(p^12) physical remainder is in the note.',
            'mean_defect_coefficients':dict(sorted(mean.items())),
            'covariance_by_offset':offsets,
            'integrated_query_covariance':dict(sorted(spectra.items())),
            'J_coefficients':dict(sorted(residual.items())),
            'thermal_projection_first_possible_order':15}

def independent_block_control() -> list[dict]:
    """Exact rational control of the rare-clock limit; explicitly a toy."""
    ans=[]
    for p in (F(1,2),F(1,4),F(1,8)):
        rho=p**8
        for u,z in ((F(0),F(0)),(F(1),F(1,2)),(F(2),F(1))):
            a=(1-rho+rho*z)**8; s=1/(1+u*rho)
            val=rho*s/(1-(1-rho)*s*a)
            lim=1/(1+u+8*(1-z))
            ans.append({'p':packed(p),'u':packed(u),'z':packed(z),
                        'rational_discrete_transform':packed(val),
                        'limiting_transform':packed(lim),
                        'absolute_difference':packed(abs(val-lim))})
    return ans

def report() -> dict:
    return {'scope':'finite exact controls, independent of historical engines',
            'shield':[shield_control(q) for q in (F(1,3),F(1,2),F(3,4),F(9,10))],
            'positive_polynomial_coefficients_ascending':P_COEFF,
            'animals':animal_control(), 'patterns':pattern_control(),
            'limiting_model':limiting_laws(), 'defect_series':sparse_defect_series(),
            'independent_block_toy':independent_block_control()}

def main() -> None:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();result=report()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    print(json.dumps({'output':str(args.output),'shield_configurations':2048,
                      'animal_classes':result['animals']['checked_classes']}))

if __name__=='__main__': main()

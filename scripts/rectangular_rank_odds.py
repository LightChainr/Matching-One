#!/usr/bin/env python3
"""Finite controls for aspect-uniform homological-balance consistency.

The all-width theorem is in the accompanying proof, not an extrapolation of
these tiny tables. Uses exact integers/Fractions; no Monte Carlo, no p_c fit.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import json
from math import log, log1p
from pathlib import Path
import time


def steps(matching: bool):
    out=[(1,0),(-1,0),(0,1),(0,-1)]
    if matching:out += [(1,1),(1,-1),(-1,1),(-1,-1)]
    return out


def lifted_rank(mask: int, width: int, length: int, matching: bool=False) -> int:
    """Physical integer lifts, independent of local-arm / slab classification."""
    if min(width,length)<2 or mask<0 or mask>>(width*length):
        raise ValueError('invalid periodic graph')
    pos={};direction=None
    for root in range(width*length):
        if root in pos or not (mask>>root)&1:continue
        pos[root]=(0,0);stack=[root]
        while stack:
            v=stack.pop();x,y=v%width,v//width;px,py=pos[v]
            for dx,dy in steps(matching):
                u=((y+dy)%length)*width+(x+dx)%width
                if not (mask>>u)&1:continue
                candidate=(px+dx,py+dy)
                if u not in pos:pos[u]=candidate;stack.append(u)
                else:
                    a,b=candidate[0]-pos[u][0],candidate[1]-pos[u][1]
                    assert a%width==0 and b%length==0
                    a,b=a//width,b//length
                    if a or b:
                        if direction is None:direction=(a,b)
                        elif direction[0]*b-direction[1]*a:return 2
    return int(direction is not None)


def local_arm(mask: int,w: int,m: int,v: int,radius: int,matching: bool) -> bool:
    if radius<1 or 2*radius+2>min(w,m):raise ValueError('local box is not embedded')
    x0,y0=v%w,v//w
    if not (mask>>v)&1:return False
    seen={(0,0)};stack=[(0,0)]
    while stack:
        x,y=stack.pop()
        if max(abs(x),abs(y))==radius:return True
        for dx,dy in steps(matching):
            a,b=x+dx,y+dy
            if max(abs(a),abs(b))>radius or (a,b) in seen:continue
            u=((y0+b)%m)*w+(x0+a)%w
            if (mask>>u)&1:seen.add((a,b));stack.append((a,b))
    return False


def slab_crossing(mask: int,w: int,m: int,y0: int,radius: int,matching: bool) -> bool:
    if not 0<=y0<=m-radius-1:raise ValueError('slab crosses the closing seam')
    seen={(x,0) for x in range(w) if (mask>>(y0*w+x))&1};stack=list(seen)
    while stack:
        x,y=stack.pop()
        if y==radius:return True
        for dx,dy in steps(matching):
            a,b=(x+dx)%w,y+dy
            if not 0<=b<=radius or (a,b) in seen:continue
            if (mask>>((y0+b)*w+a))&1:seen.add((a,b));stack.append((a,b))
    return False


def path_arm_bound(degree: int,p: F,radius: int) -> F:
    if degree<2 or not 0<p<1 or radius<1:raise ValueError('invalid path bound inputs')
    return min(F(1),degree*(degree-1)**(radius-1)*p**(radius+1))


def sign_certificate(width: int,radius: int,alpha: F) -> dict:
    """An upper bound alpha on the unconditional local one-arm probability.

    Exact positive margin implies P2/P0 <= exp(-gamma*m)<1 for EVERY m>=w.
    No claim is made that any arbitrary supplied alpha bounds a physical law.
    """
    if type(width) is not int or type(radius) is not int or radius<1 or width<2*radius+2:
        raise ValueError('require 1<=R and 2R+2<=width')
    if not 0<alpha<1:raise ValueError('alpha must lie strictly between 0 and 1')
    q=width*alpha
    exponent=2*width*(radius+1)
    margin=(1-alpha)**exponent-q
    holds=(q<1 and margin>0)
    return {'width':width,'radius':radius,'one_arm_upper_bound':str(alpha),
            'slab_upper_bound':str(q),'comparison_exponent':exponent,
            'exact_margin_positive':margin>0,'uniform_sign_certified':holds,
            'exact_margin_expression':'(1-alpha)^(2*w*(R+1))-w*alpha',
            'margin_numerator_bits':abs(margin.numerator).bit_length(),
            'margin_denominator_bits':margin.denominator.bit_length(),
            'gamma_diagnostic':(-log(float(q))/(2*(radius+1))+width*log1p(-float(alpha))) if q<1 else None,
            'scope':'all integer lengths m>=width, if the supplied one-arm upper bound is valid'}


def evaluate_counts(counts: list[int],p: F) -> tuple[F,F]:
    n=len(counts)-1;value=F(0);derivative=F(0)
    for k,c in enumerate(counts):
        value+=c*p**k*(1-p)**(n-k)
        if k:derivative+=c*k*p**(k-1)*(1-p)**(n-k)
        if k<n:derivative-=c*(n-k)*p**k*(1-p)**(n-k-1)
    return value,derivative


def fisher_decomposition(p0: F,p2: F,d0: F,d2: F) -> dict:
    e=p0+p2
    if not 0<p0<1 or not 0<p2<1 or not e<1:raise ValueError('positive probabilities needed')
    h=p2/e;de=d0+d2;dh=(d2*e-p2*de)/e**2
    direct=d0*d0/p0+d2*d2/p2+de*de/(1-e)
    activity=de*de/(e*(1-e));conditional=e*dh*dh/(h*(1-h))
    assert direct==activity+conditional
    assert dh>0 or d0==d2==0
    return {'E':str(e),'conditional_H':str(h),'H_derivative':str(dh),
            'rank_Fisher_information':str(direct),'activity_information':str(activity),
            'conditional_information_per_original_snapshot':str(conditional),
            'exact_decomposition':True}


def tiny_controls() -> dict:
    w=m=4;radius=1;n=w*m;allmask=(1<<n)-1;rows=[];ranks={}
    for matching in (False,True):
        counts=[[0]*(n+1) for _ in range(3)];noarms=[0]*(n+1)
        slab_counts=[0]*(2*w+1);stored=[]
        for mask in range(1<<n):
            rank=lifted_rank(mask,w,m,matching);stored.append(rank)
            k=mask.bit_count();counts[rank][k]+=1
            arms=any(local_arm(mask,w,m,v,radius,matching) for v in range(n))
            if not arms:
                assert rank==0;noarms[k]+=1
            if rank==2:
                assert all(slab_crossing(mask,w,m,y,radius,matching) for y in (0,2))
        for mask in range(1<<(2*w)):
            if slab_crossing(mask,w,m,0,1,matching):slab_counts[mask.bit_count()]+=1
        comparisons=[]
        for p in (F(1,10),F(1,3),F(1,2)):
            p0,d0=evaluate_counts(counts[0],p);p2,d2=evaluate_counts(counts[2],p)
            independent_set,_=evaluate_counts(noarms,p)
            exact_slab,_=evaluate_counts(slab_counts,p)
            degree=8 if matching else 4
            arm=p*(1-(1-p)**degree)  # exact R=1 local arm
            lower=(1-arm)**n;upper=min(F(1),w*arm)**2
            assert p0>=independent_set>=lower
            assert p2<=exact_slab**2<=upper
            comparisons.append({'p':str(p),'P0':str(p0),'P2':str(p2),
                 'all_no_arm_probability':str(independent_set),'FKG_P0_lower':str(lower),
                 'exact_slab_crossing_probability':str(exact_slab),'independent_slab_P2_upper':str(exact_slab**2),
                 'one_arm_slab_P2_upper':str(upper),'inequalities_exact':True,
                 'information':fisher_decomposition(p0,p2,d0,d2)})
        ranks[matching]=stored
        rows.append({'graph':'NN+NNN' if matching else 'NN','width':w,'length':m,
                     'configurations_checked':1<<n,'rank_totals':list(map(sum,counts)),
                     'rank_counts_by_occupation':counts,'comparisons':comparisons})
    assert all(ranks[False][mask]+ranks[True][allmask^mask]==2 for mask in range(1<<n))
    # Larger supports, but a fixed structural set, not another census or Monte Carlo.
    structural=0
    for w,m,r in ((6,8,2),(8,12,3),(16,32,4)):
        masks=[0,(1<<(w*m))-1]
        for y in range(m):masks.append(((1<<w)-1)<<(w*y))
        for x in range(w):masks.append(sum(1<<(w*y+x) for y in range(m)))
        for shift in range(w):masks.append(sum(1<<(w*y+(y+shift)%w) for y in range(m)))
        for mask in masks:
            for matching in (False,True):
                rank=lifted_rank(mask,w,m,matching)
                if rank>0:assert any(local_arm(mask,w,m,v,r,matching) for v in range(w*m))
                if rank==2:assert all(slab_crossing(mask,w,m,y,r,matching)
                                      for y in range(0,m-r,r+1))
                structural+=1
    return {'exhaustive_graph_configurations':2*(1<<16),'models':rows,
            'complementary_rank_sum_checked':1<<16,'additional_structural_graph_cases':structural}


def report() -> dict:
    start=time.perf_counter();certificates=[]
    for degree,p,w,r in ((4,F(1,10),16,2),(8,F(1,20),16,2),(4,F(1,4),256,32)):
        cert=sign_certificate(w,r,path_arm_bound(degree,p,r))
        assert cert['uniform_sign_certified']
        cert.update({'p':str(p),'graph_degree':degree,'input_certificate':'simple-path union bound; no unknown empirical constant'})
        certificates.append(cert)
    return {'schema':'matching-one.rectangle-rank-odds-controls.v1',
       'status':'finite controls for a proof using Harris association and subcritical one-arm decay',
       'sign_certificates':certificates,'tiny_controls':tiny_controls(),
       'limits':['No numerical p_c interval or rate near criticality is inferred.',
                 'The all-aspect root theorem is proved in the accompanying note, not from these sizes.',
                 'Fisher decomposition concerns independent one-p snapshots, not full permutation paths.'],
       'elapsed_seconds':time.perf_counter()-start}


def main() -> None:
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--out',type=Path);args=ap.parse_args()
    text=json.dumps(report(),indent=2,allow_nan=False)+'\n'
    if args.out:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        with args.out.open('x',encoding='utf-8') as f:f.write(text)
    else:print(text,end='')

if __name__=='__main__':main()

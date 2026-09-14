#!/usr/bin/env python3
"""Finite motifs and exact controls for a monotone percolation cut/mark limit.

No Monte Carlo and no exponential-height simulation.  The script verifies
finite interfaces; it is not a numerical proof of a scaling limit.
"""
from __future__ import annotations

import argparse
from collections import Counter, deque
from fractions import Fraction as F
from itertools import combinations, product
import json
from math import comb, factorial, prod
from pathlib import Path
from typing import Iterable, NamedTuple

Point = tuple[int, int]

class Motif(NamedTuple):
    kind: str
    sites: frozenset[Point]


def steps(matching: bool) -> tuple[Point, ...]:
    if matching:
        return tuple((x, y) for x in (-1, 0, 1) for y in (-1, 0, 1) if x or y)
    return ((1, 0), (-1, 0), (0, 1), (0, -1))


def normalized(sites: Iterable[Point], width: int) -> frozenset[Point]:
    points = frozenset((x % width, y) for x, y in sites)
    low = min(y for _, y in points)
    return frozenset((x, y-low) for x, y in points)


def motifs(width: int, black_matching: bool) -> tuple[Motif, ...]:
    """The two resonant pairs: (w=4, black king), (w=8, black NN)."""
    if (width, black_matching) not in ((4, True), (8, False)):
        raise ValueError('Only the two stated resonance pairs are supported')
    barrier = set()
    for word in product((-1, 0, 1) if black_matching else (0,), repeat=width):
        if sum(word):
            continue
        y = 0
        points = []
        for x, jump in enumerate(word):
            points.append((x, y))
            y += jump
        barrier.add(normalized(points, width))
    result = [Motif('barrier', s) for s in sorted(barrier, key=lambda a: sorted(a))]
    for x in range(width):
        s = normalized(((x+dx, dy) for dx, dy in steps(not black_matching)), width)
        result.append(Motif('hole', s))
    assert len({m.sites for m in result}) == len(result)
    return tuple(result)


def components(sites: Iterable[Point], width: int, matching: bool):
    """Lifted BFS. Stores physical x-displacements, not wrap booleans."""
    remaining = set(sites)
    answer = []
    while remaining:
        root = min(remaining)
        seen = {root: 0}
        queue = deque([root])
        wound = False
        while queue:
            x, y = queue.popleft()
            for dx, dy in steps(matching):
                z = ((x+dx) % width, y+dy)
                if z not in remaining:
                    continue
                value = seen[(x, y)] + dx
                if z in seen:
                    wound |= seen[z] != value
                else:
                    seen[z] = value
                    queue.append(z)
        remaining.difference_update(seen)
        answer.append((frozenset(seen), wound))
    return answer


def dsu_winds(sites: Iterable[Point], width: int, matching: bool) -> bool:
    """Independent potential-DSU check of the physical winding flag."""
    vertices = tuple(sorted(sites)); index = {v:i for i,v in enumerate(vertices)}
    par = list(range(len(vertices))); delta = [0]*len(vertices)
    def find(i):
        if par[i] != i:
            root, d = find(par[i]); delta[i] += d; par[i] = root
        return par[i], delta[i]
    flag = False
    for a, (x,y) in enumerate(vertices):
        for dx,dy in steps(matching):
            z = ((x+dx) % width, y+dy)
            if z not in index:
                continue
            b = index[z]; ra,da=find(a); rb,db=find(b)
            if ra==rb:
                flag |= db-da != dx
            else:
                par[rb]=ra; delta[rb]=dx+da-db
    return flag


def physical_minimal_control() -> dict:
    width, height = 4, 5
    vertices = tuple(product(range(width), range(height)))
    universe = set(product(range(width), range(-1,height+1)))
    catalog = motifs(4, True)
    count=0; winding_sets=0; hole_sets=0; empty_hole_count=0
    for n in range(5):
        for raw in combinations(vertices,n):
            black = frozenset(raw)
            cc = components(black,width,True)
            winding = any(w for _,w in cc)
            assert winding == dsu_winds(black,width,True)
            normalized_black = normalized(black,width) if black else frozenset()
            predicted = any(m.kind=='barrier' and m.sites==normalized_black for m in catalog)
            assert winding == predicted
            white_cc = components(universe-black,width,False)
            essential = set().union(*(set(s) for s,w in white_cc if w))
            query = set(essential)
            for x,y in essential:
                for dx,dy in steps(False):
                    z=((x+dx)%width,y+dy)
                    if z in black:
                        query.add(z)
            found = set(vertices)-query
            expected = {v for v in vertices if all(((v[0]+dx)%width,v[1]+dy) in black
                                                    for dx,dy in steps(False))}
            assert found == expected
            count += 1
            winding_sets += winding
            hole_sets += bool(found)
            empty_hole_count += len(found)
    return dict(configurations=count, independent_black_winding_checks=count,
                complete_white_explorations=count, winding_sets=winding_sets,
                hole_sets=hole_sets, hole_centres=empty_hole_count)


def dependency_polynomials(width: int, matching: bool) -> dict:
    ms=motifs(width,matching)
    b1=0; b2=Counter(); typed=Counter()
    for i,a in enumerate(ms):
        for j,b in enumerate(ms):
            amin=min(y for _,y in a.sites); amax=max(y for _,y in a.sites)
            bmin=min(y for _,y in b.sites); bmax=max(y for _,y in b.sites)
            for shift in range(amin-bmax,amax-bmin+1):
                shifted=frozenset((x,y+shift) for x,y in b.sites)
                if not a.sites & shifted:
                    continue
                b1+=1
                if i==j and shift==0:
                    continue
                assert a.sites != shifted
                u=len(a.sites|shifted)
                b2[u]+=1; typed[(a.kind,b.kind,u)]+=1
    return dict(width=width, black_matching=matching,
                barrier_templates=sum(m.kind=='barrier' for m in ms),
                hole_templates=sum(m.kind=='hole' for m in ms),
                size=width, b1_coefficient=b1,
                b1_power=2*width,
                b2={str(k):v for k,v in sorted(b2.items())},
                smallest_union=min(b2),
                typed=[dict(first=a,second=b,union=u,count=n) for (a,b,u),n in sorted(typed.items())])


def intensities(width: int, matching: bool, a: tuple[F,...]):
    if len(a)!=width or any(x<=0 for x in a):
        raise ValueError('Need one positive amplitude per column')
    vals={k:sum((prod((a[x] for x,y in m.sites),start=F(1))
                 for m in motifs(width,matching) if m.kind==k),F(0)) for k in ('barrier','hole')}
    vals['ratio']=vals['hole']/vals['barrier']
    return vals


def two_time(k:F, r:F, z:F, v:F, s:F=F(0), t:F=F(0)) -> F:
    """E exp(-s Y_1 -t Y_k) z^H_1 v^H_k for a fixed-location tag."""
    if k<1 or r<0 or not 0<=z<=1 or not 0<=v<=1 or s<0 or t<0:
        raise ValueError('Need k>=1,r,s,t>=0 and z,v in [0,1]')
    a=1+s+r*(1-z); b=k*t+r*(1-v)*(z+k-1)
    return ((a+k-1)/(a*(a+k-1+b)))**2


def two_time_by_integrals(k,r,z,v,s=F(0),t=F(0)):
    """Independent no-new-record / new-record integral decomposition."""
    a=1+s+r*(1-z); b=k*t+r*(1-v)*(z+k-1); c=k-1
    # min(x, Z), x~Exp(1), Z~Exp(k-1). Keep the atom Z>=x.
    no_jump=1/(a+b+c)
    jump=F(0) if c==0 else c/(a*(a+b+c))
    return (no_jump+jump)**2


# Two-variable truncated Taylor algebra. Coefficients are not derivatives.
DEGREES=tuple((i,j) for i in range(3) for j in range(3))
def jadd(a,b):
    return {ij:a.get(ij,F(0))+b.get(ij,F(0)) for ij in DEGREES}
def jscale(a,c): return {ij:v*c for ij,v in a.items()}
def jmul(a,b):
    out={ij:F(0) for ij in DEGREES}
    for (i,j),x in a.items():
        for (k,l),y in b.items():
            if i+k<=2 and j+l<=2: out[(i+k,j+l)]+=x*y
    return out

def jinv(a):
    a0=a[(0,0)]
    if not a0: raise ZeroDivisionError('Zero constant Taylor coefficient')
    q={ij:F(0) for ij in DEGREES}; q[(0,0)]=1/a0
    for n in range(1,5):
        for i,j in DEGREES:
            if i+j!=n: continue
            q[(i,j)]=-sum((a.get((k,l),0)*q.get((i-k,j-l),0)
                            for k in range(i+1) for l in range(j+1) if k or l),F(0))/a0
    return q


def mark_moments(k:F,r:F) -> dict:
    one={(0,0):F(1)}; a={(0,0):F(1),(1,0):-r}
    c={(0,0):k,(1,0):-r}
    d={(0,0):k,(1,0):-r,(0,1):-r*k,(1,1):-r}
    side=jmul(c,jinv(jmul(a,d))); total=jmul(side,side)
    mu1=total[(1,0)]; mu2=total[(0,1)]
    var=2*total[(2,0)]+mu1-mu1**2
    cov=total[(1,1)]-mu1*mu2
    assert mu1==mu2==2*r
    assert var==2*r*(1+r)
    assert cov==2*r*(1+r)/k
    return dict(mean=mu1,variance=var,covariance=cov,correlation=cov/var if var else None)


def generator_moments(n:int) -> F:
    # Stationarity against Exp(1): A x^n = n*x^n-n*x^(n+1)/(n+1).
    return F(n*factorial(n))-F(n, n+1)*factorial(n+1)


def adjoint_moment(m:int,n:int) -> tuple[F,F]:
    """<x^m,A x^n>_Exp = <A_reverse x^m,x^n>_Exp."""
    forward=F(n*factorial(m+n))-F(n,n+1)*factorial(m+n+1)
    reverse=-F(m*factorial(m+n))+sum((F(comb(m,j)*factorial(j)*factorial(m+n-j))
                                      for j in range(1,m+1)),F(0))
    return forward,reverse


def motif_multitime_check() -> int:
    count=0
    # Exact common-label threshold probability for two overlapping black motifs.
    # P(A completed at t1, B completed at t2), t1<=t2.
    for width,matching in ((4,True),(8,False)):
        ms=motifs(width,matching)
        eps=F(1,20); t1=F(1,2); t2=F(3,2)
        for a in ms:
            for b in ms:
                shifted=frozenset((x,y+1) for x,y in b.sites)
                union=a.sites|shifted
                direct=prod((eps*min([t for motif,t in ((a.sites,t1),(shifted,t2)) if v in motif])
                             for v in union), start=F(1))
                separated=(eps*t1)**len(a.sites)*(eps*t2)**len(shifted-a.sites)
                assert direct==separated
                count+=1
    return count


def report() -> dict:
    profiles=[]
    for w,m in ((4,True),(8,False)):
        for changed in (False,True):
            a=[F(1)]*w
            if changed: a[0]=F(2);a[w//2]=F(1,2)
            v=intensities(w,m,tuple(a))
            assert v['barrier']==(19 if m else 1)
            profiles.append(dict(width=w,black_matching=m,amplitudes=a,**v,
                                 component_no_hole=1/(1+v['ratio']),
                                 tagged_no_hole=1/(1+v['ratio'])**2))
    integral_checks=0
    for k,r,z,v,s,t in product(map(F,[1,2,4]),(F(4,19),F(8)),
                              (F(0),F(1,2),F(1)),(F(0),F(1,3),F(1)),
                              (F(0),F(1,2)),(F(0),F(1,4))):
        assert two_time(k,r,z,v,s,t)==two_time_by_integrals(k,r,z,v,s,t)
        integral_checks+=1
    moments=[]
    for r,k in product((F(4,19),F(8)),(F(1),F(2),F(4))):
        moments.append(dict(ratio=r,time_ratio=k,**mark_moments(k,r)))
    assert all(generator_moments(n)==0 for n in range(1,9))
    assert all(adjoint_moment(m,n)[0]==adjoint_moment(m,n)[1] for m,n in product(range(7),repeat=2))
    return dict(scope='Finite exact interfaces; limits are proved in the companion notes.',
                physical_control=physical_minimal_control(),
                motif_dependencies=[dependency_polynomials(4,True),dependency_polynomials(8,False)],
                profiles=profiles, two_time_integral_equalities=integral_checks,
                two_time_moments=moments, generator_stationarity_orders=list(range(1,9)),
                common_label_threshold_equalities=motif_multitime_check(),
                time_reverse_polynomial_identities=49)


def jsonable(x):
    if isinstance(x,F): return dict(fraction=str(x),decimal=float(x))
    if isinstance(x,dict): return {str(k):jsonable(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)): return [jsonable(v) for v in x]
    return x


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    text=json.dumps(jsonable(report()),ensure_ascii=False,indent=2,sort_keys=True)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(text,encoding='utf-8')
    else: print(text,end='')

if __name__=='__main__': main()

#!/usr/bin/env python3
"""Exact geometry controls for balance roots on arbitrary square-lattice tori.

The asymptotic theorem is proved in the accompanying note; this program checks
finite geometry, rank/crossing implications, and rational probability bounds.
It does not estimate a critical probability or fit an exponent.
"""
from __future__ import annotations
import argparse
from dataclasses import dataclass
from fractions import Fraction as F
from math import comb, isqrt
from pathlib import Path
import json
import time

NN = ((1,0),(-1,0),(0,1),(0,-1))
MATCHING = NN + ((1,1),(1,-1),(-1,1),(-1,-1))

def det(u, v): return u[0]*v[1]-u[1]*v[0]
def dot(u, v): return u[0]*v[0]+u[1]*v[1]
def norm2(u): return dot(u,u)
def nearest_integer(q: F) -> int:
    """Nearest integer, with a fixed tie convention (towards +infinity)."""
    return (2*q.numerator+q.denominator)//(2*q.denominator)

def reduced_basis(a, b):
    """Lagrange/Gauss reduction, using integer arithmetic only.

    Output u is shortest; v is a reduced completion. The proof is in the note.
    Neither u nor v need be primitive in ambient Z^2.
    """
    u, v = tuple(a), tuple(b)
    if det(u,v)==0: raise ValueError('periods must be independent')
    for _ in range(10000):
        if norm2(v)<norm2(u): u,v=v,u
        k=nearest_integer(F(dot(u,v),norm2(u)))
        if k==0: break
        v=(v[0]-k*u[0],v[1]-k*u[1])
    else: raise RuntimeError('lattice reduction failed to terminate')
    if det(u,v)<0: v=(-v[0],-v[1])
    assert norm2(u)<=norm2(v) and 2*abs(dot(u,v))<=norm2(u)
    assert 4*det(u,v)**2>=3*norm2(u)**2
    return u,v

@dataclass(frozen=True)
class Torus:
    """HNF columns (a,0),(b,c), 0<=b<a; exactly a*c vertices."""
    a: int
    b: int
    c: int
    def __post_init__(self):
        if type(self.a) is not int or type(self.b) is not int or type(self.c) is not int:
            raise TypeError('HNF entries must be integers')
        if self.a<1 or self.c<1 or not 0<=self.b<self.a:
            raise ValueError('require a,c>=1 and 0<=b<a')
    @property
    def n(self): return self.a*self.c
    @property
    def points(self): return [(x,y) for y in range(self.c) for x in range(self.a)]
    def reduce(self,x,y):
        q,r=divmod(y,self.c)
        return (x-q*self.b)%self.a,r
    def index(self,x,y):
        x,y=self.reduce(x,y); return y*self.a+x
    def basis(self): return reduced_basis((self.a,0),(self.b,self.c))
    def edges(self,steps=NN):
        return [[(self.index(x+dx,y+dy),dx,dy) for dx,dy in steps]
                for x,y in self.points]


def rank_lift(mask, edges):
    """Physical Euclidean-lift traversal, independent of reduced coordinates."""
    positions={}; first=None
    for root in range(len(edges)):
        if not (mask>>root)&1 or root in positions: continue
        positions[root]=(0,0); stack=[root]
        while stack:
            a=stack.pop(); x,y=positions[a]
            for b,dx,dy in edges[a]:
                if not (mask>>b)&1: continue
                proposed=(x+dx,y+dy)
                if b not in positions:
                    positions[b]=proposed; stack.append(b)
                else:
                    delta=(proposed[0]-positions[b][0],proposed[1]-positions[b][1])
                    if delta!=(0,0):
                        if first is None: first=delta
                        elif det(first,delta): return 2
    return int(first is not None)


def slab_geometry(torus: Torus, steps=NN, *, coarse=False):
    """Exact transverse bands in q(x)=det(u,x) modulo N.

    The theorem uses B=|u|^2/8. Tiny controls use B=2D+|u|^2/32,
    still satisfying (B-2D)/|u| > r, with r=|u|/64.
    The offset 1/17 avoids integer vertex levels and 1/8,1/32 boundaries.
    """
    u,v=torus.basis(); s=norm2(u); n=torus.n
    D=max(abs(det(u,e)) for e in steps)
    B=F(s,8) if coarse else 2*D+F(s,32)
    k=int(F(n)/B)
    if k<1: raise ValueError('chosen bands do not fit')
    if B-2*D<=F(s,64): raise ValueError('slab is too thin for the local arm radius')
    # Local first-exit support fits in radius r+sqrt(2). Use a rational
    # upper bound 3/2 instead of sqrt(2) for the exact injectivity check.
    if F(s)*(F(1,2)-F(1,64))**2<=F(9,4):
        raise ValueError('control torus too short for this local support')
    qs=[(F(det(u,x))-F(1,17))%n for x in torus.points]
    edges=torus.edges(steps)
    slabs=[]
    for j in range(k):
        low=j*B; high=(j+1)*B
        nodes=[i for i,q in enumerate(qs) if low<q<high]
        bottom=[i for i in nodes if qs[i]<low+D]
        top=[i for i in nodes if qs[i]>high-D]
        inside=set(nodes); adj={}
        for i in nodes:
            adj[i]=[z for z,dx,dy in edges[i] if z in inside
                    and qs[z]-qs[i]==det(u,(dx,dy))]
        assert len(bottom)<=4*(isqrt(s)+1)
        slabs.append({'nodes':nodes,'bottom':bottom,'top':top,'adj':adj})
    joined=[i for slab in slabs for i in slab['nodes']]
    assert len(joined)==len(set(joined))
    return {'u':u,'v':v,'systole_squared':s,'D':D,'B':B,'count':k,'slabs':slabs}


def has_crossing(mask, slab):
    target=set(slab['top']); seen={i for i in slab['bottom'] if (mask>>i)&1}
    todo=list(seen)
    while todo:
        i=todo.pop()
        if i in target: return True
        for j in slab['adj'][i]:
            if (mask>>j)&1 and j not in seen: seen.add(j);todo.append(j)
    return False


def probability(counts, p):
    n=len(counts)-1
    return sum((F(c)*p**k*(1-p)**(n-k) for k,c in enumerate(counts)),F(0))

def sign_certificate(torus, degree, p):
    """A rational sign test using elementary self-avoiding-path domination.

    For r=ell/64 and max step<=sqrt(2), at least L edges are needed.
    Use a slightly smaller integer L with L^2*8192 <= ell^2.
    a_r <= d (d-1)^(L-1) p^(L+1), valid for L>=1.
    """
    u,v=torus.basis(); s=norm2(u); n=torus.n
    if s<64**2: raise ValueError('asymptotic geometry certificate needs ell>=64')
    L=isqrt(s//8192)
    if L<1: raise ValueError('path certificate needs at least one edge')
    alpha=degree*(degree-1)**(L-1)*p**(L+1)
    B=F(s,8); k=int(F(n)/B)
    M=4*(isqrt(s)+int(isqrt(s)**2<s))
    q=M*alpha
    assert 4*n<=k*s
    r=(n+k-1)//k
    # Bernoulli: (1-alpha)^r >= 1-r*alpha; avoids enormous powers.
    good=(0<alpha<1 and q<1-r*alpha)
    return {'hnf':[torus.a,torus.b,torus.c],'shortest_u':u,'completion_v':v,
            'systole_squared':s,'slabs':k,'entry_count_upper':M,'path_length_lower':L,
            'degree':degree,'p':str(p),'arm_upper':str(alpha),
            'all_lengths_not_claimed_from_one_instance':True,
            'ceil_N_over_slabs':r,'bernoulli_lower':str(1-r*alpha),
            'strict_sign_certified':good}


def check_small(torus):
    n=torus.n; allmask=(1<<n)-1
    models={}
    ranks=[]
    for name,steps in [('NN',NN),('matching',MATCHING)]:
        edges=torus.edges(steps); geometry=slab_geometry(torus,steps)
        bins=[[0]*(n+1) for _ in range(3)]
        crossing_bins=[[0]*(n+1) for _ in geometry['slabs']]
        rs=[]
        for mask in range(1<<n):
            r=rank_lift(mask,edges); k=mask.bit_count();bins[r][k]+=1;rs.append(r)
            for idx,slab in enumerate(geometry['slabs']):
                cross=has_crossing(mask,slab)
                crossing_bins[idx][k]+=int(cross)
                if r==2 and not cross: raise AssertionError(('rank2 without slab',torus,mask))
        checks=[]
        for p in [F(1,32),F(1,16),F(1,8)]:
            a=p*(1-(1-p)**len(steps)) # r<1, exactly the occupied-neighbour arm.
            assert geometry['systole_squared']<64**2
            P0=probability(bins[0],p);P2=probability(bins[2],p)
            probs=[probability(c,p) for c in crossing_bins]
            assert P0>=(1-a)**n
            product=F(1)
            for e,slab in zip(probs,geometry['slabs']):
                assert e<=len(slab['bottom'])*a
                product*=e
            assert P2<=product
            checks.append({'p':str(p),'P0':str(P0),'P2':str(P2),
                           'harris_lower':str((1-a)**n),'crossing_product':str(product)})
        models[name]={'degree':len(steps),'rank_totals':list(map(sum,bins)),
                     'rank_counts_by_occupation':bins,'slabs':geometry['count'],
                     'entry_counts':[len(t['bottom']) for t in geometry['slabs']],
                     'exact_probability_checks':checks}
        ranks.append(rs)
    assert all(ranks[0][m]+ranks[1][allmask^m]==2 for m in range(1<<n))
    u,v=torus.basis()
    return {'hnf':[torus.a,torus.b,torus.c],'N':n,'u':u,'v':v,
            'systole_squared':norm2(u),'configurations':1<<n,'models':models,
            'complementary_rank_identity':True}


def check_reductions(maxdet=80):
    count=0
    for a in range(1,maxdet+1):
      for c in range(1,maxdet//a+1):
       for b in range(a):
        t=Torus(a,b,c);u,v=t.basis();s=norm2(u)
        assert det(u,v)==t.n
        # Exhaust all ambient integer vectors inside the candidate shortest box.
        R=isqrt(s)
        for x in range(-R,R+1):
         for y in range(-R,R+1):
          if (x,y)!=(0,0) and norm2((x,y))<s:
           assert t.reduce(x,y)!=(0,0)
        count+=1
    return count


def check_large():
    cases=[Torus(64,0,80),Torus(65,17,130),Torus(4225,268,1),Torus(4273,1200,1)]
    # Last two are Gaussian ideals (63+16i), (32+57i); reduction checks the systole.
    out=[]
    for t in cases:
        u,v=t.basis();s=norm2(u)
        if s<4096: continue
        masks=[0,(1<<t.n)-1]
        pts=t.points
        for q in range(8):
            masks.append(sum(1<<i for i,(x,y) in enumerate(pts)
                             if ((17*x+23*y+q)%11)<q+2))
        for steps in (NN,MATCHING):
            geo=slab_geometry(t,steps,coarse=True);edges=t.edges(steps)
            for m in masks:
                if rank_lift(m,edges)==2:
                    assert all(has_crossing(m,b) for b in geo['slabs'])
            out.append({'hnf':[t.a,t.b,t.c],'N':t.n,'systole_squared':s,
                        'degree':len(steps),'slabs':geo['count'],'fixed_masks':len(masks),
                        'all_slabs_vertex_disjoint':True})
    return out


def report(small=True):
    start=time.perf_counter()
    reductions=check_reductions(50)
    cases=[Torus(10,3,1),Torus(13,5,1),Torus(5,2,3),Torus(4,1,4)] if small else []
    smallout=[check_small(t) for t in cases]
    certs=[sign_certificate(Torus(4096,113,4096),d,p)
           for d,p in [(4,F(1,10)),(8,F(1,20))]]
    assert all(c['strict_sign_certified'] for c in certs)
    return {'schema':'matching-one.oblique-balance-controls.v1',
            'scope':'arbitrary HNF square-cell tori, root odds not full birth CDF',
            'reduced_hnf_bases_checked':reductions,
            'small_controls':smallout,'graph_configuration_checks':2*sum(x['configurations'] for x in smallout),
            'complement_pair_checks':sum(x['configurations'] for x in smallout),
            'large_fixed_structure_controls':check_large(),
            'finite_subcritical_sign_certificates':certs,
            'elapsed_seconds':time.perf_counter()-start}

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',type=Path);ap.add_argument('--skip-small',action='store_true')
    args=ap.parse_args();text=json.dumps(report(not args.skip_small),indent=2,allow_nan=False)+'\n'
    if args.out:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        with args.out.open('x') as f:f.write(text)
    else:print(text,end='')

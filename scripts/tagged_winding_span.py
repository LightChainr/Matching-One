#!/usr/bin/env python3
"""All-height winding-component span by a single tagged lineage.

Only colours, connectivity and winding are stored. No age/depth cutoff.
The source is a finite two-row Bernoulli experiment, not a stationary solve.
For a finite closed state set, d_h = alpha R**(h-1) b is exact.
"""
from __future__ import annotations
import argparse
from collections import Counter
from fractions import Fraction
import json
from pathlib import Path
from time import perf_counter
from typing import NamedTuple

NEUTRAL, TAG, FORBIDDEN = 0, 1, 2

class State(NamedTuple):
    labels: tuple[int, ...]
    gains: tuple[int, ...]
    winding: tuple[int, ...]
    colours: tuple[int, ...]

class DSU:
    def __init__(self, n: int):
        self.par=list(range(n)); self.delta=[0]*n
        self.wind=[False]*n; self.col=[0]*n
    def find(self, a):
        if self.par[a] != a:
            r,g=self.find(self.par[a]); self.delta[a]+=g; self.par[a]=r
        return self.par[a],self.delta[a]
    def join(self,a,b,gain):
        ra,da=self.find(a); rb,db=self.find(b)
        if ra==rb: self.wind[ra] |= (db-da!=gain)
        else:
            self.par[rb]=ra; self.delta[rb]=gain+da-db
            self.wind[ra] |= self.wind[rb]; self.col[ra] |= self.col[rb]

def empty(width):
    if width < 2: raise ValueError('width must be >=2; lifted parallel edges retained')
    return State((-1,)*width,(0,)*width,(),())

def step(state: State, mask: int, matching=False, tagged=True, connected=False):
    """Return (next state, outcome). outcome=0 continue,1 accept,-1 reject.

    In source mode tagged=False, colours propagate but retirement is ignored.
    """
    w=len(state.labels)
    if not 0 <= mask < 1<<w: raise ValueError('invalid mask')
    d=DSU(2*w); reps={}
    old=[i for i,k in enumerate(state.labels) if k>=0]
    new=[i for i in range(w) if mask>>i&1]
    for i in old:
        k=state.labels[i]
        if k in reps: d.join(reps[k],i,state.gains[i])
        else: reps[k]=i
    for k,i in reps.items():
        root,_=d.find(i); d.wind[root]=bool(state.winding[k]); d.col[root]=state.colours[k]
    for i in new:
        j=(i+1)%w
        if mask>>j&1: d.join(w+i,w+j,(i+1)//w)
        for dx in ((-1,0,1) if matching else (0,)):
            j=(i+dx)%w
            if state.labels[j]>=0: d.join(w+i,j,(i+dx)//w)
    roots={d.find(i)[0] for i in old+[w+i for i in new]}
    kept={d.find(w+i)[0] for i in new}
    if connected and (roots-kept or not kept):
        return None, 1 if (not kept and len(roots)==1 and d.wind[next(iter(roots))]) else -1
    if tagged:
        marked=[r for r in roots if d.col[r]&TAG]
        if len(marked)!=1: raise AssertionError('exactly one live tag required')
        r=marked[0]
        if d.col[r]&FORBIDDEN: return None,-1
        if r not in kept: return None,1 if d.wind[r] else -1
    labels=[-1]*w; gains=[0]*w; wind=[]; colours=[]; canon={}
    for i in new:
        r,g=d.find(w+i)
        if r not in canon:
            canon[r]=(len(wind),g); wind.append(int(d.wind[r])); colours.append(d.col[r])
        k,g0=canon[r]; labels[i]=k; gains[i]=0 if d.wind[r] else g-g0
    return State(tuple(labels),tuple(gains),tuple(wind),tuple(colours)),0

def source_entries(width, matching=False):
    """Integer multiplicities indexed by state and occupancy of two source rows.

    All old row (-1) components are forbidden. A row-zero candidate is tagged;
    smaller row-zero candidate components are also forbidden (tie breaking).
    """
    out=Counter()
    for prev in range(1<<width):
        old,_=step(empty(width),prev,matching,False)
        old=old._replace(colours=(FORBIDDEN,)*len(old.winding))
        for mask in range(1<<width):
            cur,_=step(old,mask,matching,False)
            candidates=[k for k,c in enumerate(cur.colours) if not c]
            for k in candidates:
                col=tuple(TAG if j==k else FORBIDDEN if (c or j<k) else NEUTRAL
                          for j,c in enumerate(cur.colours))
                out[(cur._replace(colours=col),prev.bit_count()+mask.bit_count())]+=1
    return out

def build(width, matching=False, state_cap=100000):
    """Complete reachable tagged automaton; abort rather than truncate."""
    if width>8: raise ValueError('reference source enumeration limited to width <=8')
    sources=source_entries(width,matching)
    states=[]; index={}
    def add(s):
        if s not in index:
            if len(states)>=state_cap: raise RuntimeError('state cap reached')
            index[s]=len(states); states.append(s)
        return index[s]
    for s,k in sources: add(s)
    transitions=[]
    for s in states:
        row=[]
        for mask in range(1<<width):
            nxt,out=step(s,mask,matching)
            row.append((add(nxt) if out==0 else (-1 if out==1 else -2),out))
        transitions.append(row)
    alpha=Counter({(index[s],k):c for (s,k),c in sources.items()})
    return states,transitions,alpha

def lump(transitions, source):
    """All-p exit-labelled probabilistic bisimulation, not a minimality claim."""
    blocks=[0]*len(transitions)
    while True:
        lookup={}; nxt=[]
        for row in transitions:
            sig=tuple(sorted(Counter((mask.bit_count(), blocks[j] if j>=0 else j)
                                     for mask,(j,_) in enumerate(row)).items()))
            if sig not in lookup: lookup[sig]=len(lookup)
            nxt.append(lookup[sig])
        if nxt==blocks: break
        blocks=nxt
    reps=[blocks.index(i) for i in range(max(blocks)+1)]
    reduced=[[(blocks[j] if j>=0 else j,o) for j,o in transitions[i]] for i in reps]
    src=Counter()
    for (i,k),c in source.items(): src[(blocks[i],k)]+=c
    return reduced,src,blocks

def numeric_system(transitions, source, p, exact=True):
    n=len(transitions); w=(len(transitions[0])-1).bit_length()
    zero=Fraction(0) if exact else 0.0
    pp=Fraction(p) if exact else float(p)
    if not 0 < pp < 1: raise ValueError("requires 0<p<1")
    weights=[pp**mask.bit_count()*(1-pp)**(w-mask.bit_count()) for mask in range(1<<w)]
    R=[[zero for _ in range(n)] for _ in range(n)]; b=[zero]*n; alpha=[zero]*n
    for i,row in enumerate(transitions):
        for wt,(j,_) in zip(weights,row):
            if j>=0: R[i][j]+=wt
            elif j==-1: b[i]+=wt
    for (i,k),count in source.items(): alpha[i]+=count*pp**k*(1-pp)**(2*w-k)
    delta=(1-pp)**w
    if exact:
        assert all(sum(row)<=1-delta for row in R)
        assert all(sum(row)+bb<=1 for row,bb in zip(R,b))
    return alpha,R,b,delta

def solve(A, b):
    n=len(b); M=[[Fraction(x) for x in row]+[Fraction(y)] for row,y in zip(A,b)]
    for j in range(n):
        pivot=next((i for i in range(j,n) if M[i][j]),None)
        if pivot is None: raise ValueError('singular matrix')
        M[j],M[pivot]=M[pivot],M[j]; v=M[j][j]; M[j]=[x/v for x in M[j]]
        for i in range(n):
            if i!=j and M[i][j]:
                v=M[i][j]; M[i]=[a-v*b for a,b in zip(M[i],M[j])]
    return [row[-1] for row in M]

def moments(transitions, source, p, bins=12, exact=True):
    alpha,R,b,delta=numeric_system(transitions,source,p,exact)
    n=len(b)
    if exact:
        A=[[int(i==j)-R[i][j] for j in range(n)] for i in range(n)]
        z1=solve(A,b); z2=solve(A,z1); z3=solve(A,z2)
        dot=lambda v:sum(a*x for a,x in zip(alpha,v))
        nu=dot(z1); mean=dot(z2)/nu; second=dot([2*x-y for x,y in zip(z3,z2)])/nu
        v=alpha[:]; dh=[]
        for _ in range(bins):
            dh.append(sum(x*y for x,y in zip(v,b)))
            v=[sum(v[i]*R[i][j] for i in range(n)) for j in range(n)]
        tail=sum(x*y for x,y in zip(v,z1))
        t1=sum(x*(bins*y+z) for x,y,z in zip(v,z1,z2))
        t2=sum(x*(bins*bins*y+2*bins*z+2*u-z) for x,y,z,u in zip(v,z1,z2,z3))
        assert sum(dh)+tail==nu
        assert sum((i+1)*x for i,x in enumerate(dh))+t1==nu*mean
        assert sum((i+1)**2*x for i,x in enumerate(dh))+t2==nu*second
    else:
        import numpy as np
        from scipy.linalg import lu_factor, lu_solve
        rr=np.asarray(R); aa=np.asarray(alpha); bb=np.asarray(b)
        lu=lu_factor(np.eye(n)-rr)
        z1=lu_solve(lu,bb); z2=lu_solve(lu,z1); z3=lu_solve(lu,z2)
        nu=float(aa@z1); mean=float(aa@z2)/nu; second=float(aa@(2*z3-z2))/nu
        v=aa.copy(); dh=[]
        for _ in range(bins): dh.append(float(v@bb)); v=v@rr
        tail=float(v@z1); t1=float(v@(bins*z1+z2)); t2=float(v@(bins*bins*z1+(2*bins-1)*z2+2*z3))
    var=second-mean*mean
    return dict(nu=nu,mean=mean,second=second,variance=var,cv2=var/mean**2,
                d_h=dh,tail=tail,tail_first=t1,tail_second=t2,delta=delta)

def jsonable(x):
    if isinstance(x,Fraction): return str(x)
    if isinstance(x,dict):return {str(k):jsonable(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [jsonable(v) for v in x]
    return x

def certified_moments(transitions, source, p, bins=12):
    """Exact rational residual enclosures after one floating correction per solve.

    Returned centres are rational corrected solutions, not the pre-rounding
    float vectors. The inverse infinity norm is at most 1/delta.
    """
    import numpy as np
    from scipy.linalg import lu_factor, lu_solve
    alpha,R,b,delta=numeric_system(transitions,source,Fraction(p),True)
    n=len(b); A=[[Fraction(int(i==j))-R[i][j] for j in range(n)] for i in range(n)]
    af=np.array([[float(x) for x in row] for row in A]); lu=lu_factor(af)
    def exact_residual(x,rhs):
        return [rhs[i]-sum(v*y for v,y in zip(A[i],x)) for i in range(n)]
    def corrected(rhs):
        xf=lu_solve(lu,np.array([float(x) for x in rhs]))
        x=[Fraction.from_float(float(v)) for v in xf]
        r=exact_residual(x,rhs)
        dx=lu_solve(lu,np.array([float(v) for v in r]))
        x=[v+Fraction.from_float(float(d)) for v,d in zip(x,dx)]
        r=exact_residual(x,rhs)
        return x,max(map(abs,r))
    x1,r1=corrected(b); e1=r1/delta
    x2,r2=corrected(x1); e2=(r2+e1)/delta
    x3,r3=corrected(x2); e3=(r3+e2)/delta
    mass=sum(alpha); dot=lambda x:sum(a*b for a,b in zip(alpha,x))
    centres=[dot(x1),dot(x2),dot([2*x-y for x,y in zip(x3,x2)])]
    errors=[mass*e1,mass*e2,mass*(2*e3+e2)]
    raw=[(c-e,c+e) for c,e in zip(centres,errors)]
    if raw[0][0]<=0:raise ArithmeticError('density not separated from zero')
    mean=(raw[1][0]/raw[0][1],raw[1][1]/raw[0][0])
    second=(raw[2][0]/raw[0][1],raw[2][1]/raw[0][0])
    var=(second[0]-mean[1]**2,second[1]-mean[0]**2)
    cv=(var[0]/mean[1]**2,var[1]/mean[0]**2)
    return dict(density_interval=raw[0],mean_interval=mean,second_interval=second,
                variance_interval=var,cv2_interval=cv,raw_moment_intervals=raw,
                residuals=[r1,r2,r3],inverse_bound=1/delta,
                centres=dict(nu=centres[0],mean=centres[1]/centres[0],
                    second=centres[2]/centres[0],
                    cv2=centres[2]*centres[0]/centres[1]**2-1))


def weighted_choice(weights, rng):
    from math import lcm
    weights=[Fraction(x) for x in weights]
    if min(weights)<0 or sum(weights)<=0: raise ValueError('invalid exact weights')
    denominator=1
    for x in weights: denominator=lcm(denominator,x.denominator)
    ints=[x.numerator*(denominator//x.denominator) for x in weights]
    pick=rng.randrange(sum(ints))
    for i,v in enumerate(ints):
        if pick<v:return i
        pick-=v
    raise AssertionError('choice')


def source_paths(width,matching=False):
    for prev in range(1<<width):
        old,_=step(empty(width),prev,matching,False)
        old=old._replace(colours=(FORBIDDEN,)*len(old.winding))
        for mask in range(1<<width):
            cur,_=step(old,mask,matching,False)
            for k,c in enumerate(cur.colours):
                if c:continue
                col=tuple(TAG if j==k else FORBIDDEN if (c0 or j<k) else NEUTRAL
                          for j,c0 in enumerate(cur.colours))
                anchor=cur.labels.index(k)
                yield prev,mask,anchor,cur._replace(colours=col)


def sample_component(width,matching,p,seed=0,max_rows=10000):
    """Exact complete-component Palm sampler (small widths, rational p).

    The row cap raises, never returns a censored sample. Exact weights are
    used; the final path telescoping identity is checked as a Fraction.
    """
    import random
    p=Fraction(p); rng=random.Random(seed)
    states,tr,src=build(width,matching); red,ss,blocks=lump(tr,src)
    a,R,b,delta=numeric_system(red,ss,p,True); n=len(b)
    A=[[int(i==j)-R[i][j] for j in range(n)] for i in range(n)]
    h=solve(A,b); nu=sum(x*y for x,y in zip(a,h)); index={s:i for i,s in enumerate(states)}
    wt=[p**m.bit_count()*(1-p)**(width-m.bit_count()) for m in range(1<<width)]
    paths=list(source_paths(width,matching))
    weights=[wt[prev]*wt[mask]*h[blocks[index[s]]] for prev,mask,anchor,s in paths]
    assert sum(weights)==nu
    choice=weighted_choice(weights,rng);prev,mask,anchor,s=paths[choice]
    qpath=weights[choice]/nu; prior=wt[prev]*wt[mask]; rows=[prev,mask]
    for _ in range(max_rows):
        i=index[s]; hi=h[blocks[i]]; choices=[]
        for m,(j,out) in enumerate(tr[i]):
            hm=h[blocks[j]] if j>=0 else int(j==-1)
            choices.append(wt[m]*hm/hi)
        assert sum(choices)==1
        m=weighted_choice(choices,rng); qpath*=choices[m];prior*=wt[m];rows.append(m)
        nxt,out=step(s,m,matching)
        if out:
            assert out==1 and qpath==prior/nu
            return dict(rows=rows,anchor=anchor,span=len(rows)-2,
                        path_probability=qpath,unconditioned_row_weight=prior,nu=nu)
        s=nxt
    raise RuntimeError('sampling row cap reached; no sample returned')


def expand_mask(mask,width,matching=False):
    if not matching:return mask
    full=(1<<width)-1
    return mask | ((mask<<1)&full) | (mask>>(width-1)) | (mask>>1) | ((mask&1)<<(width-1))


def activity_transfer(width,matching=False,state_cap=100000):
    """Direct complete-component activity; track two rows plus connectivity.

    A transition label (k,b) means u**k v**b for occupied and DISTINCT
    boundary sites finalized in the current row. Connected pieces may not
    retire before the whole selected component. No random environment needed.
    """
    if not 2<=width<=6:raise ValueError('reference activity builder width 2..6')
    states=[];idx={};src=Counter()
    def add(s):
        if s not in idx:
            if len(states)>=state_cap:raise RuntimeError('activity state cap')
            idx[s]=len(states);states.append(s)
        return idx[s]
    for mask in range(1,1<<width):
        first,_=step(empty(width),mask,matching,False)
        src[(add((0,first)),expand_mask(mask,width,matching).bit_count())]+=1
    tr=[]
    full=(1<<width)-1
    for old,s in states:
        current=sum(1<<i for i,k in enumerate(s.labels) if k>=0)
        horizontal=expand_mask(current,width,True)
        k=current.bit_count();row=[]
        for mask in range(1<<width):
            nxt,out=step(s,mask,matching,False,True)
            boundary=(horizontal|expand_mask(old,width,matching)|expand_mask(mask,width,matching))&(~current)&full
            b=boundary.bit_count()
            if out==1:row.append((-1,k,b+expand_mask(current,width,matching).bit_count()))
            elif out==0:row.append((add((current,nxt)),k,b))
        tr.append(row)
    blocks=[0]*len(tr)
    while True:
        lookup={};nxt=[]
        for row in tr:
            sig=tuple(sorted(Counter((blocks[j] if j>=0 else j,k,b) for j,k,b in row).items()))
            if sig not in lookup:lookup[sig]=len(lookup)
            nxt.append(lookup[sig])
        if nxt==blocks:break
        blocks=nxt
    reduced=[[(blocks[j] if j>=0 else j,k,b) for j,k,b in tr[blocks.index(i)]] for i in range(max(blocks)+1)]
    source=Counter()
    for (i,b),c in src.items():source[(blocks[i],b)]+=c
    return states,reduced,source


def activity_joint_moments(tr,src,p):
    """Exact all-height complete-Palm moments of span L and occupation K."""
    p=Fraction(p);n=len(tr)
    R=[[[Fraction(0) for _ in range(n)] for _ in range(n)] for order in range(3)]
    b=[[Fraction(0)]*n for _ in range(3)];alpha=[Fraction(0)]*n
    for i,row in enumerate(tr):
        for j,k,nb in row:
            weight=p**k*(1-p)**nb
            factors=(1,k,k*(k-1))
            for q,f in enumerate(factors):
                if j<0:b[q][i]+=f*weight
                else:R[q][i][j]+=f*weight
    for (i,nb),c in src.items():alpha[i]+=c*(1-p)**nb
    A=[[int(i==j)-R[0][i][j] for j in range(n)] for i in range(n)]
    mv=lambda mat,x:[sum(v*y for v,y in zip(row,x)) for row in mat]
    add=lambda *xs:[sum(z) for z in zip(*xs)]
    dot=lambda x:sum(a*b for a,b in zip(alpha,x))
    y=solve(A,b[0]);yL=solve(A,y);z3=solve(A,yL)
    yK=solve(A,add(b[1],mv(R[1],y)))
    yKK=solve(A,add(b[2],mv(R[2],y),[2*x for x in mv(R[1],yK)]))
    yLK=solve(A,add(mv(R[1],yL),yK))
    nu=dot(y);EL=dot(yL)/nu;EK=dot(yK)/nu
    varL=dot([2*x-y for x,y in zip(z3,yL)])/nu-EL**2
    varK=dot(add(yKK,yK))/nu-EK**2;cov=dot(yLK)/nu-EL*EK
    assert varL>=0 and varK>=0 and cov**2<=varL*varK
    return dict(nu=nu,mean_span=EL,mean_occupation=EK,variance_span=varL,
                variance_occupation=varK,covariance_span_occupation=cov,
                correlation_squared=cov**2/(varL*varK) if varL*varK else Fraction(0))


if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--width',type=int,default=4)
    ap.add_argument('--matching',action='store_true'); ap.add_argument('--p',default='1/4')
    ap.add_argument('--float',action='store_true'); ap.add_argument('--output',type=Path)
    args=ap.parse_args(); start=perf_counter()
    states,tr,src=build(args.width,args.matching); tr,src,blocks=lump(tr,src)
    out={'width':args.width,'matching':args.matching,'p':args.p,'states':len(states),
         'blocks':len(tr),'all_height_moments':moments(tr,src,Fraction(args.p),exact=not args.float),
         'seconds':perf_counter()-start,'mode':'floating diagnostic' if args.float else 'exact rational'}
    text=json.dumps(jsonable(out),indent=2)
    if args.output:
        if args.output.exists():raise FileExistsError(args.output)
        args.output.parent.mkdir(parents=True,exist_ok=True); args.output.write_text(text+'\n')
    else:print(text)

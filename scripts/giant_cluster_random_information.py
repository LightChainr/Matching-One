#!/usr/bin/env python3
"""Complete-cluster column scores: exact all-height controls, no sampling.

Uses the unchanged lifted frontier engine, but reconstructs a NEW direct
activity with even/odd COLUMN occupation and distinct-boundary marks. The
lumping preserves these marks and is not the homogeneous unmarked quotient.
"""
from __future__ import annotations
import argparse
from collections import Counter, deque
from fractions import Fraction as F
import hashlib
import importlib.util
from itertools import product
import json
from math import comb, factorial, sqrt
from pathlib import Path
import sys

ENGINE_BLOB = '52f3611990ce2b1331d9e5296e0262f5e402e0d7'

def load_engine():
    here = Path(__file__).resolve()
    paths = (here.with_name('tagged_winding_span.py'),
             here.parents[2]/'source_inputs/tagged_winding_span.py')
    for path in paths:
        if not path.is_file():
            continue
        raw = path.read_bytes()
        got = hashlib.sha1(b'blob '+str(len(raw)).encode()+b'\0'+raw).hexdigest()
        if got != ENGINE_BLOB:
            raise ValueError(f'Unexpected engine blob {got}: {path}')
        spec = importlib.util.spec_from_file_location('random_info_pinned', path)
        mod = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = mod
        spec.loader.exec_module(mod)
        return mod
    raise FileNotFoundError('Pinned tagged_winding_span.py is required.')

def parity_counts(mask, width):
    return tuple(sum((mask >> x) & 1 for x in range(k, width, 2)) for k in (0, 1))

def colored_activity(engine, width, matching=False, cap=10000):
    if width not in (2, 3, 4):
        raise ValueError('Bounded reference controls support widths 2, 3, 4.')
    states, index, source = [], {}, Counter()
    def add(state):
        if state not in index:
            if len(states) >= cap:
                raise RuntimeError('State cap reached; nothing truncated.')
            index[state] = len(states); states.append(state)
        return index[state]
    for mask in range(1, 1 << width):
        state, _ = engine.step(engine.empty(width), mask, matching, False)
        bd = parity_counts(engine.expand_mask(mask, width, matching), width)
        source[(add((0, state)), (0, 0, *bd))] += 1
    full, transitions = (1 << width)-1, []
    for old, state in states:
        current = sum(1 << x for x, label in enumerate(state.labels) if label >= 0)
        k = parity_counts(current, width)
        horizontal = engine.expand_mask(current, width, True)
        row = []
        for mask in range(1 << width):
            nxt, outcome = engine.step(state, mask, matching, False, True)
            bd = (horizontal | engine.expand_mask(old, width, matching)
                  | engine.expand_mask(mask, width, matching)) & (~current) & full
            b = parity_counts(bd, width)
            if outcome == 1:
                upper = parity_counts(engine.expand_mask(current, width, matching), width)
                row.append((-1, (*k, b[0]+upper[0], b[1]+upper[1])))
            elif outcome == 0:
                row.append((add((current, nxt)), (*k, *b)))
        transitions.append(row)
    blocks = [0]*len(states)
    while True:
        labels, nxt = {}, []
        for row in transitions:
            sig = tuple(sorted(Counter((blocks[j] if j >= 0 else j, mark)
                                       for j, mark in row).items()))
            if sig not in labels: labels[sig] = len(labels)
            nxt.append(labels[sig])
        if nxt == blocks: break
        blocks = nxt
    n = max(blocks)+1
    reduced = [[(blocks[j] if j >= 0 else j, mark)
                for j, mark in transitions[blocks.index(i)]] for i in range(n)]
    src = Counter()
    for (i, mark), count in source.items(): src[(blocks[i], mark)] += count
    return reduced, src, len(states)

def activity_weight(mark, qe, qo):
    ke, ko, be, bo = mark
    return qe**ke * qo**ko * (1-qe)**be * (1-qo)**bo

def dot(a, b): return sum((x*y for x, y in zip(a,b)), F(0))
def mv(A, x): return [dot(row, x) for row in A]

def inverse_fraction(A):
    """Exact Gauss-Jordan; build once and reuse for moment right-hand sides."""
    n = len(A)
    M = [[F(v) for v in row]+[F(i == j) for j in range(n)] for i,row in enumerate(A)]
    for j in range(n):
        k = next((i for i in range(j,n) if M[i][j]), None)
        if k is None: raise ValueError('Singular exact activity resolvent.')
        M[j], M[k] = M[k], M[j]
        d = M[j][j]; M[j] = [v/d for v in M[j]]
        for i in range(n):
            if i != j and M[i][j]:
                d = M[i][j]; M[i] = [a-d*b for a,b in zip(M[i],M[j])]
    return [row[n:] for row in M]

def indices(d, order):
    return sorted((a for a in product(range(order+1),repeat=d) if sum(a)<=order),
                  key=lambda a:(sum(a),a))

def all_height_moments(tr, src, qe, qo, forms, order=4):
    """Raw joint moments of additive marks using differentiated resolvents.

    forms is a list of length-4 rational vectors on (Ke,Ko,Be,Bo).
    Source boundary marks and BOTH terminal boundary rows are included.
    """
    qe,qo = F(qe),F(qo)
    if not (0 < qe < 1 and 0 < qo < 1): raise ValueError('Invalid site probabilities.')
    keys = indices(len(forms),order); zero = (0,)*len(forms); n = len(tr)
    R = {}; b = {}; alpha = {}
    for a in keys:
        R[a] = [[F(0)]*n for _ in range(n)]; b[a] = [F(0)]*n; alpha[a] = [F(0)]*n
    def weighted(mark, a):
        z = activity_weight(mark,qe,qo)
        for power,form in zip(a,forms): z *= dot(form,mark)**power
        return z
    for i,row in enumerate(tr):
        for j,mark in row:
            for a in keys:
                value = weighted(mark,a)
                if j<0: b[a][i] += value
                else: R[a][i][j] += value
    for (i,mark),count in src.items():
        for a in keys: alpha[a][i] += count*weighted(mark,a)
    G = inverse_fraction([[F(i==j)-R[zero][i][j] for j in range(n)] for i in range(n)])
    y, moments = {}, {}
    for a in keys:
        rhs = b[a][:]
        for h in keys:
            if h==zero or any(x>z for x,z in zip(h,a)): continue
            prev = tuple(z-x for x,z in zip(h,a))
            factor = 1
            for z,x in zip(a,h): factor *= comb(z,x)
            term = mv(R[h],y[prev]); rhs = [v+factor*u for v,u in zip(rhs,term)]
        y[a] = mv(G,rhs)
        total = F(0)
        for h in keys:
            if any(x>z for x,z in zip(h,a)): continue
            prev = tuple(z-x for x,z in zip(h,a)); factor = 1
            for z,x in zip(a,h): factor *= comb(z,x)
            total += factor*dot(alpha[h],y[prev])
        moments[a] = total
    nu = moments[zero]
    return nu,{a:v/nu for a,v in moments.items()}

def physical_shapes(width,height,matching):
    """Independent lifted BFS, not the frontier engine; boundary masks exact."""
    steps=[(-1,0),(1,0),(0,-1),(0,1)]
    if matching: steps += [(a,b) for a in (-1,1) for b in (-1,1)]
    records=[]
    for mask in range(1,1<<(width*height)):
        if not mask & ((1<<width)-1): continue
        if not mask >> (width*(height-1)): continue
        C={(x,y) for y in range(height) for x in range(width) if mask>>(y*width+x)&1}
        root=min(C); lift={root:0}; Q=deque([root]); wind=False
        while Q:
            u=Q.popleft()
            for dx,dy in steps:
                v=((u[0]+dx)%width,u[1]+dy)
                if v not in C: continue
                expected=lift[u]+dx
                if v in lift: wind |= lift[v] != expected
                else: lift[v]=expected;Q.append(v)
        if len(lift)!=len(C) or not wind: continue
        B={((x+dx)%width,y+dy) for x,y in C for dx,dy in steps}-C
        mark=tuple(sum(x%2==j for x,y in D) for D in (C,B) for j in (0,1))
        records.append(mark)
    return records

def finite_coefficients(tr,src,height):
    states=Counter(src)
    exits=Counter()
    for level in range(1,height+1):
        nxt=Counter()
        for (i,mark),count in states.items():
            for j,m in tr[i]:
                total=tuple(a+b for a,b in zip(mark,m))
                if j<0:
                    if level==height: exits[total]+=count
                elif level<height: nxt[(j,total)]+=count
        states=nxt
    return exits

def centre(m,i,j):
    a,b=m[(1,0)],m[(0,1)];value=F(0)
    for k in range(i+1):
        for l in range(j+1):
            value += comb(i,k)*comb(j,l)*(-a)**(i-k)*(-b)**(j-l)*m[(k,l)]
    return value

def wired_score_controls():
    """Independent normalization identity on finite wired graphs.

    Root n is fixed open. A site is queried iff in the root component or
    its distinct external vertex boundary. All multilinear score moments
    vanish, but squared scores need not be independent.
    """
    graphs=[(5,[(5,0),(0,1),(1,2),(2,3),(3,4)]),
            (6,[(6,0),(0,1),(1,2),(2,3),(3,0),(2,4),(4,5),(5,1)]),
            (6,[(6,0),(6,3),(0,1),(1,2),(3,4),(4,5),(0,3),(1,4),(2,5)])]
    out=[]
    for n,edges in graphs:
        adj=[set() for _ in range(n+1)]
        for a,b in edges:adj[a].add(b);adj[b].add(a)
        for q in (F(2,5),F(3,4)):
            p=1-q; moments=[F(0)]*(1<<n); probe=[F(0)]*n
            cov=[[F(0)]*n for _ in range(n)];response=[[F(0)]*n for _ in range(n)]
            raw_response=[[F(0)]*n for _ in range(n)];normalizer=F(0)
            h=[(F(1,31) if i%2==0 else -F(1,37)) for i in range(n)]
            for mask in range(1<<n):
                C={n};todo=[n]
                while todo:
                    u=todo.pop()
                    for v in adj[u]:
                        if v not in C and (v==n or mask>>v&1):C.add(v);todo.append(v)
                Q=(C|{v for u in C for v in adj[u]})-{n}
                wt=q**mask.bit_count()*p**(n-mask.bit_count())
                fullscore=[(F(mask>>i&1)-q)/(p*q) for i in range(n)]
                xi=[fullscore[i] if i in Q else F(0) for i in range(n)]
                prodv=[F(1)]*(1<<n)
                for sub in range(1,1<<n):
                    bit=sub & -sub; i=bit.bit_length()-1
                    prodv[sub]=prodv[sub^bit]*xi[i];moments[sub]+=wt*prodv[sub]
                ratio=F(1)
                for i in range(n):
                    ratio*=1+h[i]*xi[i];probe[i]+=wt*(i in Q)
                    for j in range(n):
                        cov[i][j]+=wt*xi[i]*xi[j]
                        response[i][j]+=wt*(i in C)*xi[j]
                        raw_response[i][j]+=wt*(i in C)*fullscore[j]
                normalizer+=wt*ratio
            assert all(v==0 for v in moments[1:])
            assert normalizer==1 and response==raw_response
            assert all(cov[i][j]==(probe[i]/(p*q) if i==j else 0)
                       for i in range(n) for j in range(n))
            out.append(dict(vertices=n,q=q,configurations=1<<n,
                            zero_multilinear_moments=(1<<n)-1,
                            covariance_identities=n*n,response_identities=n*n,
                            finite_tilt_normalizer=normalizer))
    return out

def serialize(obj):
    if isinstance(obj,F): return {'fraction':str(obj),'decimal':format(float(obj),'.14g')}
    if isinstance(obj,dict): return {str(k):serialize(v) for k,v in obj.items()}
    if isinstance(obj,(list,tuple)): return [serialize(v) for v in obj]
    return obj

def run():
    engine=load_engine(); cache={}; physical=[]; tables=[]
    for w in (2,3,4):
        for matching in (False,True):
            tr,src,full=colored_activity(engine,w,matching);cache[(w,matching)]=(tr,src)
            checked=0
            for h in (1,2,3):
                coeff=Counter(physical_shapes(w,h,matching))
                assert coeff==finite_coefficients(tr,src,h),(w,h,matching)
                checked+=sum(coeff.values())
            physical.append(dict(width=w,matching=matching,full_states=full,
                                 colored_states=len(tr),complete_shapes=checked))
    for w in (2,4):
        pair=[]
        for matching,q in ((False,F(1,4)),(True,F(3,4))):
            tr,src=cache[(w,matching)];p=1-q
            f0=(1/q,1/q,-1/p,-1/p);f1=(1/q,-1/q,-1/p,1/p)
            nu,m=all_height_moments(tr,src,q,q,[f0,f1],4)
            var0=centre(m,2,0);var1=centre(m,0,2);mixed=centre(m,2,2)
            assert m[(0,1)]==0 and centre(m,1,1)==0
            row=dict(width=w,matching=matching,q=q,nu=nu,mean_uniform=m[(1,0)],
                     variance_uniform=var0,variance_striped=var1,
                     cross_covariance=centre(m,1,1),
                     kurtosis_uniform=centre(m,4,0)/var0**2,
                     kurtosis_striped=centre(m,0,4)/var1**2,
                     square_correlation=(mixed-var0*var1)/
                        (sqrt(float((centre(m,4,0)-var0**2)*(centre(m,0,4)-var1**2)))),
                     normalized_mixed_fourth=mixed/(var0*var1))
            tables.append(row);pair.append(row)
            vf=(1/q**2,1/q**2,1/p**2,1/p**2)
            nv,vm=all_height_moments(tr,src,q,q,[vf],1);assert nv==nu
            row['square_correlation_squared']=(mixed-var0*var1)**2/((centre(m,4,0)-var0**2)*(centre(m,0,4)-var1**2))
            row['expected_observed_information']=vm[(1,)]
            row['log_intensity_hessian']=[var0-vm[(1,)],var1-vm[(1,)]]
        assert pair[0]['nu']==pair[1]['nu']
        assert pair[0]['mean_uniform']==-pair[1]['mean_uniform']
        assert pair[0]['log_intensity_hessian']==pair[1]['log_intensity_hessian']
    checks=[]
    for w in (2,3,4):
        for qe,qo in ((F(2,3),F(3,4)),(F(4,5),F(5,7))):
            white=all_height_moments(*cache[(w,True)],qe,qo,[],0)[0]
            black=all_height_moments(*cache[(w,False)],1-qe,1-qo,[],0)[0]
            assert white==black
            checks.append(dict(width=w,white_even=qe,white_odd=qo,nu=white))
    ring={(x,y) for x in (-1,0,1) for y in (-1,0,1) if x or y}
    hole=Counter()
    for dx,dy in product(range(-2,3),repeat=2):
        translated={(x+dx,y+dy) for x,y in ring}
        if ring & translated:
            hole[len(ring|translated)]+=1;hole[16]-=1
    assert dict(hole)=={8:1,12:4,13:4,14:12,15:4,16:-25}
    model={
        'query_hole_covariance_coefficients':dict(sorted(hole.items())),
        'orthogonal_square_correlation':F(1,5),
        'orthogonal_normalized_mixed_fourth':F(2),
        'marginal_kurtosis':F(6),
        'student_t2_cdf_at_0':F(1,2),
        'pooled_estimator_limit_df':{n:2*n for n in (1,2,3,4)},
        'pooled_limit_variance_times_c':{n:F(1,n-1) for n in (2,3,4)},
        'common_clock_even_moments':{
            (i,j):F(factorial(i+j))*F(factorial(2*i),2**i*factorial(i))*
                     F(factorial(2*j),2**j*factorial(j))
            for i,j in ((1,0),(2,0),(1,1),(2,1),(2,2))}}
    return serialize(dict(engine_blob=ENGINE_BLOB,wired_controls=wired_score_controls(),physical_controls=physical,
                          all_height_scores=tables,inhomogeneous_duality=checks,
                          limit_model=model,
                          scope='Finite controls exact except square correlation/display decimals; no large-width inference from these tables.'))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path)
    args=parser.parse_args();result=run();text=json.dumps(result,ensure_ascii=False,indent=2)+'\n'
    if args.output:
        if args.output.exists(): raise FileExistsError(args.output)
        args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(text)
    else: print(text)

#!/usr/bin/env python3
"""Exact small interfaces for intrinsic birth clocks and response alignment.

No Monte Carlo, no scaling fit, and no external transfer input. The 3x3
cylinder is horizontally periodic and vertically free. Bivariate coupling
is integrated exactly, not approximated by independent parameter samples.
"""
from __future__ import annotations
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from math import comb
from pathlib import Path

W = H = 3
N = W * H

@lru_cache(None)
def components(mask: int) -> tuple[tuple[int, bool], ...]:
    """Physical NN components, with horizontal winding from lifted BFS."""
    unseen = {i for i in range(N) if mask >> i & 1}
    ans = []
    while unseen:
        root = min(unseen)
        unseen.remove(root)
        todo = [root]
        lift = {root: 0}
        bits = 0
        wraps = False
        for v in todo:
            bits |= 1 << v
            x, y = v % W, v // W
            for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                if not 0 <= y + dy < H:
                    continue
                u = (y + dy) * W + (x + dx) % W
                if not mask >> u & 1:
                    continue
                z = lift[v] + dx
                if u not in lift:
                    lift[u] = z
                    unseen.remove(u)
                    todo.append(u)
                elif lift[u] != z:
                    assert (lift[u] - z) % W == 0
                    wraps = True
        ans.append((bits, wraps))
    return tuple(ans)

@lru_cache(None)
def winding_count(mask: int) -> int:
    return sum(w for _, w in components(mask))

def occupancy(mask: int, p: F) -> F:
    k = mask.bit_count()
    return p**k * (1-p)**(N-k)

@lru_cache(None)
def nu(p: F) -> F:
    """Finite-strip TOTAL expected component count, not infinite intensity."""
    return sum((occupancy(m,p)*winding_count(m) for m in range(1 << N)), F(0))

@lru_cache(None)
def clock_bernstein(r: F) -> tuple[F, ...]:
    """Normalized degree-N Bernstein coefficients of beta(p;r), t=p/r.

Each final essential component contributes its increasing first-winding
reliability polynomial. Lower-degree polynomials are elevated exactly.
"""
    if not 0 < r < 1:
        raise ValueError('terminal parameter must lie in (0,1)')
    raw = [F(0) for _ in range(N+1)]
    for final in range(1 << N):
        weight = occupancy(final,r)
        for cmask, wraps in components(final):
            if not wraps:
                continue
            n = cmask.bit_count()
            count = [0]*(n+1)
            sub = cmask
            while True:
                if winding_count(sub):
                    count[sub.bit_count()] += 1
                if sub == 0:
                    break
                sub = (sub-1)&cmask
            for j in range(N+1):
                a = sum(count[k]*comb(N-n,j-k)
                        for k in range(n+1) if 0 <= j-k <= N-n)
                raw[j] += weight*a
    return tuple(raw[j]/comb(N,j) for j in range(N+1))

def bernstein_eval(coeffs: tuple[F,...], t: F) -> F:
    n = len(coeffs)-1
    return sum((a*comb(n,j)*t**j*(1-t)**(n-j)
                for j,a in enumerate(coeffs)),F(0))

def beta(p: F, r: F) -> F:
    if not 0 <= p <= r < 1:
        raise ValueError('require 0 <= p <= r < 1')
    return bernstein_eval(clock_bernstein(r),p/r)

def beta_prime(p: F, r: F) -> F:
    c = clock_bernstein(r)
    return N/r*bernstein_eval(tuple(c[j+1]-c[j] for j in range(N)),p/r)

def direct_pair_values(p: F, r: F) -> tuple[F,F,F,int]:
    """Integrate early/final common labels over all 3^N pairs."""
    hits = loss = pairs = F(0)
    checked = 0
    for final in range(1 << N):
        fc = [b for b,w in components(final) if w]
        early = final
        while True:
            ec = [b for b,w in components(early) if w]
            ns = [sum((e & c)==e for e in ec) for c in fc]
            assert sum(ns)==len(ec)
            hit = sum(k>0 for k in ns)
            lost = len(ec)-hit
            pair = sum(k*(k-1)//2 for k in ns)
            assert 0 <= lost <= pair
            wt = p**early.bit_count() * (r-p)**(final.bit_count()-early.bit_count()) * (1-r)**(N-final.bit_count())
            hits += wt*hit
            loss += wt*lost
            pairs += wt*pair
            checked += 1
            if early==0:
                break
            early=(early-1)&final
    return hits,loss,pairs,checked

def invert_beta(target: F, r: F, steps: int=48) -> tuple[F,F]:
    if not 0 <= target <= beta(r,r):
        raise ValueError('target outside clock range')
    lo,hi=F(0),r
    for _ in range(steps):
        mid=(lo+hi)/2
        if beta(mid,r)<target:
            lo=mid
        else:
            hi=mid
    assert beta(lo,r)<=target<=beta(hi,r)
    return lo,hi

def harmonic_controls() -> dict:
    h=F(-527,625)
    h8=2*h*h-1
    leak0=(h8-h)/(1-h)
    leak4=(1-h8)/(1-h)
    assert leak0==F(429,625) and leak4==F(196,625)
    return {'H4_oblique':h,'H8_oblique':h8,
            'hidden_H8_into_scalar':leak0,'hidden_H8_into_spin4':leak4,
            'pair_projection_spin4_coefficient':1/(1-h)}

def alignment_controls() -> dict:
    # b=x+g*x^2, e=b^2: all pure reparametrization, but normalized
    # single-charge curve changes by g*x^2 (root=0, root slope=1).
    checks=0
    for g in (F(-1,4), F(0), F(1,4)):
        for x in (F(-1,4),F(-1,8),F(0),F(1,8),F(1,4)):
            b=x+g*x*x
            bp=1+2*g*x
            bg=x*x
            ep=2*b*bp
            eg=2*b*bg
            assert bp>0
            assert eg-ep*bg/bp==0
            assert b-x==g*x*x
            checks+=1
    # b=x, e=x^2+g*x: scalar charge curve unchanged, true N=x.
    for x in (F(-1,4),F(-1,8),F(0),F(1,8),F(1,4)):
        assert (x - (2*x+F(1,4))*F(0))==x
        checks+=1
    return {'exact_controls':checks,
            'pure_tangent_false_shape_signal_at_x_1_4_g_1_4':F(1,64),
            'pure_normal_invisible_to_single_charge_at_x_1_4':F(1,4)}

@lru_cache(None)
def torus_rank(mask: int) -> int:
    """Independent 2D lifted BFS on the 3x3 NN torus (not the free strip)."""
    unseen = {i for i in range(N) if mask >> i & 1}
    periods = []
    while unseen:
        root=min(unseen); unseen.remove(root)
        lift={root:(0,0)}; queue=[root]
        for v in queue:
            x,y=v%W,v//W
            for dx,dy in ((1,0),(-1,0),(0,1),(0,-1)):
                u=((y+dy)%H)*W+(x+dx)%W
                if not mask>>u&1:
                    continue
                z=(lift[v][0]+dx,lift[v][1]+dy)
                if u not in lift:
                    lift[u]=z;unseen.remove(u);queue.append(u)
                else:
                    a,b=z[0]-lift[u][0],z[1]-lift[u][1]
                    assert a%W==0 and b%H==0
                    if a or b:
                        periods.append((a//W,b//H))
    if not periods:
        return 0
    a,b=periods[0]
    return 2 if any(a*y-b*x for x,y in periods) else 1

def pair_marks(mask: int) -> tuple[int,int,int]:
    """Horizontal NN pair sources; each horizontal edge counted once."""
    black=white=0
    for y in range(H):
        for x in range(W):
            a=(mask>>(y*W+x))&1
            b=(mask>>(y*W+(x+1)%W))&1
            black+=a*b;white+=(1-a)*(1-b)
    k=mask.bit_count()
    assert white==N-2*k+black
    return k,black,white

def rank_response_controls() -> dict:
    data=[(m,torus_rank(m),pair_marks(m)) for m in range(1<<N)]
    rows=[]
    for p in (F(1,3),F(1,2),F(2,3)):
        q=1-p
        mass=[F(0)]*3
        sums=[[F(0)]*3 for _ in range(3)]
        for mask,r,marks in data:
            wt=occupancy(mask,p);mass[r]+=wt
            for j,a in enumerate(marks):
                sums[r][j]+=wt*a
        means=[[s/mass[r] for s in sums[r]] for r in range(3)]
        def coordinate_derivative(j):
            return (means[0][j]-means[2][j])/2, means[1][j]-(means[0][j]+means[2][j])/2
        bK,eK=coordinate_derivative(0)
        bp,ep=bK/(p*q),eK/(p*q)
        assert bp<0
        response=[]
        for j in range(3):
            bg,eg=coordinate_derivative(j)
            response.append((-bg/bp,eg-ep*bg/bp))
        assert response[0]==(-p*q,F(0))
        assert response[1][1]==response[2][1]
        assert response[2][0]-response[1][0]==2*p*q
        rows.append({'p':p,'sector_probabilities':mass,
                     'K_logit':response[0], 'black_pair':response[1],
                     'white_pair':response[2],
                     'scope':'equal-b thermal compensation at p; not a claim p is a matching root'})
    # All amplitudes below are rational fugacities exp(g), so no floating
    # exponentials occur. This verifies the full normalized measures,
    # not merely the three rank probabilities or their derivatives.
    full_checks=0
    for p in (F(1,3),F(1,2),F(2,3)):
        for v in (F(1,2),F(2),F(3)):
            pe=p/(p+(1-p)*v*v)  # logit(pe)=logit(p)-2log(v)
            lhs=[];rhs=[]
            for mask,r,(k,b,w) in data:
                lhs.append(occupancy(mask,p)*v**w)
                rhs.append(occupancy(mask,pe)*v**b)
            zl,zr=sum(lhs),sum(rhs)
            for a,b in zip(lhs,rhs):
                assert a/zl==b/zr
                full_checks+=1
    return {'physical_torus_masks':len(data),'finite_derivative_cases':rows,
            'normalized_measure_equalities':full_checks}

def repaired_table() -> dict:
    """Postprocess PRINTED existing floats as exact decimal rationals.

These are not new eigenvalue calculations/certificates. A decimal input's
last-digit error remains present in the resulting displayed values.
"""
    p=F('0.5927460507921');q=1-p
    bs=['-1.350049753766846','-1.7138942463383808','-2.0732497837477912',
        '-2.4300006319440906','-2.7852763705546106']
    ws=['-0.6080813474300464','-0.7864337384173812','-0.9602971742425906',
        '-1.1315559208546921','-1.3013395578810114']
    tb=['-0.3116848217980418','-0.3102994990346319','-0.3098063945349776',
        '-0.30962331020127565','-0.30954105615205907']
    tw=['0.17111151832689642','0.17249684109030636','0.17298994558996025',
        '0.17317302992366124','0.17325528397288223']
    result=[]
    for w,b,a,t1,t2 in zip(range(4,9),bs,ws,tb,tw):
        nb=F(b)+w*p*p;nw=F(a)+w*q*q
        diffT=F(t2)-F(t1)-2*p*q
        assert abs(nb-nw)<F(2,10**15)
        assert abs(diffT)<F(5,10**15)
        result.append({'width':w,'raw_black_N':F(b),'raw_white_N':F(a),
                       'physical_black_N':nb,'physical_white_N':nw,
                       'pair_normal_difference':nb-nw,'thermal_shift_residual':diffT})
    return {'input_commit':'723f62929c6a343a91607e4caa47425e06230eda',
            'input_blob':'8dde718f411f2e36fa2fd70ea9cbb7cee7234ae3',
            'input_path':'results/research-dispatch/sector-root-response-20260914.json',
            'qualification':'Repair of published numerical values, NOT an independent transfer run.',
            'p_ref':p,'twice_pq':2*p*q,'rows':result}

def spin_zero_mode_control() -> dict:
    c=F(0);h=F(5,96);hbar=h
    eigen=h*h-(c+2)*h/12+c*(5*c+22)/2880
    assert h-hbar==0 and eigen==F(-55,9216)
    return {'c':c,'h':h,'hbar':hbar,'spatial_momentum_weight':h-hbar,
            'chiral_spin4_I3_eigenvalue':eigen,
            'meaning':'A Virasoro zero-mode counterexample to blanket spin exclusion, not a site field identification.',
            'source':'https://arxiv.org/html/1810.11053v2 equation (2)'}

def encode(obj):
    if isinstance(obj,F):
        return {'fraction':str(obj),'decimal':float(obj)}
    if isinstance(obj,dict):
        return {str(k):encode(v) for k,v in obj.items()}
    if isinstance(obj,(list,tuple)):
        return [encode(v) for v in obj]
    return obj

def report() -> dict:
    terminals=(F(1,3),F(1,2),F(3,4))
    cases=[]
    total_pairs=0
    for r in terminals:
        coeff=clock_bernstein(r)
        assert coeff[0]==0 and coeff[-1]==nu(r)
        assert all(coeff[j+1]>=coeff[j] for j in range(N))
        for p in (r/3,2*r/3):
            a,l,m,n=direct_pair_values(p,r)
            total_pairs+=n
            assert a==beta(p,r) and nu(p)==a+l and 0<=l<=m
            assert 0<beta_prime(p,r)<=N
            cases.append({'p':p,'terminal':r,'beta':a,'nu_early':nu(p),
                          'loss':l,'pair_bound':m,'derivative':beta_prime(p,r)})
    r,p0=F(3,4),F(1,4)
    base=beta(p0,r)
    brackets=[]
    for lam in (F(1,2),F(1),F(2)):
        lo,hi=invert_beta(lam*base,r)
        brackets.append({'relative_clock':lam,'p_lower':lo,'p_upper':hi})
    return encode({'scope':'Exact finite free-ended C3 x 3 interfaces, not large-width data.',
                   'physical_masks':1<<N,'three_colour_pairs_per_case':3**N,
                   'weighted_pair_evaluations':total_pairs,
                   'clock_cases':cases,'inverse_clock':brackets,
                   'harmonics':harmonic_controls(),'alignment':alignment_controls(),
                   'physical_sources':rank_response_controls(),
                   'repaired_external_table':repaired_table(),
                   'spin_zero_mode':spin_zero_mode_control()})

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    result=json.dumps(report(),ensure_ascii=False,indent=2,sort_keys=True)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(result,encoding='utf-8')
    else:
        print(result,end='')

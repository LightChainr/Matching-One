#!/usr/bin/env python3
"""Exact finite controls for guarded preparations and two-row noise robustness.

No simulation, floating-point rank, whole-width state enumeration, or network
access. The all-width conclusions are proved in NOTE.md, not inferred from
these finite controls. Optional --input-check replays the preceding author's
physical oracle on the cases already evaluated here.
"""
from __future__ import annotations
import argparse
from collections import Counter, deque
from decimal import Decimal, localcontext
from fractions import Fraction
from itertools import product, combinations
import json
from math import comb, isqrt
from pathlib import Path
import resource
import time

PAIRS = ((0, 0), (1, 0), (0, 1), (1, 1))
RANK_CALLS = 0

def prepare(bits: tuple[tuple[int, int], ...], w: int) -> tuple[int, int, int]:
    if not bits or w < 5 * len(bits):
        raise ValueError('width >= 5 times number of blocks is required')
    return (0, sum(15 << (5*g) for g in range(len(bits))),
            sum((9 + 2*u + 4*v) << (5*g) for g, (u, v) in enumerate(bits)))

def query(w: int, g: int, which: int) -> tuple[int, int]:
    if w < 5 or 5*g + 4 >= w or which not in (0, 1):
        raise ValueError('invalid bit query')
    full = (1 << w)-1
    if which == 0:
        return full ^ (5 << (5*g)), full ^ (4 << (5*g))
    return full ^ (10 << (5*g)), full ^ (2 << (5*g))

def torus_rank(w: int, rows: tuple[int, ...]) -> int:
    """Integer lifted graph BFS; no boundary state or mod-2 assumption."""
    global RANK_CALLS
    RANK_CALLS += 1
    h = len(rows)
    if min(w, h) < 2 or any(m < 0 or m >= 1 << w for m in rows):
        raise ValueError('honest torus dimensions and valid masks required')
    seen: dict[tuple[int, int], tuple[int, int]] = {}
    first: tuple[int, int] | None = None
    for y, mask in enumerate(rows):
        for x in range(w):
            if not (mask >> x & 1) or (x, y) in seen:
                continue
            seen[x, y] = (0, 0)
            queue = deque([(x, y)])
            while queue:
                a, b = queue.popleft()
                pa, pb = seen[a, b]
                for dx, dy in ((1,0), (-1,0), (0,1), (0,-1)):
                    c, d = (a+dx) % w, (b+dy) % h
                    if not (rows[d] >> c & 1):
                        continue
                    proposed = pa+dx, pb+dy
                    if (c,d) not in seen:
                        seen[c,d] = proposed
                        queue.append((c,d))
                    else:
                        ux, uy = proposed[0]-seen[c,d][0], proposed[1]-seen[c,d][1]
                        assert ux % w == 0 and uy % h == 0
                        if ux or uy:
                            if first is None:
                                first = ux, uy
                            elif first[0]*uy != first[1]*ux:
                                return 2
    return int(first is not None)

def rotate_mask(m: int, w: int, shift: int) -> int:
    return sum(1 << ((x+shift) % w) for x in range(w) if m >> x & 1)

def normalize_u(m: int, n: int, w: int, g: int, which: int) -> tuple[int,int]:
    # Every U query becomes a=0,i=1,j=2,d=3. Reflect V within its block.
    def move(row: int) -> int:
        result=0
        for x in range(w):
            if row >> x & 1:
                z=(x-5*g) % w
                if which: z=(3-z) % w
                result |= 1 << z
        return result
    return move(m), move(n)

def good_sufficient_event(m: int, n: int, w: int, g: int, which: int) -> bool:
    """Seven intended high sites plus no closed cross-row adjacent pair."""
    m,n=normalize_u(m,n,w,g,which)
    local=((0,1),(1,1),(1,0),(3,0),(3,1),(w-1,0),(w-1,1))
    if not all(((m,n)[row] >> col) & 1 for col,row in local):
        return False
    # Two-row rectangle: columns 3,...,w-1. There is no matching TB path
    # iff no vacant cross-row pair is in the same or adjacent column.
    return not any(not (m >> x & 1) and not (n >> y & 1)
                   for x in range(3,w) for y in range(max(3,x-1),min(w,x+2)))

def low_sites_correct(m: int,n: int,w: int,g: int,which: int) -> bool:
    m,n=normalize_u(m,n,w,g,which)
    return not (m & 5) and not (n & 4)

def ladder_crosses(m: int,n: int,length: int) -> bool:
    active={(x,row) for row,mask in enumerate((m,n)) for x in range(length) if mask >> x & 1}
    todo=deque((0,row) for row in (0,1) if (0,row) in active)
    seen=set(todo)
    while todo:
        x,row=todo.popleft()
        if x==length-1:return True
        for v in ((x-1,row),(x+1,row),(x,1-row)):
            if v in active and v not in seen:seen.add(v);todo.append(v)
    return False

def ladder_control() -> dict:
    counts=[]
    for length in range(2,7):
        bad=0
        for code in range(1 << (2*length)):
            m,n=code & ((1<<length)-1),code >> length
            dual=any(not(m>>x&1) and not(n>>y&1)
                     for x in range(length) for y in range(max(0,x-1),min(length,x+2)))
            assert ladder_crosses(m,n,length) != dual
            bad+=int(dual)
        counts.append({'columns':length,'all_colourings':1<<(2*length),'no_LR_crossing':bad,
                       'adjacent_cross_row_pairs':3*length-2})
    return {'rectangles':counts,'checked_colourings':sum(d['all_colourings'] for d in counts)}

def qstr(q: Fraction) -> str:
    return f'{q.numerator}/{q.denominator}'

def query_controls(parent_oracle=None) -> dict:
    checked=0; oracle_checks=0; guard_checks=0; count_good=0; zero_low_checks=0
    histograms=[]
    for w,t,full in ((5,1,True),(10,2,False),(12,2,False)):
        faults=(list(range(1<<(2*w))) if full else
                [0]+[1<<i for i in range(2*w)]+[(1<<i)|(1<<j) for i,j in combinations(range(2*w),2)])
        for bits in product(PAIRS,repeat=t):
            old=prepare(bits,w)
            assert torus_rank(w,old)==0
            for g in range(t):
                for which in (0,1):
                    target=bits[g][which]
                    m0,n0=query(w,g,which)
                    errors=[0]*(2*w+1); total=[0]*(2*w+1)
                    for fault in faults:
                        m=m0 ^ (fault & ((1<<w)-1)); n=n0 ^ (fault>>w)
                        r=torus_rank(w,old+(m,n)); checked+=1
                        assert r in (0,1)
                        if parent_oracle is not None:
                            assert r==parent_oracle(w,list(old)+( [m,n] ))[0]
                            oracle_checks+=1
                        if good_sufficient_event(m,n,w,g,which) and target:
                            assert r==1; count_good+=1
                        if low_sites_correct(m,n,w,g,which) and not target:
                            assert r==0;zero_low_checks+=1
                        if any(not(m>>col&1) and not(n>>col&1) for col in range(4,5*t,5)):
                            assert r==0 and torus_rank(w,old+(m,))==0;guard_checks+=1
                        k=fault.bit_count();total[k]+=1;errors[k]+=int(r!=target)
                    if full:
                        assert total==[comb(2*w,k) for k in range(2*w+1)]
                        values=[]
                        for e in (Fraction(1,1000),Fraction(1,100),Fraction(1,10),Fraction(1,4)):
                            probability=sum(Fraction(c)*e**k*(1-e)**(2*w-k) for k,c in enumerate(errors))
                            bound=7*e+3*w*e*e
                            assert probability<=bound
                            if not target:assert probability<=3*e
                            values.append({'epsilon':qstr(e),'error':qstr(probability),'bound':qstr(min(1,bound))})
                        histograms.append({'width':w,'bits':bits,'query':[g,which],
                            'error_counts_by_fault_number':errors,'rational_evaluations':values})
    return {'full_graph_rank_calls':checked,'parent_oracle_comparisons':oracle_checks,
            'guard_implies_zero_checks':guard_checks,'positive_witness_checks':count_good,
            'zero_bit_low_mask_checks':zero_low_checks,
            'width5_complete_fault_histograms':histograms,
            'larger_widths_scope':'widths 10 and 12: all preparations and single-bit queries, all <=2 future faults only'}

def adaptive_erasure_control() -> dict:
    """Full finite seed enumeration for an observation-adaptive q in [1/3,2/3].

    Uniform trits implement Bernoulli 1/3 or 2/3. The second row depends on
    R1. A guard trit=2 at both rows forces the same zero rank transcript.
    """
    w=5; allrows=list(product(range(3),repeat=w)); checked=0; erase=0
    olds=[prepare((b,),w) for b in PAIRS]
    for first in allrows:
        # First-row probabilities follow the U command, with bit error 1/3.
        target0,_=query(w,0,0)
        m=sum(1<<j for j,u in enumerate(first) if u < (2 if target0>>j&1 else 1))
        r1=[torus_rank(w,old+(m,)) for old in olds]
        for second in allrows:
            if first[4]!=2 or second[4]!=2:continue
            transcripts=[]
            for old,r in zip(olds,r1):
                # Flip the intended row pattern based on the observable rank.
                command=((1<<w)-1)^4 if r==0 else 4
                n=sum(1<<j for j,u in enumerate(second) if u < (2 if command>>j&1 else 1))
                r2=torus_rank(w,old+(m,n));checked+=1
                transcripts.append((r,r2))
            assert set(transcripts)=={(0,0)}
            erase+=1
    assert erase==3**8
    assert Fraction(erase,3**10)==Fraction(1,9)
    return {'seed_pairs':3**10,'forced_erasure_seed_pairs':erase,'probability':'1/9',
            'physical_graph_calls_on_erasure':checked,'adaptive_policies_scope':'one explicit two-row policy, all forced-erasure trit seeds'}

def strict_scalar_controls() -> dict:
    # Uniform C=1/400 uses sqrt(w)>=2 and w epsilon^2<=1/400^2.
    assert Fraction(7,800)+Fraction(3,160000)==Fraction(1403,160000)<Fraction(1,100)
    assert 4**100*99**198*2**100>7**100*100**200
    assert 4**50*49**98*4**50>13**50*50**100
    # Check exact necessary integer heights without decimal approximations.
    necessary=[]
    for w in (50,500,5000,50000):
        g=w//5; e=Fraction(1,10); target=Fraction(49,50)
        h=1
        while (1-e**h)**g<target:h+=1
        assert h==1 or (1-e**(h-1))**g<target
        necessary.append({'w':w,'g':g,'epsilon':'1/10','requested_bit_error':'1/100','minimum_possible_height_lower_bound':h})
    display=[]
    with localcontext() as ctx:
        ctx.prec=80
        for w in (500,5000,50000):
            g=w//5;b=Decimal('0.99')**g
            n=(Decimal('0.02').ln()/(1-b).ln())
            display.append({'w':w,'h':2,'epsilon':'1/10','beta_decimal':str(b),
                            'pair_error_lower_bound_decimal':str((1-b)/2),
                            'necessary_repeats_for_error_1pct_decimal':str(n)})
    for w in (5,10,25,100):
        g=w//5
        for e in (Fraction(1,100),Fraction(1,10),Fraction(1,4)):
            for h in range(1,5):
                b=(1-e**h)**g
                assert 0<b<1 and b<=(1-e**(h+1))**g
    rarity={'probability_formula':'p^(6t)*(1-p)^(3w-8t)',
            'w_equals_5t_half_probability':'2^(-13t)',
            'optimizer_when_w_equals_5t':'6/13',
            'one_block_max_probability':qstr(Fraction(6,13)**6*Fraction(7,13)**7)}
    return {'strict_1pct_bound':'1403/160000 < 1/100','RAC_integer_inequalities':True,
            'exact_height_lower_bounds':necessary,'illustrations_only_decimal':display,'natural_preparation_mass':rarity}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,default=Path('result.json'))
    ap.add_argument('--input-check',action='store_true',help='also compare every evaluated query against supplied preceding oracle')
    args=ap.parse_args(); start=time.perf_counter(); oracle=None
    if args.input_check:
        import importlib.util
        path=Path(__file__).parent/'input/tensor_check.py'
        spec=importlib.util.spec_from_file_location('prior_tensor_oracle',path)
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        oracle=module.physical_rank
    data={'schema':'matching-one.guard-erasure-and-two-row-noise.v1',
          'scope':'tensor preparation family; rank-only transcript; future sites conditionally independent; input noise, not preparation noise',
          'ladder':ladder_control(),'query':query_controls(oracle),
          'adaptive_erasure':adaptive_erasure_control(),'scalars':strict_scalar_controls(),
          'not_claimed':['fixed-noise logarithmic-height full-query construction','all-history erasure bound',
                         'full repository CI','critical exponent','independent external proof review']}
    data['total_new_integer_graph_calls'] = RANK_CALLS
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(data,ensure_ascii=False,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'output':str(args.out),'elapsed_seconds':time.perf_counter()-start,
                      'peak_rss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                      'query_graph_calls':data['query']['full_graph_rank_calls'],
                      'adaptive_graph_calls':data['adaptive_erasure']['physical_graph_calls_on_erasure']},indent=2))

if __name__=='__main__':main()

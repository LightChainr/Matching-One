#!/usr/bin/env python3
"""Exact site-cluster sewing, Palm debiasing, and readout of #741 values.

This is a small finite control, not a re-run of the width-12 builder.
Coordinates are (column,row); horizontal edges retain their integer lift.
Only complete components with minimum row zero are counted.  Boundary sites
outside the occupied window are present and carry their correct void weight.
Python standard library only; source logarithms are pinned decimal inputs.
"""
from __future__ import annotations
import argparse
from collections import Counter
from decimal import Decimal, localcontext
from fractions import Fraction
from itertools import product
import json
from pathlib import Path

Point = tuple[int, int]
STEPS4 = ((1,0),(-1,0),(0,1),(0,-1))
STEPS8 = tuple((dx,dy) for dx in (-1,0,1) for dy in (-1,0,1) if dx or dy)
SOURCE = '3745b13b8e1126017a1e88567a6af44881b8a1ff'
# Reported values, not recomputed stationary distributions.  The original
# certificate bounds are imported, not independently validated here.
RETURNED = {
 'NN': {'p':'1/4', 'logs': {'2':'-2.868379786881457','4':'-5.419513505387646',
   '8':'-10.023930900250674','12':'-14.400335955876589'}, 'err12':'5.546825221358631e-10'},
 'matching': {'p':'1/8', 'logs': {'2':'-3.3390984215201938','4':'-5.742669143987384',
   '8':'-10.180101524506828','12':'-14.471261214034028'}, 'err12':'1.4529177860822529e-10'}
}

def steps(matching: bool):
    return STEPS8 if matching else STEPS4

def vertices(mask: int, width: int, height: int) -> frozenset[Point]:
    return frozenset((i % width, i // width) for i in range(width*height) if mask >> i & 1)

def components(sites: frozenset[Point], width: int, matching: bool):
    """Independent BFS with full horizontal displacement, no vertical seam."""
    unseen = set(sites)
    out = []
    while unseen:
        root = min(unseen)
        unseen.remove(root)
        potential = {root:0}
        todo = [root]
        winding = False
        for a in todo:
            for dx,dy in steps(matching):
                b = ((a[0]+dx) % width, a[1]+dy)
                if b not in sites:
                    continue
                value = potential[a]+dx
                if b in potential:
                    residual = value-potential[b]
                    assert residual % width == 0
                    winding |= residual != 0
                else:
                    potential[b] = value
                    unseen.remove(b)
                    todo.append(b)
        out.append((frozenset(potential),winding))
    return out

def dsu_winding(sites: frozenset[Point], width: int, matching: bool) -> bool:
    """Independent weighted union/find, with seam gains rather than dx."""
    items = sorted(sites)
    index = {v:i for i,v in enumerate(items)}
    parent=list(range(len(items))); delta=[0]*len(items); wind=[False]*len(items)
    def find(i):
        if parent[i]!=i:
            old=parent[i]; r,d=find(old)
            delta[i]+=d; parent[i]=r
        return parent[i],delta[i]
    for a in items:
        for dx,dy in steps(matching):
            b=((a[0]+dx)%width,a[1]+dy)
            if b not in index:
                continue
            i,j=index[a],index[b]
            ri,di=find(i); rj,dj=find(j)
            gain=(a[0]+dx)//width
            if ri==rj:
                wind[ri] |= dj-di != gain
            else:
                parent[rj]=ri; delta[rj]=gain+di-dj
                wind[ri] |= wind[rj]
    return any(wind[find(i)[0]] for i in range(len(items)))

def boundary_direct(sites: frozenset[Point], width: int, matching: bool) -> frozenset[Point]:
    return frozenset(((x+dx)%width,y+dy) for x,y in sites for dx,dy in steps(matching))-sites

def boundary_columns(sites: frozenset[Point], width: int, matching: bool):
    cols=[{y for x,y in sites if x==i} for i in range(width)]
    out=[]
    for i,s in enumerate(cols):
        b={y+d for y in s for d in (-1,1)}
        for neighbor in (cols[(i-1)%width],cols[(i+1)%width]):
            b.update(y+d for y in neighbor for d in ((-1,0,1) if matching else (0,)))
        out.append(frozenset(b-s))
    return cols,out

def seam_marks(sites: frozenset[Point], width: int, matching: bool):
    # Mark every oriented edge crossing the *fixed* seam w-1 -> 0.
    # The row convention for a diagonal is its tail's row.
    edges=tuple(((width-1,y),(0,y+dy)) for x,y in sorted(sites) if x==width-1
                for dy in ((-1,0,1) if matching else (0,)) if (0,y+dy) in sites)
    return edges

def shape_table(width: int, height: int, matching: bool):
    """All complete candidate shapes, min row 0, contained in rows [0,H)."""
    table=Counter(); masks=0; connected_winding=0
    for mask in range(1,1 << (width*height)):
        c=vertices(mask,width,height)
        cols,bs=boundary_columns(c,width,matching)
        direct=boundary_direct(c,width,matching)
        assert direct==frozenset((i,y) for i,b in enumerate(bs) for y in b)
        parts=components(c,width,matching)
        assert any(w for _,w in parts)==dsu_winding(c,width,matching)
        masks+=1
        if len(parts)!=1 or not parts[0][1] or min(y for _,y in c)!=0:
            continue
        marks=seam_marks(c,width,matching)
        assert marks
        table[(len(c),len(direct),len(marks),len({a[1] for a,b in marks}))]+=1
        connected_winding+=1
    return table,{'nonempty_masks':masks,'anchored_winding_shapes':connected_winding}

def evaluate_table(table, p: Fraction, mode: str='edges'):
    """A component-Palm law and its c-size-biased marked law."""
    weights=Counter()
    for (n,b,ce,cr),multiplicity in table.items():
        c=ce if mode=='edges' else cr
        weights[c]+=multiplicity*p**n*(1-p)**b
    nu=sum(weights.values(),Fraction())
    mu=sum((c*a for c,a in weights.items()),Fraction())
    if not nu:
        raise ValueError('no winding shape in this window')
    e_c=mu/nu
    e_inv_component=sum((a/c for c,a in weights.items()),Fraction())/nu
    # Mark-Palm probability is c*a/mu, NOT a/nu.
    e_inv_mark=sum((c*a/mu/Fraction(c) for c,a in weights.items()),Fraction())
    assert mu*e_inv_mark==nu and e_inv_mark==1/e_c
    return {'nu_truncated':str(nu),'marked_intensity':str(mu),
            'component_mean_c':str(e_c),
            'component_mean_inverse_c':str(e_inv_component),
            'mark_mean_inverse_c':str(e_inv_mark),
            'wrong_component_reciprocal_estimate':str(mu*e_inv_component),
            'c_weights':{str(c):str(a) for c,a in sorted(weights.items())}}

def guard_enumeration(width: int, height: int, matching: bool, p: Fraction):
    """Enumerate H+2 actual rows; require target component min row 1,
    max row <= H. Guards are random; unrelated guard sites are NOT forced shut.
    Returns a direct product-measure expectation and the number enumerated.
    """
    n=width*(height+2)
    by_k=Counter()
    for mask in range(1 << n):
        s=vertices(mask,width,height+2)
        count=sum(w and min(y for _,y in c)==1 and max(y for _,y in c)<=height
                  for c,w in components(s,width,matching))
        by_k[mask.bit_count()]+=count
    val=sum((Fraction(c)*p**k*(1-p)**(n-k) for k,c in by_k.items()),Fraction())
    return val,1 << n

def exact_marking_counterexample():
    # C=[full bottom ring]+two separated teeth in the upper row.
    c=frozenset([(x,0) for x in range(4)]+[(0,1),(2,1)])
    cols,bs=boundary_columns(c,4,False)
    sum_contacts=0
    for x,y in c:
        sum_contacts+=sum(((x+dx)%4,y+dy) not in c for dx,dy in STEPS4)
    return {'sites':sorted(c),'occupied':len(c),'external_void_sites':len(boundary_direct(c,4,False)),
            'open_to_void_incidences':sum_contacts,'void_rows_by_column':[sorted(b) for b in bs],
            'overcount_if_weighted_per_incidence':sum_contacts-len(boundary_direct(c,4,False))}


def matmul(a, b):
    n=len(a)
    return [[sum((a[i][k]*b[k][j] for k in range(n)),Fraction())
             for j in range(n)] for i in range(n)]

def two_row_transfer_trace(width: int, matching: bool, p: Fraction) -> Fraction:
    """Exact complete winding-component activity inside a two-row strip.

    Column symbols 1,2,3 are nonempty subsets of the two rows. For NN,
    consecutive symbols must overlap; for matching they all communicate.
    Pair states retain the external-boundary overlap memory exactly.
    """
    alphabet=(1,2,3)
    def compatible(a,b): return matching or bool(a&b)
    states=[(a,b) for a in alphabet for b in alphabet if compatible(a,b)]
    ix={s:i for i,s in enumerate(states)}
    k=[[Fraction() for _ in states] for _ in states]
    def rows(a): return {r for r in (0,1) if a>>r&1}
    for i,(a,b) in enumerate(states):
        sb=rows(b)
        for c in alphabet:
            if not compatible(b,c): continue
            void={y+d for y in sb for d in (-1,1)}
            for t in (rows(a),rows(c)):
                void.update(y+d for y in t for d in ((-1,0,1) if matching else (0,)))
            void-=sb
            k[i][ix[b,c]]=p**len(sb)*(1-p)**len(void)
    result=[[Fraction(i==j) for j in range(len(states))] for i in range(len(states))]
    power=k
    t=width
    while t:
        if t&1: result=matmul(result,power)
        t//=2
        if t: power=matmul(power,power)
    return sum((result[i][i] for i in range(len(states))),Fraction())

def range_cdf(h: Decimal) -> Decimal:
    """Brownian-bridge RANGE CDF, not a proven site-cluster limit.

    Spectral series for small h, Poisson-dual series for large h.
    65-digit Decimal with an 60-digit internal truncation threshold.
    """
    if h<=0: return Decimal(0)
    with localcontext() as ctx:
        ctx.prec=65
        pi=Decimal('3.1415926535897932384626433832795028841971693993751058209749445923078')
        eps=Decimal('1e-60')
        if h<=Decimal('1.4'):
            answer=Decimal(0)
            for n in range(1,10000):
                term=pi*pi*n*n/h**3 * (-pi*pi*n*n/(2*h*h)).exp()
                answer+=term
                if term<eps: break
            return +(2*pi).sqrt()*answer
        answer=Decimal(1)
        for n in range(1,10000):
            term=2*(1-4*n*n*h*h)*(-2*n*n*h*h).exp()
            answer+=term
            if abs(term)<eps: break
        return +answer

def contrast_weights(widths):
    if len(widths)!=3 or len(set(widths))!=3 or min(widths)<=0:
        raise ValueError('three distinct positive widths required')
    x,y,z=widths
    return (z-y,x-z,y-x)

def contrast(logs, widths):
    with localcontext() as ctx:
        ctx.prec=65
        c=contrast_weights(widths)
        den=-sum(Decimal(a)*Decimal(w).ln() for a,w in zip(c,widths))
        return sum(Decimal(a)*logs[w] for a,w in zip(c,widths))/den

def contrast_report():
    with localcontext() as ctx:
        ctx.prec=65
        out={'source_commit':SOURCE,'source_path':'results/geometric-consistency/winding-prefactor-contrast.json',
             'source_note':'Input density values and solver error bounds imported; width-12 certificate not regenerated.',
             'models':{}}
        for name,raw in RETURNED.items():
            y={int(k):Decimal(v) for k,v in raw['logs'].items()}
            result={'source_log_nu':raw['logs'],'correct_beta_2_4_8':str(contrast(y,(2,4,8))),
                    'correct_beta_4_8_12':str(contrast(y,(4,8,12))),
                    'invalid_2_4_8_same_spacing_formula':str((y[2]-2*y[4]+y[8])/(Decimal(4)/3).ln()),
                    'correct_2_4_8_weights':[2,-3,1]}
            # .5 + c1/w interpolates the returned three rows exactly; no DOF.
            denominator=(Decimal(4)/3).ln()
            c1=Decimal(12)*(y[4]-2*y[8]+y[12]-Decimal('.5')*denominator)
            b=[y[w]+Decimal('.5')*Decimal(w).ln()-c1/Decimal(w) for w in (4,8,12)]
            kappa=-(b[1]-b[0])/4
            log_a=b[0]+4*kappa
            result['conditional_beta_half_plus_c_over_w']={
                'beta':'0.5','c1':str(c1),'kappa':str(kappa),'log_A':str(log_a),
                'scope':'post-return 3-parameter interpolation of 3 values, zero residual DOF; NOT evidence for beta=1/2'}
            # A four-width contrast cancels 1, w AND 1/w. Higher errors remain.
            c=(-4,15,-20,9); ws=(2,4,8,12)
            result['four_width_eliminate_one_over_w']=str(sum(Decimal(a)*y[w] for a,w in zip(c,ws))/
                (-sum(Decimal(a)*Decimal(w).ln() for a,w in zip(c,ws))))
            err=Decimal(raw['err12'])+Decimal('3e-15') # conservative displayed-input rounding
            result['propagated_4_8_12_log_input_bound']=str(err/denominator)
            out['models'][name]=result
        # Exact ansatz check: the wrong unequal-width statistic retains -2*kappa.
        return out

def run():
    output={'schema':'matching-one.sewing-identities.v1','date':'2026-09-13',
            'scope':'Exact finite-window algebra + return-value reanalysis; conjectures live in note, not numerical conclusions.',
            'factorization_controls':[], 'guard_controls':[], 'small_Palm_laws':[], 'two_row_trace_controls':[],
            'boundary_counterexample':exact_marking_counterexample(),
            'prefactor_return_reanalysis':contrast_report()}
    cache={}
    for matching in (False,True):
        for w,h in ((3,1),(3,2),(3,3),(4,1),(4,2),(4,3),(5,2)):
            table,counts=shape_table(w,h,matching); cache[w,h,matching]=table
            output['factorization_controls'].append({'width':w,'height':h,'matching':matching,**counts})
        for w in (3,4,5):
            for p in (Fraction(1,4),Fraction(1,2)):
                xi2=two_row_transfer_trace(w,matching,p)
                xi1=(p*(1-p)**2)**w
                anchored=Fraction(evaluate_table(cache[w,2,matching],p)['nu_truncated'])
                assert xi2-xi1==anchored
                output['two_row_trace_controls'].append({'width':w,'matching':matching,'p':str(p),
                    'pair_states':9 if matching else 7,'Xi2':str(xi2),'Xi1':str(xi1),
                    'anchored_height_at_most_two':str(anchored)})
        for p in (Fraction(1,4),Fraction(1,2)):
            for mode in ('edges','rows'):
                output['small_Palm_laws'].append({'width':4,'height':3,'matching':matching,'p':str(p),
                     'mark':mode,**evaluate_table(cache[4,3,matching],p,mode)})
        for w,h in ((3,1),(3,2),(4,1)):
            p=Fraction(1,3)
            actual,count=guard_enumeration(w,h,matching,p)
            expected=Fraction(evaluate_table(cache[w,h,matching],p)['nu_truncated'])
            assert actual==expected
            output['guard_controls'].append({'width':w,'height':h,'matching':matching,'p':str(p),
                    'enumerated_configurations':count,'direct_expectation':str(actual),
                    'component_activity_expectation':str(expected)})
    output['counts']={'factorization_nonempty_masks':sum(r['nonempty_masks'] for r in output['factorization_controls']),
                      'guard_configurations':sum(r['enumerated_configurations'] for r in output['guard_controls'])}
    output['conjectural_Brownian_range_targets']={str(h):str(range_cdf(Decimal(str(h)))) for h in (0.5,1,1.5,2,3)}
    return output

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    result=run()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    # Exclusive create: do not silently replace a published result.
    with args.output.open('x',encoding='utf-8') as f:
        json.dump(result,f,indent=2,ensure_ascii=False); f.write('\n')
    print(json.dumps(result['counts']))

if __name__=='__main__':
    main()

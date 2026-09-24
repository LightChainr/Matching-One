#!/usr/bin/env python3
"""Exact controls for a fixed-noise wedge query. Python standard library only.

This checker is not a Monte Carlo calculation. Its finite tests check the
geometric implications used in NOTE.md; the all-width bound is proved there.
"""
from __future__ import annotations
import argparse
from fractions import Fraction
from itertools import combinations, product
import importlib.util
import json
from pathlib import Path
import time

NN=((1,0),(-1,0),(0,1),(0,-1))
STAR=NN+((1,1),(1,-1),(-1,1),(-1,-1))
CALLS=0

def reaches(active, starts, targets, star=False):
    targets=set(targets); seen=set(starts)&active; todo=list(seen)
    steps=STAR if star else NN
    while todo:
        x,y=todo.pop()
        if (x,y) in targets:return True
        for dx,dy in steps:
            z=(x+dx,y+dy)
            if z in active and z not in seen:
                seen.add(z);todo.append(z)
    return False

def physical_rank(w, rows):
    """Lifted graph traversal on the physical torus, not a row-state solver."""
    global CALLS
    CALLS+=1
    h=len(rows); coords={}; first=None
    for y in range(h):
        for x in range(w):
            root=(x,y)
            if not rows[y]>>x&1 or root in coords:continue
            coords[root]=(0,0);todo=[root]
            while todo:
                a,b=todo.pop();px,py=coords[(a,b)]
                for dx,dy in NN:
                    z=((a+dx)%w,(b+dy)%h)
                    if not rows[z[1]]>>z[0]&1:continue
                    q=(px+dx,py+dy)
                    if z not in coords:coords[z]=q;todo.append(z)
                    else:
                        v=(q[0]-coords[z][0],q[1]-coords[z][1])
                        assert v[0]%w==0 and v[1]%h==0
                        if v!=(0,0):
                            if first is None:first=v
                            elif first[0]*v[1]-first[1]*v[0]:return 2
    return int(first is not None)

def prepare(w,bits):
    assert w>=5*len(bits)
    return [0,sum(15<<(5*g) for g in range(len(bits))),
            sum((9+2*u+4*v)<<(5*g) for g,(u,v) in enumerate(bits))]

def orient(row,w,g,bit):
    return sum(1<<((5*g+(x if bit==0 else 3-x))%w)
               for x in range(w) if row>>x&1)

def query(w,H,g=0,bit=0):
    """Two initial rows followed by H growing wedge layers; h=H+2."""
    assert H>=0 and 2*(H//4)+3<w
    full=(1<<w)-1;rows=[full^5]
    for y in range(H+1):
        r=y//4
        low=sum(1<<((2+x)%w) for x in range(-r,r+1))
        rows.append(full^low)
    return [orient(row,w,g,bit) for row in rows]

def cone(H):
    """Truncated cone after the four-site narrow stem, with NN boundary arcs."""
    assert H>=4
    D={(x,y) for y in range(4,H+1) for x in range(-(y//4),y//4+1)}
    left={(0,4),(-1,4)}
    for y in range(5,H+1):
        old=-((y-1)//4);new=-(y//4)
        left.add((old,y));left.add((new,y))
    right={(-x,y) for x,y in left}
    return D,left,right,{(x,H) for x in range(-(H//4),H//4+1)}

def trapezoid(W,H):
    assert W>2*(H//4)
    D={(x,y) for y in range(H+1) for x in range(y//4,W-y//4+1)}
    left={(0,0)}
    for y in range(1,H+1):
        f0=(y-1)//4;f1=y//4
        if f0!=f1:left.add((f1,y-1))
        left.add((f1,y))
    right={(W-x,y) for x,y in left}
    top={(x,H) for x in range(H//4,W-H//4+1)}
    bottom={(x,0) for x in range(W+1)}
    return D,left|right|top,bottom

def configurations(n,cutoff=None):
    if cutoff is None:
        for m in range(1<<n):yield [j for j in range(n) if m>>j&1]
    else:
        for k in range(min(cutoff,n)+1):yield from combinations(range(n),k)

def duality_controls():
    out=[]
    for H in range(4,9):
        D,L,R,top=cone(H);v=sorted(D);n=0
        for badidx in configurations(len(v)):
            bad={v[j] for j in badidx}
            good=reaches(D-bad,{(0,4)},top)
            cut=reaches(bad,L,R,True)
            assert good!=cut,('cone',H,bad,good,cut)
            n+=1
        out.append({'domain':'cone','height':H,'sites':len(v),'colorings':n})
    for W,H,c in ((2,2,None),(3,3,None),(6,4,3),(8,8,2),(16,12,2)):
        D,outer,bottom=trapezoid(W,H);v=sorted(D);n=0
        for badidx in configurations(len(v),c):
            bad={v[j] for j in badidx}
            good=reaches(D-bad,{(0,0)},{(W,0)})
            cut=reaches(bad,outer,bottom,True)
            assert good!=cut,('trapezoid',W,H,bad,good,cut)
            n+=1
        out.append({'domain':'trapezoid','width':W,'height':H,'sites':len(v),
                    'maximum_faults':c,'colorings':n})
    return out

def boundary_counts():
    controls=0
    for H in range(4,65):
        D,L,R,_=cone(H)
        assert L<=D and R<=D
        for n in range(1,H+2):
            possible={a for a in L if any(max(abs(a[0]-b[0]),abs(a[1]-b[1]))<=n-1 for b in R)}
            assert all(y<=4*n+12 for x,y in possible)
            assert len(possible)<=8*n+26
            controls+=1
        W=2*H+5;D,outer,bottom=trapezoid(W,H)
        assert outer<=D and bottom<=D
        assert len(outer)<=4*(H+1)+W+1
        for n in range(1,H+1):
            possible={a for a in outer if a[1]<=n-1}
            assert len(possible)<=4*n
            controls+=1
    return controls

def good_zero(rows,w,H):
    if rows[0]&5:return False
    C={(x,y) for y in range(H+1) for x in range(-(y//4),y//4+1)}
    white={z for z in C if not rows[z[1]+1]>>((2+z[0])%w)&1}
    return reaches(white,{(0,0)},{(x,H) for x in range(-(H//4),H//4+1)})

def good_one(rows,w,H):
    if rows[0]&10!=10:return False
    W=w-2;D,_,_=trapezoid(W,H)
    black={z for z in D if rows[z[1]+1]>>((3+z[0])%w)&1}
    return reaches(black,{(0,0)},{(W,0)})

def physical_controls(peer=None):
    out=[];peer_count=0;witness0=witness1=0
    for w,H,t,c in ((5,0,1,None),(10,4,2,1),(12,8,2,1),(16,8,3,1),(12,4,1,2)):
        base=query(w,H);sites=w*len(base);cases=0
        for bs in product(range(4),repeat=t):
            bits=[(b&1,b>>1) for b in bs];old=prepare(w,bits);target=bits[0][0]
            for f in configurations(sites,c):
                rows=base.copy()
                for j in f:rows[j//w]^=1<<(j%w)
                r=physical_rank(w,old+rows)
                if peer is not None and w!=16:
                    assert r==peer(w,old+rows)[0];peer_count+=1
                if not f:assert r==target,('fault-free',w,H,bits,r)
                if not target and good_zero(rows,w,H):
                    assert r==0,('zero-witness',w,H,bits,f,r);witness0+=1
                if target and good_one(rows,w,H):
                    assert r==1,('one-witness',w,H,bits,f,r);witness1+=1
                cases+=1
        out.append({'width':w,'height_H':H,'new_rows':H+2,'blocks':t,
                    'maximum_faults':c,'configurations':cases})
    rot=0
    for w,H,t in ((10,4,2),(12,8,2),(16,8,3)):
        for bs in product(range(4),repeat=t):
            bits=[(b&1,b>>1) for b in bs];old=prepare(w,bits)
            for g in range(t):
                for bit in (0,1):
                    assert physical_rank(w,old+query(w,H,g,bit))==bits[g][bit]
                    rot+=1
    return {'cases':out,'peer_comparisons':peer_count,'zero_witnesses':witness0,
            'one_witnesses':witness1,'rotated_queries':rot}

def open_signature(w,rows):
    """Exact integer boundary elimination on an open cylinder; no seam closure."""
    h=len(rows);seen={};blocks=[];horizontal=False
    for y in range(h):
        for x in range(w):
            if not rows[y]>>x&1 or (x,y) in seen:continue
            root=(x,y);seen[root]=0;todo=[root];live=[]
            while todo:
                a,b=todo.pop();px=seen[(a,b)]
                if b==h-1:live.append((a,px))
                for dx,dy in NN:
                    if not 0<=b+dy<h:continue
                    z=((a+dx)%w,b+dy)
                    if not rows[z[1]]>>z[0]&1:continue
                    q=px+dx
                    if z not in seen:seen[z]=q;todo.append(z)
                    elif q!=seen[z]:
                        assert (q-seen[z])%w==0
                        horizontal=True
            if live:blocks.append(live)
    normalized=[]
    for b in blocks:
        b.sort();o=b[0][1]
        normalized.append(tuple((x,0 if horizontal else g-o) for x,g in b))
    return horizontal,tuple(sorted(normalized))

def bottleneck_controls():
    """All masks, not only an event subset: retain the old bit iff m&15==10."""
    out=[]
    for w,t in ((5,1),(10,2),(12,2)):
        n=live=0
        for other in product((0,1),repeat=2*t-1):
            x=[0]+list(other);bits=list(zip(x[::2],x[1::2]));a=prepare(w,bits)
            bits[0]=(1,bits[0][1]);b=prepare(w,bits)
            for m in range(1<<w):
                different=open_signature(w,a+[m])!=open_signature(w,b+[m])
                assert different==(m&15==10),(w,other,m,different)
                n+=1;live+=different
        out.append({'width':w,'blocks':t,'queried_bit':0,
                    'first_row_comparisons':n,'difference_retaining_masks':live,
                    'necessary_and_sufficient_local_mask':10})
    # Recompute the first-interface TV independently from its two boundary laws.
    # The optional other bit is fixed in turn; all 2^5 masks are weighted exactly.
    laws=[]
    for eps in (Fraction(1,10),Fraction(1,100),Fraction(1,100000)):
        for other in (0,1):
            rows=[prepare(5,[(u,other)]) for u in (0,1)]
            distributions=[{},{}]
            for m in range(32):
                distance=(m^26).bit_count()  # high 1,3,4; low 0,2
                weight=eps**distance*(1-eps)**(5-distance)
                for u in (0,1):
                    key=open_signature(5,rows[u]+[m])
                    distributions[u][key]=distributions[u].get(key,Fraction(0))+weight
            keys=set(distributions[0])|set(distributions[1])
            tv=sum(abs(distributions[0].get(k,0)-distributions[1].get(k,0)) for k in keys)/2
            assert tv==(1-eps)**4
            laws.append({'epsilon':str(eps),'other_bit':other,'boundary_TV':str(tv)})
    return {'all_mask_checks':out,'exact_boundary_channel_checks':laws}

def probability_controls():
    e=Fraction(1,100000);a=7*e
    out=[]
    for w in (5,31,32,64,1000,10**6,10**12,10**24,10**48):
        if w<32:
            h=2;bound=2*w*e
        else:
            H=8
            while w*a**(H+1)>Fraction(1,10000):H+=1
            assert 2*H+8<=w
            h=H+2
            bound=512*e/(1-a)**2+8*w*a**(H+1)/(1-a)
        assert bound<Fraction(1,100)
        out.append({'width':str(w),'new_rows':h,'error_upper_bound':str(bound)})
    universal=512*e/(1-a)**2+Fraction(8,10000)/(1-a)
    assert universal<Fraction(3,500)
    # At each parity, w*a^(floor((w-8)/2)+1) decreases as w grows by 2.
    for w in (32,33):
        assert w*a**((w-8)//2+1)<=Fraction(1,10000)
        assert Fraction(w+2,w)*a<1
    # Exact finite geometric-sum check behind the contour count.
    for n in range(1,25):
        assert sum(Fraction(j)*a**(j-1) for j in range(1,n+1))==(
            1-(n+1)*a**n+n*a**(n+1))/(1-a)**2
    return {'epsilon':str(e),'universal_error_bound':str(universal),'examples':out,
            'finite_geometric_sum_checks':24,'interface_pair_error_floor':str((1-(1-e)**4)/2),
            'interface_error_floor_at_one_percent_noise':str((1-Fraction(99,100)**4)/2),
            'interface_error_floor_at_ten_percent_noise':str((1-Fraction(9,10)**4)/2)}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--input-check',action='store_true');args=ap.parse_args()
    peer=None
    if args.input_check:
        path=Path(__file__).resolve().parent/'input'/'tensor_check.py'
        spec=importlib.util.spec_from_file_location('parent_tensor',path)
        mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
        peer=mod.physical_rank
    start=time.perf_counter()
    out={'schema':'matching-one.fixed-noise-wedge-controls.v1',
         'duality':duality_controls(),'boundary_count_inequalities':boundary_counts(),
         'physical':physical_controls(peer),'interface_erasure':bottleneck_controls(),
         'probabilities':probability_controls()}
    out['physical_graph_rank_calls']=CALLS
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'elapsed_seconds_diagnostic':round(time.perf_counter()-start,6),
                      'physical_graph_rank_calls':CALLS,'out':str(args.out)}))

if __name__=='__main__':main()

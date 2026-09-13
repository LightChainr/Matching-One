#!/usr/bin/env python3
"""P398 is the periodic identified-connectivity O(1) TL chain on 2w ends.

Constructive all-width map; bounded exhaustive checks at w=2..8. No simulation.
Row generators act on functions. A deterministic map contributes Q_map-I.
Fattening identifies detach_i with e_(2i), join_i with e_(2i+1), 0-based.
The half-step rotation exchanges the two rate families, not arbitrary readouts.
"""
from __future__ import annotations
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
import json
from math import comb, factorial, gcd, lcm
from pathlib import Path
import time
from typing import Sequence

State = tuple[int, ...]


def canonical(values: Sequence[int]) -> State:
    labels: dict[int, int] = {}
    return tuple(labels.setdefault(value, len(labels)) for value in values)


@lru_cache(None)
def link_patterns(ends: int) -> tuple[State, ...]:
    if type(ends) is not int or ends < 0 or ends % 2:
        raise ValueError('a nonnegative even endpoint count is required')
    if not ends:
        return ((),)
    result = []
    for mate in range(1, ends, 2):
        for inside in link_patterns(mate-1):
            for outside in link_patterns(ends-mate-1):
                p = [-1]*ends
                p[0], p[mate] = mate, 0
                for i,j in enumerate(inside): p[i+1] = j+1
                for i,j in enumerate(outside): p[i+mate+1] = j+mate+1
                result.append(tuple(p))
    return tuple(result)


def fatten(state: Sequence[int]) -> State:
    """Block cyclic successor i->j gives a pair (2i+1,2j)."""
    if not state:
        return ()
    state = canonical(state)
    p = [-1]*(2*len(state))
    for label in range(max(state)+1):
        vertices = [i for i,x in enumerate(state) if x == label]
        for i,j in zip(vertices, vertices[1:]+vertices[:1]):
            p[2*i+1], p[2*j] = 2*j, 2*i+1
    return tuple(p)


def unfatten(pairing: Sequence[int]) -> State:
    n = len(pairing)
    if n % 2 or any(type(j) is not int or not 0 <= j < n
                    or pairing[j] != i or j == i or (i+j)%2 != 1
                    for i,j in enumerate(pairing)):
        raise ValueError('expected a bipartite fixed-point-free pairing')
    labels = [-1]*(n//2)
    for i in range(n//2):
        if labels[i] >= 0: continue
        j = i
        while labels[j] < 0:
            labels[j] = i
            j = pairing[2*j+1]//2
    return canonical(labels)


def join(state: Sequence[int], i: int) -> State:
    if not 0 <= i < len(state): raise ValueError('point outside state')
    a,b = state[i],state[(i+1)%len(state)]
    return canonical([a if x == b else x for x in state])


def detach(state: Sequence[int], i: int) -> State:
    if not 0 <= i < len(state): raise ValueError('point outside state')
    q = list(state); q[i] = max(state)+1
    return canonical(q)


def tl_reconnect(pairing: Sequence[int], a: int) -> State:
    n = len(pairing)
    if not 0 <= a < n: raise ValueError('endpoint outside pairing')
    b = (a+1)%n
    c,d = pairing[a],pairing[b]
    if c == b: return tuple(pairing)
    q = list(pairing)
    q[a],q[b],q[c],q[d] = b,a,d,c
    return tuple(q)


def rotate(pairing: Sequence[int], shift: int = 1) -> State:
    n = len(pairing)
    if not n: return ()
    q = [-1]*n
    for i,j in enumerate(pairing): q[(i+shift)%n] = (j+shift)%n
    return tuple(q)


def half_step(state: Sequence[int]) -> State:
    return unfatten(rotate(fatten(state)))


def reflection(pairing: Sequence[int]) -> State:
    n=len(pairing); q=[-1]*n
    for i,j in enumerate(pairing): q[n-1-i] = n-1-j
    return tuple(q)


def independent_rgs(width: int) -> list[State]:
    """Independent restricted-growth enumeration, only used at widths <=6."""
    quads=list(combinations(range(width),4))
    result=[]
    def grow(prefix: State):
        if len(prefix) == width:
            if not any(prefix[a]==prefix[c] and prefix[b]==prefix[d]
                       and prefix[a]!=prefix[b] for a,b,c,d in quads):
                result.append(prefix)
            return
        for k in range(max(prefix, default=-1)+2): grow(prefix+(k,))
    grow(())
    return result


def state_space(width: int) -> list[State]:
    if type(width) is not int or not 2 <= width <= 8:
        raise ValueError('bounded executable controls require width 2..8')
    return sorted(unfatten(p) for p in link_patterns(2*width))


def generator(states: Sequence[State], eta: Fraction = Fraction(0)) -> list[list[Fraction]]:
    if abs(eta)>1: raise ValueError('Markov rates require |eta|<=1')
    n=len(states);w=len(states[0]);index={s:i for i,s in enumerate(states)}
    out=[[Fraction(0)]*n for _ in range(n)]
    for row,s in enumerate(states):
        for i in range(w):
            for fn,rate in ((join,1+eta),(detach,1-eta)):
                out[row][index[fn(s,i)]] += rate
                out[row][row] -= rate
    return out


def stationary(g: list[list[Fraction]]) -> list[Fraction]:
    n=len(g)
    a=[[g[i][j] for i in range(n)]+[Fraction(0)] for j in range(n)]
    a[-1]=[Fraction(1)]*(n+1)
    for k in range(n):
        pivot=next((j for j in range(k,n) if a[j][k]),None)
        if pivot is None: raise ValueError('stationary system singular')
        a[k],a[pivot]=a[pivot],a[k]
        div=a[k][k];a[k]=[x/div for x in a[k]]
        for j in range(n):
            if j != k and a[j][k]:
                mult=a[j][k];a[j]=[x-mult*y for x,y in zip(a[j],a[k])]
    pi=[row[-1] for row in a]
    assert sum(pi)==1 and min(pi)>0
    assert all(sum(pi[i]*g[i][j] for i in range(n))==0 for j in range(n))
    return pi


def asm_number(width: int) -> int:
    ans=Fraction(1)
    for j in range(width): ans*=Fraction(factorial(3*j+1),factorial(width+j))
    if ans.denominator!=1: raise AssertionError('ASM product not integral')
    return ans.numerator


def report(max_width: int = 8) -> dict:
    if not 2 <= max_width <= 8: raise ValueError('max_width must be 2..8')
    start=time.perf_counter();rows=[];total_moves=0
    for w in range(2,max_width+1):
        ss=state_space(w);pairset=set(link_patterns(2*w));n=len(ss)
        assert n==len(set(ss))==comb(2*w,w)//(w+1)
        if w<=6: assert ss==independent_rgs(w)
        for s in ss:
            p=fatten(s);k=half_step(s)
            assert p in pairset and unfatten(p)==s and fatten(unfatten(p))==p
            assert max(k)+max(s)+2==w+1
            assert fatten(canonical(s[::-1]))==reflection(p)
            assert half_step(k)==canonical((s[-1],)+s[:-1])
            for i in range(w):
                assert fatten(detach(s,i))==tl_reconnect(p,2*i)
                assert fatten(join(s,i))==tl_reconnect(p,2*i+1)
                assert half_step(detach(s,i))==join(k,i)
                assert half_step(join(s,i))==detach(k,(i+1)%w)
                assert int(s[i]==s[(i+1)%w])==int(k.count(k[(i+1)%w])==1)
                assert int(s.count(s[i])==1)==int(k[i]==k[(i+1)%w])
            # Named TL relations: idempotence, adjacent sandwich and distant commuting.
            for a in range(2*w):
                e=tl_reconnect(p,a)
                assert tl_reconnect(e,a)==e
                for b in ((a-1)%(2*w),(a+1)%(2*w)):
                    assert tl_reconnect(tl_reconnect(e,b),a)==e
                for b in range(a+1,2*w):
                    if (b-a)%(2*w) not in (1,2*w-1):
                        assert tl_reconnect(e,b)==tl_reconnect(tl_reconnect(p,b),a)
        item={'width':w,'states':n,'move_equalities':2*w*n,
              'bijection_and_cyclic_seam':True,'TL_relations':True,
              'half_step_square_is_site_rotation':True,
              'rate_swap_and_observer_identities':True}
        total_moves+=2*w*n
        if w<=5:
            g=generator(ss);pi=stationary(g)
            den=lcm(*(q.denominator for q in pi));ivec=[int(q*den) for q in pi]
            common=gcd(*ivec);ivec=[x//common for x in ivec]
            assert sum(ivec)==asm_number(w)
            mean=sum(q*(max(s)+1) for q,s in zip(pi,ss))
            assert mean==Fraction(w+1,2)
            index={s:i for i,s in enumerate(ss)};perm=[index[half_step(s)] for s in ss]
            gp=generator(ss,Fraction(1,4));gm=generator(ss,Fraction(-1,4))
            assert all(gp[i][j]==gm[perm[i]][perm[j]] for i in range(n) for j in range(n))
            pp,pm=stationary(gp),stationary(gm)
            assert all(pp[i]==pm[perm[i]] for i in range(n))
            item['stationary']={'states_lexicographic':[list(s) for s in ss],
                  'primitive_integer_weights':ivec,'weight_sum':sum(ivec),
                  'ASM_product':asm_number(w),'mean_blocks':str(mean),
                  'eta_quarter_law_conjugacy':True}
        rows.append(item)
    return {'schema':'matching-one.p398-tl-identification.v1',
        'status':'constructive all-width identity; finite controls not a novelty claim',
        'map':'detach_i <-> e_(2i); join_i <-> e_(2i+1), zero-based; loop weight 1',
        'boundary':'periodic identified-connectivity disk link patterns, 2w endpoints',
        'row_generator':'sum_a (E_a-I)^T in the standard column-state TL convention',
        'total_move_equalities':total_moves,'widths':rows,
        'elapsed_seconds':time.perf_counter()-start}


def main() -> None:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',type=Path);ap.add_argument('--max-width',type=int,default=8)
    args=ap.parse_args();text=json.dumps(report(args.max_width),indent=2,allow_nan=False)+'\n'
    if args.out:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        with args.out.open('x',encoding='utf-8') as f:f.write(text)
    else: print(text,end='')

if __name__=='__main__':main()

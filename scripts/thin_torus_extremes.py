#!/usr/bin/env python3
"""Exact geometric controls for thin square-site torus rank births.

This is a finite control for the proofs in thin-torus-two-birth-limits-20260912.md.
It runs no Monte Carlo and uses no fixed-width spectrum to classify a graph.
All counts and birth moments in this file are integers or Fraction strings.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, permutations, product
import json
from math import comb, factorial
from pathlib import Path
import time


def central_trinomial(width: int) -> int:
    if not isinstance(width, int) or isinstance(width, bool) or width < 2:
        raise ValueError('integer width >= 2 is required')
    return sum(factorial(width)//(factorial(k)**2*factorial(width-2*k))
               for k in range(width//2+1))


def ambient_rank(mask: int, width: int, length: int, matching: bool = False) -> int:
    """Independent physical-displacement traversal; parallel edges retained."""
    if min(width, length) < 2 or mask < 0 or mask >> (width*length):
        raise ValueError('invalid honest rectangular torus or occupation mask')
    steps = [(1,0),(-1,0),(0,1),(0,-1)]
    if matching:
        steps += [(1,1),(1,-1),(-1,1),(-1,-1)]
    positions = {}
    direction = None
    for root in range(width*length):
        if not (mask >> root) & 1 or root in positions:
            continue
        positions[root] = (0,0)
        stack = [root]
        while stack:
            v = stack.pop()
            x,y = v % width, v // width
            px,py = positions[v]
            for dx,dy in steps:
                u = ((y+dy) % length)*width+(x+dx) % width
                if not (mask >> u) & 1:
                    continue
                proposal = (px+dx,py+dy)
                if u not in positions:
                    positions[u] = proposal
                    stack.append(u)
                else:
                    hx = proposal[0]-positions[u][0]
                    hy = proposal[1]-positions[u][1]
                    if hx % width or hy % length:
                        raise AssertionError('cycle did not close modulo periods')
                    hx,hy = hx//width,hy//length
                    if hx or hy:
                        if direction is None:
                            direction = (hx,hy)
                        elif direction[0]*hy-direction[1]*hx:
                            return 2
    return int(direction is not None)


def minimal_matching_motifs(width: int, length: int) -> set[int]:
    """All size-w transverse matching cycles, encoded by their vertex sets."""
    if length <= 2*width:
        raise ValueError('motif control uses length > 2*width')
    motifs = set()
    for steps in product((-1,0,1), repeat=width):
        if sum(steps):
            continue
        for y0 in range(length):
            y = y0
            mask = 0
            for x,dy in enumerate(steps):
                mask |= 1 << ((y % length)*width+x)
                y += dy
            motifs.add(mask)
    return motifs


def motif_control(width: int) -> dict:
    length = 2*width+1
    predicted = minimal_matching_motifs(width,length)
    observed = set()
    checked = 0
    for subset in combinations(range(width*length),width):
        mask = sum(1 << i for i in subset)
        checked += 1
        if ambient_rank(mask,width,length,matching=True):
            observed.add(mask)
    if observed != predicted or len(observed) != length*central_trinomial(width):
        raise AssertionError('minimal-cycle classification failed')
    overlap = Counter()
    motifs = sorted(predicted)
    for a,b in combinations(motifs,2):
        if a & b:
            overlap[(a|b).bit_count()] += 1
    return {'width':width,'length':length,'subsets_checked':checked,
            'central_trinomial':central_trinomial(width),
            'minimal_cycle_sets':len(observed),
            'classification_exact':True,
            'unordered_overlapping_pairs_by_union_size':dict(sorted(overlap.items()))}


def barrier_control(width: int, length: int) -> dict:
    """Every configuration with both a full and an empty row has rank one."""
    full = (1 << width)-1
    checked = qualifying = 0
    for mask in range(1 << (width*length)):
        rows = [(mask >> (j*width)) & full for j in range(length)]
        checked += 1
        if full in rows and 0 in rows:
            qualifying += 1
            if ambient_rank(mask,width,length) != 1:
                raise AssertionError('full/empty barrier implication failed')
    return {'width':width,'length':length,'configurations_checked':checked,
            'full_and_empty_row_cases':qualifying,'rank_one_in_every_case':True}


def first_two_birth_counts(width: int, length: int) -> dict:
    """Subset DP for exact ordered-permutation birth ranks; not a joint-law fit."""
    n = width*length
    if n > 12:
        raise ValueError('this control is deliberately limited to 12 sites')
    ranks = [ambient_rank(mask,width,length) for mask in range(1 << n)]
    # State is (current subset, first birth rank, zero if not yet born).
    layer = {(0,0):1}
    joint = Counter()
    for k in range(n):
        new = defaultdict(int)
        for (mask,first),count in layer.items():
            for v in range(n):
                if (mask >> v) & 1:
                    continue
                nxt = mask | (1 << v)
                r = ranks[nxt]
                a = first or ((k+1) if r else 0)
                if r == 2:
                    joint[a,k+1] += count*factorial(n-k-1)
                else:
                    new[nxt,a] += count
        layer = new
    if sum(joint.values()) != factorial(n):
        raise AssertionError('permutation mass lost')
    total = factorial(n)
    def mean(func):
        return sum(Fraction(v,total)*func(a,b) for (a,b),v in joint.items())
    t1 = mean(lambda a,b: Fraction(a,n+1))
    t2 = mean(lambda a,b: Fraction(b,n+1))
    # Conditional on the permutation: E U_(a) U_(b)=a(b+1)/[(n+1)(n+2)].
    product_mean = mean(lambda a,b: Fraction(a*(b+1),(n+1)*(n+2)))
    # Check the histogram/binomial reconstruction at three rational p values.
    recon = []
    for p in (Fraction(1,3),Fraction(1,2),Fraction(2,3)):
        tail = [sum(Fraction(comb(n,k))*p**k*(1-p)**(n-k)
                    for k in range(j,n+1)) for j in range(n+1)]
        by_births = mean(lambda a,b: (tail[a]+tail[b])/2)
        by_graph = sum(Fraction(ranks[mask],2)*p**mask.bit_count()*
                       (1-p)**(n-mask.bit_count()) for mask in range(1 << n))
        if by_births != by_graph:
            raise AssertionError('birth histogram does not reconstruct the rank CDF')
        recon.append({'p':str(p),'F_exact':str(by_births)})
    return {'width':width,'length':length,'sites':n,'permutations_counted':total,
            'joint_K1_K2':[[a,b,c] for (a,b),c in sorted(joint.items())],
            'continuous_birth_mean_T1':str(t1),'continuous_birth_mean_T2':str(t2),
            'continuous_birth_covariance':str(product_mean-t1*t2),
            'finite_births_claimed_independent':False,'rank_CDF_reconstruction':recon}


def brute_permutation_control() -> dict:
    w,m,n = 2,3,6
    ranks = [ambient_rank(mask,w,m) for mask in range(1 << n)]
    counts = Counter()
    for word in permutations(range(n)):
        mask = a = b = 0
        for j,v in enumerate(word,1):
            mask |= 1 << v
            if ranks[mask] >= 1 and not a:
                a = j
            if ranks[mask] == 2:
                b = j
                break
        counts[a,b] += 1
    dp = first_two_birth_counts(w,m)
    if [[a,b,c] for (a,b),c in sorted(counts.items())] != dp['joint_K1_K2']:
        raise AssertionError('independent permutation loop disagrees')
    return {'width':w,'length':m,'permutations_explicitly_checked':factorial(n),
            'DP_matches_explicit_permutations':True}


def report() -> dict:
    started = time.perf_counter()
    return {'schema':'matching-one.thin-torus-geometric-controls.v1',
            'standing':'Exact finite controls; the all-width asymptotic proof is in the note.',
            'motifs':[motif_control(w) for w in (2,3,4)],
            'barriers':[barrier_control(2,3),barrier_control(3,4)],
            'births':[first_two_birth_counts(2,3),first_two_birth_counts(3,3)],
            'permutation_control':brute_permutation_control(),
            'elapsed_seconds':time.perf_counter()-started,
            'limits':['No new Monte Carlo.','Motif checks do not replace the all-width proof.',
                      'Finite joint birth times are not asserted independent.',
                      'Thin rectangles are not the existing Gaussian square-size lineage.']}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path)
    args = parser.parse_args()
    text = json.dumps(report(),indent=2,allow_nan=False)+'\n'
    if args.out:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        with args.out.open('x',encoding='utf-8') as f:
            f.write(text)
    else:
        print(text,end='')

if __name__ == '__main__':
    main()

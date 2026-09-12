#!/usr/bin/env python3
"""Independent exact-grid and projection checks of a saved source quotient.

This verifier does not call the polynomial-coefficient refinement routine.
A degree-complete rational tensor grid certifies each full polynomial family.
"""
from __future__ import annotations
import argparse
from collections import defaultdict
from fractions import Fraction
from itertools import product
import json
from pathlib import Path
import time
from width4_site_sources import DEFAULT_CERTIFICATE, ROOT, load_certificate

DEFAULT_RESULT = ROOT/'results/research-control-20260912/width4-site-source-quotients.json'


def verify(source, result):
    started=time.perf_counter()
    rows=source['quotient_transitions']
    ranks=source['quotient_rank_output']
    checks=0
    by_profile={}
    for name, entry in result['profiles'].items():
        groups=entry['groups'];labels=entry['labels']
        if len(labels)!=509 or len(set(labels))!=entry['classes']:
            raise AssertionError('incorrect partition dimensions')
        grid=[tuple(Fraction(k+1,len(g)+2) for k in range(len(g)+1)) for g in groups]
        point_count=0
        for values in product(*grid):
            ps=[None]*4
            for g,p in zip(groups,values):
                for j in g:ps[j]=p
            denominator=1
            for p in ps:denominator*=p.denominator
            weights=[]
            for mask in range(16):
                probability=Fraction(1)
                for j,p in enumerate(ps):probability*=p if mask&(1<<j) else 1-p
                integer=probability*denominator
                if integer.denominator!=1:raise AssertionError('bad common denominator')
                weights.append(integer.numerator)
            if sum(weights)!=denominator:raise AssertionError('row not normalized')
            representatives={}
            for i,row in enumerate(rows):
                target=defaultdict(int)
                for b,j in enumerate(row):target[labels[j]]+=weights[b]
                signature=(ranks[i],tuple(sorted(target.items())))
                if labels[i] in representatives:
                    if signature!=representatives[labels[i]]:
                        raise AssertionError('saved partition fails exact tensor-grid row equality')
                else:representatives[labels[i]]=signature
                checks+=1
            point_count+=1
        by_profile[name]={'rational_tensor_grid_points':point_count,
                          'classes':entry['classes'],'all_source_rows_checked':509}
    # Project H from both sides onto D4-invariant functions. Integer block
    # sums suffice for the zero check, so no inverse or numerical rank is used.
    labels=result['profiles']['homogeneous']['labels']
    totals=defaultdict(int)
    individual_nonzero=0
    for i,row in enumerate(rows):
        image=defaultdict(int)
        for mask,j in enumerate(row):
            image[labels[j]]+=4*((mask&1)-((mask>>2)&1))
        individual_nonzero+=sum(bool(v) for v in image.values())
        for block,value in image.items():totals[labels[i],block]+=value
    if any(totals.values()) or individual_nonzero==0:
        raise AssertionError('odd projection test fails')
    return {'schema':'matching-one.width4-site-source-independent-grid.v1',
            'profiles':by_profile,'exact_state_rows_checked':checks,
            'all_family_identities_certified_by_degree_complete_grids':True,
            'D4_averaged_odd_operator_is_zero':True,
            'unprojected_H_times_invariant_indicators_nonzero_entries':individual_nonzero,
            'standing':'finite exact common-lumping/projection verification, not generic positive minimality',
            'elapsed_seconds':time.perf_counter()-started}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--certificate',type=Path,default=DEFAULT_CERTIFICATE)
    p.add_argument('--result',type=Path,default=DEFAULT_RESULT)
    p.add_argument('--out',type=Path)
    args=p.parse_args()
    result=verify(load_certificate(args.certificate),json.loads(args.result.read_text()))
    text=json.dumps(result,indent=2)+'\n'
    if args.out:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        with args.out.open('x',encoding='utf-8') as f:f.write(text)
    else:print(text,end='')

if __name__=='__main__':main()

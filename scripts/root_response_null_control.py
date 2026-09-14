#!/usr/bin/env python3
"""Small algebraic checks used to CHOOSE a research task, not a size census.

Virasoro scope: negative-mode PBW algebra, quotient by the level-two
singular vector at c=0,h=5/8, followed by quotient by L_-1 derivatives.
This does not identify the actual percolation thermal representation.
The probability family below is deliberately a countermodel, not site data.
"""
from __future__ import annotations
from collections import defaultdict
from fractions import Fraction as F
from functools import lru_cache
import json
import math
from pathlib import Path
import argparse

@lru_cache(None)
def partitions(n: int, cap: int | None = None) -> tuple[tuple[int, ...], ...]:
    if n == 0:
        return ((),)
    cap = n if cap is None else min(cap, n)
    return tuple((k,) + rest for k in range(cap, 0, -1)
                 for rest in partitions(n-k, k))

@lru_cache(None)
def normal(word: tuple[int, ...]) -> tuple[tuple[tuple[int, ...], F], ...]:
    """L_-a L_-b = L_-b L_-a + (b-a)L_-(a+b) when a<b."""
    for i in range(len(word)-1):
        a, b = word[i:i+2]
        if a < b:
            out = defaultdict(F)
            swapped = word[:i] + (b, a) + word[i+2:]
            joined = word[:i] + (a+b,) + word[i+2:]
            for w, c in normal(swapped):
                out[w] += c
            for w, c in normal(joined):
                out[w] += (b-a)*c
            return tuple(sorted((w,c) for w,c in out.items() if c))
    return ((word, F(1)),)

def vector(word: tuple[int, ...], basis: tuple[tuple[int,...], ...]) -> list[F]:
    d = dict(normal(word))
    return [d.get(x, F(0)) for x in basis]

def rank(rows: list[list[F]], ncols: int) -> int:
    a = [row[:] for row in rows]
    r = 0
    for c in range(ncols):
        pivot = next((i for i in range(r, len(a)) if a[i][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        val = a[r][c]
        a[r] = [x/val for x in a[r]]
        for i in range(r+1, len(a)):
            val = a[i][c]
            if val:
                a[i] = [x-val*y for x,y in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r

def module_rows(max_level: int = 4) -> list[dict]:
    rows = []
    for n in range(max_level+1):
        basis = partitions(n)
        nulls = []
        if n >= 2:
            for a in partitions(n-2):
                v2 = vector(a+(2,), basis)
                v11 = vector(a+(1,1), basis)
                nulls.append([x-F(2,3)*y for x,y in zip(v2, v11)])
        derivs = [vector((1,)+a, basis) for a in partitions(n-1)] if n else []
        nr = rank(nulls, len(basis))
        both = rank(nulls+derivs, len(basis))
        rows.append(dict(level=n, verma_dimension=len(basis), null_rank=nr,
                         null_quotient_dimension=len(basis)-nr,
                         derivative_image_rank_in_quotient=both-nr,
                         nonderivative_quotient_dimension=len(basis)-both,
                         pbw_basis=[list(t) for t in basis],
                         null_rows=[[str(x) for x in v] for v in nulls],
                         derivative_rows=[[str(x) for x in v] for v in derivs]))
    return rows

def check() -> dict:
    h, central = F(5,8), F(0)
    singular_checks = {
        'L1_coefficient': str(3-F(2,3)*(4*h+2)),
        'L2_coefficient': str(4*h+central/2-F(2,3)*6*h),
    }
    assert all(v == '0' for v in singular_checks.values())
    rows = module_rows()
    assert [r['null_quotient_dimension'] for r in rows] == [1,1,1,2,3]
    assert [r['nonderivative_quotient_dimension'] for r in rows] == [1,0,0,1,1]
    assert dict(normal((1,2))) == {(2,1):F(1), (3,):F(1)}
    assert dict(normal((1,3))) == {(3,1):F(1), (4,):F(2)}
    # Fully positive laws: intrinsic c=2cosh(b), independent of root displacement.
    curves = []
    for a in (F(-1,5), F(0), F(1,5)):
        for L in (4,8,16):
            root = F(3,5) + a/F(L**4)
            for b in (F(-1,2), F(0), F(1,2)):
                cb = 2*math.cosh(float(b))
                den = 2*(cb+math.cosh(float(b)))
                probs = [math.exp(float(b))/den, 2*cb/den,
                         math.exp(-float(b))/den]
                assert min(probs)>0 and abs(sum(probs)-1)<1e-14
                assert abs(math.log(probs[1]/(2*math.sqrt(probs[0]*probs[2])))-math.log(cb)) < 1e-14
            curves.append({'L':L,'amplitude':str(a),'root':str(root),
                           'intrinsic_odd_part':'0 exactly (c=2cosh b)'})
    # Exact tangent/normal decomposition and nuisance-coordinate invariance.
    response_checks=0
    for bp in (F(-2),F(-1,3)):
        for ep in (F(0),F(3,7)):
            for bg in (F(1,5),F(-4,9)):
                for eg in (F(0),F(7,11)):
                    normal_response=eg-ep*bg/bp
                    for mixing in (F(-3),F(0),F(2,7)):
                        assert (eg+mixing*ep)-ep*(bg+mixing*bp)/bp == normal_response
                        response_checks+=1
    return {'scope':'Algebraic selection aid; no new percolation simulation or all-size proof',
            'virasoro_c':'0','highest_weight':'5/8',
            'singular_vector':'(L_-2-(2/3)L_-1^2)|h>',
            'singular_checks':singular_checks, 'levels':rows,
            'physical_module_identified':False,
            'counterfamily':curves,
            'response_mixing_equalities':response_checks,
            'check_summary':{'singular_equalities':2, 'module_dimension_sequences':2,
                             'PBW_commutators':2,'positive_law_cases':27,
                             'response_mixing_equalities':response_checks}}

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    report=check()
    text=json.dumps(report,ensure_ascii=False,indent=2,sort_keys=True)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(text,encoding='utf-8')
    else:
        print(text,end='')

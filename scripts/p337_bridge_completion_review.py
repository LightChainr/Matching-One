#!/usr/bin/env python3
"""Exact two-shared-component factorization and completion-robust pair check.

Standard library only. Colours are equality labels, not an integer-Q fit.
This script sums finite equality patterns and creates no random configurations.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import factorial
from pathlib import Path
from typing import Tuple

Pattern = Tuple[int, ...]

@lru_cache(maxsize=None)
def partitions(n: int) -> tuple:
    if type(n) is not int or n < 1:
        raise ValueError("positive integer size required")
    out = []
    def visit(a: tuple, largest: int) -> None:
        if len(a) == n:
            out.append(a)
            return
        for x in range(largest + 2):
            visit(a + (x,), max(largest, x))
    visit((0,), 0)
    return tuple(out)


def validate_pattern(p: Pattern, length: int) -> Pattern:
    p = tuple(p)
    if len(p) != length or any(type(x) is not int or x < 0 for x in p):
        raise ValueError("invalid exterior pattern")
    if set(p) != set(range(max(p) + 1)):
        raise ValueError("labels must be contiguous from zero")
    return p


def raw_kernel_at_one(c: Pattern, alpha: F = F(0)) -> F:
    a, b, c_, d = c
    if a == b or c_ == d:
        return F(0)
    crossing = int(a == c_ and b == d) + int(a == d and b == c_)
    singles = int(a == c_) + int(a == d) + int(b == c_) + int(b == d)
    return F(crossing + singles - 4, 2) + alpha


@lru_cache(maxsize=None)
def kernel_at_one(c: Pattern, alpha: F = F(0)) -> F:
    """K2+c(Q)K0 at Q1, with c(1)=1 and c'(1)=alpha."""
    return (raw_kernel_at_one(c, alpha) + raw_kernel_at_one(c[1:]+c[:1], alpha))/2


def raw_kernel_at_integer(q: int, c: Pattern, alpha: F = F(0)) -> F:
    a, b, c_, d = c
    if q < 3:
        raise ValueError("finite-colour check requires Q>=3")
    if a == b or c_ == d:
        return F(0)
    crossing = int(a == c_ and b == d) + int(a == d and b == c_)
    singles = int(a == c_) + int(a == d) + int(b == c_) + int(b == d)
    return F(1,2)*(crossing-F(singles,q-2)+F(4,q*(q-2)))+alpha/F(q)


def kernel_at_integer(q: int, c: Pattern, alpha: F = F(0)) -> F:
    return (raw_kernel_at_integer(q,c,alpha)+raw_kernel_at_integer(q,c[1:]+c[:1],alpha))/2


@lru_cache(maxsize=None)
def pair_activation(pi: Pattern, alpha_x: F = F(0), alpha_y: F = F(0)) -> F:
    pi = validate_pattern(pi, 8)
    out = F(0)
    for rho in partitions(max(pi)+1):
        k = max(rho)+1
        if k < 2:
            continue
        c = tuple(rho[x] for x in pi)
        multiplicity = (-1)**(k-2)*factorial(k-2)
        out += multiplicity*kernel_at_one(c[:4],alpha_x)*kernel_at_one(c[4:],alpha_y)
    return out


def distinguished_patterns() -> tuple:
    """62 four-port patterns with ordered bridge labels 0 and 1."""
    out = set()
    for pi in partitions(4):
        b = max(pi)+1
        for a in range(b):
            for c in range(b):
                if a == c:
                    continue
                mapping = {a: 0, c: 1}
                for x in pi:
                    if x not in mapping:
                        mapping[x] = len(mapping)
                out.add(tuple(mapping[x] for x in pi))
    return tuple(sorted(out))


@lru_cache(maxsize=None)
def endpoint_signature(pi: Pattern, alpha: F = F(0)) -> F:
    pi = validate_pattern(pi, 4)
    if 0 not in pi or 1 not in pi:
        raise ValueError("both bridge components must reach this mark")
    out = F(0)
    for rho in partitions(max(pi)+1):
        if rho[0] == rho[1]:
            continue
        fresh = max(rho)-1
        multiplicity = (-1)**fresh*factorial(fresh)
        out += multiplicity*kernel_at_one(tuple(rho[x] for x in pi),alpha)
    return out


def bridge_endpoints(pi: Pattern) -> tuple:
    pi = validate_pattern(pi,8)
    common = sorted(set(pi[:4]) & set(pi[4:]))
    if len(common) != 2:
        raise ValueError("exactly two shared exterior components required")
    result = []
    for half in (pi[:4],pi[4:]):
        mapping = {common[0]:0,common[1]:1}
        for x in half:
            if x not in mapping:
                mapping[x]=len(mapping)
        result.append(tuple(mapping[x] for x in half))
    return tuple(result)


def integer_two_bridge_check(q: int, px: Pattern, py: Pattern, alpha: F=F(0)) -> dict:
    """Independent literal-colour conditional sum at integer Q."""
    def conditional(pi: Pattern, a: int, b: int) -> F:
        n=max(pi)-1
        total=F(0)
        for rest in product(range(q),repeat=n):
            colours=(a,b)+rest
            total += kernel_at_integer(q,tuple(colours[x] for x in pi),alpha)
        return total/q**n
    ex,nx=conditional(px,0,0),conditional(px,0,1)
    ey,ny=conditional(py,0,0),conditional(py,0,1)
    bx=(ex+(q-1)*nx)/q
    by=(ey+(q-1)*ny)/q
    joint=(ex*ey+(q-1)*nx*ny)/q
    residual=joint-bx*by-F(q-1,q*q)*(ex-nx)*(ey-ny)
    return {'residual':str(residual),'joint':str(joint)}


def build_report() -> dict:
    local=distinguished_patterns()
    canonical=Counter(str(endpoint_signature(p)) for p in local)
    checked=0
    for pi in partitions(8):
        if len(set(pi[:4]) & set(pi[4:])) != 2:
            continue
        px,py=bridge_endpoints(pi)
        actual=pair_activation(pi)
        expected=endpoint_signature(px)*endpoint_signature(py)
        if actual != expected:
            raise AssertionError((pi,actual,expected))
        checked+=1
    wire=(0,1,2,3,0,1,2,3)
    gram=[[pair_activation(wire),F(-1,4)],[F(-1,4),F(1,2)]]
    for ax in (F(-3),F(0),F(1,2),F(1),F(5,3)):
        for ay in (F(-3),F(0),F(1,2),F(1),F(5,3)):
            actual=pair_activation(wire,ax,ay)
            expected=F(3,2)+(ax-F(1,2))*(ay-F(1,2))/2
            if actual != expected:
                raise AssertionError((ax,ay,actual,expected))
    literal=[]
    pairs=[((0,0,1,1),(0,0,1,1)),((0,1,2,3),(0,2,1,3)),((0,1,0,2),(0,1,1,2))]
    for q in (4,5,6):
        for px,py in pairs:
            row=integer_two_bridge_check(q,px,py)
            if F(row['residual']):
                raise AssertionError(row)
            literal.append({'Q':q,'x':px,'y':py,**row})
    return {
        'schema':'matching-one/p337-bridge-completion-review/v1',
        'priority':'P0 research recommendation; not a continuum claim grade',
        'science_base':'baa5d33b2f87b2868aa0cb9d3f6518c93dbf3bff',
        'kernel_base':'2ba8863f75e0ced211b7b5442e8cddbe2fbd3deb',
        'distinguished_endpoint_patterns':len(local),
        'canonical_endpoint_signature_counts':dict(sorted(canonical.items())),
        'canonical_endpoint_max_abs':str(max(abs(endpoint_signature(p)) for p in local)),
        'two_shared_exterior_patterns_checked':checked,
        'two_shared_factorization_failures':0,
        'four_wire_bivariate_checks':25,
        'four_wire_activation':'3/2+(alpha_x-1/2)*(alpha_y-1/2)/2',
        'same_completion_lower_bound':'3/2',
        'first_Q_Gram':[[str(x) for x in row] for row in gram],
        'Gram_determinant':str(gram[0][0]*gram[1][1]-gram[0][1]**2),
        'counterterm_quotient_norm':'3/2',
        'literal_colour_checks':literal,
        'endpoint_table':[{'pattern':p,'canonical_kappa':str(endpoint_signature(p)),
                           'alpha_coefficient':str(endpoint_signature(p,F(1))-endpoint_signature(p))} for p in local],
        'new_random_samples':0,'new_lattice_population_enumeration':False,
        'full_repository_suite_run':False,
        'scope':'Finite equality-tensor contractions; conditional coefficients and bounds only. No ensemble sign, exponent, CFT field, or production authorization.'
    }


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--verify',type=Path)
    args=parser.parse_args()
    report=build_report()
    if args.verify:
        if report != json.loads(args.verify.read_text(encoding='utf-8')):
            # JSON turns tuples into lists.
            if json.loads(json.dumps(report)) != json.loads(args.verify.read_text(encoding='utf-8')):
                raise SystemExit('certificate differs from exact regeneration')
        print('Exact certificate verified')
    if args.output:
        if args.output.exists():
            raise SystemExit('refusing to overwrite an existing result')
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    if not args.output and not args.verify:
        print(json.dumps(report,ensure_ascii=False,indent=2))

if __name__=='__main__':
    main()

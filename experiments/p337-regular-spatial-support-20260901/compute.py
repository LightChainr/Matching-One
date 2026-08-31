#!/usr/bin/env python3
"""Exact canonical completed-pair Q1 activation on exterior Bell8 states.

Independent-component colours are not required to be different. Enumerating
coarsenings accounts for coincident colours via falling factorials.
"""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import csv
import json
from math import factorial
from pathlib import Path


@lru_cache(None)
def partitions(n):
    if n == 0:
        return ((),)
    return tuple(p + (v,) for p in partitions(n-1)
                 for v in range(max(p, default=-1)+2))


def k4(port):
    """4*Kreg(1), with port order N,E,S,W."""
    def bracket(z):
        a,b,c,d = z
        if a == b or c == d:
            return 0
        return (int(a == c and b == d) + int(a == d and b == c)
                + int(a == c) + int(a == d) + int(b == c)
                + int(b == d) - 4)
    return bracket(port) + bracket(port[1:] + port[:1])


def activation16(pi):
    b = max(pi)+1
    result = 0
    for sigma in partitions(b):
        k = max(sigma)+1
        if k == 1:
            continue
        rho = tuple(sigma[v] for v in pi)
        result += ((-1)**(k-2) * factorial(k-2)
                   * k4(rho[:4]) * k4(rho[4:]))
    return result


def lattice_witness():
    """One prescribed 16x16 occupation; no configuration enumeration."""
    length = 16
    marks = [(3,7),(11,7)]
    upper = {(x,9) for x in range(3,13)} | {
        (3,8),(4,8),(4,7),(11,8),(12,8),(12,7)}
    lower = {(x,5) for x in range(2,12)} | {
        (2,6),(2,7),(3,6),(10,6),(10,7),(11,6)}
    occupied = upper | lower
    assert not occupied.intersection(marks)
    directions = [(0,1),(1,0),(0,-1),(-1,0)]
    def neighbours(v):
        return [((v[0]+dx)%length,(v[1]+dy)%length)
                for dx,dy in directions]
    components = {}
    for start in sorted(occupied):
        if start in components:
            continue
        cid = len(set(components.values()))
        stack = [start]; components[start] = cid
        while stack:
            v=stack.pop()
            for w in neighbours(v):
                if w in occupied and w not in components:
                    components[w]=cid; stack.append(w)
    assert len(set(components.values())) == 2
    assert all(w in occupied for v in marks for w in neighbours(v))
    labels = [components[w] for v in marks for w in neighbours(v)]
    renumber = {}
    pi = tuple(renumber.setdefault(k,len(renumber)) for k in labels)
    assert pi == (0,0,1,1,0,0,1,1)
    # No occupied edge crosses the fundamental square; every occupied cycle
    # therefore has zero lifted displacement, including cycles inside a band.
    assert all(1 <= x <= 14 and 1 <= y <= 14 for x,y in occupied)
    assert activation16(pi) == 1
    return {'periods': [[16,0],[0,16]], 'marks': marks,
            'occupied': sorted(occupied), 'occupied_components': 2,
            'occupied_rank': 0, 'port_order': 'N,E,S,W at x, then at y',
            'exterior_partition': ''.join(map(str,pi)),
            'shared': 2, 'activation': '1/16',
            'relative_double_coefficient': '(Q-1)*(Q-2)^2/(16*Q^3)',
            'scope': 'A single fixed occupation, not an annealed correlator.'}


def main(outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    for pi in partitions(8):
        shared = len(set(pi[:4]) & set(pi[4:]))
        a16 = activation16(pi)
        rows.append({'partition': ''.join(map(str, pi)),
                     'components': max(pi)+1, 'shared': shared,
                     'activation16': a16})
    assert len(rows) == 4140
    assert sum(len(partitions(max(pi)+1)) for pi in partitions(8)) == 167894
    assert all(r['activation16'] == 0 for r in rows if r['shared'] <= 1)
    assert next(r['activation16'] for r in rows
                if r['partition'] == '01230123') == 26
    with (outdir/'kernel.csv').open('w') as f:
        writer=csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    grouped = {}
    for s in range(5):
        subset=[r for r in rows if r['shared'] == s]
        hist=Counter(r['activation16'] for r in subset)
        grouped[str(s)]={'count':len(subset),
                        'histogram_activation': {str(Fraction(k,16)):v
                                                for k,v in sorted(hist.items())}}
    extreme=max(abs(r['activation16']) for r in rows)
    summary={'source_commit':'2ba8863f75e0ced211b7b5442e8cddbe2fbd3deb',
             'source_note':'notes/local-pair-two-insertion-algebra.md, equation (5)',
             'partitions':len(rows), 'coarsening_pairs':167894,
             'shared_groups':grouped, 'max_abs_activation':str(Fraction(extreme,16)),
             'max_abs_witnesses':[r for r in rows if abs(r['activation16'])==extreme],
             'interpretation':'All Bell8 partitions; not a census of realizable planar occupations.'}
    (outdir/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (outdir/'witness.json').write_text(json.dumps(lattice_witness(),indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    import sys
    main(Path(sys.argv[1]))

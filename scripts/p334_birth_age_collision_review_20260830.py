#!/usr/bin/env python3
"""Exact P334 birth-age and collision certificate on small Gaussian tori.

Python 3.10+, standard library only. Run:
    python scripts/p334_birth_age_collision_review_20260830.py --output exact.json

Enumerates occupied subsets, not full permutations. An independent enumeration
of all 10P5 ordered prefixes verifies the N10 history witness. The graph is the
nearest-neighbor square-site quotient with period columns (a,b), (-b,a).
This is a bounded verification oracle, not a replacement production engine.
"""
from __future__ import annotations

import argparse
import itertools
import json
from collections import Counter, defaultdict, deque
from fractions import Fraction
from math import comb, factorial, gcd
from pathlib import Path

Line = tuple[int, int] | None
State = tuple[int, Line]


def enumerate_states(a: int, b: int) -> tuple[int, list[State]]:
    n = a*a + b*b
    if not 5 <= n <= 20:
        raise ValueError('This exhaustive oracle requires 5 <= a*a+b*b <= 20.')

    def key(x: int, y: int) -> tuple[int, int]:
        return ((a*x+b*y) % n, (-b*x+a*y) % n)

    representatives = [(0, 0)]
    ids = {key(0, 0): 0}
    steps = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    for x, y in representatives:
        for dx, dy in steps:
            point = key(x+dx, y+dy)
            if point not in ids:
                ids[point] = len(representatives)
                representatives.append((x+dx, y+dy))
    assert len(representatives) == n
    neighbors = [[(ids[key(x+dx, y+dy)], dx, dy) for dx, dy in steps]
                 for x, y in representatives]
    states: list[State] = []
    for mask in range(1 << n):
        lifted: dict[int, tuple[int, int]] = {}
        first_cycle = None
        rank, line = 0, None
        for start in range(n):
            if not ((mask >> start) & 1) or start in lifted:
                continue
            lifted[start] = (0, 0)
            queue = deque([start])
            while queue:
                node = queue.popleft()
                x, y = lifted[node]
                for other, dx, dy in neighbors[node]:
                    if not ((mask >> other) & 1):
                        continue
                    if other not in lifted:
                        lifted[other] = (x+dx, y+dy)
                        queue.append(other)
                    else:
                        vx = x+dx-lifted[other][0]
                        vy = y+dy-lifted[other][1]
                        if vx or vy:
                            if first_cycle is None:
                                first_cycle, rank = (vx, vy), 1
                            elif first_cycle[0]*vy-first_cycle[1]*vx:
                                rank = 2
        if rank == 1:
            assert first_cycle is not None
            x, y = first_cycle
            u, v = a*x+b*y, -b*x+a*y
            assert u % n == 0 and v % n == 0
            u, v = u//n, v//n
            divisor = gcd(abs(u), abs(v))
            u, v = u//divisor, v//divisor
            if u < 0 or (u == 0 and v < 0):
                u, v = -u, -v
            line = (u, v)
        states.append((rank, line))
    return n, states


def census(a: int, b: int) -> tuple[dict, list[State]]:
    n, states = enumerate_states(a, b)
    edges = {name: [0]*n for name in ('01', '02', '12')}
    state_counts = [[0, 0, 0] for _ in range(n+1)]
    line_counts: Counter = Counter()
    for mask, (rank, line) in enumerate(states):
        k = mask.bit_count()
        state_counts[k][rank] += 1
        if line is not None:
            line_counts[line] += 1
        for vertex in range(n):
            if (mask >> vertex) & 1:
                continue
            after = states[mask | (1 << vertex)][0]
            assert after >= rank
            name = f'{rank}{after}'
            if name in edges:
                edges[name][k] += 1
    masses = {name: sum((Fraction(c, n*comb(n-1, k))
                        for k, c in enumerate(counts)), Fraction())
              for name, counts in edges.items()}
    d = masses['02']
    assert masses['01'] == masses['12'] == 1-d
    for k, counts in enumerate(state_counts):
        assert sum(counts) == comb(n, k)
    return {
        'a': a, 'b': b, 'N': n, 'edges_by_k': edges,
        'direct_edges': sum(edges['02']), 'collision_probability': str(d),
        'collision_probability_float': float(d),
        'integrated_flux': {name: str(value) for name, value in masses.items()},
        'expected_quadratic_variation': str(2+2*d),
        'states_by_k': state_counts,
        'line_states': [{'line': list(line), 'count': count}
                        for line, count in sorted(line_counts.items())],
    }, states


def memory_witness(n: int, states: list[State]) -> dict:
    assert n == 10
    prefixes: dict[int, dict[int, int]] = {}
    totals: Counter = Counter()
    exits: Counter = Counter()
    for mask, (rank, line) in enumerate(states):
        k = mask.bit_count()
        if rank != 1:
            continue
        counts: dict[int, int] = defaultdict(int)
        for vertex in range(n):
            if not ((mask >> vertex) & 1):
                continue
            previous = mask ^ (1 << vertex)
            previous_rank, previous_line = states[previous]
            if previous_rank == 0:
                counts[k] += factorial(k-1)
            elif previous_rank == 1:
                assert previous_line == line
                for birth, count in prefixes[previous].items():
                    counts[birth] += count
        assert sum(counts.values()) == factorial(k)
        prefixes[mask] = counts
        if k == 5 and line == (1, 0):
            exit_count = sum(states[mask | (1 << v)][0] == 2
                             for v in range(n) if not ((mask >> v) & 1))
            for birth, count in counts.items():
                totals[birth] += count
                exits[birth] += count*exit_count

    # Verification independent of the prefix dynamic-programming recurrence.
    brute_totals: Counter = Counter()
    brute_exits: Counter = Counter()
    for prefix in itertools.permutations(range(n), 5):
        mask, birth = 0, None
        for k, vertex in enumerate(prefix, start=1):
            mask |= 1 << vertex
            if birth is None and states[mask][0] >= 1:
                birth = k
        if states[mask] == (1, (1, 0)):
            brute_totals[birth] += 1
            brute_exits[birth] += sum(states[mask | (1 << v)][0] == 2
                                      for v in range(n) if not ((mask >> v) & 1))
    assert totals == brute_totals == {4: 1440, 5: 4560}
    assert exits == brute_exits == {4: 1920, 5: 6480}
    hazards = {birth: Fraction(exits[birth], 5*totals[birth]) for birth in totals}
    assert hazards[4] == Fraction(4, 15)
    assert hazards[5] == Fraction(27, 95)
    assert hazards[5]-hazards[4] == Fraction(1, 57)
    assert sum((Fraction(totals[j], 6000)*hazards[j] for j in totals), Fraction()) == Fraction(7, 25)
    return {
        'N': n, 'k': 5, 'line': [1, 0],
        'rows': [{'first_birth': j, 'prefix_count': totals[j], 'exit_weight': exits[j],
                  'hazard': str(hazards[j])} for j in sorted(totals)],
        'hazard_difference': '1/57', 'uniform_layer_hazard': '7/25',
        'brute_force_prefixes_checked': 30240,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    expected = {(2, 1): '0', (3, 0): '3/35', (3, 1): '5/63',
                (3, 2): '304/3465', (4, 0): '2809/45045', (4, 1): '3511/60060'}
    rows = []
    witness = None
    for (a, b), fraction in expected.items():
        row, states = census(a, b)
        assert row['collision_probability'] == fraction
        if (a, b) == (3, 1):
            witness = memory_witness(row['N'], states)
        if (a, b) == (4, 1):
            assert row['direct_edges'] == 8823
            lc = {tuple(item['line']): item['count'] for item in row['line_states']}
            assert lc[1, 0]+lc[0, 1] == 36516
            assert lc[1, 1]+lc[1, -1] == 2380
        rows.append(row)
    result = {'schema': 'p334-birth-age-collision-v1', 'status': 'exact_finite_volume',
              'geometries': rows, 'memory_witness': witness,
              'claim_boundary': 'No asymptotic memory law, six-arm lemma, or CFT field identification is proved.'}
    text = json.dumps(result, indent=2)+'\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding='utf-8')
    print(text, end='')


if __name__ == '__main__':
    main()

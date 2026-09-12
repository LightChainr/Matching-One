#!/usr/bin/env python3
"""Exact rank-preserving row composition with a pinned torus seam.

Boundary labels 0..w-1 pin row zero; w..2w-1 are the live row. Gains
count crossings of the horizontal cut, NOT Euclidean displacements. All
arithmetic used by transitions and rank tests is integer. This is a small
reference algebra, not an optimized pTL implementation.
"""
from __future__ import annotations
from functools import lru_cache
from math import comb
from typing import NamedTuple


class State(NamedTuple):
    horizontal: bool
    entries: tuple[tuple[int, int], ...]  # canonical component root, relative gain


class GainDSU:
    def __init__(self, size: int, horizontal: bool = False):
        self.parent = list(range(size))
        self.gain = [0] * size
        self.horizontal = horizontal

    def find(self, vertex: int) -> tuple[int, int]:
        p = self.parent[vertex]
        if p == vertex:
            return vertex, 0
        root, gain = self.find(p)
        self.gain[vertex] += gain
        self.parent[vertex] = root
        return root, self.gain[vertex]

    def edge(self, u: int, v: int, gain: int) -> None:
        ru, du = self.find(u)
        rv, dv = self.find(v)
        if ru == rv:
            self.horizontal |= du + gain != dv
        else:
            self.parent[rv] = ru
            self.gain[rv] = du + gain - dv


def _validate(width: int, mask: int) -> None:
    if not isinstance(width, int) or width < 2:
        raise ValueError('integer width >= 2 required')
    if not isinstance(mask, int) or not 0 <= mask < (1 << width):
        raise ValueError('row mask outside width')


def _row(dsu: GainDSU, offset: int, width: int, mask: int) -> None:
    # Distinct periodic edges survive w=2: 0->1 with gain 0; 1->0 with gain 1.
    for x in range(width):
        y = (x + 1) % width
        if mask >> x & 1 and mask >> y & 1:
            dsu.edge(offset+x, offset+y, int(x == width-1))


def _canonical(dsu: GainDSU, keep: list[tuple[int, int]], width: int,
               shear: bool) -> State:
    entries = [(-1, 0)] * (2*width)
    groups: dict[int, list[tuple[int, int]]] = {}
    for old, new in keep:
        r, d = dsu.find(old)
        groups.setdefault(r, []).append((new, d))
    for group in groups.values():
        root, d0 = min(group)
        for v, d in group:
            entries[v] = (root, 0 if dsu.horizontal else d-d0)
    if shear and not dsu.horizontal:
        # A global Dehn shear: add one common integer to ALL live endpoints.
        # This changes only the horizontal part of eventual vertical cycles,
        # via (hx,hy)->(hx+c*hy,hy), so total rational rank is preserved.
        mixed = [(j, r, d) for j, (r, d) in enumerate(entries)
                 if j >= width and 0 <= r < width]
        if mixed:
            _, _, d = min(mixed)
            c = -d
            entries = [(r, g + c*(int(j >= width)-int(r >= width)))
                       if r >= 0 else (r, g)
                       for j, (r, g) in enumerate(entries)]
    return State(dsu.horizontal, tuple(entries))


def initial_state(width: int, mask: int, *, shear: bool = True) -> State:
    _validate(width, mask)
    dsu = GainDSU(2*width)
    keep = []
    for x in range(width):
        if mask >> x & 1:
            dsu.edge(x, width+x, 0)  # virtual alias; no site is counted twice
            keep.extend(((x, x), (width+x, width+x)))
    _row(dsu, width, width, mask)
    return _canonical(dsu, keep, width, shear)


@lru_cache(maxsize=None)
def advance(state: State, mask: int, shear: bool = True) -> State:
    width = len(state.entries)//2
    _validate(width, mask)
    dsu = GainDSU(3*width, state.horizontal)
    for v, (root, gain) in enumerate(state.entries):
        if root >= 0 and v != root:
            dsu.edge(root, v, gain)
    keep = [(x, x) for x in range(width) if state.entries[x][0] >= 0]
    for x in range(width):
        if mask >> x & 1:
            keep.append((2*width+x, width+x))
            if state.entries[width+x][0] >= 0:
                dsu.edge(width+x, 2*width+x, 0)
    _row(dsu, 2*width, width, mask)
    return _canonical(dsu, keep, width, shear)


def close_rank(state: State) -> int:
    """Close live row onto pinned first row by edges with vertical cut gain 1."""
    width = len(state.entries)//2
    adjacency: dict[int, list[tuple[int, int, int]]] = {
        i: [] for i, (r, _) in enumerate(state.entries) if r >= 0}
    def edge(u, v, x, y):
        adjacency[u].append((v, x, y))
        adjacency[v].append((u, -x, -y))
    for v, (root, gain) in enumerate(state.entries):
        if root >= 0 and v != root:
            edge(root, v, gain, 0)
    for x in range(width):
        if x in adjacency and width+x in adjacency:
            edge(width+x, x, 0, 1)
    first = (1, 0) if state.horizontal else None
    positions = {}
    for root in adjacency:
        if root in positions:
            continue
        positions[root] = (0, 0)
        stack = [root]
        while stack:
            u = stack.pop()
            xu, yu = positions[u]
            for v, x, y in adjacency[u]:
                proposed = (xu+x, yu+y)
                if v not in positions:
                    positions[v] = proposed
                    stack.append(v)
                else:
                    dx, dy = proposed[0]-positions[v][0], proposed[1]-positions[v][1]
                    if dx or dy:
                        if first is None:
                            first = (dx, dy)
                        elif first[0]*dy-first[1]*dx:
                            return 2
    return int(first is not None)


def state_for_rows(width: int, rows: tuple[int, ...] | list[int], *,
                   shear: bool = True) -> State:
    if not rows:
        raise ValueError('at least one row required')
    state = initial_state(width, rows[0], shear=shear)
    for mask in rows[1:]:
        state = advance(state, mask, shear)
    return state


def torus_rank(width: int, rows, *, shear: bool = True) -> int:
    if len(rows) < 2:
        raise ValueError('length >= 2 required to retain both interfaces')
    return close_rank(state_for_rows(width, rows, shear=shear))


def reachable_states(width: int, *, max_states: int = 100000, shear: bool = True):
    """Exhaust all row transitions; raises rather than calling a cutoff complete."""
    initial = [initial_state(width, mask, shear=shear) for mask in range(1<<width)]
    states = list(dict.fromkeys(initial))
    index = {s: i for i, s in enumerate(states)}
    transitions = []
    cursor = 0
    while cursor < len(states):
        row = []
        state = states[cursor]
        for mask in range(1<<width):
            new = advance(state, mask, shear)
            if new not in index:
                if len(states) >= max_states:
                    raise RuntimeError('state cap reached; closure NOT established')
                index[new] = len(states)
                states.append(new)
            row.append(index[new])
        transitions.append(row)
        cursor += 1
    return states, transitions, [index[s] for s in initial]


def rank_count_polynomials(width: int, length: int, *, shear: bool = True):
    """Exact occupation-count distributions, no p-grid and no random sampling."""
    if length < 2:
        raise ValueError('length >= 2 required')
    current = {}
    for mask in range(1<<width):
        s = initial_state(width, mask, shear=shear)
        a = current.setdefault(s, [0]*(width+1))
        a[mask.bit_count()] += 1
    state_counts = [len(current)]
    for height in range(2, length+1):
        nxt = {}
        for state, coefficients in current.items():
            for mask in range(1<<width):
                new = advance(state, mask, shear)
                output = nxt.setdefault(new, [0]*(width*height+1))
                k = mask.bit_count()
                for j, value in enumerate(coefficients):
                    if value:
                        output[j+k] += value
        current = nxt
        state_counts.append(len(current))
    counts = [[0]*(width*length+1) for _ in range(3)]
    for state, coefficients in current.items():
        output = counts[close_rank(state)]
        for k, c in enumerate(coefficients):
            output[k] += c
    for k in range(width*length+1):
        if sum(row[k] for row in counts) != comb(width*length, k):
            raise AssertionError('occupation-count conservation failed')
    return counts, state_counts




def minimize_rank_machine(states, transitions):
    """Continuation-equivalence partition; no claim about linear/Hankel order."""
    outputs = [close_rank(s) for s in states]
    labels = outputs[:]
    refinement_sizes = [len(set(labels))]
    while True:
        index = {}
        new = []
        for i, row in enumerate(transitions):
            key = (labels[i], tuple(labels[j] for j in row))
            if key not in index:
                index[key] = len(index)
            new.append(index[key])
        refinement_sizes.append(len(index))
        if len(index) == len(set(labels)):
            labels = new
            break
        labels = new
    representatives = [labels.index(i) for i in range(len(set(labels)))]
    quotient = [[labels[j] for j in transitions[i]] for i in representatives]
    rank = [outputs[i] for i in representatives]
    # Congruence, checked for EVERY state, not only the representative.
    for i, row in enumerate(transitions):
        if outputs[i] != rank[labels[i]]:
            raise AssertionError('output not constant on a quotient class')
        if [labels[j] for j in row] != quotient[labels[i]]:
            raise AssertionError('quotient is not a transition congruence')
    return labels, representatives, quotient, rank, refinement_sizes


if __name__ == '__main__':
    import argparse, json, time
    ap = argparse.ArgumentParser()
    ap.add_argument('--width', type=int, default=4)
    ap.add_argument('--length', type=int, default=4)
    ap.add_argument('--closure', action='store_true')
    args = ap.parse_args()
    t = time.perf_counter()
    if args.closure:
        states, transitions, initial = reachable_states(args.width)
        report = {'states': len(states), 'transitions': sum(map(len, transitions)),
                  'initial_states': initial,
                  'max_abs_gain': max(abs(g) for s in states for _,g in s.entries)}
    else:
        counts, sizes = rank_count_polynomials(args.width, args.length)
        report = {'width': args.width, 'length': args.length, 'rank_by_k': counts,
                  'matching_bernstein_counts': [c-a for a,c in zip(counts[0],counts[2])],
                  'rank_totals': list(map(sum,counts)), 'prefix_state_counts': sizes}
    report['elapsed_seconds_diagnostic'] = time.perf_counter()-t
    print(json.dumps(report, indent=2))

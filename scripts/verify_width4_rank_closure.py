#!/usr/bin/env python3
"""Build and independently check the bounded width-four closure certificate.

The graph oracle traverses physical lifted edges on the full torus, without
using a boundary summary, row compatibility, or gauge normalization. Exact
artifacts contain integers/booleans/strings only; timing is a separate log.
"""
from __future__ import annotations
import argparse
from collections import Counter, deque
from itertools import product
import json
from pathlib import Path
import time
from lifted_boundary_rank import (State, initial_state, advance, close_rank,
    state_for_rows, torus_rank, reachable_states, rank_count_polynomials,
    minimize_rank_machine)


def graph_rank(width, rows):
    height = len(rows)
    positions = {}
    first = None
    for y in range(height):
        for x in range(width):
            v = y*width+x
            if not (rows[y] >> x & 1) or v in positions:
                continue
            positions[v] = (0, 0)
            stack = [v]
            while stack:
                a = stack.pop()
                ax, ay = a%width, a//width
                px, py = positions[a]
                for dx, dy in ((-1,0), (1,0), (0,-1), (0,1)):
                    bx, by = (ax+dx)%width, (ay+dy)%height
                    if not (rows[by] >> bx & 1):
                        continue
                    b = by*width+bx
                    q = (px+dx, py+dy)
                    if b not in positions:
                        positions[b] = q
                        stack.append(b)
                    else:
                        gx, gy = q[0]-positions[b][0], q[1]-positions[b][1]
                        if gx%width or gy%height:
                            raise AssertionError('cycle displacement not a period')
                        if gx or gy:
                            if first is None:
                                first = (gx, gy)
                            elif first[0]*gy-first[1]*gx:
                                return 2
    return int(first is not None)


def independent_advance(state, mask):
    """Graph traversal of a boundary star, independent of weighted union/find."""
    w = len(state.entries)//2
    active = {i for i,(r,_) in enumerate(state.entries) if r>=0}
    active.update(2*w+x for x in range(w) if mask>>x & 1)
    adj = {i:[] for i in active}
    def edge(u,v,g):
        adj[u].append((v,g));adj[v].append((u,-g))
    for v,(r,g) in enumerate(state.entries):
        if r>=0 and v!=r:edge(r,v,g)
    for x in range(w):
        if mask>>x & 1:
            if state.entries[w+x][0]>=0:edge(w+x,2*w+x,0)
            y=(x+1)%w
            if mask>>y & 1:edge(2*w+x,2*w+y,int(x==w-1))
    seen={};groups=[];horizontal=state.horizontal
    for root in sorted(active):
        if root in seen:continue
        seen[root]=0;stack=[root];group=[]
        while stack:
            u=stack.pop()
            if u<w:group.append((u,seen[u]))
            elif u>=2*w:group.append((u-w,seen[u]))
            for v,g in adj[u]:
                proposed=seen[u]+g
                if v not in seen:seen[v]=proposed;stack.append(v)
                elif proposed!=seen[v]:horizontal=True
        if group:groups.append(group)
    entries=[(-1,0)]*(2*w)
    for group in groups:
        r,d0=min(group)
        for v,d in group:entries[v]=(r,0 if horizontal else d-d0)
    mixed=[(j,d) for j,(r,d) in enumerate(entries) if j>=w and 0<=r<w]
    if mixed and not horizontal:
        c=-min(mixed)[1]
        entries=[(r,d+c*(int(j>=w)-int(r>=w))) if r>=0 else (r,d)
                 for j,(r,d) in enumerate(entries)]
    return State(horizontal,tuple(entries))


def minimality_signature(q, rank, state):
    """All suffixes of lengths 0,1,2 in fixed lexicographic order."""
    return (rank[state], *[rank[x] for x in q[state]],
            *[rank[q[x][b]] for x in q[state] for b in range(len(q[0]))])


def certificate():
    states, transitions, initial = reachable_states(4)
    labels, reps, quotient, outputs, refinements = minimize_rank_machine(states, transitions)
    for i,row in enumerate(transitions):
        for mask,target in enumerate(row):
            if independent_advance(states[i],mask) != states[target]:
                raise AssertionError('DFS and DSU disagree on a certified transition')
    # Independently express separation without reusing the refinement step.
    signatures = [minimality_signature(quotient, outputs, s) for s in range(len(quotient))]
    if len(set(signatures)) != len(quotient):
        raise AssertionError('the announced two-row distinguishing depth fails')
    # Every state must occur at a physically valid length >=2, not only in
    # the convenient aliased first-row initializer.
    reached = {transitions[s][mask] for s in initial for mask in range(16)}
    todo = deque(reached)
    while todo:
        s = todo.popleft()
        for t in transitions[s]:
            if t not in reached:
                reached.add(t); todo.append(t)
    if reached != set(range(len(states))):
        raise AssertionError('one state is not reachable at length >=2')
    return {
        'schema': 'matching-one.width4-rank-closure.v1',
        'scope': 'square-site NN, periods (4,0),(0,m), m>=2; rank only',
        'state_order': 'BFS: initial masks 0..15, successors 0..15',
        'gain_convention': 'oriented horizontal cut crossing; internal vertical gain zero; closure seam gain (0,1)',
        'states': [[int(s.horizontal), [list(x) for x in s.entries]] for s in states],
        'transitions': transitions,
        'initial': initial,
        'rank_output': [close_rank(s) for s in states],
        'quotient_map': labels,
        'quotient_representatives': reps,
        'quotient_transitions': quotient,
        'quotient_rank_output': outputs,
        'quotient_initial': [labels[s] for s in initial],
        'refinement_counts': refinements,
        'distinguishing_suffix_length_bound': 2,
        'all_states_reachable_after_at_least_two_rows': True,
        'all_transitions_checked_by_independent_graph_traversal': True,
    }


def report():
    start = time.perf_counter()
    cert = certificate()
    print('complete closure',len(cert['states']),'quotient',len(cert['quotient_transitions']), flush=True)
    checks = {}
    for width, height in ((2,2),(2,3),(3,2),(3,3),(3,4),(4,2),(4,3),(4,4)):
        count = 0
        tally = [[0]*(width*height+1) for _ in range(3)]
        for rows in product(range(1<<width), repeat=height):
            expected = graph_rank(width, rows)
            got = torus_rank(width, rows)
            raw = torus_rank(width, rows, shear=False)
            if expected != got or expected != raw:
                raise AssertionError((width, rows, expected, got, raw))
            tally[expected][sum(x.bit_count() for x in rows)] += 1
            count += 1
        weighted, prefixes = rank_count_polynomials(width,height)
        if weighted != tally:
            raise AssertionError('per-occupation polynomial mismatch')
        checks[f'{width}x{height}'] = {'configurations':count,
            'rank_totals':list(map(sum,tally)), 'rank_by_occupation':tally,
            'oracle_vs_raw_vs_shear':True, 'all_coefficients_match':True}
        print('checked',width,height,count,flush=True)
    extra = {}
    for m in (5,6):
        counts, sizes = rank_count_polynomials(4,m)
        extra[str(m)] = {'rank_by_occupation':counts, 'rank_totals':list(map(sum,counts)),
                        'prefix_state_counts':sizes, 'new_enumeration':False}
    words = [[13,5,13,0],[7,5,13,0],[1,5,4,5],[11,14],
             [0,7,5,7],[0,13,5,7]]
    witnesses = [{'rows':w,'physical_rank':graph_rank(4,w),
                  'compressed_rank':torus_rank(4,w)} for w in words]
    summary = {'schema':'matching-one.width4-rank-closure-summary.v1',
        'geometric_boundary_states':len(cert['states']),
        'mask_transitions':sum(map(len,cert['transitions'])),
        'all_transitions_independently_reconstructed':True,
        'largest_stored_abs_gain':max(abs(d) for h,e in cert['states'] for _,d in e),
        'deterministic_rank_future_classes':len(cert['quotient_transitions']),
        'refinement_counts':cert['refinement_counts'],
        'distinct_suffix_signatures_length_at_most_two':len(cert['quotient_transitions']),
        'all_class_representatives_reachable_at_physical_lengths':True,
        'independent_graph_checks':checks,
        'total_independent_graph_configurations':sum(v['configurations'] for v in checks.values()),
        'width4_additional_DP_only':extra,'compulsory_witnesses':witnesses,
        'not_claimed':['minimal linear realization','pTL intertwiner','continuum state count',
                       'width-independent bounds','new infinite-lattice critical probability'],
    }
    print('elapsed_seconds_diagnostic',time.perf_counter()-start,flush=True)
    return cert, summary


def write_new(path, data):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x', encoding='utf-8') as f:
        if 'certificate' in path.name:
            items=list(data.items()); f.write('{\n')
            for j,(key,value) in enumerate(items):
                f.write('  '+json.dumps(key)+': ')
                if isinstance(value,list) and value and isinstance(value[0],list):
                    f.write('[\n')
                    for i,row in enumerate(value):
                        f.write('    '+json.dumps(row,separators=(',',':'),allow_nan=False))
                        f.write(',\n' if i+1<len(value) else '\n')
                    f.write('  ]')
                else:
                    f.write(json.dumps(value,separators=(',',':'),ensure_ascii=False,allow_nan=False))
                f.write(',\n' if j+1<len(items) else '\n')
            f.write('}\n')
        else:
            json.dump(data,f,indent=2,ensure_ascii=False,allow_nan=False)
            f.write('\n')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--outdir',type=Path,required=True);args=ap.parse_args()
    cert, summary = report()
    write_new(args.outdir/'width4-rank-closure-certificate.json',cert)
    write_new(args.outdir/'width4-rank-closure-summary.json',summary)

#!/usr/bin/env python3
"""Bounded feasibility probe: exact tiny controls; capped N50 frontier only."""
import json
import math
from pathlib import Path
import resource
import signal
import sys
import time

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'scripts'))
from integer_period_torus import integer_torus_geometry, classify_configuration
from planar_transition_table import build_width_table

resource.setrlimit(resource.RLIMIT_CPU, (55, 55))
START = time.process_time()
RSS_LIMIT_MIB = 500
PROGRESS = {}


class ProbeBudgetExceeded(RuntimeError):
    pass


def rss_mib():
    raw = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return raw / (1024**2 if sys.platform == 'darwin' else 1024)


def check_budget(signum, frame):
    if rss_mib() > RSS_LIMIT_MIB:
        raise ProbeBudgetExceeded('500 MiB process RSS stop')
    if time.process_time()-START > 45:
        raise ProbeBudgetExceeded('45 second shared CPU stop')


def prepare(a, b):
    geo = integer_torus_geometry(((a, -b), (b, a)), name=f'{a},{b}')
    n = geo.n
    # One fixed, deterministic repository vertex order; no physical parameter scan.
    previous = [[] for _ in range(n)]
    adjacency = [set() for _ in range(n)]
    primal = []
    for edge in geo.primal_edges:
        u, v = sorted((edge.i, edge.j))
        assert u != v
        primal.append((u, v))
    pset = set(primal)
    assert len(pset) == 2*n, 'probe scope: simple NN graph'
    matching = []
    for edge in geo.matching_edges:
        u, v = sorted((edge.i, edge.j))
        assert u != v
        matching.append((u, v))
    assert len(set(matching)) == 4*n, 'probe scope: simple matching graph'
    for u, v in matching:
        previous[v].append((u, (u, v) in pset))
        adjacency[u].add(v)
        adjacency[v].add(u)
    faces = []
    closing = [[] for _ in range(n)]
    for x, y in geo.coordinates:
        face = tuple(geo.vertex((xx, yy)) for xx, yy in
                     ((x, y), (x+1, y), (x, y+1), (x+1, y+1)))
        assert len(set(face)) == 4
        faces.append(face)
        closing[max(face)].append(face)
    last = [max(adjacency[v] | {v}) for v in range(n)]
    frontiers = [tuple(u for u in range(v) if last[u] >= v) for v in range(n+1)]
    for v in range(n):
        allowed = set(frontiers[v]) | {v}
        assert all(u in allowed for u, _ in previous[v])
        assert all(all(u in allowed for u in face) for face in closing[v])
    return geo, previous, faces, closing, frontiers


def canonical(labels):
    positive, negative = {}, {}
    out = []
    for label in labels:
        remap = positive if label > 0 else negative
        if label not in remap:
            remap[label] = (len(remap)+1)*(1 if label > 0 else -1)
        out.append(remap[label])
    return tuple(out)


def frontier(a, b, cap):
    start = time.process_time()
    geo, previous, faces, closing, fs = prepare(a, b)
    n = geo.n
    # key=(signed canonical connectivity labels, K, partial q)
    # value=(count,sum partial Sstar).  Sstar starts with its constant 2N.
    states = {((), 0, 0): (1, 2*n)}
    rows = []
    PROGRESS['current_frontier'] = {'geometry': [a, b], 'N': n, 'layers': rows}
    stopped = None
    for v in range(n):
        vertices = fs[v] + (v,)
        positions = {u: i for i, u in enumerate(vertices)}
        retained = [positions[u] for u in fs[v+1]]
        edge_positions = [(positions[u], nn) for u, nn in previous[v]]
        face_positions = [tuple(positions[u] for u in face) for face in closing[v]]
        following = {}
        for (oldlabels, k, partialq), (count, sum_s) in states.items():
            for occupied in (0, 1):
                newlabel = (max((x for x in oldlabels if x > 0), default=0)+1
                            if occupied else min((x for x in oldlabels if x < 0), default=0)-1)
                labels = list(oldlabels) + [newlabel]
                dq, ds = (0, -3) if occupied else (-1, 1)
                for pos, nearest in edge_positions:
                    same_colour = (labels[pos] > 0) == bool(occupied)
                    if not same_colour or (occupied and not nearest):
                        continue
                    if occupied:
                        dq += 1  # occupied NN edge contributes to Euler term
                        ds += 1
                    first, second = labels[-1], labels[pos]
                    if first != second:
                        labels = [first if x == second else x for x in labels]
                        dq += -1 if occupied else 1
                        ds -= 1
                if occupied:
                    complete_faces = sum(all(labels[p] > 0 for p in face) for face in face_positions)
                    dq -= complete_faces
                    ds += complete_faces
                target = (canonical(labels[p] for p in retained), k+occupied, partialq+dq)
                c0, s0 = following.get(target, (0, 0))
                following[target] = (c0+count, s0+sum_s+ds*count)
        states = following
        rows.append({'processed_vertices': v+1, 'frontier_vertices': len(fs[v+1]),
                     'states_including_K_and_partial_q': len(states),
                     'represented_prefix_assignments': sum(c for c, _ in states.values()),
                     'partial_q_min': min(key[2] for key in states),
                     'partial_q_max': max(key[2] for key in states),
                     'cpu_seconds': time.process_time()-start,
                     'cumulative_process_peak_rss_mib': rss_mib()})
        assert rows[-1]['represented_prefix_assignments'] == 2**(v+1)
        if len(states) > cap:
            stopped = 'declared_state_cap_exceeded_after_complete_layer'
            break
        if time.process_time()-START > 45:
            stopped = 'shared_45_second_soft_CPU_budget_reached'
            break
    histogram = None
    if stopped is None:
        assert all(not labels and q in (-1, 0, 1) for labels, k, q in states)
        histogram = {(k, q): value for (labels, k, q), value in states.items()}
        for k in range(n+1):
            assert sum(c for (kk, q), (c, s) in histogram.items() if kk == k) == math.comb(n, k)
    return {'geometry': [a, b], 'N': n, 'status': stopped or 'complete',
            'state_cap': cap, 'fixed_order': 'integer_period_torus quotient-key sorted representatives',
            'max_frontier_vertices': max(map(len, fs)),
            'frontier_widths_all_layers': list(map(len, fs)),
            'max_states': max(x['states_including_K_and_partial_q'] for x in rows),
            'layers': rows, 'cpu_seconds': time.process_time()-start,
            'cumulative_process_peak_rss_mib': rss_mib()}, histogram


def existing_oracle(a, b):
    start = time.process_time()
    geo, previous, faces, closing, fs = prepare(a, b)
    n = geo.n
    hist = {}
    for mask in range(1 << n):
        active = [bool(mask & (1 << i)) for i in range(n)]
        channel, black_components = classify_configuration(geo, active)
        white_channel, white_components = classify_configuration(geo, [not v for v in active], matching=True)
        k = sum(active)
        edges = sum(active[e.i] and active[e.j] for e in geo.primal_edges)
        full_faces = sum(all(active[u] for u in face) for face in faces)
        q = channel.max_rank-1
        euler_q = len(black_components)-len(white_components)-k+edges-full_faces
        assert q == euler_q
        source = len(black_components)+len(white_components)+full_faces+2*n-4*k+edges
        assert source == 2*len(black_components)+2*edges-5*k-channel.max_rank+2*n+1
        old_c, old_s = hist.get((k, q), (0, 0))
        hist[(k, q)] = old_c+1, old_s+source
    return hist, {'existing_entry': 'integer_period_torus.classify_configuration',
                  'full_configurations_checked': 2**n, 'cpu_seconds': time.process_time()-start,
                  'cumulative_process_peak_rss_mib': rss_mib()}


def main():
    result = {'scope': 'feasibility only; no N50 scientific score; fixed t=0 first source moments',
              'CPU_hard_limit_seconds': 55, 'CPU_soft_limit_seconds': 45,
              'RSS_stop_mib': RSS_LIMIT_MIB, 'RSS_check_interval_seconds': 0.05,
              'tiny_controls': [], 'N50_prefix_probes': []}
    PROGRESS['result'] = result
    start = time.process_time()
    summary, serialized = build_width_table(8)
    result['existing_planar_entry'] = {'entry': 'planar_transition_table.build_width_table(8)',
        'summary': summary, 'serialized_bytes': len(serialized),
        'cpu_seconds': time.process_time()-start,
        'cumulative_process_peak_rss_mib': rss_mib(),
        'boundary': 'planar detach/join only; not weighted torus propagation'}
    del serialized
    for a, b in ((3, 0), (3, 2)):
        stage, histogram = frontier(a, b, 50000)
        reference, oracle = existing_oracle(a, b)
        assert histogram == reference
        stage['independent_reference'] = oracle
        stage['exact_K_q_count_and_sumS_equal'] = True
        stage['histogram'] = [{'k': k, 'q': q, 'count': c, 'sum_sstar': s}
                              for (k, q), (c, s) in sorted(histogram.items())]
        result['tiny_controls'].append(stage)
    for a, b in ((5, 5), (1, 7)):
        stage, histogram = frontier(a, b, 50000)
        result['N50_prefix_probes'].append(stage)
    result['total_cpu_seconds'] = time.process_time()-START
    result['process_peak_rss_mib'] = rss_mib()
    (ROOT/'probe_result.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'total_cpu_seconds': result['total_cpu_seconds'],
                      'peak_rss_mib': result['process_peak_rss_mib'],
                      'tiny': [(r['N'], r['max_states'], r['exact_K_q_count_and_sumS_equal']) for r in result['tiny_controls']],
                      'N50': [(r['geometry'], r['status'], r['max_frontier_vertices'], r['max_states'], len(r['layers'])) for r in result['N50_prefix_probes']]}, indent=2))


if __name__ == '__main__':
    signal.signal(signal.SIGALRM, check_budget)
    signal.setitimer(signal.ITIMER_REAL, 0.05, 0.05)
    try:
        main()
    except ProbeBudgetExceeded as error:
        signal.setitimer(signal.ITIMER_REAL, 0)
        PROGRESS['budget_stop'] = str(error)
        PROGRESS['total_cpu_seconds'] = time.process_time()-START
        PROGRESS['process_peak_rss_mib'] = rss_mib()
        (ROOT/'probe_partial.json').write_text(json.dumps(PROGRESS, indent=2)+'\n')
        print(json.dumps({'budget_stop': str(error), 'partial_saved': True}))

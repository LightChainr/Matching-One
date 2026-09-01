#!/usr/bin/env python3
"""Second bounded probe: move partial q into six exact mixed moments."""
import json
import math
from pathlib import Path
import resource
import signal
import time
import probe as base

ROOT = Path(__file__).resolve().parent
START = time.process_time()
resource.setrlimit(resource.RLIMIT_CPU, (40, 40))
PROGRESS = {'scope': 'resource/tiny validation only; no N50 scoring',
            'RSS_stop_mib': 400, 'shared_previous_probe_CPU_seconds': 1.644584,
            'CPU_soft_limit_seconds_this_process': 35,
            'tiny': [], 'N50': []}


def watchdog(signum, frame):
    if base.rss_mib() > 400 or time.process_time()-START > 35:
        raise base.ProbeBudgetExceeded('400 MiB RSS / 35 second CPU guard')


def translate(v, dq, ds):
    n, q, q2, s, qs, q2s = v
    return (n, q+dq*n, q2+2*dq*q+dq*dq*n, s+ds*n,
            qs+dq*s+ds*q+dq*ds*n,
            q2s+2*dq*qs+dq*dq*s+ds*q2+2*dq*ds*q+dq*dq*ds*n)


def run(a, b, cap=400000):
    start = time.process_time()
    geo, previous, faces, closing, fs = base.prepare(a, b)
    n = geo.n
    states = {((), 0): (1, 0, 0, 2*n, 0, 0)}
    rows = []
    record = {'geometry': [a, b], 'N': n, 'max_frontier_vertices': max(map(len, fs)),
              'state_cap_after_complete_layer': cap, 'layers': rows,
              'value_fields': ['count', 'sum_q', 'sum_q2', 'sum_S', 'sum_qS', 'sum_q2S'],
              'key': ['signed_canonical_boundary_partition', 'K'],
              'status': 'running'}
    PROGRESS['current'] = record
    for vertex in range(n):
        positions = {u: i for i, u in enumerate(fs[vertex]+(vertex,))}
        retained = [positions[u] for u in fs[vertex+1]]
        edges = [(positions[u], nn) for u, nn in previous[vertex]]
        face_positions = [tuple(positions[u] for u in face) for face in closing[vertex]]
        following = {}
        for (oldlabels, k), value in states.items():
            for occupied in (0, 1):
                newlabel = (max((x for x in oldlabels if x > 0), default=0)+1
                            if occupied else min((x for x in oldlabels if x < 0), default=0)-1)
                labels = list(oldlabels)+( [newlabel] )
                dq, ds = (0, -3) if occupied else (-1, 1)
                for pos, nearest in edges:
                    if (labels[pos] > 0) != bool(occupied) or (occupied and not nearest):
                        continue
                    if occupied:
                        dq += 1
                        ds += 1
                    first, second = labels[-1], labels[pos]
                    if first != second:
                        labels = [first if x == second else x for x in labels]
                        dq += -1 if occupied else 1
                        ds -= 1
                if occupied:
                    full = sum(all(labels[p] > 0 for p in face) for face in face_positions)
                    dq -= full
                    ds += full
                target = base.canonical(labels[p] for p in retained), k+occupied
                shifted = translate(value, dq, ds)
                if target in following:
                    following[target] = tuple(x+y for x, y in zip(following[target], shifted))
                else:
                    following[target] = shifted
        states = following
        rows.append({'processed_vertices': vertex+1, 'frontier_vertices': len(fs[vertex+1]),
                     'states': len(states), 'represented_prefixes': sum(x[0] for x in states.values()),
                     'cpu_seconds': time.process_time()-start, 'peak_rss_mib': base.rss_mib()})
        assert rows[-1]['represented_prefixes'] == 2**(vertex+1)
        if len(states) > cap:
            record['status'] = 'declared_state_cap_exceeded_after_complete_layer'
            break
    else:
        record['status'] = 'complete'
    record['cpu_seconds'] = time.process_time()-start
    record['peak_rss_mib'] = base.rss_mib()
    record['max_states'] = max(row['states'] for row in rows)
    return record, states if record['status'] == 'complete' else None


def main():
    old = json.loads((ROOT/'probe_result.json').read_text())
    for a, b in ((3, 0), (3, 2)):
        row, states = run(a, b)
        reference = next(x for x in old['tiny_controls'] if x['geometry'] == [a, b])
        expected = {}
        for item in reference['histogram']:
            k, q, c, s = (item[key] for key in ('k', 'q', 'count', 'sum_sstar'))
            packet = (c, q*c, q*q*c, s, q*s, q*q*s)
            expected[k] = tuple(x+y for x, y in zip(expected.get(k, (0,)*6), packet))
        assert {k: packet for (labels, k), packet in states.items() if labels == ()} == expected
        row['all_six_final_moments_match_independent_tiny_oracle'] = True
        row['per_layer_q_key_state_counts'] = [x['states_including_K_and_partial_q'] for x in reference['layers']]
        row['moments_by_K'] = {str(k): packet for k, packet in sorted(expected.items())}
        PROGRESS['tiny'].append(row)
    for a, b in ((5, 5), (1, 7)):
        row, states = run(a, b)
        PROGRESS['N50'].append(row)
    PROGRESS.pop('current', None)


if __name__ == '__main__':
    signal.signal(signal.SIGALRM, watchdog)
    signal.setitimer(signal.ITIMER_REAL, 0.05, 0.05)
    try:
        main()
    except base.ProbeBudgetExceeded as error:
        PROGRESS['guard_stop'] = str(error)
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        PROGRESS['total_cpu_seconds_this_process'] = time.process_time()-START
        PROGRESS['process_peak_rss_mib'] = base.rss_mib()
        (ROOT/'moment_probe_result.json').write_text(json.dumps(PROGRESS, indent=2)+'\n')
        print(json.dumps({k: v for k, v in PROGRESS.items() if k not in ('tiny', 'N50', 'current')}, indent=2))

#!/usr/bin/env python3
"""Fixed finite Euler-preserving next-label policies, on the saved paired tails."""
import argparse
import gzip
import hashlib
import json
import platform
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.stats import binom

ROOT = Path(__file__).resolve().parent
P_REF = 0.59274605079
DELTA = {325: -0.7634556213017751, 425: -0.8928996539792388}
POLICIES = ['plus.-1', 'plus.+1', 'minus.-1', 'minus.+1', 'dot.plus', 'dot.minus']


def read_csv(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as f:
        names = f.readline().strip().split(',')
        x = np.loadtxt(f, delimiter=',', dtype=np.int64)
    return {s: i for i, s in enumerate(names)}, x


def policy_tables(n, census_path):
    h, x = read_csv(census_path)
    counters, ids = np.unique(x[:, h['counter']], return_inverse=True)
    assert len(counters) == 20000
    counts = np.zeros((20000, 5, 5, 4, 4), dtype=np.int16)
    idx = (ids, x[:, h['first_e']], x[:, h['second_e']],
           x[:, h['L_first']], x[:, h['L_second']])
    counts[idx] = x[:, h['count']]
    k0 = np.zeros(20000, dtype=int)
    k0[ids] = x[:, h['k0']]
    d = n-k0
    class_count = counts.sum(axis=(-2, -1), keepdims=True)
    pi = class_count / d[:, None, None, None, None]
    lf, ls = np.arange(4)[:, None], np.arange(4)[None, :]
    gs = [(lf+ls)/2, (lf-ls)/2]
    weights, dots, diagnostics = [], [], []
    outside = d-counts.sum(axis=(1, 2, 3, 4))
    for g in gs:
        mu = np.divide((counts*g).sum(axis=(-2, -1), keepdims=True), class_count,
                       out=np.zeros_like(pi), where=class_count > 0)
        dot = pi*(g-mu)
        dots.append(dot)
        for t in (-1., 1.):
            unnormalized = np.exp(t*pi*g)
            z = np.divide((counts*unnormalized).sum(axis=(-2, -1), keepdims=True),
                          class_count, out=np.ones_like(pi), where=class_count > 0)
            w = unnormalized/z
            class_error = float(np.max(np.abs((counts*w).sum(axis=(-2, -1), keepdims=True)-class_count)))
            second_moment = (outside+(counts*w*w).sum(axis=(1, 2, 3, 4)))/d
            ess_fraction = 1/second_moment
            present = w[counts > 0]
            diagnostics.append({'class_mass_error_count_units': class_error,
                                'weight_min': float(min(1, present.min())),
                                'weight_max': float(max(1, present.max())),
                                'ess_fraction_min': float(ess_fraction.min()),
                                'ess_fraction_mean': float(ess_fraction.mean()),
                                'ess_fraction_p05': float(np.quantile(ess_fraction, .05))})
            weights.append(w)
    assert max(z['class_mass_error_count_units'] for z in diagnostics) < 1e-10
    return counters, np.array(weights+dots), diagnostics, counts


def score(n, output):
    counters, tables, diagnostics, counts = policy_tables(n, ROOT/f'census/N{n}/census.csv.gz')
    batch_values, batch_labels = [], None
    source_hashes = {}
    cells = ['all', '00', '01', '02', '10', '20']
    for b in range(20):
        raw_path = ROOT/f'results/p334-nested-next-label-forks/N{n}/N{n}.batch{b:02}.csv.gz'
        mark_path = ROOT/f'results/p334-next-label-contact-coordinates/N{n}/N{n}.batch{b:02}.csv.gz'
        h, raw = read_csv(raw_path)
        c, marks = read_csv(mark_path)
        raw = raw[np.lexsort(tuple(raw[:, h[k]] for k in ('replica', 'group', 'quartet', 'counter')))].reshape(1000, 8, 2, 2, -1)
        marks = marks[np.lexsort(tuple(marks[:, c[k]] for k in ('group', 'quartet', 'counter')))].reshape(1000, 8, 2, -1)
        for k in ('counter', 'quartet', 'group', 'next_label'):
            assert np.array_equal(raw[..., 0, h[k]], marks[..., c[k]])
        ids = np.searchsorted(counters, marks[..., c['counter']])
        assert np.array_equal(counters[ids], marks[..., c['counter']])
        ranks = [raw[:, 0, 0, 0, h[f'{o}_rank']] for o in ('first', 'second')]
        cell_id = 3*ranks[0]+ranks[1]
        safe = np.ones((1000, 8, 2), dtype=bool)
        es, loops = [], []
        for o, r in zip(('first', 'second'), ranks):
            safe &= marks[..., c[f'{o}_rank_after']] == r[:, None, None]
            e = marks[..., c[f'{o}_e']]
            es.append(e)
            loops.append((r == 0)[:, None, None]*(e-marks[..., c[f'{o}_c']]))
        ix = (ids[safe], es[0][safe], es[1][safe], loops[0][safe], loops[1][safe])
        assert np.all(counts[ix] > 0), 'sampled mark missing from exact census'
        sampled = np.ones((6, 1000, 8, 2))
        sampled[4:] = 0
        sampled[:, safe] = tables[(slice(None),)+ix]
        dw = sampled[..., 0]-sampled[..., 1]
        # Fixed-time CDF values conditional on the two actual birth clocks.
        ys = []
        for endpoint in ('p_ref', 'p_integral'):
            yo = []
            for o in ('first', 'second'):
                k1, k2 = (raw[..., h[f'{o}_{k}']] for k in ('k1', 'k2'))
                if endpoint == 'p_ref':
                    f1, f2 = (binom.sf(k-1, n, P_REF) for k in (k1, k2))
                else:
                    f1, f2 = 1-k1/(n+1), 1-k2/(n+1)
                yo.append(np.stack((f1+f2-1, 1-f1+f2), axis=-1).mean(axis=3))
            for output_name, obs in [('S', (yo[0]+yo[1])/2), ('D', (yo[0]-yo[1])/DELTA[n])]:
                for j, name in enumerate(('A', 'E')):
                    ys.append((f'{endpoint}.{output_name}.{name}', obs[..., j]))
        values, labels = [], []
        for cell in cells:
            chosen = np.ones(1000, dtype=bool) if cell == 'all' else cell_id == 3*int(cell[0])+int(cell[1])
            for response, y in ys:
                dy = y[..., 0]-y[..., 1]
                # E[(w(U)-w(V))*(Y_U-Y_V)/2] = E_qt[Y]-E_uniform[Y].
                v = (dw[:, chosen]*dy[chosen]).sum(axis=(1, 2))/(2*8000)
                for policy, value in zip(POLICIES, v):
                    labels.append(f'{cell}.{response}.{policy}')
                    values.append(float(value))
                for i, mark in ((0, 'plus'), (2, 'minus')):
                    dp = v[4+i//2]
                    for label, value in [('odd', (v[i+1]-v[i])/2),
                                         ('even', (v[i+1]+v[i])/2),
                                         ('odd_minus_tangent', (v[i+1]-v[i])/2-dp),
                                         ('plus_remainder', v[i+1]-dp),
                                         ('minus_remainder', v[i]+dp)]:
                        labels.append(f'{cell}.{response}.{mark}.{label}')
                        values.append(float(value))
        if batch_labels is None:
            batch_labels = labels
        else:
            assert labels == batch_labels
        batch_values.append(values)
        for path in (raw_path, mark_path):
            source_hashes[str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
        print('scored', n, b, flush=True)
    x = np.array(batch_values)
    mean = x.mean(axis=0)
    factor = (x-mean)/np.sqrt(20*19)
    result = {'N': n, 'batch_ids': list(range(20)), 'population_per_batch': 1000,
              'labels': batch_labels, 'batch_values': x.tolist(), 'estimate': mean.tolist(),
              'se': np.linalg.norm(factor, axis=0).tolist(), 'factor': factor.tolist(),
              'policy_diagnostics_order': POLICIES[:4], 'policy_diagnostics': diagnostics,
              'input_sha256': source_hashes}
    (output/f'N{n}.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=ROOT/'output')
    args = parser.parse_args()
    args.output.mkdir(exist_ok=False)
    started = time.perf_counter()
    result = {'schema': 'p334.finite-source.v1', 'finite_t': [-1, 1], 'new_tails': 0,
              'new_independent_prefixes': 0, 'p_ref': P_REF,
              'interpretation': 'Existing-data paired importance response for the exactly normalized finite policy; not direct intervention resampling.',
              'environment': {'hostname': platform.node(), 'machine': platform.machine(),
                              'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__},
              'sizes': {}}
    for n in (325, 425):
        result['sizes'][str(n)] = score(n, args.output)
    result['elapsed_seconds'] = time.perf_counter()-started
    (args.output/'latest.json').write_text(json.dumps(result, separators=(',', ':'), allow_nan=False)+'\n')
    print('complete', result['elapsed_seconds'], flush=True)


if __name__ == '__main__':
    main()

#!/usr/bin/env python3
"""Existing64+8 conditional birth-shape response, no sampling or p scan."""
import os
for variable in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[variable] = '1'
import argparse
from concurrent.futures import ProcessPoolExecutor
import gzip
import hashlib
import json
from pathlib import Path
import platform
import time
import numpy as np
import scipy
from ustat import shape_estimators, check_algebra

ROOT = Path(__file__).resolve().parent
DELTA = {325: -0.7634556213017751, 425: -0.8928996539792388}
MODES = ('old8', 'new64', 'combined72')
SHAPES = ('VarC', 'CovCW', 'VarW')
MARKS = ('plus', 'minus')
ORIENTATIONS = ('first', 'second')
BASE_LABELS = [f'{shape}.{output}.{mark}' for shape in SHAPES
               for output in ('S', 'D') for mark in MARKS]+['shape_Frobenius_energy']


def array_csv(path):
    with gzip.open(path, 'rt') as stream:
        names = stream.readline().strip().split(',')
        values = np.loadtxt(stream, delimiter=',', dtype=np.int64)
    return {name: i for i, name in enumerate(names)}, values


def score_from_marks(e, c, after, counts, totals, d):
    p, q = e.shape[:2]
    ids = np.broadcast_to(np.arange(p)[:, None, None], (p, q, 2))
    safe = np.all(after == 0, axis=-1)
    ix = (ids, e[..., 0], e[..., 1])
    mass, means_numerator = counts[ix], totals[ix]
    assert np.all(mass[safe] > 0)
    numerator = mass[..., None]*(e-c)-means_numerator
    numerator *= safe[..., None]
    combined = np.stack((numerator.sum(-1), numerator[..., 0]-numerator[..., 1]), -1)
    return combined/(2*d[:, None, None, None])


def moments(birth, score, k0, n):
    # Birth axes P,Q,label,suffix,physical-orientation,K1/K2.
    c = ((birth[..., 0]+birth[..., 1])/2-k0[:, None, None, None, None])/(n+1)
    w = (birth[..., 1]-birth[..., 0])/(n+1)
    label = np.stack((c, w, c*c, c*w, w*w), -1).mean(axis=3)
    b = label.mean(axis=2)
    h = 0.5*(score[:, :, 0]-score[:, :, 1])[:, :, :, None, None]*(label[:, :, 0]-label[:, :, 1])[:, :, None]
    return b, h


def transform_saved_moments(values, a, baseline):
    x, y, xx, xy, yy = [values[..., i] for i in range(5)]
    a = a.reshape((len(a),)+(1,)*(x.ndim-1))
    return np.stack(((x+y)/2-a if baseline else (x+y)/2,
                     y-x, (xx+2*xy+yy)/4-a*(x+y)+(a*a if baseline else 0),
                     (yy-xx)/2-a*(y-x), xx-2*xy+yy), -1)


def batch(task):
    n, batch_id, old, input_root = task
    input_root = Path(input_root)
    header, raw = array_csv(input_root/'new64'/f'N{n}.batch{batch_id:02}.csv.gz')
    raw = raw[np.lexsort(tuple(raw[:, header[k]] for k in ('replica', 'group', 'quartet', 'counter')))]
    p = len(old['counter'])
    raw = raw.reshape(p, 64, 2, 2, -1)
    assert np.array_equal(raw[:, 0, 0, 0, header['counter']], old['counter'])
    assert np.all(raw[..., header['first_rank']] == 0) and np.all(raw[..., header['second_rank']] == 0)
    assert np.all(raw[..., header['k0']] == old['k0'][:, None, None, None])
    assert np.all(raw[:, :, 0, 0, header['quartet']] == np.arange(8, 72)[None])
    assert np.all(raw[:, :, :, 0, header['group']] == np.arange(2)[None, None])
    assert np.all(raw[..., header['replica']] == np.arange(2)[None, None, None])
    for key in ('next_label', 'first_next_rank', 'second_next_rank', 'first_e', 'first_c', 'second_e', 'second_c'):
        assert np.array_equal(raw[..., 0, header[key]], raw[..., 1, header[key]])
    birth = np.stack([np.stack([raw[..., header[f'{o}_{k}']] for k in ('k1', 'k2')], -1)
                      for o in ORIENTATIONS], -2)
    assert np.all((birth[..., 0] > old['k0'][:, None, None, None, None])
                  & (birth[..., 0] <= birth[..., 1]) & (birth[..., 1] <= n))
    es, cs, after = [np.stack([raw[:, :, :, 0, header[f'{o}_{field}']] for o in ORIENTATIONS], -1)
                      for field in ('e', 'c', 'next_rank')]
    new_score = score_from_marks(es, cs, after, old['census_class_count'], old['census_class_loop_sum'], old['d'])
    own_old_score = score_from_marks(old['contact_e'], old['contact_c'], old['rank_after'],
                                    old['census_class_count'], old['census_class_loop_sum'], old['d'])
    score_error = float(np.max(np.abs(own_old_score-old['score'])))
    assert score_error == 0
    old_b, old_h = moments(old['birth_k'].astype(np.int64), old['score'], old['k0'], n)
    new_b, new_h = moments(birth, new_score, old['k0'], n)
    # Reuse/check the existing archive; no old hierarchy is recomputed.
    old_b_expected = transform_saved_moments(old['b'], old['k0']/(n+1), True)
    old_h_expected = transform_saved_moments(old['h'], old['k0']/(n+1), False)
    moment_error = max(float(np.max(abs(old_b-old_b_expected))), float(np.max(abs(old_h-old_h_expected))))
    assert moment_error < 2e-13
    result = {'N': n, 'batch': batch_id, 'counter': old['counter'],
              'score_archive_max_error': score_error, 'old_moment_transform_max_error': moment_error}
    for mode, b, h in (('old8', old_b, old_h), ('new64', new_b, new_h),
                        ('combined72', np.concatenate((old_b, new_b), 1), np.concatenate((old_h, new_h), 1))):
        cov, energy, squares = shape_estimators(b, h)
        values = []
        for j in range(3):
            for output in ('S', 'D'):
                physical = cov[..., j]  # P,source,orientation
                v = physical.mean(-1) if output == 'S' else (physical[..., 0]-physical[..., 1])/DELTA[n]
                values.extend(v[:, s] for s in range(2))
        values.append(energy)
        result[mode] = np.column_stack(values)
        result[mode+'_physical'] = cov
        result[mode+'_shape_squared'] = squares
    return result


def write_report(result, target):
    lines = ['# P334 conditional shape: fixed cell00, existing 8+64 quartets', '',
             'All values use the original 20000-prefix denominator. Cell00 is the selected contribution; no claim about unsampled non00 shape derivatives.', '',
             '| N | mode | target | estimate | original20-batch SE | LOO minimum | LOO maximum |',
             '|---:|---|---|---:|---:|---:|---:|']
    for ns, row in result['sizes'].items():
        for mode in ('new64', 'combined72'):
            for key in BASE_LABELS:
                i = row['labels'].index(mode+'.cell00.'+key)
                loo = np.array(row['LOO'])[:, i]
                lines.append(f"| {ns} | {mode} | {key} | {row['estimate'][i]:.12g} | {row['se'][i]:.6g} | {loo.min():.12g} | {loo.max():.12g} |")
    lines += ['', 'Existing local-mean Jacobian determinants are appended to the same factor without rerunning their analysis.',
              'Shape energy is an unbiased signed finite-sample estimator of a nonnegative population target; it is not clipped.',
              'Old8, new64 and combined72 share the original prefixes and are not independent confirmations.', '']
    target.write_text('\n'.join(lines))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers', type=int, default=16)
    parser.add_argument('--output', type=Path, default=ROOT/'results')
    args = parser.parse_args()
    args.output.mkdir(exist_ok=False)
    started = time.time()
    algebra = check_algebra()
    manifest = json.loads((ROOT/'inputs/manifest.json').read_text())
    assert hashlib.sha256((ROOT/'PLAN.md').read_bytes()).hexdigest() == manifest['plan_sha256_before_scoring']
    for entry in manifest['files']:
        assert hashlib.sha256((ROOT/'inputs'/entry['local_path']).read_bytes()).hexdigest() == entry['sha256']
    keys = ('b', 'h', 'score', 'birth_k', 'contact_e', 'contact_c', 'rank_after', 'counter', 'k0', 'd', 'census_class_count', 'census_class_loop_sum')
    tasks = []
    retained_counts = {}
    for n in (325, 425):
        with np.load(ROOT/'inputs/old8'/f'N{n}.npz') as archive:
            selected = archive['rankcell'] == 0
            batches = archive['batch'][selected]
            old = {key: archive[key][selected] for key in keys}
        assert len(old['counter']) == {325: 1502, 425: 1551}[n]
        retained_counts[str(n)] = np.bincount(batches, minlength=20).tolist()
        for b in range(20):
            tasks.append((n, b, {key: value[batches == b] for key, value in old.items()}, str(ROOT/'inputs')))
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        blocks = list(pool.map(batch, tasks))
    old_mean = json.loads((ROOT/'inputs/existing_local_mean_rank.json').read_text())
    result = {'schema': 'p334.conditional-shape.v1', 'sizes': {}, 'algebra_checks': algebra,
              'population': 'Cell00 contribution to original20000 prefixes per N; 20 batches x1000 denominator. Other cells are not declared zero.',
              'modes': 'old8, existing new64 and combined72; all share original prefixes and original batch IDs.',
              'new_sampling': 0, 'cell00_counts': retained_counts,
              'target_coordinates': 'C=(K1+K2-2*k0)/(2*(N+1)), W=(K2-K1)/(N+1); the prefix-fixed centering does not affect conditional covariances.',
              'fourth_readout': 'E[1_cell00 sum_(physical geometry,plus/minus source) {HVarC^2+2*HCovCW^2+HVarW^2}]; unbiased U2/U3/U4, never directly squared or clipped.',
              'null': 'For every prefix and physical geometry, arbitrary source-dependent deterministic translations of each birth preserve all three conditional covariance entries. Energy>0 excludes this class somewhere in cell00; zero pooled means do not establish it.',
              'source_hashes': manifest}
    for n in (325, 425):
        bb = [block for block in blocks if block['N'] == n]
        labels, vectors = [], []
        for mode in MODES:
            for j, key in enumerate(BASE_LABELS):
                labels.append(mode+'.cell00.'+key)
                vectors.append(np.array([block[mode][:, j].sum()/1000 for block in bb]))
        ref = old_mean['sizes'][str(n)]
        ref_matrix = np.array(ref['joint_20_batch_means'])
        for mode in MODES:
            for endpoint in ('p_ref.A', 'p_integral.A'):
                key = mode+'.cell00.'+endpoint+'.E_det_JZ'
                labels.append('existing_mean_rank.'+key)
                vectors.append(ref_matrix[:, ref['base_labels'].index(key)])
        x = np.column_stack(vectors)
        mean = x.mean(0)
        loo = (20*mean[None]-x)/19
        factor = np.sqrt(19/20)*(loo-loo.mean(0))
        comparisons = {}
        for key in BASE_LABELS:
            differences = {}
            for mode in ('new64', 'combined72'):
                v = x[:, labels.index(mode+'.cell00.'+key)]-x[:, labels.index('old8.cell00.'+key)]
                differences[mode+'_minus_old8'] = {'estimate': float(v.mean()), 'se': float(v.std(ddof=1)/np.sqrt(20)), 'batch_values': v.tolist()}
            comparisons[key] = differences
        result['sizes'][str(n)] = {'batch_ids': list(range(20)), 'labels': labels, 'estimate': mean.tolist(),
            'se': np.linalg.norm(factor, axis=0).tolist(), 'joint_20_batch_means': x.tolist(), 'LOO': loo.tolist(), 'factor': factor.tolist(),
            'paired_mode_differences': comparisons,
            'checks': {'score_archive_max_error': max(b['score_archive_max_error'] for b in bb),
                       'old_moment_transform_max_error': max(b['old_moment_transform_max_error'] for b in bb)},
            'cell00_count': sum(retained_counts[str(n)])}
        np.savez_compressed(args.output/f'prefix_shape_N{n}.npz', counter=np.concatenate([b['counter'] for b in bb]),
            batch=np.concatenate([np.full(len(b['counter']), b['batch']) for b in bb]), labels=np.array(BASE_LABELS),
            physical_axis_order=np.array(['prefix', 'source_plus_minus', 'physical_first_second', 'VarC_CovCW_VarW']),
            **{mode: np.concatenate([b[mode] for b in bb]) for mode in MODES},
            **{mode+'_physical': np.concatenate([b[mode+'_physical'] for b in bb]) for mode in MODES},
            **{mode+'_shape_squared': np.concatenate([b[mode+'_shape_squared'] for b in bb]) for mode in MODES})
        for key in BASE_LABELS:
            j = labels.index('combined72.cell00.'+key)
            print(n, key, f'{mean[j]:.12g} +/- {np.linalg.norm(factor[:,j]):.6g}', flush=True)
    result['run_receipt'] = {'hostname': platform.node(), 'machine': platform.machine(), 'python': platform.python_version(),
        'numpy': np.__version__, 'scipy': scipy.__version__, 'workers': args.workers,
        'started_unix': started, 'finished_unix': time.time(), 'elapsed_seconds': time.time()-started,
        'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'ustat_sha256': hashlib.sha256((ROOT/'ustat.py').read_bytes()).hexdigest(), 'new_sampling': 0}
    (args.output/'latest.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    write_report(result, args.output/'REPORT.md')
    (args.output/'run_receipt.json').write_text(json.dumps(result['run_receipt'], indent=2)+'\n')
    print(json.dumps(result['run_receipt']), flush=True)


if __name__ == '__main__':
    main()

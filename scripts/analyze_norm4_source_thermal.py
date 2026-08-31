#!/usr/bin/env python3
"""Actual source tangents on both norm-4 lineages, from old NZ permutations."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
import subprocess
import time
from fractions import Fraction
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import brentq
from scipy.special import gammaln

from analyze_p40_source_thermal import direction_values

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / 'analysis/p40_source_thermal_chain_candidates.json'
HYPOTHESES = ROOT / 'analysis/norm4_source_tangent_hypotheses.json'
OUTPUT = ROOT / 'results/norm4-source-thermal'
NS = (65, 85, 130, 170, 260, 340)
LINEAGES = ((65, 130, 260), (85, 170, 340))
FIELDS = ('U', 'Udot_density', 'Udot_fugacity', 'rootdot_fugacity', 'p0', 'D',
          'Udot_direct_fugacity', 'Udot_rootmotion_fugacity',
          'Udot_slope_source_fugacity', 'Udot_slope_root_fugacity')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_raw(path, n):
    counts = np.zeros((100, 2, n + 1), dtype=np.int64)
    sums = np.zeros((100, 2, n + 1, 5), dtype=np.float64)
    seen = set()
    with path.open() as handle:
        for row in csv.DictReader(handle):
            b, k = int(row['batch']), int(row['k'])
            g = ('first', 'second').index(row['orientation'])
            key = (b, g, k)
            if key in seen or int(row['n']) != n or not 0 <= b < 100 or not 0 <= k <= n:
                raise ValueError(f'invalid or duplicate NZ row: {path} {key}')
            seen.add(key)
            counts[b, g, k] = int(row['samples'])
            sums[b, g, k] = [int(row[f]) for f in ('sum_q', 'sum_e', 'sum_s', 'sum_qs', 'sum_es')]
    if len(seen) != 100 * 2 * (n + 1) or not np.all(counts == 1000):
        raise ValueError(f'{path}: expected 100 paired batches of 1000 complete permutations')
    sums[:, :, :, 2:] /= n
    return sums


def binomial_moments(sums, samples, p, n):
    """Every permutation supplies all K; derivatives integrate Binomial weights."""
    k = np.arange(n + 1, dtype=float)
    log_b = gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)
    log_b += k * math.log(p) + (n - k) * math.log1p(-p)
    weights = np.exp(log_b)
    score = k / p - (n - k) / (1 - p)
    wp = weights * score
    wpp = weights * (score * score - k / p**2 - (n - k) / (1 - p)**2)
    z, zp, zpp = weights.sum(), wp.sum(), wpp.sum()
    mean = weights @ sums / (samples * z)
    first = (wp @ sums / samples - mean * zp) / z
    second = (wpp @ sums / samples - mean * zpp - 2 * first * zp) / z
    return mean, first, second, {'binomial_mass': float(z), 'permutations': samples}


def at_root(sums, samples, n, delta, bracket):
    def packets(p):
        return [binomial_moments(sums[g], samples, p, n) for g in range(2)]

    def root_function(p):
        return sum(packet[0][0] for packet in packets(p)) / 2

    p0 = float(brentq(root_function, *bracket, xtol=5e-14, rtol=5e-14))
    packet = packets(p0)
    rows = [direction_values(p) for p in packet]
    d = (rows[0]['q_p'] + rows[1]['q_p']) / 2
    if not math.isfinite(d) or d <= 0:
        raise ValueError(f'N{n}: nonpositive local matching slope')
    jq = (rows[0]['Jq'] + rows[1]['Jq']) / 2
    jqp = (rows[0]['Jq_p'] + rows[1]['Jq_p']) / 2
    qpp = (rows[0]['q_pp'] + rows[1]['q_pp']) / 2
    b = (rows[0]['E_p'] - rows[1]['E_p']) / delta
    je_p = (rows[0]['JE_p'] - rows[1]['JE_p']) / delta
    e_pp = (rows[0]['E_pp'] - rows[1]['E_pp']) / delta
    prefactor = n**(13 / 8) / 2
    rootdot = -jq / d
    pieces = [prefactor * je_p / d, prefactor * e_pp * rootdot / d,
              -prefactor * b * jqp / d**2, -prefactor * b * qpp * rootdot / d**2]
    udot = math.fsum(pieces)
    result = dict(zip(FIELDS, [prefactor * b / d, udot, n * udot, n * rootdot, p0, d,
                               *(n * term for term in pieces)]))
    diagnostic = {'p0': p0, 'pooled_q': float(root_function(p0)), 'D': float(d),
                  'permutations': samples, 'direction': dict(zip(('first', 'second'), rows)),
                  'integration': packet[0][3], 'delta_cos4': delta,
                  'Udot_density_terms': dict(zip(('direct', 'rootmotion', 'slope_source', 'slope_root'), pieces))}
    return {key: float(value) for key, value in result.items()}, diagnostic


def vectorize(by_n):
    values = {f'N{n}.{field}': by_n[n][field] for n in NS for field in FIELDS}
    for base, middle, end in LINEAGES:
        for model, coefficients in (('q2', (1, -3, 2)), ('Jordan', (1, -2, 1))):
            for field in ('U', 'Udot_density', 'Udot_fugacity'):
                values[f'{model}.{base}.{field}'] = sum(c * by_n[n][field]
                                                       for c, n in zip(coefficients, (base, middle, end)))
        values[f'drift_shape.{base}'] = by_n[middle]['U'] - by_n[base]['U']
    for model, factor in (('q2', .5), ('Jordan', 1.0)):
        values[f'{model}.generator_drift_determinant'] = factor * (
            values[f'{model}.65.Udot_fugacity'] * values['drift_shape.85']
            - values[f'{model}.85.Udot_fugacity'] * values['drift_shape.65'])
    return values


def joint_zero(labels, central, covariance, selected):
    indices = [labels.index(label) for label in selected]
    block = covariance[np.ix_(indices, indices)]
    means = central[indices]
    eigenvalues = np.linalg.eigvalsh(block)
    if eigenvalues[0] <= 0 or not np.isfinite(block).all():
        return {'status': 'unresolved_covariance', 'labels': selected, 'covariance': block.tolist()}
    statistic = float(means @ np.linalg.solve(block, means))
    return {'status': 'computed', 'labels': selected, 'residuals': means.tolist(),
            'covariance': block.tolist(), 'chi_square': statistic, 'df': 2,
            'nominal_p': math.exp(-statistic / 2),
            'interpretation': 'Finite archived-subset diagnostic of the explicitly stated extension, not model certification'}


def main():
    started = time.perf_counter()
    destination = OUTPUT / 'latest.json'
    if destination.exists():
        raise ValueError('saved result exists; do not overwrite or repeat')
    manifest = json.loads(MANIFEST.read_text())
    hypotheses = json.loads(HYPOTHESES.read_text())
    runs = {run['N']: run for run in manifest['runs']}
    bracket = hypotheses['root_bracket']
    raw = {n: read_raw(OUTPUT / 'raw' / f'n{n}.csv', n) for n in NS}
    totals = {n: raw[n].sum(axis=0) for n in NS}
    deltas = {n: float(Fraction(runs[n]['delta_cos4'])) for n in NS}
    central_by_n, central_diagnostics = {}, {}
    for n in NS:
        central_by_n[n], central_diagnostics[n] = at_root(totals[n], 100000, n, deltas[n], bracket)
    central_map = vectorize(central_by_n)
    labels = list(central_map)
    central = np.array(list(central_map.values()))
    covariance = np.zeros((len(labels), len(labels)))
    groups = {}
    for n in NS:
        groups.setdefault(runs[n]['dependency_group'], []).append(n)
    group_results = {}
    for group, ns in groups.items():
        vectors, roots = [], {n: [] for n in ns}
        for batch in range(100):
            omitted = dict(central_by_n)
            for n in ns:
                omitted[n], diagnostic = at_root(totals[n] - raw[n][batch], 99000, n, deltas[n], bracket)
                roots[n].append(diagnostic['p0'])
            vectors.append(list(vectorize(omitted).values()))
        vectors = np.asarray(vectors)
        deviations = vectors - vectors.mean(axis=0)
        component = 99 / 100 * deviations.T @ deviations
        covariance += component
        group_results[group] = {'Ns': ns, 'delete_one_batch_ids': list(range(100)),
                                 'delete_one_vectors': vectors.tolist(), 'covariance_contribution': component.tolist(),
                                 'root_ranges': {n: [min(roots[n]), max(roots[n])] for n in ns},
                                 'operation': 'omit the same configuration batch in every N of this group; other groups stay at their central values'}
    errors = np.sqrt(np.maximum(0, np.diag(covariance)))
    estimates = {label: {'value': float(value), 'se': float(error),
                          'z': float(value / error) if error > 0 else None}
                 for label, value, error in zip(labels, central, errors)}
    summaries = {}
    for model in ('q2', 'Jordan'):
        summaries[model] = {
            'unperturbed_recurrence': joint_zero(labels, central, covariance, [f'{model}.{base}.U' for base in (65, 85)]),
            'common_bulk_source_rigidity': joint_zero(labels, central, covariance, [f'{model}.{base}.Udot_fugacity' for base in (65, 85)]),
            'one_generator_drift': estimates[f'{model}.generator_drift_determinant']}
    result = {
        'schema': 'matching-one.norm4-source-thermal.v1', 'status': 'computed_existing_production_subset',
        'execution_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        'hypotheses': hypotheses, 'replay_manifest': manifest, 'labels': labels,
        'estimates': estimates, 'covariance': covariance.tolist(), 'by_N': central_diagnostics,
        'dependency_groups': group_results, 'model_extension_diagnostics': summaries,
        'source_coordinates': {'density': 'S=(CB+CW)/N', 'common_bulk_fugacity': 's=CB+CW; v_N=N*Udot_density'},
        'estimator': 'Every archived permutation supplies all occupation prefixes. Integrate q,E,S,qS,ES with exact Binomial weights; analytic first and second p derivatives; refind the pooled matching root inside every delete-one. The K prefixes are not independent samples.',
        'finite_subset_scope': '100000 old permutations per N, drawn as the declared initial contiguous production subset; these source observations are new, the random configurations are old. Original full-production errors and model-exclusion claims do not transfer to this subset.',
        'source_extension_scope': 'Rigidity and one-generator drift are explicit hypotheses under the same paired-cluster fugacity, not consequences of the old t=0 recurrence. Check the reported baseline jointly when interpreting them; a compatible tangent does not rescue an incompatible baseline.',
        'inputs': [{'path': str((OUTPUT / 'raw' / f'n{n}.csv').relative_to(ROOT)), 'sha256': sha(OUTPUT / 'raw' / f'n{n}.csv')} for n in NS],
        'code': [{'path': str(path.relative_to(ROOT)), 'sha256': sha(path)} for path in
                 (Path(__file__), ROOT / 'scripts/analyze_p40_source_thermal.py', MANIFEST, HYPOTHESES)],
        'environment': {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__, 'machine': platform.machine()},
        'elapsed_seconds': time.perf_counter() - started, 'new_samples': 0, 'server_actions': 0, 'test_suites': []}
    with destination.open('x') as handle:
        handle.write(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({'elapsed_seconds': result['elapsed_seconds'], 'model_extension_diagnostics': summaries,
                       'source_response': {n: estimates[f'N{n}.Udot_fugacity'] for n in NS}}, indent=2))


if __name__ == '__main__':
    main()

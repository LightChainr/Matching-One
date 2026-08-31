#!/usr/bin/env python3
"""Use the full-production complement to estimate one existing cluster source."""
from __future__ import annotations

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

import analyze_norm4_source_thermal as old
from norm4_source_two_phase_core import baseline, estimate
from norm4_source_two_phase_inputs import load_complement

ROOT = old.ROOT
CONTRACT = ROOT / 'analysis/norm4_source_two_phase_contract.json'
DESTINATION = ROOT / 'results/norm4-source-two-phase/latest.json'
NS = old.NS
VIEWS = ('raw_bulk', 'anchor_only_bulk', 'two_phase_bulk')
FIELDS = (*VIEWS, 'topology_bulk', 'residual_bulk', 'U_source', 'U_anchor', 'p0_anchor')
MODELS = {'q2': (1, -3, 2), 'Jordan': (1, -2, 1)}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def point_record(n, fitted, anchor, raw_state):
    return {**{field: fitted[field] for field in FIELDS if field in fitted},
            'raw_bulk': raw_state['raw_bulk'], 'U_source': raw_state['U_source'],
            'U_anchor': anchor['U'], 'p0_anchor': anchor['p0']}


def vectorize(points):
    values = {f'N{n}.{field}': float(points[n][field]) for n in NS for field in FIELDS}
    for start, middle, end in old.LINEAGES:
        for model, coeff in MODELS.items():
            for view in (*VIEWS, 'U_anchor', 'U_source'):
                values[f'{model}.{start}.{view}'] = sum(c * points[n][view]
                                                      for c, n in zip(coeff, (start, middle, end)))
    increments = {start: points[middle]['U_anchor'] - points[start]['U_anchor']
                  for start, middle, _ in old.LINEAGES}
    for model, factor in (('q2', .5), ('Jordan', 1)):
        for view in VIEWS:
            values[f'{model}.{view}.generator_drift'] = factor * (
                values[f'{model}.65.{view}'] * increments[85]
                - values[f'{model}.85.{view}'] * increments[65])
    return values


def covariance_component(central, vectors, counts):
    """First-order group-mean influence covariance with true block weights."""
    counts = np.asarray(counts, dtype=float)
    h = counts / counts.sum()
    z = ((1 - h) / h)[:, None] * (central[None, :] - vectors)
    z -= h @ z
    return (z.T * h) @ z / (len(counts) - 1)


def joint_zero(labels, central, covariance, wanted):
    idx = [labels.index(label) for label in wanted]
    c = covariance[np.ix_(idx, idx)]
    value = central[idx]
    if np.linalg.eigvalsh(c)[0] <= 0:
        return {'status': 'unresolved_covariance', 'labels': wanted, 'values': value.tolist(), 'covariance': c.tolist()}
    statistic = float(value @ np.linalg.solve(c, value))
    return {'status': 'computed', 'labels': wanted, 'values': value.tolist(), 'covariance': c.tolist(),
            'chi_square': statistic, 'df': 2, 'nominal_p': math.exp(-statistic / 2)}


def main():
    started = time.perf_counter()
    if DESTINATION.exists():
        raise ValueError('saved two-phase result exists; do not overwrite or repeat')
    contract = json.loads(CONTRACT.read_text())
    manifest = json.loads(old.MANIFEST.read_text())
    hypothesis = json.loads(old.HYPOTHESES.read_text())
    previous_path = ROOT / contract['prior_result']
    previous = json.loads(previous_path.read_text())
    previous_indices = {label: i for i, label in enumerate(previous['labels'])}
    runs = {run['N']: run for run in manifest['runs']}
    deltas = {n: float(Fraction(runs[n]['delta_cos4'])) for n in NS}
    bracket = hypothesis['root_bracket']
    profiles, complements, anchors, fits, raw_states, points = {}, {}, {}, {}, {}, {}
    groups = {}
    for n in NS:
        path = old.OUTPUT / 'raw' / f'n{n}.csv'
        profiles[n] = old.read_raw(path, n)
        complements[n] = load_complement(n, profiles[n])
        complement = complements[n]
        anchors[n] = baseline(complement['sums'].sum(axis=0), int(complement['counts'].sum()), n, deltas[n], bracket)
        fits[n] = estimate(profiles[n], anchors[n], n, deltas[n])
        raw_states[n] = {'raw_bulk': previous['estimates'][f'N{n}.Udot_fugacity']['value'],
                         'U_source': previous['estimates'][f'N{n}.U']['value']}
        points[n] = point_record(n, fits[n], anchors[n], raw_states[n])
        groups.setdefault(runs[n]['dependency_group'], []).append(n)
    central_map = vectorize(points)
    labels = list(central_map)
    central = np.array(list(central_map.values()))
    covariance = np.zeros((len(labels), len(labels)))
    contributions = {}
    for group, sizes in groups.items():
        vectors = []
        prior_vectors = previous['dependency_groups'][group]['delete_one_vectors']
        for batch in range(100):
            changed = dict(points)
            for n in sizes:
                fitted = estimate(profiles[n], anchors[n], n, deltas[n], omitted_batch=batch)
                raw_state = {'raw_bulk': prior_vectors[batch][previous_indices[f'N{n}.Udot_fugacity']],
                             'U_source': prior_vectors[batch][previous_indices[f'N{n}.U']]}
                changed[n] = point_record(n, fitted, anchors[n], raw_state)
            vectors.append(list(vectorize(changed).values()))
        vectors = np.asarray(vectors)
        counts = np.repeat(1000, 100)
        component = covariance_component(central, vectors, counts)
        covariance += component
        contributions['source:' + group] = {
            'stage': 'marked_source', 'Ns': sizes, 'batch_counts': counts.tolist(),
            'delete_one_vectors': vectors.tolist(), 'covariance': component.tolist(),
            'operation': 'same marked batch omitted across this dependency group; affected training folds refitted; complement held fixed'}
        print(f'Completed marked-source covariance: {sizes}', flush=True)

    for group, sizes in groups.items():
        vectors, roots = [], {n: [] for n in sizes}
        counts = complements[sizes[0]]['counts']
        if not all(np.array_equal(counts, complements[n]['counts']) for n in sizes):
            raise ValueError('joint complement block sizes differ within a declared dependency group')
        for batch in range(100):
            changed = dict(points)
            for n in sizes:
                complement = complements[n]
                anchor = baseline(complement['sums'].sum(axis=0) - complement['sums'][batch],
                                  int(complement['counts'].sum() - complement['counts'][batch]),
                                  n, deltas[n], bracket)
                fitted = estimate(profiles[n], anchor, n, deltas[n])
                changed[n] = point_record(n, fitted, anchor, raw_states[n])
                roots[n].append(anchor['p0'])
            vectors.append(list(vectorize(changed).values()))
        vectors = np.asarray(vectors)
        component = covariance_component(central, vectors, counts)
        covariance += component
        contributions['complement:' + group] = {
            'stage': 'unmarked_complement', 'Ns': sizes, 'batch_counts': counts.tolist(),
            'delete_one_vectors': vectors.tolist(), 'covariance': component.tolist(), 'delete_one_roots': roots,
            'operation': 'same original production block omitted with its actual remaining sample count; source marks and fixed-p_ref learning data unchanged; all baseline-dependent responses updated'}
        print(f'Completed complement covariance: {sizes}', flush=True)

    errors = np.sqrt(np.maximum(0, np.diag(covariance)))
    estimates = {label: {'value': float(value), 'se': float(error), 'z': float(value / error) if error else None}
                 for label, value, error in zip(labels, central, errors)}
    precision, differences = {}, {}
    for n in NS:
        indices = {view: labels.index(f'N{n}.{view}') for view in VIEWS}
        precision[n] = {'standard_errors': {view: float(errors[i]) for view, i in indices.items()},
                        'raw_to_two_phase_SE_ratio': float(errors[indices['raw_bulk']] / errors[indices['two_phase_bulk']]),
                        'variance_by_stage': {view: {stage: float(sum(item['covariance'][i][i] for item in contributions.values() if item['stage'] == stage))
                                                    for stage in ('marked_source', 'unmarked_complement')} for view, i in indices.items()}}
        for view in ('anchor_only_bulk', 'two_phase_bulk'):
            i, j = indices[view], indices['raw_bulk']
            se = math.sqrt(max(0.0, covariance[i, i] + covariance[j, j] - 2 * covariance[i, j]))
            value = float(central[i] - central[j])
            differences[f'N{n}.{view}_minus_raw'] = {'value': value, 'se': se,
                                                    'z': value / se if se else None}
    model_results = {}
    for model in MODELS:
        model_results[model] = {
            'unperturbed_complement': joint_zero(labels, central, covariance, [f'{model}.{n}.U_anchor' for n in (65, 85)]),
            'source_rigidity': {view: joint_zero(labels, central, covariance, [f'{model}.{n}.{view}' for n in (65, 85)]) for view in VIEWS},
            'one_generator_drift': {view: estimates[f'{model}.{view}.generator_drift'] for view in VIEWS}}
    paths = [Path(__file__), ROOT / 'scripts/norm4_source_two_phase_core.py', ROOT / 'scripts/norm4_source_two_phase_inputs.py',
             ROOT / 'scripts/analyze_norm4_source_thermal.py', CONTRACT, old.MANIFEST, old.HYPOTHESES]
    result = {
        'schema': 'matching-one.norm4-source-two-phase.v1', 'status': 'computed_zero_replay_two_phase_source_analysis',
        'execution_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        'contract': contract, 'labels': labels, 'estimates': estimates, 'covariance': covariance.tolist(),
        'covariance_contributions': contributions, 'by_N': {n: {'anchor': anchors[n], 'source_fit': fits[n], 'points': points[n]} for n in NS},
        'precision_comparison': precision, 'paired_estimator_differences': differences, 'model_extension_diagnostics': model_results,
        'unmarked_inputs': {n: complements[n]['provenance'] for n in NS},
        'marked_inputs': [{'N': n, 'path': str((old.OUTPUT / 'raw' / f'n{n}.csv').relative_to(ROOT)), 'sha256': digest(old.OUTPUT / 'raw' / f'n{n}.csv')} for n in NS],
        'saved_raw_estimator': {'path': contract['prior_result'], 'sha256': digest(previous_path), 'raw_delete_one_vectors_reused': True},
        'code': [{'path': str(path.relative_to(ROOT)), 'sha256': digest(path)} for path in paths],
        'interpretation': 'The primary two-phase view, baseline-only comparator and original raw view estimate one source response; their precision and differences use joint covariance. No automatic improvement, model certification, field identity or independent evidence count follows from using a larger unmarked archive.',
        'uncertainty_boundary': 'First-order influence covariance, independent marked/complement counter blocks under the declared PRNG model; learning coefficients refitted in source omissions, baseline uncertainty propagated with actual unequal block weights. No exact finite-sample coverage or arbitrary-cluster robustness claim.',
        'new_samples': 0, 'configuration_replays': 0, 'server_actions': 0, 'test_suites': [],
        'environment': {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__, 'machine': platform.machine()},
        'elapsed_seconds': time.perf_counter() - started}
    DESTINATION.parent.mkdir(parents=True, exist_ok=True)
    with DESTINATION.open('x') as handle:
        handle.write(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({'elapsed_seconds': result['elapsed_seconds'], 'precision': precision,
                      'source_responses': {n: estimates[f'N{n}.two_phase_bulk'] for n in NS},
                      'model_extension_diagnostics': model_results}, indent=2))


if __name__ == '__main__':
    main()

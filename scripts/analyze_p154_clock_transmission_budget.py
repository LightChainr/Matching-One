#!/usr/bin/env python3
"""Paired calibration uncertainty for a specified birth-clock transmission map."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import platform
import subprocess
import time

import numpy as np

from analyze_norm4_source_endpoint_1m import load_profile
from analyze_norm4_source_thermal import binomial_moments, direction_values

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / 'results/p154-clock-transmission-budget/latest.json'
CONTRACT = ROOT / 'analysis/p154_clock_transmission_budget.json'
FIELDS = ('p0', 'D', 'c', 'K_rel', 'delta', 'R', 'v_clock', 'v_measured', 'v_minus_clock')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def point(sums, samples, n, p, r, measured_v, angular_delta):
    rows = [direction_values(binomial_moments(sums[g], samples, p, n)) for g in range(2)]
    d = np.mean([z['q_p'] for z in rows])
    c = np.mean([z['E_p'] for z in rows]) / d
    b = (rows[0]['E_p'] - rows[1]['E_p']) / angular_delta
    h = (rows[0]['E_pp'] - rows[1]['E_pp']) / angular_delta
    q2 = (rows[0]['q_pp'] - rows[1]['q_pp']) / angular_delta
    t = np.mean([z['q_pp'] for z in rows])
    e2 = np.mean([z['E_pp'] for z in rows])
    if d <= 0 or abs(c) >= 1:
        raise ValueError('nonpositive birth slopes in the declared clock model')
    k_rel = n**(13/8) / (4*d) * (c*h-q2-(b/d)*(c*t-e2))
    delta = 2*r/(d*(1-c*c))
    v_clock = delta*k_rel
    return dict(zip(FIELDS, map(float, (p, d, c, k_rel, delta, r, v_clock,
                                       measured_v, measured_v-v_clock))))


def main():
    start = time.perf_counter()
    if DEST.exists():
        raise ValueError('saved result exists; reproduce in a separate checkout')
    contract = json.loads(CONTRACT.read_text())
    ns = contract['sizes']
    source_path = ROOT / 'results/norm4-source-endpoint-1m/latest.json'
    lag_path = ROOT / 'results/norm4-lagged-source/latest.json'
    source = json.loads(source_path.read_text())
    lag = json.loads(lag_path.read_text())
    if lag['source_result_sha256'] != sha(source_path):
        raise ValueError('lag and source baseline archives no longer align')
    runs = {r['N']: r for r in json.loads((ROOT/'analysis/p40_source_thermal_chain_candidates.json').read_text())['runs']}
    si = {v:i for i,v in enumerate(source['labels'])}
    li = {v:i for i,v in enumerate(lag['labels'])}
    profiles, central, provenance = {}, {}, []
    for n in ns:
        old = ROOT/f'results/norm4-source-thermal/raw/n{n}.csv'
        inc = ROOT/f'results/norm4-source-endpoint-1m/increment/raw/n{n}.csv'
        profiles[n] = load_profile(old, n, 1000, runs[n]) + load_profile(inc, n, 9000, runs[n])
        for path in (old, inc):
            provenance.append({'path': str(path.relative_to(ROOT)), 'sha256': sha(path)})
        saved = lag['by_N'][str(n)]['points']
        central[n] = point(profiles[n].sum(axis=0), 1000000, n, saved['p0'],
                           saved['total.rank1_rootdot'], saved['total.v'],
                           source['by_N'][str(n)]['source']['delta_cos4'])
    labels = [f'N{n}.{field}' for n in ns for field in FIELDS]
    vector = lambda rows: np.array([rows[n][f] for n in ns for f in FIELDS])
    mean = vector(central)
    covariance = np.zeros((len(labels), len(labels)))
    groups = {}
    for name, sg in source['covariance_contributions'].items():
        selected = [n for n in sg['Ns'] if n in ns]
        if not name.startswith('source:') or not selected:
            continue
        lg = lag['covariance_contributions'][name]
        if sg['delete_one_batch_ids'] != list(range(100)) or lg['delete_one_batch_ids'] != list(range(100)):
            raise ValueError('changed source/lag paired deletion order')
        sv, lv = np.asarray(sg['delete_one_vectors']), np.asarray(lg['delete_one_vectors'])
        deleted = []
        for batch in range(100):
            changed = dict(central)
            for n in selected:
                p = lv[batch, li[f'N{n}.p0']]
                if p != sv[batch, si[f'N{n}.p0']]:
                    raise ValueError('source and lag roots are from different omissions')
                changed[n] = point(profiles[n].sum(axis=0)-profiles[n][batch], 990000, n, p,
                                   lv[batch, li[f'N{n}.total.rank1_rootdot']],
                                   lv[batch, li[f'N{n}.total.v']],
                                   source['by_N'][str(n)]['source']['delta_cos4'])
            deleted.append(vector(changed))
        deleted = np.asarray(deleted)
        factor = math.sqrt(.99)*(deleted-deleted.mean(axis=0))
        covariance += factor.T@factor
        groups[name] = {'Ns': selected, 'batch_ids': list(range(100)),
                        'delete_one_vectors': deleted.tolist(), 'factor': factor.tolist()}
    errors = np.sqrt(np.maximum(0, covariance.diagonal()))
    estimates = {name: {'value': float(value), 'se': float(error),
                        'z': float(value/error) if error else None}
                 for name, value, error in zip(labels, mean, errors)}
    planning = {}
    z = contract['planning']['separation_SE']
    budget = contract['planning']['reference_new_samples_per_N']
    for n in ns:
        forecast = estimates[f'N{n}.v_clock']
        measured = estimates[f'N{n}.v_measured']
        gap, calibration_se = forecast['value'], forecast['se']
        target_variance = measured['se']**2
        available_variance = gap**2/z**2-calibration_se**2
        planning[str(n)] = {
            'point_forecast_only_samples_for_3SE': 1e6*z*z*target_variance/gap**2 if gap else None,
            'calibration_aware_samples_for_3SE': 1e6*target_variance/available_variance if available_variance > 0 else None,
            'null_budget_reason': None if available_variance > 0 else 'fixed_discovery_calibration_alone_overlaps_screening_at_3SE',
            'prediction_interval_3SE': [gap-z*calibration_se, gap+z*calibration_se],
            'new_block_samples': budget,
            'new_block_3SE_resolution': z*math.sqrt(target_variance*1e6/budget),
            'training_plus_new_block_separation_in_SE': abs(gap)/math.sqrt(calibration_se**2+target_variance*1e6/budget),
            'boundary': 'independent_future_block_unchanged_efficiency_Gaussian_planning_not_power_certificate',
        }
    DEST.parent.mkdir(parents=True, exist_ok=True)
    result = {'schema': contract['schema'], 'status': 'computed_discovery_calibration_for_prospective_budget',
              'execution_commit': subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
              'contract': contract, 'contract_sha256': sha(CONTRACT), 'script_sha256': sha(Path(__file__)),
              'inputs': provenance+[{'path':str(p.relative_to(ROOT)), 'sha256':sha(p)} for p in (source_path,lag_path)],
              'labels': labels, 'estimates': estimates, 'covariance': covariance.tolist(),
              'dependency_groups': groups, 'planning': planning,
              'environment': {'python':platform.python_version(),'machine':platform.machine(),'numpy':np.__version__},
              'elapsed_seconds':time.perf_counter()-start, 'new_random_samples':0,
              'configuration_replays':0, 'root_solvers':0, 'test_suites':0, 'server_actions':0}
    DEST.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'planning':planning,'clock_predictions':{n:estimates[f'N{n}.v_clock'] for n in ns},
                      'elapsed_seconds':result['elapsed_seconds']},indent=2))


if __name__ == '__main__':
    main()

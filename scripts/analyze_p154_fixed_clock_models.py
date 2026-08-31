#!/usr/bin/env python3
"""Two fixed transmission laws: pooled calibration, angular/global predictions."""
from __future__ import annotations

from fractions import Fraction
import json
import math
from pathlib import Path
import platform
import subprocess
import time

import numpy as np

from analyze_norm4_lagged_source import load_events, point as lag_point
from analyze_norm4_source_endpoint_1m import load_profile
from analyze_norm4_source_thermal import binomial_moments, direction_values
from analyze_p154_clock_transmission_budget import sha

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / 'analysis/p154_fixed_clock_models.json'
DEST = ROOT / 'results/p154-fixed-clock-models/latest.json'


def cos4(period):
    a, b = period
    return Fraction(a**4-6*a*a*b*b+b**4, (a*a+b*b)**2)


def point(profile, events, samples, n, p, weights):
    response, kernels = lag_point(profile, events, samples, n, p)
    jets = [direction_values(binomial_moments(profile[g], samples, p, n)) for g in range(2)]
    qp, ep, qpp, epp = (np.array([z[key] for z in jets])
                        for key in ('q_p', 'E_p', 'q_pp', 'E_pp'))
    f1, f2 = (qp-ep)/2, (qp+ep)/2
    births = np.stack((kernels[..., 0]+kernels[..., 1],
                       kernels[..., 2]+kernels[..., 1]), axis=-1)
    observed = np.array([binomial_moments(births[g], 1, p, n)[0] for g in range(2)])
    angular = float(weights[0]-weights[1])
    p4 = lambda values: (values[0]-values[1])/angular
    d, b, h, t = qp.mean(), p4(ep), p4(epp), qpp.mean()
    r, rootdot, v = (response[key] for key in ('total.rank1_rootdot', 'total.rootdot', 'total.v'))
    out = {'p0': p, 'R': r, 'rootdot': rootdot, 'V': v,
           'entry_contrast': observed[0, 0]-observed[1, 0],
           'completion_contrast': observed[0, 1]-observed[1, 1]}
    for name, w in (('scalar', np.ones(2)), ('fourfold', np.array(list(map(float, weights))))):
        cw = (w*ep).mean()/d
        sw = (w*qp).mean()-cw*ep.mean()
        if sw == 0:
            raise ValueError(f'{name} relative shift is unidentifiable at N{n}')
        delta = 2*r/sw
        m = rootdot-delta*cw/2
        predicted = np.stack((-(m-delta*w/2)*f1, -(m+delta*w/2)*f2), axis=-1)
        k = n**(13/8)/(4*d)*(cw*h-p4(w*qpp)-(b/d)*(cw*t-(w*epp).mean()))
        pred_v = delta*k
        out.update({f'{name}.m': m, f'{name}.delta': delta,
                    f'{name}.K': k, f'{name}.predicted_V': pred_v,
                    f'{name}.residual_V': v-pred_v})
        for j, label in enumerate(('entry', 'completion')):
            contrast = predicted[0, j]-predicted[1, j]
            out[f'{name}.predicted_{label}_contrast'] = contrast
            out[f'{name}.residual_{label}_contrast'] = out[f'{label}_contrast']-contrast
        # Pooled calibration identities, retained as deterministic residuals.
        predicted_qh = predicted.sum(axis=1).mean()
        predicted_eh = (predicted[:, 1]-predicted[:, 0]).mean()
        out[f'{name}.calibration_root_residual'] = -predicted_qh/d-rootdot
        out[f'{name}.calibration_rank1_residual'] = -predicted_eh-ep.mean()*rootdot-r
    return {key: float(value) for key, value in out.items()}


def main():
    start = time.perf_counter()
    if DEST.exists():
        raise ValueError('saved result exists; reproduce in a separate checkout')
    contract = json.loads(CONTRACT.read_text())
    ns = contract['sizes']
    source_path = ROOT/'results/norm4-source-endpoint-1m/latest.json'
    lag_path = ROOT/'results/norm4-lagged-source/latest.json'
    source, lag = json.loads(source_path.read_text()), json.loads(lag_path.read_text())
    if lag['source_result_sha256'] != sha(source_path):
        raise ValueError('source and lag archives differ')
    runs = {r['N']: r for r in json.loads((ROOT/'analysis/p40_source_thermal_chain_candidates.json').read_text())['runs']}
    li = {name: i for i, name in enumerate(lag['labels'])}
    profiles, events, weights, central, inputs = {}, {}, {}, {}, []
    for n in ns:
        old = ROOT/f'results/norm4-source-thermal/raw/n{n}.csv'
        inc = ROOT/f'results/norm4-source-endpoint-1m/increment/raw/n{n}.csv'
        evt = ROOT/f'results/norm4-lagged-source/raw/n{n}.csv.gz'
        profiles[n] = load_profile(old, n, 1000, runs[n])+load_profile(inc, n, 9000, runs[n])
        events[n] = load_events(evt, n, 10000, runs[n])
        weights[n] = [cos4(runs[n][g]) for g in ('first', 'second')]
        central[n] = point(profiles[n].sum(axis=0), events[n].sum(axis=0), 1000000, n,
                           lag['by_N'][str(n)]['points']['p0'], weights[n])
        inputs.extend({'path': str(path.relative_to(ROOT)), 'sha256': sha(path)} for path in (old, inc, evt))
    fields = list(central[ns[0]])
    labels = [f'N{n}.{field}' for n in ns for field in fields]
    vector = lambda rows: np.array([rows[n][field] for n in ns for field in fields])
    mean = vector(central)
    covariance = np.zeros((len(labels), len(labels)))
    groups = {}
    for name, group in lag['covariance_contributions'].items():
        selected = [n for n in group['Ns'] if n in ns]
        if not selected:
            continue
        if group['delete_one_batch_ids'] != list(range(100)):
            raise ValueError('unaligned deletion order')
        saved, deleted = np.asarray(group['delete_one_vectors']), []
        for batch in range(100):
            changed = dict(central)
            for n in selected:
                changed[n] = point(profiles[n].sum(axis=0)-profiles[n][batch],
                                   events[n].sum(axis=0)-events[n][batch], 990000, n,
                                   saved[batch, li[f'N{n}.p0']], weights[n])
            deleted.append(vector(changed))
        deleted = np.asarray(deleted)
        factor = math.sqrt(.99)*(deleted-deleted.mean(axis=0))
        covariance += factor.T@factor
        groups[name] = {'Ns': selected, 'batch_ids': list(range(100)),
                        'delete_one_vectors': deleted.tolist(), 'factor': factor.tolist()}
    se = np.sqrt(np.maximum(covariance.diagonal(), 0))
    estimates = {name: {'value': float(value), 'se': float(error),
                        'residual_over_SE': float(value/error) if error else None}
                 for name, value, error in zip(labels, mean, se)}
    result = {'schema': contract['schema'], 'contract': contract,
              'execution_commit': subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
              'contract_sha256': sha(CONTRACT), 'script_sha256': sha(Path(__file__)),
              'inputs': inputs+[{'path':str(path.relative_to(ROOT)), 'sha256':sha(path)} for path in (source_path,lag_path)],
              'cos4_exact': {n: list(map(str, weights[n])) for n in ns},
              'labels': labels, 'estimates': estimates, 'covariance': covariance.tolist(),
              'dependency_groups': groups,
              'environment': {'python':platform.python_version(),'machine':platform.machine(),'numpy':np.__version__},
              'new_random_samples':0, 'configuration_replays':0, 'root_solvers':0, 'test_suites':0,
              'server_actions':0, 'elapsed_seconds':time.perf_counter()-start}
    DEST.parent.mkdir(parents=True, exist_ok=True)
    DEST.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({name: value for name,value in estimates.items()
                      if any(key in name for key in ('residual_entry', 'residual_completion', 'predicted_V', 'residual_V'))},indent=2))
    print('elapsed_seconds', result['elapsed_seconds'])


if __name__ == '__main__':
    main()

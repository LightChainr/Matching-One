#!/usr/bin/env python3
"""Exact cost of a specified oracle-centred iid estimator at the existing m64.

This is a post-result feasibility calculation, not a prospective validation
or an extension of the frozen multiplier grid. No target data are generated.
"""
from pathlib import Path
from fractions import Fraction as F
import csv
import hashlib
import json
import math
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
INPUT = ROOT.parent/'p337-finite-law-window-20260831'
BASE = 'cae9c8997b5994c218bfe060f75656137f745755'
N, M = 25, 64


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def encode(x):
    scale = 10**50
    lo, hi = F(math.floor(x.lo*scale), scale), F(math.ceil(x.hi*scale), scale)
    return {'lower_fraction': str(lo), 'upper_fraction': str(hi),
            'midpoint_approx': float((x.lo+x.hi)/2),
            'serialization': 'outward_1e_minus50_rational_grid'}


def main():
    tic = time.perf_counter()
    repo = Path(subprocess.check_output(['git', 'rev-parse', '--show-toplevel'], text=True).strip())
    checked = {}
    for name in ('inputs/axis.csv', 'inputs/tilted.csv', 'vendor/interval_backend.py',
                 'results/latest.json'):
        path = INPUT/name
        data = path.read_bytes()
        assert data == subprocess.check_output(['git', 'show', BASE+':'+path.relative_to(repo).as_posix()])
        checked[name] = digest(path)
    sys.path.insert(0, str(INPUT/'vendor'))
    from interval_backend import Interval as I
    original = json.loads((INPUT/'results/latest.json').read_text())
    rows = {}
    for geometry in ('axis', 'tilted'):
        with (INPUT/'inputs'/f'{geometry}.csv').open() as fh:
            reader = csv.DictReader(fh)
            assert reader.fieldnames == ['k', 'g', 'q', 'count']
            rows[geometry] = [tuple(int(row[x]) for x in reader.fieldnames) for row in reader]
        assert [sum(c for k, _, _, c in rows[geometry] if k == j) for j in range(N+1)] == [math.comb(N, j) for j in range(N+1)]
    laws = {}
    for law, delta in (('star', 0), ('drop', 1)):
        old = original['laws'][law]
        root = old['root_certificate']
        h = I(F(root['h_root_lower_fraction']), F(root['h_root_upper_fraction']))
        powers = [I.of(1)]
        for _ in range(N):
            powers.append(powers[-1]*h)
        cells, covariances, variances = {}, {}, {}
        for geometry in ('axis', 'tilted'):
            sums = {name: I.of(0) for name in ('Z', 'K', 'q', 'Kq', 'I1', 'KI1', 'K2I1')}
            for k, g, q, count in rows[geometry]:
                w = powers[k]*(F(count)*F(M)**(-g+delta*(q+1)))
                values = (1, k, q, k*q, int(q == 0), k*int(q == 0), k*k*int(q == 0))
                for name, value in zip(sums, values):
                    sums[name] += w*value
            moments = {name: value/sums['Z'] for name, value in sums.items() if name != 'Z'}
            mu, p1, k1 = moments['K'], moments['I1'], moments['KI1']
            c = k1-mu*p1
            second = moments['K2I1']-2*mu*k1+mu**2*p1
            variance = second-c**2
            conditional_mu = k1/p1
            conditional_variance = moments['K2I1']/p1-conditional_mu**2
            assert variance.lo > 0 and conditional_variance.lo > 0
            covariances[geometry], variances[geometry] = c, variance
            cells[geometry] = {
                'mean_K': encode(mu), 'P1': encode(p1),
                'rank1_mean_K': encode(conditional_mu),
                'rank1_variance_K': encode(conditional_variance),
                'rank1_mean_K_minus_mean_K': encode(conditional_mu-mu),
                'Cov_K_I1': encode(c), 'E_X_squared': encode(second),
                'Var_X': encode(variance),
                'Cov_K_q': encode(moments['Kq']-mu*moments['q'])
            }
        theta = covariances['axis']-covariances['tilted']
        assert theta.lo > 0 if law == 'star' else theta.hi < 0
        variance_sum = variances['axis']+variances['tilted']
        n_unit = variance_sum/(theta**2)
        n_three = 9*n_unit
        laws[law] = {
            'geometry': cells, 'theta_difference': encode(theta),
            'variance_sum': encode(variance_sum),
            'iid_draws_per_geometry_for_SNR1_interval': encode(n_unit),
            'necessary_draws_per_geometry_for_SNR1': math.ceil(n_unit.lo),
            'iid_draws_per_geometry_for_SNR3_interval': encode(n_three),
            'necessary_draws_per_geometry_for_SNR3': math.ceil(n_three.lo),
            'root_interval_reused_without_refinement': {
                'lo': root['h_root_lower_fraction'], 'hi': root['h_root_upper_fraction']}
        }
    out = {
        'status': 'EXACT_VARIANCE_FOR_SPECIFIED_ORACLE_IID_ESTIMATOR',
        'input_commit': BASE, 'input_sha256': checked,
        'N': N, 'm': M, 'laws': laws,
        'estimator': 'X_g=(K-E_g K)*1_rank1; independent iid samples with n per geometry; theta=E X_axis-E X_tilted',
        'oracle_inputs': ['true pooled root', 'each exact mean K', 'true pooled thermal denominator'],
        'formula': 'Var(theta_hat)=(Var X_axis+Var X_tilted)/n; SNR=abs(theta)/sqrt(Var(theta_hat))',
        'scope': 'Estimator-specific second-moment budget, not a 95 percent interval or all-algorithm complexity lower bound; autocorrelation, root and mean estimation excluded.',
        'new_samples': 0, 'new_configurations': 0, 'cloud_jobs': 0,
        'elapsed_seconds': time.perf_counter()-tic,
        'code_sha256': digest(Path(__file__))
    }
    (ROOT/'result.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({'status': out['status'], 'elapsed_seconds': out['elapsed_seconds'],
                     'SNR3_necessary_n_per_geometry': {k: v['necessary_draws_per_geometry_for_SNR3'] for k, v in laws.items()},
                     'rank1_delta_K': {k: {g: c['rank1_mean_K_minus_mean_K']['midpoint_approx'] for g, c in v['geometry'].items()} for k, v in laws.items()}}, indent=2))


if __name__ == '__main__':
    main()

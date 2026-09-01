#!/usr/bin/env python3
"""Independent Decimal row-square review; never imports/runs the original scorer."""
from datetime import datetime, timezone
from decimal import Decimal as D, localcontext, ROUND_CEILING
from fractions import Fraction as F
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
INPUT = ROOT.parent / 'p337-finite-law-window-20260831'
BASE = 'cae9c8997b5994c218bfe060f75656137f745755'
CHECKS = []


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def decimal(fraction):
    x = F(fraction)
    return D(x.numerator)/D(x.denominator)


def contained(value, interval, name):
    # Decimal result and saved rational endpoints are compared exactly here.
    lo, hi = F(interval['lower_fraction']), F(interval['upper_fraction'])
    assert lo <= F(value) <= hi, (name, str(value), interval)
    CHECKS.append(name)


def cell(rows, h, delta):
    weighted = [(k, q, D(count)*h**k*D(64)**(-g+delta*(q+1)))
                for k, g, q, count in rows]
    Z = sum(w for k, q, w in weighted)
    mu = sum(D(k)*w for k, q, w in weighted)/Z
    mean_q = sum(D(q)*w for k, q, w in weighted)/Z
    rank1 = [(k, w) for k, q, w in weighted if q == 0]
    z1 = sum(w for k, w in rank1)
    mean1 = sum(D(k)*w for k, w in rank1)/z1
    covariance = sum((D(k)-mu)*w for k, w in rank1)/Z
    # The defining pointwise square, not a K^2 - 2 mu K + mu^2 expansion.
    second = sum((D(k)-mu)**2*w for k, w in rank1)/Z
    # Also compute variance directly over every state, including rank 0/2.
    variance = sum(((D(k)-mu if q == 0 else D(0))-covariance)**2*w
                   for k, q, w in weighted)/Z
    cov_q = sum((D(k)-mu)*(D(q)-mean_q)*w for k, q, w in weighted)/Z
    return {
        'mean_K': mu, 'P1': z1/Z, 'rank1_mean_K': mean1,
        'rank1_variance_K': sum((D(k)-mean1)**2*w for k, w in rank1)/z1,
        'rank1_mean_K_minus_mean_K': mean1-mu, 'Cov_K_I1': covariance,
        'E_X_squared': second, 'Var_X': variance, 'Cov_K_q': cov_q,
    }, mean_q


def main():
    started = time.perf_counter()
    started_utc = datetime.now(timezone.utc).isoformat()
    expected = json.loads((ROOT/'result.json').read_text())
    original = json.loads((INPUT/'results/latest.json').read_text())
    assert expected['input_commit'] == BASE
    assert expected['code_sha256'] == sha(ROOT/'score.py')
    assert original['fixed_parameters']['N'] == expected['N'] == 25
    assert original['fixed_parameters']['m'] == expected['m'] == 64
    assert F(original['fixed_parameters']['DeltaCos4_fraction']) == F(1152, 625)
    repo = Path(subprocess.check_output(['git', 'rev-parse', '--show-toplevel'],
                                      cwd=ROOT, text=True).strip())
    for name, digest in expected['input_sha256'].items():
        path = INPUT/name
        assert sha(path) == digest
        assert path.read_bytes() == subprocess.check_output(
            ['git', 'show', BASE+':'+path.relative_to(repo).as_posix()], cwd=ROOT)
    rows = {}
    for geometry in ('axis', 'tilted'):
        with (INPUT/'inputs'/f'{geometry}.csv').open() as handle:
            reader = csv.DictReader(handle)
            assert reader.fieldnames == ['k', 'g', 'q', 'count']
            rows[geometry] = [tuple(int(row[k]) for k in reader.fieldnames) for row in reader]
        assert all(c > 0 and q in (-1, 0, 1) for k, g, q, c in rows[geometry])
        assert len({(k, g, q) for k, g, q, c in rows[geometry]}) == len(rows[geometry])
        assert [sum(c for k, g, q, c in rows[geometry] if k == j) for j in range(26)] == [math.comb(25, j) for j in range(26)]

    evaluated = {}
    for precision in (120, 160):
        with localcontext() as context:
            context.prec = precision
            laws = {}
            for law, delta in (('star', 0), ('drop', 1)):
                target = expected['laws'][law]
                old = original['laws'][law]
                root = target['root_interval_reused_without_refinement']
                assert root == {'lo': old['root_certificate']['h_root_lower_fraction'],
                                'hi': old['root_certificate']['h_root_upper_fraction']}
                h_lo, h_hi = F(root['lo']), F(root['hi'])
                points = {}
                for point, h_fraction in [('lower', h_lo), ('midpoint', (h_lo+h_hi)/2), ('upper', h_hi)]:
                    h = decimal(h_fraction)
                    cells, means_q = {}, []
                    for geometry in ('axis', 'tilted'):
                        values, mean_q = cell(rows[geometry], h, delta)
                        for name, value in values.items():
                            contained(value, target['geometry'][geometry][name],
                                      f'{precision}/{law}/{point}/{geometry}/{name}')
                        for name, value in [('q', mean_q), ('q_h', values['Cov_K_q']/h),
                                            ('E', 1-values['P1']), ('E_h', -values['Cov_K_I1']/h),
                                            ('P_rank1', values['P1'])]:
                            contained(value, old['geometry'][geometry][name],
                                      f'{precision}/{law}/{point}/{geometry}/original_{name}')
                        cells[geometry] = values
                        means_q.append(mean_q)
                    theta = cells['axis']['Cov_K_I1']-cells['tilted']['Cov_K_I1']
                    variance_sum = cells['axis']['Var_X']+cells['tilted']['Var_X']
                    n1 = variance_sum/(theta*theta)
                    n3 = 9*n1
                    for name, value in [('theta_difference', theta), ('variance_sum', variance_sum),
                                        ('iid_draws_per_geometry_for_SNR1_interval', n1),
                                        ('iid_draws_per_geometry_for_SNR3_interval', n3)]:
                        contained(value, target[name], f'{precision}/{law}/{point}/{name}')
                    thermal_cov = (cells['axis']['Cov_K_q']+cells['tilted']['Cov_K_q'])/2
                    U_over_A = -theta/(decimal(F(1152, 625))*thermal_cov)
                    contained(U_over_A, old['U_over_A25'], f'{precision}/{law}/{point}/original_U_over_A')
                    contained(thermal_cov/h, old['positive_slope_D'], f'{precision}/{law}/{point}/original_D_h')
                    contained(-theta/(h*decimal(F(1152, 625))), old['thermal_difference_over_Delta'],
                              f'{precision}/{law}/{point}/original_Y_h')
                    for snr, value in [(1, n1), (3, n3)]:
                        interval = target[f'iid_draws_per_geometry_for_SNR{snr}_interval']
                        lo, hi = F(interval['lower_fraction']), F(interval['upper_fraction'])
                        recorded = target[f'necessary_draws_per_geometry_for_SNR{snr}']
                        # Fraction ceil and Decimal ceil: never convert a large count to float.
                        assert math.ceil(lo) == math.ceil(hi) == recorded
                        assert int(value.to_integral_value(rounding=ROUND_CEILING)) == recorded
                    points[point] = {'geometry': {g: {k: str(v) for k, v in values.items()}
                                                  for g, values in cells.items()},
                        'theta': str(theta), 'variance_sum': str(variance_sum),
                        'SNR1_n_per_geometry': str(n1), 'SNR3_n_per_geometry': str(n3),
                        'U_over_A_from_covariances': str(U_over_A),
                        'pooled_q': str(sum(means_q)/2)}
                assert D(points['lower']['pooled_q']) < 0 < D(points['upper']['pooled_q'])
                laws[law] = points
            evaluated[str(precision)] = laws
    # Raising precision leaves the reported thresholds and all tested enclosures unchanged.
    review = {'status': 'PASS_NO_REQUIRED_CORRECTION', 'input_commit': BASE,
        'method': 'Decimal 120 and 160 digits; direct per-row centered squares; lower/mid/upper saved roots only',
        'interval_checks_passed': len(CHECKS), 'checks': CHECKS,
        'input_sha256': expected['input_sha256'], 'target_score_sha256': sha(ROOT/'score.py'),
        'target_result_sha256': sha(ROOT/'result.json'),
        'SNR3_exact_integer_per_geometry': {law: expected['laws'][law]['necessary_draws_per_geometry_for_SNR3']
                                          for law in ('star', 'drop')},
        'evaluated': evaluated,
        'scope': 'Independent high-precision numerical cross-check of reviewed exact interval arithmetic; estimator-specific SNR budget, not a confidence guarantee or universal lower bound.',
        'new_samples': 0, 'new_configurations': 0, 'new_roots': 0, 'cloud_jobs': 0,
        'original_scorer_imported_or_executed': False,
        'started_utc': started_utc, 'finished_utc': datetime.now(timezone.utc).isoformat(),
        'elapsed_seconds': time.perf_counter()-started, 'command': sys.argv,
        'review_code_sha256': sha(Path(__file__))}
    (ROOT/'review.json').write_text(json.dumps(review, indent=2)+'\n')
    print(json.dumps({k:review[k] for k in ['status', 'interval_checks_passed', 'SNR3_exact_integer_per_geometry', 'elapsed_seconds']}))


if __name__ == '__main__':
    main()

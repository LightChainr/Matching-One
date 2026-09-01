#!/usr/bin/env python3
"""Independent raw-histogram covariance reproduction at the sole frozen m64.

No import of score.py functions; only the already verified exact Interval
arithmetic is reused. No enumeration, simulation, new law or multiplier.
"""
from pathlib import Path
from fractions import Fraction as F
from datetime import datetime, timezone
import csv
import hashlib
import json
import math
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
FREEZE = '375a6f0ce67d46871ec97ea338fdf1342ed33e30'
N, M = 25, 64
DELTA = F(1152, 625)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def git_blob(rev, path):
    return subprocess.check_output(['git', 'show', rev+':'+path], cwd=ROOT)


def read_rows(path):
    reader = csv.DictReader(path.open())
    assert reader.fieldnames == ['k', 'g', 'q', 'count']
    rows = [tuple(int(row[key]) for key in reader.fieldnames) for row in reader]
    assert len({row[:3] for row in rows}) == len(rows)
    counts = [0]*(N+1)
    for k, g, q, count in rows:
        assert 0 <= k <= N and g >= 0 and q in (-1, 0, 1) and count > 0
        counts[k] += count
    assert counts == [math.comb(N, k) for k in range(N+1)]
    return rows


def evaluate_integer_polynomial(c, x):
    ans = F(0)
    for y in reversed(c):
        ans = ans*x+y
    return ans


def root_from_raw_pairs(rows, rank_coefficient, original):
    # A common positive integer rescaling of all row weights, followed by
    # direct row-pair accumulation: coefficient weight*(q_axis+q_tilt).
    shift = max(g-rank_coefficient*(q+1) for tab in rows.values() for _, g, q, _ in tab)
    weighted = {geometry: [(k, q, n*M**(shift-g+rank_coefficient*(q+1)))
                           for k, g, q, n in tab] for geometry, tab in rows.items()}
    coeff = [0]*(2*N+1)
    for ka, qa, wa in weighted['axis']:
        for kt, qt, wt in weighted['tilted']:
            coeff[ka+kt] += wa*wt*(qa+qt)
    divisor = math.gcd(*coeff)
    coeff = [x//divisor for x in coeff]
    assert coeff == list(map(int, original['primitive_integer_coefficients_ascending']))
    signs = [1 if x > 0 else -1 for x in coeff if x]
    assert signs == original['nonzero_coefficient_signs_ascending']
    variations = sum(a != b for a, b in zip(signs, signs[1:]))
    assert variations == 1 and coeff[0] < 0 < coeff[-1]
    lo = F(original['h_root_lower_fraction'])
    hi = F(original['h_root_upper_fraction'])
    assert 0 < lo < hi
    assert evaluate_integer_polynomial(coeff, lo) < 0 < evaluate_integer_polynomial(coeff, hi)
    initial_hi = F(original['initial_upper_fraction'])
    width = initial_hi/F(2)**160
    assert original['bisections_completed'] == 160 and hi-lo == width
    assert (lo/width).denominator == 1
    assert evaluate_integer_polynomial(coeff, 0) < 0 < evaluate_integer_polynomial(coeff, initial_hi)
    return Interval(lo, hi), {
        'primitive_integer_polynomial_matches_from_raw_row_pairs': True,
        'degree': len(coeff)-1, 'descartes_variations': variations,
        'exact_endpoint_signs': [-1, 1],
        'unique_positive_root': True,
        'same_160_bisection_bracket_verified': True,
        'common_integer_row_weight_rescaling_power': shift,
        'lower_fraction': str(lo), 'upper_fraction': str(hi)
    }


def raw_moments(rows, rank_coefficient, h):
    # Every original histogram row contributes directly; no grouped moment
    # polynomial or algebraically cancelled quotient derivative is imported.
    sums = {key: Interval.of(0) for key in ('Z', 'K', 'q', 'qK', 'I1', 'KI1', 'E')}
    powers = [Interval.of(1)]
    for _ in range(N):
        powers.append(powers[-1]*h)
    for k, g, q, count in rows:
        weight = powers[k]*(F(count)*F(M)**(-g+rank_coefficient*(q+1)))
        values = {'Z': 1, 'K': k, 'q': q, 'qK': q*k,
                  'I1': int(q == 0), 'KI1': k*int(q == 0), 'E': q*q}
        for key, value in values.items():
            sums[key] += weight*value
    assert sums['Z'].lo > 0
    moments = {key: value/sums['Z'] for key, value in sums.items() if key != 'Z'}
    covariance_kq = moments['qK']-moments['K']*moments['q']
    covariance_ki1 = moments['KI1']-moments['K']*moments['I1']
    return sums, moments, covariance_kq, covariance_ki1


def encode(x):
    grid = 10**45
    lo, hi = F(math.floor(x.lo*grid), grid), F(math.ceil(x.hi*grid), grid)
    return {'lower_fraction': str(lo), 'upper_fraction': str(hi),
            'midpoint_approx': float((x.lo+x.hi)/2), 'width_approx': float(hi-lo),
            'serialization': 'outward_1e_minus45_rational_grid'}


def compare(x, published):
    lo, hi = F(published['lower_fraction']), F(published['upper_fraction'])
    assert lo <= x.lo <= x.hi <= hi
    return {'independent_enclosure': encode(x), 'inside_original_published_enclosure': True}


def main():
    started = time.perf_counter()
    original_bytes = (ROOT/'results/latest.json').read_bytes()
    original = json.loads(original_bytes)
    receipt = json.loads((ROOT/'results/run.json').read_text())
    assert receipt['freeze_commit'] == original['provenance']['freeze_commit'] == FREEZE
    assert sha(ROOT/'results/latest.json') == receipt['latest_json_sha256']
    assert sha(ROOT/'results/REPORT.md') == receipt['REPORT_md_sha256']
    repo = Path(subprocess.check_output(['git', 'rev-parse', '--show-toplevel'], cwd=ROOT, text=True).strip())
    prefix = ROOT.relative_to(repo).as_posix()
    sources = json.loads((ROOT/'inputs/SOURCES.json').read_text())['entries']
    names = ['CONTRACT.md', 'score.py', 'inputs/SOURCES.json']+[x['path'] for x in sources]
    checked = {}
    for name in names:
        data = (ROOT/name).read_bytes()
        assert git_blob(FREEZE, prefix+'/'+name) == data
        assert sha(ROOT/name) == receipt['input_and_frozen_file_checks'][name]['sha256']
        checked[name] = sha(ROOT/name)
    for item in sources:
        data = (ROOT/item['path']).read_bytes()
        assert hashlib.sha256(data).hexdigest() == item['sha256']
        assert git_blob(item['source_commit'], item['source_path']) == data
    # The exact arithmetic class has now been verified against F0 and source.
    sys.path.insert(0, str(ROOT/'vendor'))
    global Interval
    from interval_backend import Interval
    assert original['fixed_parameters']['N'] == N and original['fixed_parameters']['m'] == M
    assert original['fixed_parameters']['DeltaCos4_fraction'] == str(DELTA)
    rows = {g: read_rows(ROOT/'inputs'/f'{g}.csv') for g in ('axis', 'tilted')}
    laws = {}
    for law, coefficient in (('star', 0), ('drop', 1)):
        old = original['laws'][law]
        h, root_checks = root_from_raw_pairs(rows, coefficient, old['root_certificate'])
        cells, c_kq, c_ki1 = {}, {}, {}
        for geometry in ('axis', 'tilted'):
            sums, moments, c_kq[geometry], c_ki1[geometry] = raw_moments(rows[geometry], coefficient, h)
            cell_old = old['geometry'][geometry]
            q_h, e_h = c_kq[geometry]/h, -c_ki1[geometry]/h
            checks = {key: compare(value, cell_old[published]) for key, value, published in [
                ('q', moments['q'], 'q'), ('E', moments['E'], 'E'),
                ('P1', moments['I1'], 'P_rank1'), ('q_h_from_covariance', q_h, 'q_h'),
                ('E_h_from_negative_rank1_covariance', e_h, 'E_h')]}
            probability_upper = F(cell_old['P_rank1']['upper_fraction'])
            required = math.ceil(F(19, 20)/probability_upper)
            bound_old = cell_old['ordinary_sampling_necessary_bound']
            assert probability_upper == F(bound_old['P1_upper_fraction_used'])
            assert required == bound_old['necessary_draws_lower_bound'] and required > 10**9
            cells[geometry] = {
                'raw_unnormalized_sums': {key: encode(value) for key, value in sums.items()},
                'normalized_moments': {key: encode(value) for key, value in moments.items()},
                'Cov_K_q': encode(c_kq[geometry]), 'Cov_K_I1': encode(c_ki1[geometry]),
                'comparisons': checks,
                'necessary_draws_using_published_upper': required,
                'published_P1_upper_fraction': str(probability_upper),
                'union_bound_ceil_reproduced': True
            }
        denominator = (c_kq['axis']+c_kq['tilted'])/2
        assert denominator.lo > 0
        # Cancel the same 1/h BEFORE calculating original U/A25.
        numerator = -(c_ki1['axis']-c_ki1['tilted'])/DELTA
        observer = numerator/denominator
        assert observer.hi < 0 if law == 'star' else observer.lo > 0
        laws[law] = {
            'root_checks': root_checks,
            'normalization': 'each full raw Z separately, at this law\'s own root',
            'cells': cells, 'mean_Cov_K_q': encode(denominator),
            'slope_covariance_reproduction': compare(denominator/h, old['positive_slope_D']),
            'thermal_numerator_covariance_reproduction': compare(numerator/h, old['thermal_difference_over_Delta']),
            'U_over_A25_covariance_reproduction': compare(observer, old['U_over_A25']),
            'sign': 'negative' if law == 'star' else 'positive'
        }
    assert laws['star']['root_checks']['lower_fraction'] != laws['drop']['root_checks']['lower_fraction']
    assert (ROOT/'results/latest.json').read_bytes() == original_bytes
    out = {
        'status': 'PASS_INDEPENDENT_RAW_ROW_COVARIANCE_REPRODUCTION',
        'reviewed_at_utc': datetime.now(timezone.utc).isoformat(),
        'freeze_commit': FREEZE, 'frozen_file_hashes': checked,
        'original_result_sha256': hashlib.sha256(original_bytes).hexdigest(),
        'raw_rows_read': {g: len(rows[g]) for g in rows},
        'counts_each_K_equal_binomial25': True,
        'all_inputs_match_frozen_and_original_git_blobs': True,
        'method': 'Full raw (K,g,q) row sums for Z,K,q,qK,I1,KI1,E; covariance identity; common 1/h cancels. Only Interval reused, no score functions imported.',
        'formula': 'U/A25=-(Cov(K,I1)_axis-Cov(K,I1)_tilted)/(DeltaCos4*mean_g Cov(K,q))',
        'laws': laws,
        'primary_decision_reproduced': original['primary_decision'],
        'resource_decision_reproduced': original['resource_decision'],
        'evaluated_multipliers': [M], 'new_samples': 0, 'new_enumerations': 0, 'cloud_jobs': 0,
        'scope': 'Correctness reproduction at the same fixed point, not new evidence or an efficient-estimator lower bound.',
        'elapsed_seconds': time.perf_counter()-started,
        'review_code_sha256': sha(ROOT/'review.py')
    }
    (ROOT/'review.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({'status': out['status'], 'elapsed_seconds': out['elapsed_seconds'],
                      'U_over_A25': {law: data['U_over_A25_covariance_reproduction']['independent_enclosure'] for law, data in laws.items()},
                      'necessary_draws': {law: {g: x['necessary_draws_using_published_upper'] for g, x in data['cells'].items()} for law, data in laws.items()}}, indent=2))


if __name__ == '__main__':
    main()

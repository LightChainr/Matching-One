#!/usr/bin/env python3
"""Exact information-only continuation witness; no graph enumeration or sampling.

The completions are witnesses in the stated moment/support relaxation, not
claims about the true graph's multi-hole rank function.  S and the physical
chart are unchanged; the obstruction already occurs at t=0.
"""
from fractions import Fraction as F
from pathlib import Path
import csv
import hashlib
import json
import math

ROOT = Path(__file__).resolve().parent
INPUT = ROOT.parent / 'p337-two-coupling-closure-20260831'
M = 25


def trim(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def add(a, b):
    z = [F(0)] * max(len(a), len(b))
    for i, x in enumerate(a):
        z[i] += x
    for i, x in enumerate(b):
        z[i] += x
    return trim(z)


def scale(a, s):
    return trim([s*x for x in a])


def sub(a, b):
    return add(a, scale(b, -1))


def mul(a, b):
    z = [F(0)] * (len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            z[i+j] += x*y
    return trim(z)


def der(a):
    return trim([i*a[i] for i in range(1, len(a))] or [F(0)])


def ev(a, p):
    z = F(0)
    for x in reversed(a):
        z = z*p+x
    return z


def bernstein_sums_to_power(c):
    z = [F(0)] * (M+1)
    for k, x in enumerate(c):
        for j in range(M-k+1):
            z[k+j] += x*(-1)**j*math.comb(M-k, j)
    return trim(z)


def bounds(a, lo, hi):
    mid, rad = (lo+hi)/2, (hi-lo)/2
    d = der(a)
    dd_abs = [abs(x) for x in der(d)]
    remainder = abs(ev(d, mid))*rad + ev(dd_abs, max(abs(lo), abs(hi)))*rad*rad/2
    val = ev(a, mid)
    return val-remainder, val+remainder


def iprod(a, b):
    values = [x*y for x in a for y in b]
    return min(values), max(values)


def idiv(a, b):
    assert not b[0] <= 0 <= b[1]
    return iprod(a, (1/b[1], 1/b[0]))


def encode(interval):
    lo, hi = interval
    grid = 10**35
    low = F((lo*grid).numerator//(lo*grid).denominator, grid)
    high = F(-((-hi*grid).numerator//(-hi*grid).denominator), grid)
    assert low <= lo <= hi <= high
    return {'lower': str(low), 'upper': str(high),
            'midpoint_approx': float((lo+hi)/2),
            'width_approx': float(high-low)}


def exact_scalar(x):
    return {'fraction': str(x), 'approx': float(x)}


def comb(n, k):
    return math.comb(n, k) if 0 <= k <= n else 0


def table(name, complement=False):
    path = INPUT / 'inputs' / name
    rows = list(csv.DictReader(path.open()))
    assert len(rows) == M+1
    for k, row in enumerate(rows):
        n, q, e = int(row['count']), int(row['sum_q']), int(row['sum_e'])
        assert int(row['k']) == k and n == comb(M, k)
        assert 0 <= e <= n and -e <= q <= e
        assert (e+q) % 2 == 0 and (e-q) % 2 == 0
    out = {}
    for field in ('q', 'e'):
        c = [int(row['sum_'+field]) for row in rows]
        if complement:
            c = list(reversed(c))
            if field == 'q':
                c = [-x for x in c]
        if field == 'q':
            out['q_bernstein_means'] = [F(x, comb(M, k)) for k, x in enumerate(c)]
        out[field] = bernstein_sums_to_power(c)
    return out


def score_tail_budget(p):
    """Sum only the two binomial counts, never topological configurations."""
    tau, positive, negative = F(0), F(0), F(0)
    for h in range(2, M+1):
        for k in range(M+1):
            weight = comb(M, h)*comb(M, k)*p**(M-h+k)*(1-p)**(M+h-k)
            score = F(M-h+k, 1)-2*M*p
            score /= p*(1-p)
            tau += weight
            positive += weight*max(score, 0)
            negative += weight*max(-score, 0)
    return tau, positive, negative


def main():
    provenance = json.loads((INPUT/'SOURCES.json').read_text())
    used = []
    for item in provenance['inputs']:
        if item['local_path'] not in ['inputs/endpoint-axis.csv', 'inputs/endpoint-tilted.csv',
                                     'inputs/defect-first.csv', 'inputs/defect-second.csv']:
            continue
        digest = hashlib.sha256((INPUT/item['local_path']).read_bytes()).hexdigest()
        assert digest == item['sha256']
        used.append(item)
    assert len(used) == 4
    intact = [table('endpoint-'+name+'.csv', True) for name in ('axis', 'tilted')]
    hole = [table('defect-'+name+'.csv') for name in ('first', 'second')]
    p = [F(0), F(1)]
    w0 = [F(0)]*M+[F(1)]
    w1 = [F(0)]*(M-1)+[F(M), F(-M)]
    kappa, tau = add(w0, w1), sub([F(1)], add(w0, w1))
    known = [{f: add(mul(w0, x[f]), mul(w1, y[f])) for f in ('q', 'e')}
             for x, y in zip(intact, hole)]

    # Both deterministic B-bit rules have exactly this SAME K-stratum q sum.
    # E_A = 1, E_B = 1-B1*B2-B1*B3+2*B1*B2*B3.
    for k in range(M+1):
        n = comb(M, k)
        q = 2*comb(M-2, k-2)-n
        e_a = n
        e_b = n-2*comb(M-2, k-2)+2*comb(M-3, k-3)
        for e in (e_a, e_b):
            assert abs(q) <= e <= n and (e+q) % 2 == 0
    gq = [F(-1), F(0), F(2)]
    e_gap = [F(0), F(0), F(2), F(-2)]
    q_complete = [add(x['q'], mul(tau, gq)) for x in known]
    qbar = scale(add(*q_complete), F(1, 2))
    # Swapping A/B between geometries leaves EACH complete q profile unchanged.
    e_known_gap = sub(known[0]['e'], known[1]['e'])
    e_contrasts = [add(e_known_gap, mul(tau, e_gap)),
                   sub(e_known_gap, mul(tau, e_gap))]

    bracket = (F(7, 10), F(71, 100))
    assert ev(qbar, bracket[0]) < 0 < ev(qbar, bracket[1])
    # Exclude all other roots on (0,1), without claiming cross-hole monotonicity.
    # For p<=.70: Q <= (1-kappa)*gq+kappa; both inputs of this upper bound
    # increase with p, so its value at .70 is a uniform strictly negative bound.
    left_upper = (1-ev(kappa, bracket[0]))*ev(gq, bracket[0])+ev(kappa, bracket[0])
    assert left_upper < 0
    # For p>=.71 every original 0/1-hole q profile is positive: verify its
    # normalized Bernstein coefficients are nondecreasing and its value at .71.
    original_profiles = intact+hole
    for profile in original_profiles:
        coefficients = profile['q_bernstein_means']
        assert all(a <= b for a, b in zip(coefficients, coefficients[1:]))
    original_q_at_right = [ev(profile['q'], bracket[1]) for profile in original_profiles]
    assert min(original_q_at_right) > 0 and ev(gq, bracket[1]) > 0
    # A coarse whole-bracket slope proof, independent of narrow-root evaluation.
    # |f'_0|,|f'_1| <= 2M by the product-measure derivative formula and |q|<=1.
    lo0, hi0 = bracket
    kapmax = ev(kappa, hi0)
    kapprime_max = M*(M-1)*(1-hi0)*hi0**(M-2)
    gmax = max(abs(ev(gq, lo0)), abs(ev(gq, hi0)))
    slope_coarse_lower = 4*lo0*(1-kapmax)-kapprime_max*(1+gmax)-2*M*kapmax
    assert slope_coarse_lower > 0
    lo, hi = bracket
    for _ in range(110):
        mid = (lo+hi)/2
        if ev(qbar, mid) > 0:
            hi = mid
        else:
            lo = mid
    slope = bounds(der(qbar), lo, hi)
    assert slope[0] > 0
    raw_e = [bounds(der(x), lo, hi) for x in e_contrasts]
    dc4 = F(-1152, 625)  # cos4(5,5)-cos4(1,7) = -1-527/625.
    normalized = [idiv(idiv(x, (dc4, dc4)), slope) for x in raw_e]
    # Rational enclosure of a=50^(13/8)/2: integer nested square roots are exact.
    grid = 10**40
    eighth_floor = math.isqrt(math.isqrt(math.isqrt(50**13*grid**8)))
    area = (F(eighth_floor, 2*grid), F(eighth_floor+1, 2*grid))
    assert (2*area[0])**8 <= 50**13 < (2*area[1])**8
    us = [iprod(area, x) for x in normalized]
    assert us[0][0] > 0 and us[1][1] < 0
    gap = (us[0][0]-us[1][1], us[0][1]-us[1][0])

    # Exact derivative-budget example at a rational reference, not a graph readout.
    pref = F(1, 2)
    tail, bp, bm = score_tail_budget(pref)
    assert tail == ev(tau, pref) and bp-bm == ev(der(tau), pref)
    b1 = bp+bm
    low_qprime = ev(der(scale(add(known[0]['q'], known[1]['q']), F(1, 2))), pref)
    r = {
        'decision': 'endpoint_0_1_plus_support_and_hole_tail_do_not_determine_homogeneous_U',
        'scope': 'Information-only moment/Bernstein relaxation witness, not two physical graph laws; no new source.',
        'fixed_chart': 'p_A=1-epsilon*(1-p), p_B=p; evaluation at epsilon=1,t=0',
        'source': 'S*=C+F4+Bvac unchanged; t=0 already obstructs a continuation guarantee for a family containing t=0',
        'N': 50, 'delta_cos4': str(dc4), 'inputs': used,
        'readout_rules_on_unknown_hole_layers': {'A': 'q=2*B1*B2-1; E=q^2=1',
                                                'B': 'q=B1*B2+B1*B3-1; E=q^2'},
        'shared_q_on_every_B_count': True,
        'all_original_0_1_hole_moments_unchanged': True,
        'new_hole_layer_counts_enumerated': 0,
        'new_random_samples': 0,
        'cloud_jobs': 0,
        'same_entire_q_profiles_and_common_root_in_both_completions': True,
        'root_initial_bracket': [str(x) for x in bracket],
        'unique_common_root_on_entire_open_unit_interval': True,
        'root_uniqueness_checks': {
            'uniform_Q_upper_for_0_p_le_070': exact_scalar(left_upper),
            'four_original_01_q_Bernstein_mean_sequences_nondecreasing': True,
            'four_original_01_q_at_071': [exact_scalar(x) for x in original_q_at_right],
            'gq_at_071': exact_scalar(ev(gq, bracket[1])),
            'reason': 'Q<0 on (0,.70]; Q strictly increasing on [.70,.71]; Q>0 on [.71,1). No claim of cross-hole monotonicity.'
        },
        'whole_initial_bracket_slope_lower_bound': exact_scalar(slope_coarse_lower),
        'common_root': encode((lo, hi)), 'common_positive_slope': encode(slope),
        'epsilon1_tail_at_root': encode(bounds(tau, lo, hi)),
        'raw_Eprime_orientation_differences': [encode(x) for x in raw_e],
        'U_A_first_B_second': encode(us[0]), 'U_B_first_A_second': encode(us[1]),
        'necessary_width_of_any_information_only_U_interval': encode(gap),
        'derivative_budget_at_p_one_half': {
            'tau': exact_scalar(tail), 'positive_score_mass': exact_scalar(bp),
            'negative_score_mass': exact_scalar(bm), 'B1_absolute_score_mass': exact_scalar(b1),
            'known_01_pooled_q_derivative': exact_scalar(low_qprime),
            'permitted_pooled_slope_envelope': encode((low_qprime-b1, low_qprime+b1)),
            'unscaled_P4_Eprime_remainder_bound': exact_scalar(b1/abs(dc4)),
            'normalization_denominator_positive_certified_by_this_budget': low_qprime-b1 > 0
        },
        'readme': 'See THEORY.md for the derivative/root/U error budget and the exact limits of the witness.'
    }
    (ROOT/'results.json').write_text(json.dumps(r, indent=2, ensure_ascii=False)+'\n')
    print(json.dumps({key: r[key] for key in ['decision', 'common_root', 'common_positive_slope',
                                             'U_A_first_B_second', 'U_B_first_A_second',
                                             'necessary_width_of_any_information_only_U_interval']}, indent=2))


if __name__ == '__main__':
    main()

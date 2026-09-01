#!/usr/bin/env python3
"""Independent Decimal reproduction at the fixed N25 calibration/N50 target.

Directly differentiates p**K * (1-p)**(N-K), without importing the frozen
scorer, its interval backend, or its logit-moment formulas. Decimal output is
a numerical cross-check, not a strict Fraction interval certificate.
"""
import argparse
from decimal import Decimal as D, localcontext
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import subprocess


PACKAGE = Path(__file__).resolve().parent
REPO = PACKAGE.parents[1]
FREEZE = '10c666b65566b25ddb8eaa02219947a9c5a261f2'
FIELDS = ('one', 'q', 'e', 's', 'qs', 'es')
BACKEND_SHA = '001a4ec8d85934c11690c1948f47ea8bdc892ad854ee422882cbaa4053fd09db'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_frozen():
    checks = {}
    for name in ('CONTRACT.md', 'score.py'):
        path = PACKAGE / name
        blob = subprocess.check_output(['git', '-C', str(REPO), 'show',
                                       FREEZE + ':' + str(path.relative_to(REPO))])
        if blob != path.read_bytes():
            raise ValueError('Frozen file differs: ' + name)
        checks[name] = sha(path)
    backend = PACKAGE.parent / 'p337-finite-law-window-20260831/vendor/interval_backend.py'
    if sha(backend) != BACKEND_SHA:
        raise ValueError('Pinned backend differs')
    return checks


def read_table(path, n, geometry):
    data = json.loads(path.read_text())
    if data.get('complete') is not True or data.get('N') != n or data.get('geometry') != geometry:
        raise ValueError('Wrong/incomplete fixed population')
    arrays = {name: [0] * (n + 1) for name in FIELDS}
    seen = set()
    for row in data['histogram']:
        k, q, count, sum_s = [row[name] for name in ('K', 'q', 'count', 'sum_S')]
        if any(type(x) is not int for x in (k, q, count, sum_s)):
            raise ValueError('Noninteger sufficient statistics')
        if not 0 <= k <= n or q not in (-1, 0, 1) or count <= 0 or (k, q) in seen:
            raise ValueError('Invalid/duplicate row')
        seen.add((k, q))
        values = (count, q * count, q * q * count, sum_s, q * sum_s, q * q * sum_s)
        for name, value in zip(FIELDS, values):
            arrays[name][k] += value
    if arrays['one'] != [math.comb(n, k) for k in range(n + 1)]:
        raise ValueError('Population differs from the complete binomial law')
    return arrays, {'path': str(path.resolve()), 'sha256': sha(path), 'N': n,
                    'geometry': geometry, 'rows': len(seen), 'count': sum(arrays['one'])}


def powers(p, n):
    pows, compl = [D(1)], [D(1)]
    for _ in range(n):
        pows.append(pows[-1] * p)
        compl.append(compl[-1] * (1 - p))
    return pows, compl


def direct_bernstein_weights(p, n):
    """Three derivative orders, expanded products; no score/cumulant identity."""
    pp, cc = powers(p, n)
    out = []
    for k in range(n + 1):
        l = n - k
        w = pp[k] * cc[l]
        wp = D(0)
        wpp = D(0)
        if k:
            wp += k * pp[k - 1] * cc[l]
        if l:
            wp -= l * pp[k] * cc[l - 1]
        if k >= 2:
            wpp += k * (k - 1) * pp[k - 2] * cc[l]
        if k and l:
            wpp -= 2 * k * l * pp[k - 1] * cc[l - 1]
        if l >= 2:
            wpp += l * (l - 1) * pp[k] * cc[l - 2]
        out.append((w, wp, wpp))
    return out


def root_value(pair, p, n):
    pp, cc = powers(p, n)
    return sum(D(a + b) * pp[k] * cc[n - k]
               for k, (a, b) in enumerate(zip(pair[0]['q'], pair[1]['q']))) / 2


def find_root(pair, n, precision):
    # A numerical reproduction of the same unique positive root; no scan/model selection.
    coefficients = [a + b for a, b in zip(pair[0]['q'], pair[1]['q'])]
    signs = [1 if x > 0 else -1 for x in coefficients if x]
    variations = sum(a != b for a, b in zip(signs, signs[1:]))
    if variations != 1:
        raise ValueError('Pooled positive-root uniqueness does not hold')
    lo, hi = D(0), D(4) / 5
    if not root_value(pair, lo, n) < 0 < root_value(pair, hi, n):
        raise ValueError('Fixed p=[0,4/5] root bracket failed')
    for iteration in range(4 * precision):
        mid = (lo + hi) / 2
        if mid == lo or mid == hi:
            break
        value = root_value(pair, mid, n)
        if value == 0:
            lo = hi = mid
            break
        if value < 0:
            lo = mid
        else:
            hi = mid
    root = (lo + hi) / 2
    return root, {'method': 'Decimal p-Bernstein bisection; not a rigorous interval',
                  'iterations': iteration + 1, 'descartes_sign_variations': variations,
                  'p_low': str(lo), 'p_high': str(hi),
                  'pooled_q_residual': str(root_value(pair, root, n))}


def jets(arrays, weights):
    raw = {name: tuple(sum(D(c) * w[order] for c, w in zip(values, weights))
                       for order in range(3)) for name, values in arrays.items()}
    z0, z1, z2 = raw['one']
    normalized = {}
    for name, (r0, r1, r2) in raw.items():
        f0 = r0 / z0
        f1 = (r1 - f0 * z1) / z0
        f2 = (r2 - 2 * f1 * z1 - f0 * z2) / z0
        normalized[name] = (f0, f1, f2)
    return normalized, raw['one']


def response(f, s, fs):
    # exp(t S) normalization in this geometry, including its p derivative.
    return fs[0] - f[0] * s[0], fs[1] - f[1] * s[0] - f[0] * s[1]


def calculate(pair, n, delta, precision):
    with localcontext() as ctx:
        ctx.prec = precision
        p, root_receipt = find_root(pair, n, precision)
        weights = direct_bernstein_weights(p, n)
        packets, norms, jq, je = [], [], [], []
        for arrays in pair:
            packet, norm = jets(arrays, weights)
            packets.append(packet)
            norms.append(norm)
            jq.append(response(packet['q'], packet['s'], packet['qs']))
            je.append(response(packet['e'], packet['s'], packet['es']))
        mean = lambda values: sum(values) / 2
        d = mean([x['q'][1] for x in packets])
        d_p = mean([x['q'][2] for x in packets])
        d_t = mean([x[1] for x in jq])
        if d <= 0:
            raise ValueError('Nonpositive common-root slope')
        yp = (packets[0]['e'][1] - packets[1]['e'][1]) / delta
        ypp = (packets[0]['e'][2] - packets[1]['e'][2]) / delta
        ypt = (je[0][1] - je[1][1]) / delta
        pt = -mean([x[0] for x in jq]) / d
        u = yp / d
        p_terms = {'direct': ypt / d, 'root_motion': pt * ypp / d,
                   'slope_source': -yp * d_t / (d * d),
                   'slope_root': -yp * pt * d_p / (d * d)}
        v = sum(p_terms.values())
        # Conversion to the frozen z chart uses chain rule ONLY, not logit K moments.
        g, gp = p * (1 - p), 1 - 2 * p
        zt = pt / g
        mz, mzz, mzt = g * d, g * g * d_p + g * gp * d, g * d_t
        yz, yzz, yzt = g * yp, g * g * ypp + g * gp * yp, g * ypt
        z_terms = {'direct': yzt / mz, 'root_motion': zt * yzz / mz,
                   'slope_source': -yz * mzt / (mz * mz),
                   'slope_root': -yz * zt * mzz / (mz * mz)}
        a_n = ctx.power(D(n), D(13) / 8) / 2
        core = {'p_root': p, 'h_root': p / (1 - p), 'D_p': d, 'D_z': mz,
                'root_p_t': pt, 'root_z_t': zt, 'reduced_U': u, 'reduced_V': v,
                'A': a_n, 'U': a_n * u, 'V': a_n * v,
                'p_vs_z_total_difference': v - sum(z_terms.values())}
        norm_error = max(abs(z[0] - 1) for z in norms)
        norm_p_error = max(abs(z[1]) for z in norms)
        norm_pp_error = max(abs(z[2]) for z in norms)
        if max(norm_error, norm_p_error, norm_pp_error) > D(10) ** (-precision + 12):
            raise ValueError('Bernstein normalization/derivative consistency failed')
        return {'precision': precision, 'core': {k: str(x) for k, x in core.items()},
                'p_terms': {k: str(x) for k, x in p_terms.items()},
                'z_terms_from_p_chain_rule': {k: str(x) for k, x in z_terms.items()},
                'root_receipt': root_receipt,
                'normalization_errors': {'Z_minus_1': str(norm_error), 'Z_p': str(norm_p_error),
                                         'Z_pp': str(norm_pp_error)},
                'geometry_jets_p_0_1_2': [{k: [str(y) for y in x] for k, x in packet.items()}
                                        for packet in packets],
                'geometry_source_q_0_1': [[str(x) for x in j] for j in jq],
                'geometry_source_E_0_1': [[str(x) for x in j] for j in je]}


def compare_certificate(result_path, reproduction, input_receipts, frozen_hashes):
    result = json.loads(result_path.read_text())
    if result.get('freeze_commit') != FREEZE or result.get('package_hashes') != frozen_hashes:
        raise ValueError('Published result/freeze identity differs')
    if result.get('input_sha256') != [x['sha256'] for x in input_receipts]:
        raise ValueError('Published result/input identity differs')
    if result.get('backend_sha256') != BACKEND_SHA:
        raise ValueError('Published result/backend identity differs')
    checks = {}

    def check(name, decimal, interval):
        value = Fraction(D(decimal))
        lo, hi = Fraction(interval['lower_fraction']), Fraction(interval['upper_fraction'])
        checks[name] = {'decimal_point_inside_fraction_enclosure': lo <= value <= hi,
                        'interval_width': str(hi - lo)}
        if not checks[name]['decimal_point_inside_fraction_enclosure']:
            raise ValueError('Numerical reproduction outside frozen enclosure: ' + name)

    for name in ('p_root', 'h_root', 'D_z', 'reduced_U', 'reduced_V', 'root_z_t'):
        check(name, reproduction['core'][name], result[name])
    for name, value in reproduction['z_terms_from_p_chain_rule'].items():
        check('term.' + name, value, result['reduced_terms'][name])
    for index, packet in enumerate(reproduction['geometry_jets_p_0_1_2']):
        check('mean_S.' + str(index), packet['s'][0], result['geometry_mean_S'][index])
    vlo = Fraction(result['reduced_V']['lower_fraction'])
    vhi = Fraction(result['reduced_V']['upper_fraction'])
    decisions = {'V_positive': vlo > 0, 'V_negative': vhi < 0,
                 'finite_zero_excluded': vlo > 0 or vhi < 0,
                 'endpoint_positive_sign_rejected': vhi < 0}
    if any(result.get(name) != value for name, value in decisions.items()):
        raise ValueError('Frozen result decision inconsistent with rational enclosure')
    return {'path': str(result_path.resolve()), 'sha256': sha(result_path),
            'checks': checks, 'rational_interval_decisions': decisions,
            'scope': 'Decimal reproduction agrees with existing Fraction certificate; no new certificate'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=('n25-calibration', 'n50-target'), required=True)
    parser.add_argument('--first', type=Path, required=True)
    parser.add_argument('--second', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--score-result', type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError('Will not overwrite an earlier independent review')
    if args.case == 'n25-calibration':
        n, geometries, delta_sign = 25, ([5, 0], [4, 3]), 1
        if args.score_result:
            raise ValueError('N50 certificate cannot be used for the N25 calibration')
    else:
        n, geometries, delta_sign = 50, ([5, 5], [1, 7]), -1
        if not args.score_result:
            raise ValueError('Target review requires the completed frozen primary score')
    frozen_hashes = verify_frozen()
    pair, receipts = [], []
    for path, geometry in zip((args.first, args.second), geometries):
        arrays, receipt = read_table(path, n, geometry)
        pair.append(arrays)
        receipts.append(receipt)
    reproductions = []
    for precision in (120, 160):
        with localcontext() as ctx:
            ctx.prec = precision
            delta = D(delta_sign * 1152) / 625
            reproductions.append(calculate(pair, n, delta, precision))
    with localcontext() as ctx:
        ctx.prec = 160
        differences = {name: abs(D(reproductions[0]['core'][name]) - D(reproductions[1]['core'][name]))
                       for name in reproductions[0]['core']}
        max_difference = max(differences.values())
        if max_difference > D('1e-100'):
            raise ValueError('120/160-digit reproduction is not stable')
        known = None
        if n == 25:
            expected = {'U': D('0.8804661569633677'), 'V': D('0.12616536341416915')}
            known = {name: str(D(reproductions[1]['core'][name]) - value)
                     for name, value in expected.items()}
            if max(abs(D(x)) for x in known.values()) > D('5e-15'):
                raise ValueError('Known N25 calibration failed')
    output = {'case': args.case, 'method': 'Direct Bernstein p first/second derivatives, per-geometry source normalization',
              'status': 'PASS', 'evidence_type': 'Decimal120/160 numerical cross-check; NOT a strict Fraction certificate',
              'script_sha256': sha(Path(__file__)), 'freeze_commit': FREEZE,
              'frozen_file_sha256': frozen_hashes, 'input_receipts': receipts,
              'delta': str(Fraction(delta_sign * 1152, 625)),
              'precision_stability': {'core_max_absolute_difference': str(max_difference),
                                      'differences': {k: str(v) for k, v in differences.items()}},
              'known_N25_display_differences': known, 'reproductions': reproductions,
              'coordinate_note': 'Individual root terms depend on p versus z chart. Chain-rule z terms and total V reproduce the frozen chart.'}
    if args.score_result:
        output['primary_certificate_comparison'] = compare_certificate(
            args.score_result, reproductions[1], receipts, frozen_hashes)
    args.output.write_text(json.dumps(output, indent=2, ensure_ascii=False) + '\n')
    print(json.dumps({'case': args.case, 'status': output['status'],
                      'p_root': reproductions[1]['core']['p_root'],
                      'U': reproductions[1]['core']['U'], 'V': reproductions[1]['core']['V'],
                      'precision_max_difference': str(max_difference),
                      'output': str(args.output.resolve()), 'evidence_type': output['evidence_type']}, indent=2))


if __name__ == '__main__':
    main()

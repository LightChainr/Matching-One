#!/usr/bin/env python3
"""Analytic eta=0 response and zero-frequency susceptibility of the frozen P398 model."""
from __future__ import annotations

import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
import hashlib
import json
from pathlib import Path
import platform
import time

import numpy as np
import scipy
from scipy import linalg, sparse
from scipy.optimize import brentq
from scipy.sparse.linalg import spsolve

from frozen_rate_model import construct, stationary, encode, decode, maxabs, LAGS
from frozen_model import features

ROOT = Path(__file__).resolve().parent
COUNTS = (1, 2, 4, 8, 16, 32, 93)


def covariance(source, metric, propagator):
    return (source.conj().T @ metric @ propagator @ source).conj()


def stationary_tangent(forward, perturbation, pi):
    a = forward.tolil().astype(float)
    rhs = -perturbation @ pi
    a[-1, :] = 1
    rhs[-1] = 0
    dpi = spsolve(a.tocsc(), rhs)
    residual = maxabs(forward @ dpi + perturbation @ pi)
    assert residual < 1e-10 and abs(dpi.sum()) < 1e-10
    return dpi, residual


def evaluator(h, dh, source, metric, dmetric):
    c0 = covariance(source, metric, np.eye(len(h)))
    dc0 = covariance(source, dmetric, np.eye(len(h)))

    def at(t):
        exp_h, dexp = linalg.expm_frechet(-h*t, -dh*t)
        c = covariance(source, metric, exp_h)
        stat = covariance(source, dmetric, exp_h)
        dyn = covariance(source, metric, dexp)
        u = linalg.solve(c0, c)
        du_stat = linalg.solve(c0, stat-dc0 @ u)
        du_dyn = linalg.solve(c0, dyn)
        return {'C': c, 'Cprime': stat+dyn, 'U': u,
                'Uprime': du_stat+du_dyn, 'Uprime_pi': du_stat,
                'Uprime_generator': du_dyn, 'Cprime_pi': stat, 'Cprime_generator': dyn}

    inv_h = linalg.inv(h)
    integrated_c = covariance(source, metric, inv_h)
    integrated_u = linalg.solve(c0, integrated_c)
    integrated_stat = covariance(source, dmetric, inv_h)
    integrated_dyn = covariance(source, metric, -inv_h @ dh @ inv_h)
    integrated_du_pi = linalg.solve(c0, integrated_stat-dc0 @ integrated_u)
    integrated_du_dyn = linalg.solve(c0, integrated_dyn)
    first_moment = linalg.solve(c0, covariance(source, metric, inv_h @ inv_h))
    first_derivative = linalg.solve(c0,
        covariance(source, dmetric, inv_h @ inv_h)
        + covariance(source, metric, -inv_h @ dh @ inv_h @ inv_h-inv_h @ inv_h @ dh @ inv_h)
        -dc0 @ first_moment)
    zero = {'C0': c0, 'C0prime': dc0, 'integrated_U': integrated_u,
            'integrated_Uprime': integrated_du_pi+integrated_du_dyn,
            'integrated_Uprime_pi': integrated_du_pi,
            'integrated_Uprime_generator': integrated_du_dyn,
            'integrated_Cprime': integrated_stat+integrated_dyn,
            'first_time_moment_U': first_moment,
            'first_time_moment_Uprime': first_derivative,
            'poisson_source_residual': maxabs(h @ (inv_h @ source)-source)}
    return at, zero


def encode_dict(d):
    return {k: float(v) if np.ndim(v) == 0 else encode(v) for k, v in d.items()}


def fixed_lag_crossings(at):
    roots = {}
    for i, j, label in ((0, 1, 'minus_plus'), (1, 0, 'plus_minus')):
        found = []
        for left, right in zip(LAGS[1:-1], LAGS[2:]):
            f = lambda t: float(at(t)['Uprime'][i, j].real)
            if f(left)*f(right) < 0:
                found.append({'old_lag_bracket': [left, right],
                              'root': float(brentq(f, left, right, xtol=2e-12))})
        roots[label] = found
    return roots


def spectral_response(h, dh, source, metric, dmetric, kreweras_matrix):
    values, vectors = linalg.eig(h)
    inverse = linalg.inv(vectors)
    assert np.min(values.real) > 0
    ray_phases = [np.vdot(source[:, i], kreweras_matrix @ source[:, i])
                  /np.vdot(source[:, i], source[:, i]) for i in range(2)]
    labels = []
    errors = []
    for i in range(len(values)):
        vector = vectors[:, i]
        distances = [linalg.norm(kreweras_matrix @ vector-phase*vector)
                     /linalg.norm(vector) for phase in ray_phases]
        labels.append(int(np.argmin(distances))); errors.append(min(distances))
    assert max(errors) < 1e-7
    assert labels.count(0) == labels.count(1) == 93
    c0 = covariance(source, metric, np.eye(len(h)))
    dc0 = covariance(source, dmetric, np.eye(len(h)))

    def selected(k):
        result = []
        for ray in (0, 1):
            indices = [i for i, label in enumerate(labels) if label == ray]
            indices.sort(key=lambda i: (values[i].real, abs(values[i].imag)))
            boundary = values[indices[k-1]].real
            # Preserve a complex-conjugate pair when the requested count cuts it.
            result.extend(i for i in indices if values[i].real <= boundary+1e-9)
        return np.array(result, dtype=int)

    def make_at(indices):
        lam = values[indices]
        v, w = vectors[:, indices], inverse[indices, :]
        right = w @ source
        left = source.conj().T @ metric @ v
        left_prime = source.conj().T @ dmetric @ v
        bridge = w @ dh @ v
        differences = lam[:, None]-lam[None, :]
        equal = abs(differences) < 1e-9

        def at(t):
            exponents = np.exp(-lam*t)
            divided = np.zeros(differences.shape, dtype=complex)
            np.divide(exponents[:, None]-exponents[None, :], differences,
                      out=divided, where=~equal)
            repeated = -t*np.exp(-(lam[:, None]+lam[None, :])*t/2)
            divided[equal] = repeated[equal]
            c = ((left*exponents) @ right).conj()
            dc = ((left_prime*exponents) @ right+left @ (divided*bridge) @ right).conj()
            u = linalg.solve(c0, c)
            return {'Uprime': linalg.solve(c0, dc-dc0 @ u)}

        resolvent = v @ np.diag(1/lam) @ w
        integral_u = linalg.solve(c0, covariance(source, metric, resolvent))
        derivative = linalg.solve(c0,
            covariance(source, dmetric, resolvent)
            +covariance(source, metric, -resolvent @ dh @ resolvent)-dc0 @ integral_u)
        return at, derivative

    rows = []
    for k in COUNTS:
        indices = selected(k)
        at, derivative = make_at(indices)
        rows.append({'requested_modes_per_ray': k,
                     'retained_per_ray': [sum(labels[i] == ray for i in indices) for ray in (0, 1)],
                     'retained_indices': indices.tolist(),
                     'integrated_Uprime_re_im': encode(derivative),
                     'Uprime_at_old_lags_re_im': encode(np.array([at(t)['Uprime'] for t in LAGS])),
                     'crossings': fixed_lag_crossings(at)})
    return {'eigenvalues_re_im': encode(values), 'ray_labels': labels,
            'eigenvector_condition_number': float(np.linalg.cond(vectors)),
            'eigen_equation_max_residual': maxabs(h @ vectors-vectors*values),
            'ray_max_residual': max(errors), 'rows': rows,
            'definition': 'Truncate the exact eta0 biorthogonal resolvent/Duhamel expansion on both propagator legs to the slowest modes of each protected ray; retain complete equal-real-part/conjugate groups. This is a signed modal budget, not positive mode probabilities or a new eta-dependent spectral model.'}


def analyze():
    start = time.perf_counter()
    states, parts, q, source, f, t2, complement = construct(8)
    forward = parts[0]+parts[1]
    perturbation = parts[0]-parts[1]
    pi, stationary_residual = stationary(forward)
    dpi, tangent_residual = stationary_tangent(forward, perturbation, pi)
    mass, dmass = -forward.T.tocsr(), -perturbation.T.tocsr()
    h = (q.conj().T @ mass @ q).toarray()
    dh = (q.conj().T @ dmass @ q).toarray()
    metric = (q.conj().T @ q.multiply(pi[:, None])).toarray()
    dmetric = (q.conj().T @ q.multiply(dpi[:, None])).toarray()
    kmatrix = (q.conj().T @ q[complement, :]).toarray()
    assert maxabs(dpi[complement]+dpi) < 1e-10
    assert maxabs(kmatrix @ dh+dh @ kmatrix) < 1e-10
    old = np.load(ROOT/'previous_character_i_zero.npz')
    source_diff = maxabs(source-old['source'])
    baseline_diff = max(maxabs(h-old['mass']), maxabs(metric-old['metric']), source_diff)
    assert baseline_diff < 1e-10
    full_at, zero = evaluator(h, dh, source, metric, dmetric)
    full_packets = [full_at(t) for t in LAGS]
    derivatives = np.array([row['Uprime'] for row in full_packets])
    assert maxabs(derivatives[:, (0, 1), (0, 1)]) < 1e-8
    assert maxabs(derivatives[0]) < 1e-9
    finite = json.loads((ROOT/'previous_finite_intervention.json').read_text())
    secant = decode(finite['contrasts']['finite_central_secant_U_re_im'])
    baseline_c = decode(finite['rows'][0]['C_re_im'])
    assert maxabs(np.array([r['C'] for r in full_packets])-baseline_c) < 1e-9
    named = features(states, f, t2)
    old_stages = json.loads((ROOT/'old_t4_compression.json').read_text())['stages']
    stages = [('instantaneous_two_source', source, ['psi_minus', 'psi_plus'])]
    for stage in old_stages:
        names = stage['exact_named_columns']['selected_columns_in_declaration_order']
        cols = [named[name[2:]][complement] if name.startswith('K_') else named[name] for name in names]
        z = np.asarray(q.conj().T @ np.column_stack(cols))
        stages.append((stage['stage'], z, names))
    compressed = []
    for name, z, columns in stages:
        gram = z.conj().T @ metric @ z
        dgram = z.conj().T @ dmetric @ z
        projected = linalg.solve(gram, z.conj().T @ metric @ h @ z)
        dprojected = linalg.solve(gram,
            z.conj().T @ (dmetric @ h+metric @ dh) @ z-dgram @ projected)
        coefficients = linalg.lstsq(z, source)[0]
        assert maxabs(z @ coefficients-source) < 1e-10
        at, integral = evaluator(projected, dprojected, coefficients, gram, dgram)
        packets = [at(t) for t in LAGS]
        predicted = np.array([p['Uprime'] for p in packets])
        compressed.append({'name': name, 'dimension': z.shape[1], 'columns': columns,
            'integrals': encode_dict(integral),
            'Uprime_at_old_lags_re_im': encode(predicted),
            'integrated_Uprime_error_re_im': encode(integral['integrated_Uprime']-zero['integrated_Uprime']),
            'Uprime_max_absolute_error': maxabs(predicted-derivatives),
            'crossings': fixed_lag_crossings(at),
            'metric_derivative_included': True,
            'boundary': 'Fixed old configuration span; consumes full pi0 and its exact first derivative, not a blind model for stationary weights.'})
    spectrum = spectral_response(h, dh, source, metric, dmetric, kmatrix)
    spectral_full = spectrum['rows'][-1]
    spectral_error = maxabs(decode(spectral_full['Uprime_at_old_lags_re_im'])-derivatives)
    integral_error = maxabs(decode(spectral_full['integrated_Uprime_re_im'])-zero['integrated_Uprime'])
    assert spectral_error < 1e-7 and integral_error < 1e-7
    for row in spectrum['rows']:
        row['integrated_Uprime_error_re_im'] = encode(decode(row['integrated_Uprime_re_im'])-zero['integrated_Uprime'])
        row['Uprime_max_absolute_error'] = maxabs(decode(row['Uprime_at_old_lags_re_im'])-derivatives)
    return {'schema': 'matching-one/p398-eta0-linear-response/v1',
            'scientific_scope': 'Analytic first derivative at eta0 of the same fixed width8 join/detach intervention, old attempt clock and old two sources; no finite difference, new eta point, MC or width scan.',
            'source_manifest': json.loads((ROOT/'SOURCE_MANIFEST.json').read_text()),
            'states': len(states), 'character_sector_dimension': len(h), 'lags': list(LAGS),
            'definitions': {
                'generator': 'F_eta=(1+eta)J+(1-eta)D; Fprime=J-D exactly',
                'stationary_derivative': 'F0*pi_prime=-Fprime*pi0, sum(pi_prime)=0',
                'correlation': 'C(t)=conj(Z* B exp(-H*t) Z); B=Q*diag(pi)Q; H=-Q*F.T Q',
                'normalized_derivative': 'Uprime=C0^-1*(Cprime-C0prime*U)',
                'zero_frequency': 'R=int_0^infinity U(t)dt; Rprime=C0^-1*(conj[Z*(Bprime H^-1-B H^-1 Hprime H^-1)Z]-C0prime R)',
                'split': 'Stationary reweighting pi_prime plus physical generator derivative Hprime; both share one deterministic model, not independent evidence.',
                'precision': 'Exact derivative equations evaluated in float64; not rational-arithmetic certification.'},
            'checks': {'stationary_residual': stationary_residual, 'stationary_tangent_residual': tangent_residual,
                       'stationary_derivative_sum': float(dpi.sum()), 'stationary_K_odd_residual': maxabs(dpi[complement]+dpi),
                       'old_baseline_max_difference': baseline_diff,
                       'same_ray_first_derivative_max': maxabs(derivatives[:, (0, 1), (0, 1)]),
                       'Uprime_zero_lag_max': maxabs(derivatives[0]),
                       'spectral_vs_Frechet_max_difference': spectral_error,
                       'spectral_vs_resolvent_integral_max_difference': integral_error},
            'stationary_probability': pi.tolist(), 'stationary_derivative': dpi.tolist(),
            'integrals': encode_dict(zero),
            'lag_results': [{ 'lag': lag, **encode_dict(packet)} for lag, packet in zip(LAGS, full_packets)],
            'crossings': fixed_lag_crossings(full_at),
            'finite_quarter_secant_comparison': {
                'saved_secant_re_im': encode(secant),
                'secant_minus_exact_derivative_re_im': encode(secant-derivatives),
                'max_absolute_difference': maxabs(secant-derivatives),
                'boundary': 'Reuses old revealed finite outputs solely to measure finite-step nonlinearity, not a new experiment.'},
            'compressed': compressed, 'spectral_truncation': spectrum,
            'elapsed_seconds': time.perf_counter()-start,
            'environment': {'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__,
                            'machine': platform.machine(), 'hostname': platform.node(),
                            'threads': {k:os.environ.get(k) for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS')}},
            'code_sha256': {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.glob('*.py')}}


def write_report(r, path):
    lines = ['# P398 eta=0 analytic response and zero-frequency susceptibility', '',
             'Analytic derivative equations, evaluated in float64 on the frozen 1430-state/186-sector finite model. No new eta points or MC.', '',
             '| lag | Uprime minus-plus | Uprime plus-minus | pi-only plus-minus | generator-only plus-minus |',
             '|---:|---:|---:|---:|---:|']
    for row in r['lag_results']:
        u, p, g = [decode(row[x]) for x in ('Uprime','Uprime_pi','Uprime_generator')]
        lines.append(f"| {row['lag']:g} | {u[0,1].real:.12g} | {u[1,0].real:.12g} | {p[1,0].real:.12g} | {g[1,0].real:.12g} |")
    lines += ['', 'Zero-frequency integrated Uprime (rows/columns minus, plus):', '', '```',
              str(decode(r['integrals']['integrated_Uprime'])), '```', '',
              'Stationary-only contribution:', '```', str(decode(r['integrals']['integrated_Uprime_pi'])), '```',
              'Generator-only contribution:', '```', str(decode(r['integrals']['integrated_Uprime_generator'])), '```', '',
              '| model | zero-freq minus-plus | zero-freq plus-minus | max Uprime error on old lags | crossing plus-minus |',
              '|---|---:|---:|---:|---|']
    for row in r['compressed']:
        v=decode(row['integrals']['integrated_Uprime'])
        lines.append(f"| {row['name']} | {v[0,1].real:.12g} | {v[1,0].real:.12g} | {row['Uprime_max_absolute_error']:.6g} | {row['crossings']['plus_minus']} |")
    lines += ['', '| slow modes requested per ray | retained counts | zero-freq minus-plus | zero-freq plus-minus | max Uprime error |',
              '|---:|---|---:|---:|---:|']
    for row in r['spectral_truncation']['rows']:
        v=decode(row['integrated_Uprime_re_im'])
        lines.append(f"| {row['requested_modes_per_ray']} | {row['retained_per_ray']} | {v[0,1].real:.12g} | {v[1,0].real:.12g} | {row['Uprime_max_absolute_error']:.6g} |")
    lines += ['', 'Full derivative crossings: '+str(r['crossings']), '',
              'Checks: '+json.dumps(r['checks']), '',
              f"Actual analysis wall time: {r['elapsed_seconds']:.6f} seconds; environment: {r['environment']}", '',
              'This is finite-process response, not a continuum-field or square-site identification. The old geometric models consume pi and pi-prime from the full model. Modal budgets are signed and need not converge monotonically.', '']
    path.write_text('\n'.join(lines))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'results')
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    if (args.output/'latest.json').exists():
        raise RuntimeError('Existing result: do not overwrite or repeat')
    result=analyze()
    (args.output/'latest.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    write_report(result,args.output/'REPORT.md')
    print(json.dumps({'completed':str(args.output),'elapsed_seconds':result['elapsed_seconds'],
                      'checks':result['checks'],'crossings':result['crossings'],
                      'integrated_Uprime':result['integrals']['integrated_Uprime']},indent=2),flush=True)


if __name__=='__main__':
    main()

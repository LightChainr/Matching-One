#!/usr/bin/env python3
"""Fixed width-8 join/detach rate intervention; no sampling or parameter scan."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import time
from pathlib import Path

import numpy as np
import scipy
from scipy import linalg, sparse
from scipy.sparse.linalg import spsolve

from frozen_model import (noncrossing_states, join_adjacent, detach_state,
                          rotate_state, kreweras, features)

ROOT = Path(__file__).resolve().parent
ETAS = (0.0, 0.25, -0.25)
LAGS = (0.0, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0)
PHASE = np.exp(-1j*np.pi/4)
RAY_CHANGE = np.array([[1, 1], [-PHASE, PHASE]])/np.sqrt(2)


def encode(x):
    a = np.asarray(x)
    return np.stack((a.real, a.imag), axis=-1).tolist()


def decode(x):
    a = np.asarray(x)
    return a[..., 0]+1j*a[..., 1]


def maxabs(x):
    return float(np.max(np.abs(x)))


def sparse_maxabs(x):
    return float(np.max(np.abs(x.data))) if x.nnz else 0.0


def construct(width):
    states = noncrossing_states(width)
    index = {s: i for i, s in enumerate(states)}
    count = len(states)
    parts = []
    for action in (join_adjacent, detach_state):
        rows, cols, values = [], [], []
        for x, state in enumerate(states):
            for site in range(width):
                rows.extend((index[action(state, site)], x))
                cols.extend((x, x))
                values.extend((1, -1))
        part = sparse.coo_matrix((values, (rows, cols)), shape=(count, count),
                                 dtype=np.int64).tocsc()
        part.eliminate_zeros()
        parts.append(part)
    weight = (1j)**np.arange(width)
    f = np.array([[sum(weight[j]*(s[j] == s[(j+1) % width]) for j in range(width)),
                   sum(weight[j]*(s.count(s[j]) == 1) for j in range(width))]
                  for s in states])
    t2 = np.array([sum(weight[j]*(s.count(s[j]) == 2) for j in range(width))
                   for s in states])
    seen, orbits = set(), []
    for state in states:
        if state in seen:
            continue
        orbit, current = [], state
        while current not in orbit:
            orbit.append(current)
            current = rotate_state(current)
        seen.update(orbit)
        if len(orbit) % 4 == 0:
            orbits.append([index[s] for s in orbit])
    rows, cols, values = [], [], []
    for j, orbit in enumerate(orbits):
        for step, i in enumerate(orbit):
            rows.append(i); cols.append(j)
            values.append((1j)**step/np.sqrt(len(orbit)))
    q = sparse.csc_matrix((values, (rows, cols)), shape=(count, len(orbits)))
    source = np.asarray(q.conj().T @ f) @ RAY_CHANGE
    comp = np.array([index[kreweras(s)] for s in states])
    assert maxabs(q @ source-f @ RAY_CHANGE) < 1e-12
    assert sparse_maxabs(parts[0][comp, :][:, comp]-parts[1]) == 0
    assert sparse_maxabs(parts[1][comp, :][:, comp]-parts[0]) == 0
    assert np.array_equal(f[comp, 0], -1j*f[:, 1])
    assert np.array_equal(f[comp, 1], f[:, 0])
    return states, parts, q, source, f, t2, comp


def stationary(forward):
    equations = forward.tolil().astype(float)
    equations[-1, :] = 1
    rhs = np.zeros(forward.shape[0]); rhs[-1] = 1
    pi = spsolve(equations.tocsc(), rhs)
    residual = maxabs(forward @ pi)
    assert pi.min() > 0 and abs(pi.sum()-1) < 1e-11 and residual < 1e-11
    return pi, residual


def correlations(h, source, metric):
    # Match archived convention C_ij(t)=E[psi_i(X0) conj(psi_j(Xt))].
    return np.array([(source.conj().T @ metric @ (linalg.expm(-h*t) @ source)).conj()
                     for t in LAGS])


def analyze(width, output):
    start = time.perf_counter()
    states, parts, q, source, f, t2, comp = construct(width)
    expected = (1430, 186) if width == 8 else (14, 2)
    assert (len(states), q.shape[1]) == expected
    state_count = len(states)
    named, stages = {}, []
    if width == 8:
        named = features(states, f, t2)
        stages = json.loads((ROOT/'old_t4_compression.json').read_text())['stages']
    records, internal = [], {}
    for eta in ETAS:
        # Integer numerator is retained: G_eta has exactly denominator four.
        integer_forward = int(4*(1+eta))*parts[0]+int(4*(1-eta))*parts[1]
        assert maxabs(np.asarray(integer_forward.sum(axis=0))) == 0
        forward = integer_forward.astype(float)/4
        assert sparse.csgraph.connected_components(
            forward, directed=True, connection='strong', return_labels=False) == 1
        pi, stationary_residual = stationary(forward)
        mass = -forward.T.tocsr()
        h = (q.conj().T @ mass @ q).toarray()
        invariant_residual = maxabs(mass @ q-q @ h)
        assert invariant_residual < 1e-11
        metric = (q.conj().T @ q.multiply(pi[:, None])).toarray()
        c = correlations(h, source, metric)
        c0 = c[0]
        u = np.array([linalg.solve(c0, ct) for ct in c])
        scale = np.sqrt(np.outer(c0.diagonal().real, c0.diagonal().real))
        assert np.min(c0.diagonal().real) > 0
        mean = pi @ (f @ RAY_CHANGE)
        assert maxabs(mean) < 1e-11
        row = {'eta': eta, 'join_rate': 1+eta, 'detach_rate': 1-eta,
               'attempts_per_unit_time': 2*width,
               'stationary': pi.tolist(), 'stationary_residual': stationary_residual,
               'stationary_min': float(pi.min()), 'stationary_max': float(pi.max()),
               'mean_readouts_re_im': encode(mean),
               'exact_character_sector_dimension': q.shape[1],
               'sector_invariance_max_residual': invariant_residual,
               'C_re_im': encode(c), 'C_normalized_re_im': encode(c/scale),
               'U_re_im': encode(u),
               'C0_eigenvalues': linalg.eigvalsh(c0).tolist(),
               'compressed': []}
        for stage in stages:
            names = stage['exact_named_columns']['selected_columns_in_declaration_order']
            cols = [named[name[2:]][comp] if name.startswith('K_') else named[name]
                    for name in names]
            z = np.asarray(q.conj().T @ np.column_stack(cols))
            gram = z.conj().T @ metric @ z
            assert linalg.eigvalsh(gram).min() > 1e-12
            projected = linalg.solve(gram, z.conj().T @ metric @ h @ z,
                                    assume_a='her')
            source_coefficients = linalg.lstsq(z, source)[0]
            source_residual = maxabs(z @ source_coefficients-source)
            assert source_residual < 1e-11
            prediction = correlations(projected, source_coefficients, gram)
            prediction_u = np.array([linalg.solve(prediction[0], ct) for ct in prediction])
            difference = prediction-c
            row['compressed'].append({
                'name': stage['stage'], 'dimension': len(names), 'columns': names,
                'configuration_span_fixed': True, 'metric': 'current stationary L2(pi_eta)',
                'source_representation_max_residual': source_residual,
                'gram_condition_number': float(np.linalg.cond(gram)),
                'mass_re_im': encode(projected), 'C_re_im': encode(prediction),
                'C_normalized_re_im': encode(prediction/scale),
                'U_re_im': encode(prediction_u),
                'U_max_absolute_error': maxabs(prediction_u-u),
                'absolute_error_C0_scaled': np.abs(difference/scale).tolist(),
                'same_ray_relative_error_re_im': encode(
                    np.diagonal(prediction, axis1=1, axis2=2)
                    /np.diagonal(c, axis1=1, axis2=2)-1),
                'max_error_C0_scaled': maxabs(difference/scale)})
        if output:
            key = {0.0: 'zero', .25: 'plus_quarter', -.25: 'minus_quarter'}[eta]
            sparse.save_npz(output/f'forward_{key}_times4.npz', integer_forward)
            np.savez_compressed(output/f'character_i_{key}.npz',
                                mass=h, source=source, metric=metric)
        records.append(row)
        internal[eta] = {'c': c, 'u': u, 'h': h, 'metric': metric, 'pi': pi,
                         'integer_forward': integer_forward, 'scale': scale}
        print(f'eta={eta:+.2f}: stationary_residual={stationary_residual:.3g}, '
              f'full_states={state_count}, sector={q.shape[1]}', flush=True)
    zero, plus, minus = [internal[x] for x in ETAS]
    parity = np.array([[1, -1], [-1, 1]])
    parity_residual = maxabs(plus['c']-minus['c']*parity)
    u_parity_residual = maxabs(plus['u']-minus['u']*parity)
    pi_conjugacy = maxabs(plus['pi'][comp]-minus['pi'])
    generator_conjugacy = sparse_maxabs(
        plus['integer_forward'][comp, :][:, comp]-minus['integer_forward'])
    assert generator_conjugacy == 0
    assert pi_conjugacy < 1e-11 and parity_residual < 1e-9 and u_parity_residual < 1e-9
    assert maxabs(zero['c'][:, 0, 1]) < 1e-10
    compressed_parity = []
    for jp, jm in zip(records[1]['compressed'], records[2]['compressed']):
        residual = maxabs(decode(jp['C_re_im'])-decode(jm['C_re_im'])*parity)
        assert residual < 1e-8
        compressed_parity.append({'name': jp['name'], 'max_residual': residual})
    reference_checks = []
    if width == 8:
        old = next(r for r in json.loads((ROOT/'old_source_spectrum.json').read_text())['rows']
                   if r['width'] == 8)
        for sample in old['kernel_samples']:
            i = LAGS.index(sample['s'])
            old_c = RAY_CHANGE.T @ decode(sample['C_re_im']) @ RAY_CHANGE.conj()
            error = maxabs(zero['c'][i]-old_c)
            assert error < 1e-8
            reference_checks.append({'lag': sample['s'], 'C_max_absolute_difference': error})
        for current, archive in zip(records[0]['compressed'], stages):
            cc = decode(current['C_re_im'])
            for ray_index, ray in enumerate(archive['rays']):
                for sample in ray['samples']:
                    value = cc[LAGS.index(sample['t']), ray_index, ray_index]/cc[0, ray_index, ray_index]
                    error = abs(value-sample['u_compressed'])
                    assert error < 1e-8
            reference_checks.append({'compressed': current['name'], 'old_lags_reproduced': True})
    # Exact additive finite-change split: dynamics at old stationary weights,
    # then the change of stationary weights. This is not an infinitesimal claim.
    fixed_metric = correlations(plus['h'], source, zero['metric'])
    scale0 = zero['scale']
    contrasts = {
        'eta_denominator': .5,
        'finite_central_secant_re_im': encode((plus['c']-minus['c'])/.5),
        'finite_central_secant_C0_scaled_re_im': encode((plus['c']-minus['c'])/.5/scale0),
        'finite_central_secant_U_re_im': encode((plus['u']-minus['u'])/.5),
        'even_finite_change_re_im': encode((plus['c']+minus['c'])/2-zero['c']),
        'plus_total_change_re_im': encode(plus['c']-zero['c']),
        'plus_dynamics_with_pi0_change_re_im': encode(fixed_metric-zero['c']),
        'plus_stationary_weight_change_re_im': encode(plus['c']-fixed_metric),
        'boundary': 'Central secant at fixed eta=1/4, not the derivative at eta=0; '
                    'the two finite-change terms sum exactly but are not independent processes.'}
    return {
        'schema': 'matching-one/p398-fixed-attempt-rate-intervention/v1',
        'frozen_input_commit': '1f19fc1a2d9fc59dce650e95268c716762725985',
        'width': width, 'states': [list(s) for s in states],
        'generator': 'forward G_eta=(1+eta)sum(J-I)+(1-eta)sum(D-I), Q=1',
        'time': 'old continuous attempt-clock units; 2*width total Poisson attempts per unit time, '
                'including state-preserving events; no division by escape rate',
        'readouts': 'same A,L and psi_minus=(A-phase*L)/sqrt(2), '
                    'psi_plus=(A+phase*L)/sqrt(2), phase=exp(-i*pi/4)',
        'correlation': 'C_ij(t)=E_pi_eta[psi_i(X0)*conj(psi_j(Xt))]; rows/columns minus,plus',
        'normalization': 'C_ij/sqrt(C_ii(0)*C_jj(0)) at the same eta; raw C retained',
        'U_definition': 'U_eta(t)=inverse(C_eta(0))*C_eta(t), matching the old normalized transfer readout; U_eta(0)=I even when cross covariance is nonzero. No semigroup closure is assumed.',
        'lags': list(LAGS), 'etas': list(ETAS), 'rows': records,
        'symmetry': {'generator_times4_conjugacy_residual': generator_conjugacy,
                     'stationary_K_conjugacy_max_residual': pi_conjugacy,
                     'same_even_cross_odd_max_residual': parity_residual,
                     'U_same_even_cross_odd_max_residual': u_parity_residual,
                     'compressed_parity': compressed_parity},
        'old_reference_checks': reference_checks, 'contrasts': contrasts,
        'interpretation': 'Full finite width-8 connectivity model, using its exact invariant character-i sector; '
            'the old two 93-dimensional Kreweras rays are allowed to mix. Named 14/16 configuration '
            'function spans remain fixed, with a freshly computed stationary Galerkin metric. '
            'No fit, parameter search, width expansion, new MC, or CFT/field identity claim.',
        'elapsed_seconds': time.perf_counter()-start,
        'environment': {'python': platform.python_version(), 'numpy': np.__version__,
                        'scipy': scipy.__version__, 'machine': platform.machine(),
                        'thread_env': {k: os.environ.get(k) for k in
                                       ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS')}}}


def write_report(r, destination):
    rows = r['rows']; lines = ['# P398 固定总速率干预', '',
        '固定 width=8，η=0,±1/4；保留旧 A/L、ψ± 和全部旧 lags。每个 η 重新求完整平稳分布；'
        '全 1430 态模型的 character-i 扇区是完整 186 维，未把 ±ray 分开传播。', '',
        '原 14/16 个命名配置函数保持固定，仅按各自 πη 重新计算 Galerkin 内积。', '',
        '| lag | η=1/4 C−+/sqrt(C−−(0)C++(0)) | U−+ | full normalized C−− | full normalized C++ | 14 dim max scaled error | 16 dim max scaled error |',
        '|---:|---:|---:|---:|---:|---:|---:|']
    p = rows[1]; corr = decode(p['C_normalized_re_im'])
    transfer = decode(p['U_re_im'])
    for i, lag in enumerate(r['lags']):
        c = corr[i]
        errors = [np.max(np.asarray(x['absolute_error_C0_scaled'])[i]) for x in p['compressed']]
        lines.append(f'| {lag:g} | {c[0,1].real:.9g}{c[0,1].imag:+.9g}i | '
                     f'{transfer[i,0,1].real:.9g}{transfer[i,0,1].imag:+.9g}i | '
                     f'{c[0,0].real:.9g} | {c[1,1].real:.9g} | {errors[0]:.6g} | {errors[1]:.6g} |')
    lines += ['', f"共同 K 对称预言 same-ray 偶、cross-ray 奇；最大绝对误差 {r['symmetry']['same_even_cross_odd_max_residual']:.3g}。",
              f"平稳分布 K 共轭误差 {r['symmetry']['stationary_K_conjugacy_max_residual']:.3g}；整系数生成元共轭误差为零。", '',
              '零干预时旧完整响应及 14/16 子空间旧 lags 均通过复现检查。', '',
              '完整复数双向响应、平稳分布、生成元、同一旧读出的有限中心差商、总变化的 dynamics/π 权重分解均已保存。'
              '另保留旧 U(t)=C(0)^{-1}C(t)，其 t=0 交叉项为零，以免把静态交叉协方差当作传播。'
              '中心差商仅指 η=±1/4 的固定有限对比，不能称 η=0 的精确导数。', '',
              '本结果只涉及该有限宽 join/detach 连续时间模型。14/16 维比较的误差是实测确定性近似误差，'
              '不能自动解释成闭合、连续极限或原 square-site norm-4 的验证。', '',
              f"运行 {r['elapsed_seconds']:.3f} 秒；无 MC、参数扫描或 width 扩展。详细结果见 latest.json。", '']
    destination.write_text('\n'.join(lines))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT/'results')
    parser.add_argument('--validate-width4', action='store_true')
    args = parser.parse_args()
    if args.validate_width4:
        result = analyze(4, None)
        print(json.dumps({'width4_validation': result['symmetry'],
                          'elapsed_seconds': result['elapsed_seconds']}, indent=2))
        return
    args.output.mkdir(parents=True, exist_ok=True)
    if (args.output/'latest.json').exists():
        raise RuntimeError('Existing result: do not repeat or overwrite this fixed experiment')
    result = analyze(8, args.output)
    result['input_sha256'] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                              for p in ROOT.iterdir() if p.is_file()}
    with (args.output/'latest.json').open('x') as handle:
        json.dump(result, handle, indent=2, allow_nan=False); handle.write('\n')
    write_report(result, args.output/'REPORT.md')
    print(json.dumps({'completed': str(args.output), 'symmetry': result['symmetry'],
                      'elapsed_seconds': result['elapsed_seconds']}, indent=2))


if __name__ == '__main__':
    main()

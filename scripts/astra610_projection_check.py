#!/usr/bin/env python3
"""Independent #610 check: a moving affine chart does not commute with differencing.

Reconstruct the eight committed productions by incomplete-beta mixtures and
Brent roots (production uses binomial convolution and bisection). Use NumPy GLS
to expose the change of covector between first and second differences. No new
samples, no physical-exponent inference, no modification of old artifacts.
"""
from pathlib import Path
import csv
import json
import math
import numpy as np
from scipy.optimize import brentq
from scipy.special import betainc

ROOT = Path(__file__).resolve().parents[1]
LEVELS = np.arange(1, 10) / 10


def reconstruct(row, n):
    with (ROOT / row['source']).open() as f:
        records = list(csv.DictReader(f))
    ids = sorted({int(r['batch']) for r in records})
    batches = np.zeros((len(ids), 2, 2, n + 1), dtype=np.int64)
    index = {b: i for i, b in enumerate(ids)}
    for r in records:
        batches[index[int(r['batch'])], ('first', 'second').index(r['orientation']),
                ('minus', 'plus').index(r['kind']), int(r['k'])] += int(r['count'])
    totals = batches.sum(axis=0)
    ranks = np.arange(1, n + 1)

    def quantiles(hist):
        result = []
        for orientation in hist:
            mass = (orientation / orientation.sum(axis=1)[:, None]).mean(axis=0)[1:]
            result.append([brentq(lambda p: float(mass @ betainc(ranks, n + 1 - ranks, p)) - u,
                                 0., 1., xtol=5e-15) for u in LEVELS])
        return np.array(result)

    full = quantiles(totals)
    deleted = np.array([quantiles(totals - b) for b in batches])
    out = {}
    for name, w in [('spin0', list(row['spin0_weights'][k] for k in ('first', 'second'))),
                    ('equal', [.5, .5])]:
        q = np.array(w) @ full
        d = np.einsum('o,boj->bj', w, deleted)
        centered = d - d.mean(axis=0)
        cov = (len(ids) - 1) / len(ids) * centered.T @ centered
        out[name] = dict(q=q, cov=cov, deleted=d,
                         bias=(len(ids)-1)*(d.mean(axis=0)-q),
                         reconstruction_error=float(np.max(np.abs(q-row['quantiles_'+name]))))
    return out


def gls(y, q, cov, g=None):
    basis = np.column_stack([np.ones(9), q] + ([] if g is None else [g]))
    w = np.linalg.pinv(cov, rcond=1e-12, hermitian=True)
    covectors = np.linalg.solve(basis.T @ w @ basis, basis.T @ w)
    coefficients = covectors @ y
    return coefficients, y-basis@coefficients, covectors


def run():
    flow = json.loads((ROOT/'results/wasserstein-shape-flow/latest.json').read_text())
    old = json.loads((ROOT/'results/p582-amplitude-law/latest.json').read_text())
    loaded = {int(n): reconstruct(row, int(n)) for n, row in flow['sizes'].items()}
    report = {'source_ref': '8b5f9d1acb610fd82fc8df6310a96d0fd0ccdf64',
              'status': 'C2 same productions; deterministic algebra check, no new evidence block',
              'reconstruction': {}, 'weightings': {}}
    for n, both in loaded.items():
        report['reconstruction'][n] = {name: {
            'max_quantile_difference': v['reconstruction_error'],
            'max_jackknife_bias_estimate': float(np.max(np.abs(v['bias']))),
            'quantiles': v['q'].tolist(), 'covariance': v['cov'].tolist()
        } for name, v in both.items()}
    for weighting in ('spin0', 'equal'):
        data = {n: both[weighting] for n, both in loaded.items()}
        transitions = []
        for row in old['amplitudes']:
            b,t = row['base'], row['target']; h=math.log(t/b)
            y=(data[t]['q']-data[b]['q'])/h
            cov=(data[t]['cov']+data[b]['cov'])/h**2
            _,r,_=gls(y,data[b]['q'],cov)
            transitions.append(dict(base=b,target=t,h=h,y=y,cov=cov,r=r))
        unit=np.array([r['r']/np.linalg.norm(r['r']) for r in transitions])
        _,eig=np.linalg.eigh(unit.T@unit); g=eig[:,-1]
        if g@np.array(old['frozen_direction'])<0: g=-g
        direction_reconstruction_error=float(np.max(np.abs(g-old['frozen_direction']))) if weighting=='spin0' else None
        if weighting=='spin0':
            g=np.array(old['frozen_direction'])
        # Frozen means conditional on the observed full-data direction, not prospective.
        for r in transitions:
            r['coef'],r['residual'],r['covectors']=gls(r['y'],data[r['base']]['q'],r['cov'],g)
        law=old['weighting_systematic']['per_weighting'][weighting]
        omega,lam=law['omega'],law['scale_amplitude']
        checks=[]
        for i in (0,2):
            a,b=transitions[i:i+2]; n0,n1,n2=a['base'],a['target'],b['target']
            h0,h1=a['h'],b['h']; factor=2/(h0+h1)
            weights=np.array([factor/h0,-factor*(1/h0+1/h1),factor/h1])
            observation=sum(w*data[n]['q'] for w,n in zip(weights,(n0,n1,n2)))
            cov=sum(w*w*data[n]['cov'] for w,n in zip(weights,(n0,n1,n2)))
            co,_,ell=gls(observation,data[n1]['q'],cov,g); ell=ell[2]
            k=1+h0*a['coef'][1]
            naive=factor*(b['coef'][2]-a['coef'][2])
            transport=factor*(b['coef'][2]-a['coef'][2]/k)
            remainder=factor*(ell@b['residual']-ell@a['residual']/k)
            predicted=lam*sum(w*n**(-omega) for w,n in zip(weights,(n0,n1,n2)))
            pa=lam*(n1**(-omega)-n0**(-omega))/h0
            pb=lam*(n2**(-omega)-n1**(-omega))/h1
            fixed_first=[float(ell@a['y']),float(ell@b['y'])]
            checks.append(dict(sizes=[n0,n1,n2],published_style_curvature=float(co[2]),
                exponential_prediction=float(predicted),raw_ratio=float(co[2]/predicted),
                first_amplitude_difference=float(naive),first_amplitude_difference_ratio=float(naive/predicted),
                curvature_relative_error_bound_at_3p8_percent=float(.038*(abs(pa)+abs(pb))/abs(pb-pa)),
                affine_width_multiplier=float(k),width_effective_exponent=float(-math.log(k)/h0),
                transported_first_amplitudes=float(transport),transport_remainder=float(remainder),
                transport_identity_error=float(co[2]-transport-remainder),
                chart_contribution=float(transport-naive),
                fraction_of_raw_excess_due_to_chart=float((transport-naive)/(co[2]-predicted)),
                common_covector_first_amplitudes=fixed_first,
                common_covector_identity_error=float(co[2]-factor*(fixed_first[1]-fixed_first[0])),
                conditional_predicted_curvature_with_chart=float(factor*(pb-pa/k)),
                ratio_to_conditional_chart_prediction=float(co[2]/(factor*(pb-pa/k)))))
        report['weightings'][weighting]={'g':g.tolist(),
              'rebuilt_direction_difference_from_published':direction_reconstruction_error,
              'first_amplitudes':[float(r['coef'][2]) for r in transitions],
              'curvature':checks}
        p50=transitions[-1]
        k=1+p50['h']*p50['coef'][1]
        pa=lam*(290**(-omega)-145**(-omega))/math.log(2)
        pb=lam*(725**(-omega)-290**(-omega))/math.log(2.5)
        report['weightings'][weighting]['n725_conditional_prediction']={
            'known_145_to_290_width_multiplier':float(k),
            'predicted_290_to_725_first_amplitude':float(pb),
            'naive_curvature_prediction':float(2/math.log(5)*(pb-pa)),
            'chart_corrected_curvature_prediction':float(2/math.log(5)*(pb-pa/k)),
            'status':'conditional plug-in prediction; target not read; uncertainty and residual transport not scored'}
    return report


if __name__ == '__main__':
    report=run()
    out=ROOT/'results/astra610-independent/projection-check.json'
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'reconstruction_max':max(v[n]['max_quantile_difference']
                       for v in report['reconstruction'].values() for n in v),
                      'jackknife_bias_max':max(v[n]['max_jackknife_bias_estimate']
                       for v in report['reconstruction'].values() for n in v),
                      'checks':report['weightings']},indent=2))

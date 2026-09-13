#!/usr/bin/env python3
"""Recheck new shape residuals by a full blockwise nonlinear jackknife.

Only existing histograms are read. Each independent size is deleted separately;
matching batch indices across different seeds do not make a paired experiment.
"""
from __future__ import annotations
import json
import time
from pathlib import Path
import numpy as np
import shape_lineage_review as R
import threshold_quantile_lineage as L

ROOT=Path(__file__).resolve().parents[1]
SPECS={145:'results/server-20260829/P50-n145-n290-fullcurve/raw/n145_100m',
       290:'results/server-20260829/P50-n145-n290-fullcurve/raw/n290_100m',
       725:'results/server-20260907/P612-n725-fullcurve/raw/n725_100m'}
SLICES={'quadratic_chart_source':slice(0,3),'quadratic_chart_targets':slice(3,11),
        'constant_A_over_W':slice(11,19),'normal_form_145':slice(19,22),
        'normal_form_290':slice(22,25),'normal_form_725':slice(25,28),
        'widths':slice(28,31),'interval_norm_changes':slice(31,33),'beta':slice(33,34)}


def function(qall):
    qs=[qall[9*i:9*(i+1)] for i in range(3)]
    skew=[R.normalized_skew(q) for q in qs]
    out=list(R.common_quadratic_residual(qall))
    out.extend(skew[1]-skew[0]);out.extend(skew[2]-skew[1])
    for q in qs:out.extend(R.normal_form_residual(q))
    for q in qs:out.append(R.chart(q)[0])
    for part,indices in ((1,R.ZI),(2,R.AI)):
        vectors=[R.chart(q)[part][indices] for q in qs]
        d1=vectors[1]-vectors[0];d2=vectors[2]-vectors[1]
        out.append(np.sqrt(np.sum(d2*d2))-np.sqrt(np.sum(d1*d1)))
    out.append(R.beta_from_source(qs[0]))
    return np.asarray(out)


def jack_covariance(deleted):
    deleted=np.asarray(deleted);b=len(deleted)
    if b<2:raise ValueError('at least two batches required')
    centered=deleted-deleted.mean(axis=0)
    return (b-1)/b*(centered.T@centered)


def covariance_comparison(delta,jack):
    lower=np.linalg.cholesky(delta)
    tmp=np.linalg.solve(lower,jack)
    whitened=np.linalg.solve(lower,tmp.T).T
    eig=np.linalg.eigvalsh((whitened+whitened.T)/2)
    return {'relative_generalized_eigenvalues_min_max':[float(eig[0]),float(eig[-1])],
            'max_relative_diagonal_difference':float(np.max(abs(np.diag(jack)/np.diag(delta)-1)))}


def main():
    start=time.perf_counter();reference=json.loads(R.INPUT.read_text())
    qs=[];deleted=[];covs=[];seeds=[];metadata=[]
    for n in R.SIZES:
        prefix=ROOT/SPECS[n];meta=json.loads(Path(str(prefix)+'.metadata.json').read_text())
        seeds.append(meta['seed']);t=time.perf_counter()
        loaded=L.load_batch_histograms(Path(str(prefix)+'.hist.csv'))
        if loaded['n']!=n:raise ValueError('site count mismatch')
        weights=L.spin_zero_weights(loaded['orientation_cos4theta'])
        j=L.jackknife_quantiles(loaded,weights,L.FROZEN_LEVELS)
        q=np.array(j['full']);dq=np.array(list(j['deleted'].values()))
        old=np.array(reference['sizes'][str(n)]['analysis']['spin0']['Q_pooled'])
        error=float(np.max(abs(q-old)))
        if error>1e-11:raise ValueError('pooled Q differs from reviewed source')
        qs.append(q);deleted.append(dq);covs.append(jack_covariance(dq))
        metadata.append({'N':n,'seed':meta['seed'],'batches':len(dq),'Q_reproduction_max_error':error,
                         'elapsed_seconds':time.perf_counter()-t})
        print(f'completed existing N={n}, batches={len(dq)}, pooled error={error:.3g}',flush=True)
    if len(set(seeds))!=3:raise ValueError('independent-size contract requires distinct recorded streams')
    qall=np.concatenate(qs);value=function(qall)
    covariance=np.zeros((len(value),len(value)));bias=np.zeros(len(value))
    for index,dq in enumerate(deleted):
        outputs=[]
        for row in dq:
            perturbed=qall.copy();perturbed[9*index:9*(index+1)]=row
            outputs.append(function(perturbed))
        covariance+=jack_covariance(outputs)
        bias+=(len(dq)-1)*(np.mean(outputs,axis=0)-value)
    _,delta,_=R.propagated(function,qall,R.block_diagonal(covs))
    result={'schema':'matching-one.shape-nonlinear-jackknife.v1',
            'standing':'C2 existing data; separate deletion by independent size, with shared source beta; no Monte Carlo',
            'inputs':metadata,'analysis':{}}
    for name,s in SLICES.items():
        vals=value[s];c=covariance[s,s];d=delta[s,s]
        row={'values':vals.tolist(),'standard_errors':np.sqrt(np.diag(c)).tolist(),
             'jackknife_bias_estimate':bias[s].tolist(),'covariance':c.tolist(),
             'delta_covariance_comparison':covariance_comparison(d,c)}
        if name not in ('widths','interval_norm_changes','beta'):
            row['nominal_zero_test']=R.zero_test(vals,c)
        result['analysis'][name]=row
    result['elapsed_seconds']=time.perf_counter()-start
    result['limits']=['Nonlinear jackknife checks propagation, not exact finite-sample tail calibration.',
                      'No higher polynomial is selected after the quadratic-model failure.',
                      'Raw histogram inputs are already committed; these are not new independent samples.']
    print(json.dumps(result,indent=2,allow_nan=False))

if __name__=='__main__':main()

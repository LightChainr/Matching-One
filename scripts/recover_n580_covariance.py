#!/usr/bin/env python3
"""Reanalysis of the existing N580 shards. No simulation or data overwrite.

Uses the aligned delete-one ordering recorded by the #577 replay, checking
its seed/offset/batch metadata and previously published covariance entries.
Outputs a new result to stdout. All model scores are conditional on the
stated angular/readout contract; this is not new independent evidence.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
from projective_inference import ray_residual,subspace_residual

CHANNEL='P4_S_prime'
RUNGS=(1,2,4)


def jackknife_covariance(columns):
    n=len(columns[0])
    if n<2 or any(len(c)!=n for c in columns):
        raise ValueError('aligned columns with at least two deletions required')
    if any(not math.isfinite(x) for c in columns for x in c):
        raise ValueError('delete-one estimates must be finite')
    means=[math.fsum(c)/n for c in columns]
    return [[(n-1)/n*math.fsum((columns[i][k]-means[i])*(columns[j][k]-means[j])
                              for k in range(n)) for j in range(len(columns))]
            for i in range(len(columns))]


def recover(root):
    directory=root/'results'/'aspect-ladder-n580'
    shards=[json.loads((directory/'shards'/f'rung_r{r}.json').read_text()) for r in RUNGS]
    for key in ('seed','replica_offset','samples','batches'):
        if len({str(s[key]) for s in shards})!=1:
            raise ValueError(f'unaligned replay metadata: {key}')
    n=shards[0]['batches']
    columns=[[row[CHANNEL] for row in s['measured']['_deleted']] for s in shards]
    if any(len(c)!=n for c in columns):
        raise ValueError('incomplete delete-one array')
    y=[s['measured']['channels'][CHANNEL]['value'] for s in shards]
    S=jackknife_covariance(columns)
    for i,s in enumerate(shards):
        expected=s['measured']['channels'][CHANNEL]['standard_error']**2
        if not math.isclose(S[i][i],expected,rel_tol=1e-8,abs_tol=1e-20):
            raise ValueError('jackknife diagonal does not reproduce saved standard error')
    old=json.loads((directory/'latest.json').read_text())
    for entry,j in (('r2_over_r1',1),('r4_over_r1',2)):
        item=old['ratios'][entry][CHANNEL]
        if not math.isclose(S[0][j],item['covariance'],rel_tol=1e-7,abs_tol=1e-20):
            raise ValueError('paired covariance does not reproduce #577')
    reference=json.loads((root/'results'/'aspect-ladder-n580-projective'/'latest.json').read_text())
    if any(not math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-15)
           for a,b in zip(y,reference['response_vector'])):
        raise ValueError('response vector changed')
    scores={}
    leak=1148/21025
    for name,item in reference['competitors'].items():
        v=item['ray']
        score=ray_residual(y,S,v)
        score['excluded_at_nominal_3_sigma']=score['p_value']<math.erfc(3/math.sqrt(2))
        # Hypothetical common-modulus H8 mixing: nuisance b=a*rho.
        w=[-leak*v[0],leak*v[1],-leak*v[2]]
        nuisance=subspace_residual(y,S,[v,w])
        a,b=nuisance['amplitudes']
        bounded={}
        for B in (.2,1.0):
            inside=abs(b)<=B*abs(a)
            edge=min(ray_residual(y,S,[v[i]+sign*B*w[i] for i in range(3)])['statistic']
                     for sign in (-1,1))
            bounded[str(B)]=nuisance['statistic'] if inside else edge
        scores[name]={'ray':v,'full_covariance_ray_test':score,
                      'unbounded_common_shape_H8_plane':nuisance,
                      'bounded_common_rho_min_D_no_chisquare_calibration':bounded}
    c=[1/3,-1/2,1/6]
    value=math.fsum(c[i]*y[i] for i in range(3))
    var=math.fsum(c[i]*S[i][j]*c[j] for i in range(3) for j in range(3))
    return {'schema':'matching-one.n580-complete-covariance.20260912',
            'standing':'existing-block reanalysis, not new evidence',
            'alignment':'#577 aligned delete-one convention; metadata and saved covariance checked',
            'rungs':list(RUNGS),'batches':n,'response_vector':y,'covariance':S,
            'correlation_r2_r4':S[1][2]/math.sqrt(S[1][1]*S[2][2]),
            'curvature':{'value':value,'standard_error':math.sqrt(var),'z':value/math.sqrt(var)},
            'competitors':scores,
            'boundaries':['Gaussian-reference p-values, estimated covariance; not exact finite-sample coverage',
                          'H8 plane assumes a common ratio and the same modulus shape; sensitivity only',
                          'bounded-cone residuals have no ordinary chi-square calibration here',
                          'nonfinite equivalent sigma is emitted as null, not infinite evidence']}


def finite_json(value):
    if isinstance(value,float) and not math.isfinite(value):
        return None
    if isinstance(value,dict):
        return {k:finite_json(v) for k,v in value.items()}
    if isinstance(value,list):
        return [finite_json(v) for v in value]
    return value


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    args=parser.parse_args()
    print(json.dumps(finite_json(recover(args.root)),indent=2,allow_nan=False))

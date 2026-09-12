#!/usr/bin/env python3
"""Bounded rung-specific H8 sensitivity, not a common-ratio fit or fresh experiment.

For fixed positive v, fit mu_i=a*v_i*(1+leak_i*rho_i), |rho_i|<=B,
with independently varying rho_i and a of either sign. Each sign is a cone
on eight box vertices. Conic Caratheodory in R^3 permits <=3 active rays;
enumerating all such supports gives the global projection onto each cone.
The reported reference is intersection with ONE 3-d mean ellipsoid, not a
chi-square calibration for the fitted cone. Estimated covariance makes that
reference nominal/asymptotic. No stochastic resampling, simulation, or refit of v.
"""
from __future__ import annotations
import argparse
import itertools
import json
import math
from pathlib import Path
import numpy as np
from mpmath import mp

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT/'results/research-control-20260912/n580-complete-covariance-summary.json'
LEAK = np.array([-1., 1., -1.])*(1148/21025)
SIGNS = np.array(list(itertools.product((-1,1), repeat=3)))
SUPPORTS = [s for k in (1,2,3) for s in itertools.combinations(range(8),k)]


def cone_projection(y, covariance, v, bound, leak=LEAK):
    """Global numerical projection; return both sign branches for KKT checking."""
    y,s,v,leak = (np.asarray(x,dtype=float) for x in (y,covariance,v,leak))
    if y.shape!=(3,) or s.shape!=(3,3) or v.shape!=(3,) or leak.shape!=(3,):
        raise ValueError('three rungs required')
    if not all(np.isfinite(x).all() for x in (y,s,v,leak)) or np.any(v<=0):
        raise ValueError('finite inputs and positive model shape required')
    if not math.isfinite(bound) or bound<0 or bound*np.max(np.abs(leak))>=1:
        raise ValueError('bound must preserve the sign of each shape entry')
    if not np.allclose(s,s.T,rtol=1e-13,atol=0):
        raise ValueError('covariance must be symmetric')
    L=np.linalg.cholesky(s)
    wy=np.linalg.solve(L,y)
    vertices=(v*(1+bound*SIGNS*leak)).T
    positive=np.linalg.solve(L,vertices)
    branches=[]
    for sign in (1,-1):
        G=sign*positive
        lengths=np.linalg.norm(G,axis=0)
        Gn=G/lengths
        best={'distance':float(wy@wy),'active':[],'sign':sign}
        for active in SUPPORTS:
            X=Gn[:,active]
            coeff,_,rank,_=np.linalg.lstsq(X,wy,rcond=1e-12)
            if rank!=len(active) or np.min(coeff)<-1e-12*max(1,np.linalg.norm(wy)):
                continue
            coeff=np.maximum(coeff,0)
            residual=wy-X@coeff
            d=float(residual@residual)
            if d<best['distance']:
                best={'distance':d,'active':list(active),'sign':sign}
        branches.append(best)
    best=min(branches,key=lambda b:b['distance'])
    return {'distance':best['distance'],'sign':best['sign'],'branches':branches}


def high_precision_check(y,covariance,v,bound,projection,leak=LEAK):
    """Re-solve chosen supports and check convex KKT on ALL eight rays, both signs.

    This is a numerical optimality check at 60 dps, not an exact certificate.
    Uses the same decimal inputs as float enumeration, not rounded fitted means.
    """
    with mp.workdps(60):
        vec=lambda z:mp.matrix([mp.mpf(str(t)) for t in z])
        y=vec(y); s=mp.matrix([[mp.mpf(str(t)) for t in row] for row in covariance]); W=s**-1
        vv=vec(v); ll=vec(leak); B=mp.mpf(str(bound)); checked=[]
        for branch in projection['branches']:
            G=mp.matrix(3,8)
            for j,signs in enumerate(SIGNS):
                for i in range(3):
                    G[i,j]=branch['sign']*vv[i]*(1+B*ll[i]*int(signs[i]))
            active=branch['active']
            if active:
                X=mp.matrix([[G[i,j] for j in active] for i in range(3)])
                coeff=mp.lu_solve(X.T*W*X,X.T*W*y)
                mu=X*coeff
                if min(coeff)<-mp.mpf('1e-40'):
                    raise ArithmeticError('selected cone coefficients are negative')
            else:
                coeff=mp.zeros(0,1); mu=mp.zeros(3,1)
            q=W*(y-mu); gradient=G.T*q
            scale=max(1,mp.norm(G)*mp.norm(q))
            kkt=max(mp.mpf(0),max(gradient))/scale
            complement=abs((mu.T*q)[0])/max(1,mp.norm(mu)*mp.norm(q))
            if max(kkt,complement)>mp.mpf('1e-35'):
                raise ArithmeticError('cone projection failed high-precision KKT check')
            distance=((y-mu).T*q)[0]
            if abs(distance-branch['distance'])>mp.mpf('1e-8')*max(1,distance):
                raise ArithmeticError('float and high-precision projections disagree')
            amplitude=branch['sign']*mp.fsum(coeff)
            rho=[(mu[i]/(amplitude*vv[i])-1)/ll[i] for i in range(3)] if amplitude else None
            checked.append({'distance':float(distance),'sign':branch['sign'],
                'amplitude':float(amplitude),'mean':[float(t) for t in mu],
                'rho':[float(t) for t in rho] if rho is not None else None,
                'kkt_violation_relative':float(kkt),'complementarity_relative':float(complement)})
        best=min(checked,key=lambda b:b['distance'])
        return {**{k:best[k] for k in ('distance','sign','amplitude','mean','rho')},
            'max_kkt_violation_relative':max(b['kkt_violation_relative'] for b in checked),
            'max_complementarity_relative':max(b['complementarity_relative'] for b in checked),
            'branch_distances':[b['distance'] for b in checked]}


def ellipsoid_cutoff():
    """3-d Gaussian mean ellipsoid with alpha equal to the two-sided 3-sigma tail."""
    with mp.workdps(50):
        alpha=mp.erfc(3/mp.sqrt(2)); lo=mp.mpf(0); hi=mp.mpf(100)
        for _ in range(180):
            mid=(lo+hi)/2
            if mp.gammainc(mp.mpf('1.5'),mid/2,mp.inf,regularized=True)>alpha:
                lo=mid
            else:
                hi=mid
        return float(alpha),float((lo+hi)/2)


def analyze(payload):
    y,s=payload['response_vector'],payload['covariance']
    alpha,cutoff=ellipsoid_cutoff(); results={}
    for name,source in payload['competitors'].items():
        v=source['ray']; rows={}
        for bound in (0.,.2,1.):
            proj=cone_projection(y,s,v,bound)
            checked=high_precision_check(y,s,v,bound,proj)
            rows[str(bound)]={**checked,'ellipsoid_intersects':checked['distance']<=cutoff}
        if rows['0.0']['distance']<=cutoff:
            interval=[0.,0.]
        else:
            lo=0.; hi=(1-1e-8)/max(abs(LEAK))
            if cone_projection(y,s,v,hi)['distance']>cutoff:
                raise ArithmeticError('upper bound does not reach the ellipsoid')
            for _ in range(38):
                mid=(lo+hi)/2
                if cone_projection(y,s,v,mid)['distance']>cutoff:
                    lo=mid
                else:
                    hi=mid
            interval=[lo,hi]
            endpoint_D=[high_precision_check(y,s,v,b,cone_projection(y,s,v,b))['distance'] for b in interval]
            if not endpoint_D[0]>cutoff>=endpoint_D[1]:
                raise ArithmeticError('threshold bracket failed high-precision recheck')
        z=np.asarray(y)/np.asarray(v)
        # Closed-form smallest uniform bound that interpolates the observed mean.
        point_bound=(max(z)-min(z))/(max(z)+min(z))/max(abs(LEAK)) if np.all(z>0) else None
        results[name]={'shape':v,'fixed_bounds':rows,
            'first_ellipsoid_intersection_B_bracket':interval,
            'point_interpolation_B':float(point_bound) if point_bound is not None else None}
    return {'schema':'matching-one.n580-rungwise-leakage.v1',
        'source_path':str(SOURCE.relative_to(ROOT)),
        'input_mean':y,'input_covariance':s,
        'standing':'C2 same-block retrospective sensitivity; not a new experiment or amended freeze',
        'model':'mu_i=a*v_i*(1+leak_i*rho_i), signed a, independent |rho_i|<=B',
        'leakage':LEAK.tolist(),
        'reference':{'alpha':alpha,'mean_dimension':3,'ellipsoid_cutoff':cutoff,
            'calibration':'Known Gaussian covariance: simultaneous mean region; estimated jackknife covariance: nominal/asymptotic. NOT a fitted-cone chi-square law.'},
        'results':results}


def summarize(result):
    """Compact generated artifact; --summary changes storage, not the fit."""
    rows=result['results']; checks=[f for row in rows.values() for f in row['fixed_bounds'].values()]
    result={k:v for k,v in result.items() if k!='results'}
    result['fixed_bounds']=[0.,.2,1.]
    result['validation']={
        'precision_dps':60,
        'max_kkt_violation_relative':max(f['max_kkt_violation_relative'] for f in checks),
        'max_complementarity_relative':max(f['max_complementarity_relative'] for f in checks),
        'boundary_endpoints':'38-step numerical bisection; both endpoints re-solved and KKT checked',
        'scope':'Numerical checks, not exact rational certificates; no local shard replay.'}
    result['results']={name:{
        'shape':row['shape'],
        'minimum_D':[row['fixed_bounds'][str(b)]['distance'] for b in (0.,.2,1.)],
        'ellipsoid_intersects':[row['fixed_bounds'][str(b)]['ellipsoid_intersects'] for b in (0.,.2,1.)],
        'first_intersection_B_bracket':row['first_ellipsoid_intersection_B_bracket'],
        'point_interpolation_B':row['point_interpolation_B']} for name,row in rows.items()}
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,default=SOURCE)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--summary',action='store_true')
    parser.add_argument('--source-revision',help='Provenance annotation supplied by caller, not inferred')
    args=parser.parse_args(); result=analyze(json.loads(args.source.read_text()))
    try:
        result['source_path']=str(args.source.resolve().relative_to(ROOT))
    except ValueError:
        result['source_path']=str(args.source)
    if args.source_revision:
        result['source_revision_declared']=args.source_revision
    if args.summary:
        result=summarize(result)
    text=json.dumps(result,indent=2,allow_nan=False)+'\n'
    if args.output:
        # Never overwrite a committed or pre-existing result by default.
        args.output.parent.mkdir(parents=True,exist_ok=True)
        with args.output.open('x') as f:
            f.write(text)
    else:
        print(text,end='')

if __name__=='__main__':
    main()

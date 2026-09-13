#!/usr/bin/env python3
"""Necessary full-profile common-coordinate test, without fitting that coordinate.

A=Q145 o F290 and B=Q725 o F290. If one increasing psi makes all three
quantile laws location-scale copies, psi A psi^-1 and psi B psi^-1 are affine.
Their compositions AB and BA therefore differ by a constant in psi space:
onzero AB-BA cannot change sign on the common domain.
"""
from __future__ import annotations
from fractions import Fraction
import json
import time
import numpy as np
from scipy.stats import norm
from mpmath import mp
import threshold_quantile_lineage as L
import shape_lineage_review as R
from shape_lineage_nonlinear_jackknife import SPECS,jack_covariance
from shape_reflection_commutator import Law

LEVELS=(.3,.5,.7)


def composition_difference(a,b,p):
    return a(b(p))-b(a(p))


def statistic(laws,high=False):
    q=[v.mp_quantile if high else v.quantile for v in laws]
    f=laws[1].mp_cdf if high else laws[1].cdf
    visited=[]
    def transport(i,p):
        u=f(p);visited.append(float(u));return q[i](u)
    a=lambda p:transport(0,p)
    b=lambda p:transport(2,p)
    levels=[mp.mpf(str(u)) for u in LEVELS] if high else LEVELS
    points=[q[1](u) for u in levels]
    values=[composition_difference(a,b,p) for p in points]
    if min(visited)<.05 or max(visited)>.95:raise ValueError('composition outside central domain')
    return points,values,[min(visited),max(visited)]


def main():
    start=time.perf_counter();inputs=[];pools=[];weights=[];laws=[]
    for n in R.SIZES:
        data=L.load_batch_histograms(R.ROOT/(SPECS[n]+'.hist.csv'))
        pool=L.pooled_histograms(data['batches']);cs={}
        for name in L.ORIENTATIONS:
            a,b=data['orientation_representative'][name]
            cs[name]=Fraction(a**4-6*a*a*b*b+b**4,(a*a+b*b)**2)
        gap=cs['first']-cs['second'];w={'first':-cs['second']/gap,'second':cs['first']/gap}
        inputs.append(data);pools.append(pool);weights.append(w);laws.append(Law(pool,n,w))
    points,full,domain=statistic(laws)
    with mp.workdps(55):
        _,high,_=statistic(laws,True);high=np.array([float(v) for v in high])
    cov=np.zeros((3,3));bias=np.zeros(3)
    for index,n in enumerate(R.SIZES):
        out=[]
        for batch in inputs[index]['batches'].values():
            pool={name:{kind:L._subtract(pools[index][name][kind],batch[name][kind])
                        for kind in ('minus','plus')} for name in L.ORIENTATIONS}
            altered=laws.copy();altered[index]=Law(pool,n,weights[index])
            _,vals,limits=statistic(altered);out.append(vals)
            domain=[min(domain[0],limits[0]),max(domain[1],limits[1])]
        cov+=jack_covariance(out)
        bias+=(len(out)-1)*(np.mean(out,axis=0)-full)
        print(f'transport order completed existing N={n}',flush=True)
    se=np.sqrt(np.diag(cov));critical=float(norm.ppf(1-.0027/6))
    lo,hi=high-critical*se,high+critical*se
    positive=np.where(lo>0)[0].tolist();negative=np.where(hi<0)[0].tolist()
    report={'schema':'matching-one.shape-transport-order.v1',
            'standing':'C2, existing effective quantile laws; necessary condition for full common-coordinate location-scale collapse',
            'reference_levels':list(LEVELS),'reference_points':points,
            'composition_difference_high_precision':high.tolist(),'standard_errors':se.tolist(),
            'covariance':cov.tolist(),'jackknife_bias_estimate':bias.tolist(),
            'simultaneous_nominal_alpha':.0027,'bonferroni_critical_z':critical,
            'lower':lo.tolist(),'upper':hi.tolist(),'resolved_positive':positive,'resolved_negative':negative,
            'resolved_crossing':bool(positive and negative),
            'max_abs_marginal_z':float(np.max(abs(high)/se)),
            'pooled_float_vs_55dps_max':float(np.max(abs(high-full))),
            'composition_quantiles_range_all_deletions':domain,
            'elapsed_seconds':time.perf_counter()-start,
            'limits':['A significant nonzero composition difference alone is NOT a rejection; two affine maps may fail to commute.',
                      'Opposite resolved signs at fixed points would rule out common increasing affine conjugacy on that domain.',
                      'No resolved crossing is not proof of a common coordinate or collapse. Estimated nonlinear-jackknife marginal intervals have nominal, not exact coverage.',
                      'Reference points Q290(.3/.5/.7) remain random under N290 deletion; no new Monte Carlo or polynomial fit.']}
    print(json.dumps(report,indent=2,allow_nan=False))

if __name__=='__main__':main()

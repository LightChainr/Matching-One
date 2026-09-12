#!/usr/bin/env python3
"""First-order necessary common-chart test: the reflection curves cannot cross.

One simultaneous set of 12 nominal marginal intervals, no inverse covariance.
The third pair is redundant but retained conservatively in the union bound.
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

PAIRS=((0,1),(0,2),(1,2));LEVELS=(.3,.5,.7);ALPHA=.0027


def statistic(laws,high=False):
    q=[law.mp_quantile if high else law.quantile for law in laws]
    reflection=[law.mp_reflect if high else law.reflect for law in laws]
    levels=[mp.mpf(str(u)) for u in LEVELS] if high else LEVELS
    points=[q[1](u) for u in levels]
    medians=[f(mp.mpf('.5') if high else .5) for f in q]
    maps=[[f(p) for p in points] for f in reflection]
    return [v for i,j in PAIRS for v in [medians[i]-medians[j]]+[maps[i][k]-maps[j][k] for k in range(3)]]


def summarize(value,cov):
    value=np.asarray(value,float);se=np.sqrt(np.diag(cov))
    if min(se)<=0:raise ValueError('positive marginal errors required')
    critical=norm.ppf(1-ALPHA/(2*len(value)))
    low,high=value-critical*se,value+critical*se
    pairs=[]
    for row,(i,j) in enumerate(PAIRS):
        s=slice(4*row,4*(row+1))
        positive=np.where(low[s]>0)[0].tolist();negative=np.where(high[s]<0)[0].tolist()
        pairs.append({'sizes':[R.SIZES[i],R.SIZES[j]],'values':value[s].tolist(),
                      'standard_errors':se[s].tolist(),'low':low[s].tolist(),'high':high[s].tolist(),
                      'resolved_positive_coordinates':positive,'resolved_negative_coordinates':negative,
                      'incompatible_sign_order':bool(positive and negative)})
    return {'alpha_nominal':ALPHA,'critical_z_bonferroni':float(critical),
            'coordinates_per_pair':['median_difference','R_difference_at_Q290(.3)',
                                    'R_difference_at_Q290(.5)','R_difference_at_Q290(.7)'],
            'pairs':pairs,'any_resolved_crossing':any(p['incompatible_sign_order'] for p in pairs),
            'max_abs_marginal_z':float(np.max(abs(value)/se))}


def main():
    start=time.perf_counter();loaded=[];pools=[];weights=[];laws=[]
    for n in R.SIZES:
        item=L.load_batch_histograms(R.ROOT/(SPECS[n]+'.hist.csv'))
        pool=L.pooled_histograms(item['batches']);cos={}
        for name in L.ORIENTATIONS:
            a,b=item['orientation_representative'][name]
            cos[name]=Fraction(a**4-6*a*a*b*b+b**4,(a*a+b*b)**2)
        gap=cos['first']-cos['second'];w={'first':-cos['second']/gap,'second':cos['first']/gap}
        loaded.append(item);pools.append(pool);weights.append(w);laws.append(Law(pool,n,w))
    full=statistic(laws)
    with mp.workdps(55):
        high=np.array([float(x) for x in statistic(laws,True)])
    cov=np.zeros((12,12));bias=np.zeros(12)
    for index,n in enumerate(R.SIZES):
        outputs=[]
        for batch in loaded[index]['batches'].values():
            pool={name:{kind:L._subtract(pools[index][name][kind],batch[name][kind])
                        for kind in ('minus','plus')} for name in L.ORIENTATIONS}
            changed=laws.copy();changed[index]=Law(pool,n,weights[index])
            for law in changed:law.visited=[]
            outputs.append(statistic(changed))
            visited=[u for law in changed for u in law.visited]
            if min(visited)<.05 or max(visited)>.95:raise ValueError('central range exceeded')
        cov+=jack_covariance(outputs)
        bias+=(len(outputs)-1)*(np.mean(outputs,axis=0)-full)
        print(f'order test completed existing N={n}',flush=True)
    report=summarize(high,cov)
    report.update(schema='matching-one.shape-reflection-order.v1',
                  standing='C2 same effective spin0 laws and existing random blocks; necessary sign-order condition only',
                  covariance=cov.tolist(),jackknife_bias_estimate=bias.tolist(),
                  beta_vs_bernstein_max=float(np.max(abs(high-full))),
                  elapsed_seconds=time.perf_counter()-start,
                  limits=['Nominal marginal Gaussian intervals with estimated nonlinear-jackknife errors; not exact finite-sample coverage.',
                          'Absence of a resolved crossing does not imply common symmetrization or location-scale collapse.',
                          'The point set is the already declared .3/.5/.7, not a search over favorable locations.'])
    print(json.dumps(report,indent=2,allow_nan=False))

if __name__=='__main__':main()

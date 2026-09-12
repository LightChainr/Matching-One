#!/usr/bin/env python3
"""Common-chart commutator of the three EXISTING effective quantile laws.

Positive quantile weights admit one-level CDF inversion: if p=w1*x+w2*y,
solve F1(x)=F2(y), rather than nesting inverse CDFs. Beta mixtures provide
fast CDFs; a separate high-precision Bernstein path checks cancellation.
"""
from __future__ import annotations
import json
from fractions import Fraction
from pathlib import Path
import time
import numpy as np
from scipy.special import betainc
from scipy.optimize import brentq
from mpmath import mp
import threshold_quantile_lineage as L
import shape_lineage_review as R
from shape_common_chart import reflection_commutator
from shape_lineage_nonlinear_jackknife import SPECS,jack_covariance

LEVELS=(.3,.5,.7)


def bracket_root(f,lo,hi):
    for _ in range(85):
        mid=(lo+hi)/2
        if f(mid)<0:lo=mid
        else:hi=mid
    return (lo+hi)/2


class Component:
    def __init__(self,minus,plus,n):
        self.n=n
        tm,tp=sum(minus.values()),sum(plus.values())
        increments=[Fraction(minus.get(k,0),2*tm)+Fraction(plus.get(k,0),2*tp) for k in range(n+1)]
        if sum(increments)!=1:raise ValueError('threshold mixture is not normalized')
        self.ks=np.array([k for k in range(1,n+1) if increments[k]])
        self.weights=np.array([float(increments[k]) for k in self.ks])
        self.zero=float(increments[0]);self.cumulative=[];s=Fraction(0)
        for v in increments:s+=v;self.cumulative.append(s)
        self.mp_coeff=None
    def cdf(self,p):
        if p<=0:return self.zero
        if p>=1:return 1.
        return float(self.zero+self.weights@betainc(self.ks,self.n-self.ks+1,p))
    def quantile(self,u):
        if not 0<u<1:raise ValueError('interior quantile required')
        return brentq(lambda p:self.cdf(p)-u,0.,1.,xtol=5e-15,rtol=1e-15)
    def mp_cdf(self,p):
        if p<=0:return mp.mpf(self.cumulative[0].numerator)/self.cumulative[0].denominator
        if p>=1:return mp.mpf(1)
        if self.mp_coeff is None:
            self.mp_coeff=[mp.mpf(x.numerator)/x.denominator for x in self.cumulative]
        n=self.n;k=int(mp.floor(n*p))
        mass=mp.binomial(n,k)*p**k*(1-p)**(n-k)
        total=mass*self.mp_coeff[k];cur=mass
        for j in range(k,0,-1):
            cur*=mp.mpf(j)/(n-j+1)*(1-p)/p;total+=cur*self.mp_coeff[j-1]
        cur=mass
        for j in range(k,n):
            cur*=mp.mpf(n-j)/(j+1)*p/(1-p);total+=cur*self.mp_coeff[j+1]
        return total
    def mp_quantile(self,u):
        return bracket_root(lambda p:self.mp_cdf(p)-u,mp.mpf(0),mp.mpf(1))


class Law:
    def __init__(self,pool,n,weights):
        self.components=[Component(pool[name]['minus'],pool[name]['plus'],n) for name in L.ORIENTATIONS]
        raw=[weights[name] if isinstance(weights[name],Fraction) else Fraction(str(weights[name])) for name in L.ORIENTATIONS]
        total=sum(raw);self.exact_weights=[v/total for v in raw]
        self.weights=[float(v) for v in self.exact_weights]
        if min(self.weights)<=0 or abs(sum(self.weights)-1)>1e-12:
            raise ValueError('one-level inversion requires positive normalized weights')
        self.visited=[]
    def quantile(self,u):
        return sum(w*c.quantile(u) for w,c in zip(self.weights,self.components))
    def cdf(self,p):
        a,b=self.weights;f,g=self.components
        lo=max(0.,(p-b)/a);hi=min(1.,p/a)
        x=brentq(lambda x:f.cdf(x)-g.cdf((p-a*x)/b),lo,hi,xtol=5e-15,rtol=1e-15)
        return f.cdf(x)
    def reflect(self,p):
        u=self.cdf(p);self.visited.append(u)
        return self.quantile(1-u)
    def mp_quantile(self,u):
        return sum((mp.mpf(w.numerator)/w.denominator)*c.mp_quantile(u) for w,c in zip(self.exact_weights,self.components))
    def mp_cdf(self,p):
        a,b=[mp.mpf(w.numerator)/w.denominator for w in self.exact_weights];f,g=self.components
        lo=max(mp.mpf(0),(p-b)/a);hi=min(mp.mpf(1),p/a)
        x=bracket_root(lambda x:f.mp_cdf(x)-g.mp_cdf((p-a*x)/b),lo,hi)
        return f.mp_cdf(x)
    def mp_reflect(self,p):
        return self.mp_quantile(1-self.mp_cdf(p))


def statistic(laws,high=False,levels=LEVELS):
    if not high:
        for law in laws:law.visited=[]
    q=laws[1].mp_quantile if high else laws[1].quantile
    maps=[v.mp_reflect if high else v.reflect for v in laws]
    points=[q(mp.mpf(str(u)) if high else u) for u in levels]
    values=[reflection_commutator(*maps,p)['difference'] for p in points]
    return points,values


def main():
    start=time.perf_counter();inputs=[];laws=[];pools=[];weights=[]
    old=json.loads(R.INPUT.read_text())
    for n in R.SIZES:
        path=R.ROOT/(SPECS[n]+'.hist.csv');loaded=L.load_batch_histograms(path)
        pool=L.pooled_histograms(loaded['batches'])
        cos={}
        for name in L.ORIENTATIONS:
            a,b=loaded['orientation_representative'][name]
            cos[name]=Fraction(a**4-6*a*a*b*b+b**4,(a*a+b*b)**2)
        gap=cos['first']-cos['second']
        w={'first':-cos['second']/gap,'second':cos['first']/gap};law=Law(pool,n,w)
        qs=[law.quantile(u) for u in L.FROZEN_LEVELS]
        err=max(abs(np.array(qs)-old['sizes'][str(n)]['analysis']['spin0']['Q_pooled']))
        if err>2e-12:raise ValueError('beta path changes pooled Q')
        inputs.append(loaded);pools.append(pool);weights.append(w);laws.append(law)
        print(f'N={n} pooled Q control max={err:.3g}',flush=True)
    points,pooled=statistic(laws)
    involution=max(abs(law.reflect(law.reflect(p))-p) for law in laws for p in points)
    visited=[u for law in laws for u in law.visited]
    if min(visited)<.05 or max(visited)>.95:raise ValueError('composition left the declared central domain')
    with mp.workdps(55):
        _,mp_values=statistic(laws,True)
        mp_error=max(abs(float(v)-f) for v,f in zip(mp_values,pooled))
        high_primary=[float(v) for v in mp_values]
    covariance=np.zeros((3,3));checks=[]
    for index,n in enumerate(R.SIZES):
        outputs=[]
        for number,batch in enumerate(inputs[index]['batches'].values()):
            pool={name:{kind:L._subtract(pools[index][name][kind],batch[name][kind])
                        for kind in ('minus','plus')} for name in L.ORIENTATIONS}
            changed=laws.copy();changed[index]=Law(pool,n,weights[index])
            _,vals=statistic(changed);outputs.append(vals)
            interior=[u for law in changed for u in law.visited]
            if min(interior)<.05 or max(interior)>.95:raise ValueError('deleted composition left central domain')
            if number==0:
                with mp.workdps(55):
                    _,hv=statistic(changed,True,levels=(.5,))
                    checks.append(abs(float(hv[0])-vals[1]))
        covariance+=jack_covariance(outputs)
        print(f'N={n} separate-size nonlinear deletion complete',flush=True)
    se=np.sqrt(np.diag(covariance));numerical=max([mp_error,involution]+checks)
    adequate=bool(numerical<.01*min(se))
    out={'schema':'matching-one.shape-reflection-commutator.v1',
         'standing':'C2, existing effective spin0 quantile laws; separate-size deletion; no fitted chart or new samples',
         'sizes':list(R.SIZES),'reference_levels':list(LEVELS),'reference_p_values':points,
         'commutator_high_precision_pooled':high_primary,'standard_errors':se.tolist(),
         'covariance':covariance.tolist(),'composition_quantiles_range':[min(visited),max(visited)],
         'precision_controls':{'dps':55,'pooled_beta_vs_bernstein_max':mp_error,
                               'selected_deleted_beta_vs_bernstein':checks,
                               'float_involution_max':involution,
                               'max_control_error_over_min_se':float(numerical/min(se)),
                               'adequate_at_one_percent_se':adequate},
         'nominal_zero_test':R.zero_test(high_primary,covariance) if adequate else None,
         'elapsed_seconds':time.perf_counter()-start,
         'limits':['Nonzero rules out exact common symmetrization of these effective finite laws, not asymptotic symmetry or all underlying orientation-wise scaling fields.',
                   'Zero or unresolved does not prove a common chart; this is only a necessary condition.',
                   'High-precision spot checks are numerical controls, not interval certificates or exact finite-sample Gaussian coverage.']}
    print(json.dumps(out,indent=2,allow_nan=False))

if __name__=='__main__':main()

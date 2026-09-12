#!/usr/bin/env python3
"""Deterministic spectral controls for the thin-torus probability theorem.

No sampling and no exponent fit. Widths 2/3 use the already proved row formulas.
Width 4 reads the preceding delivery's width4-parametric-definition.json.
Numerical diagnostics use mpmath; exact endpoint coefficients use integer
polynomial arithmetic and do not depend on a numerical eigenvalue threshold.
"""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
import json
from math import comb
from pathlib import Path
import time
from mpmath import mp

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DEFINITION = ROOT/'results/research-control-20260912/width4-parametric-definition.json'


def p_add(*args):
    out = [0]*max(map(len,args))
    for a in args:
        for i,v in enumerate(a): out[i] += v
    while len(out)>1 and not out[-1]: out.pop()
    return out


def p_scale(a,c):
    return [c*x for x in a]


def p_mul(a,b):
    out = [0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): out[i+j] += x*y
    return p_add(out)


def p_eval(a,p):
    out = 0
    for c in reversed(a): out = out*p+c
    return out


def one_minus(a):
    out = [0]*len(a)
    for k,c in enumerate(a):
        for j in range(k+1): out[j] += c*comb(k,j)*(-1)**j
    return p_add(out)


def normalize_fugacity(co, width):
    """f((1+t)^w lambda,t)/(1+t)^(w deg f), t=p/(1-p)."""
    result = []
    for i,row in enumerate(co):
        out = [0]*(width*i+1)
        for k,c in enumerate(row):
            if not c: continue
            if k > width*i: raise ValueError('factor coefficient exceeds row-weight degree')
            for j in range(width*i-k+1):
                out[k+j] += c*comb(width*i-k,j)*(-1)**j
        result.append(p_add(out))
    return result


class WidthSpectrum:
    def __init__(self, width: int, definition: dict | None = None):
        if width not in (2,3,4): raise ValueError('controls support widths 2, 3, 4')
        self.width = width
        if width == 2:
            self.factors = {'open':[[1],[0,-1],[0,0,0,-1,1]],
                            'side':[[1],[0,-1,1]],'closed':[[1],[-1,0,1]]}
            self.weights = {'M':{'open':1,'side':1,'closed':-1},
                            'P0':{'closed':1,'side':-2},
                            'P2':{'open':1,'side':-1}}
            self.leading = ('closed','open')
        elif width == 3:
            u,v,w = [0,1,-2,1],[0,0,1,-1],[0,0,0,1]
            uv,uw = p_mul(u,v),p_mul(u,w)
            self.factors = {
                'open':[[1],p_scale(p_add(u,p_scale(v,3),w),-1),
                        p_scale(p_add(uv,p_scale(uw,2)),-1),p_mul(uv,w)],
                'standard':[[1],p_scale(u,-1),p_scale(uv,-1)],
                'proper':[[1],p_scale(p_add(u,p_scale(v,3)),-1),p_scale(uv,-1)],
                'closed':[[1],[-1,0,0,1]]}
            self.weights = {'M':{'open':1,'standard':2,'closed':-1},
                            'P0':{'closed':1,'proper':-1,'standard':-2},
                            'P2':{'open':1,'proper':-1}}
            self.leading = ('closed','open')
        else:
            if definition is None: raise ValueError('width 4 requires its existing definition')
            self.factors = {name:normalize_fugacity(value['coefficients_desc_x_ascending_t'],4)
                            for name,value in definition['factors'].items()}
            self.factors['positive_square'] = [[1],[0,0,-1,2,-1]]
            self.weights = {}
            for which in ('M','P0','P2'):
                weights = Counter()
                if which in ('M','P0'):
                    for name,n in definition['blocks']['5']['factorization']:
                        weights[name] += (-1 if which=='M' else 1)*n
                if which in ('M','P2'):
                    for name,n in definition['blocks']['15']['factorization']: weights[name] += n
                if which != 'M':
                    for name,n in definition['blocks']['16']['factorization']: weights[name] -= n
                    weights['positive_square'] += 2
                self.weights[which] = {name:c for name,c in weights.items() if c}
            self.leading = ('closed_2','open_5')
            self.crossing_poly = definition['crossing_polynomial_ascending_t']

    def coefficients(self,name,p):
        return [p_eval(a,p) for a in self.factors[name]]

    @staticmethod
    def trace_power(co,m):
        """Small companion-matrix power; independent of large transfer iteration."""
        n = len(co)-1
        if n == 1: return (-co[1])**m
        a = mp.matrix(n)
        for j in range(n): a[0,j] = -co[j+1]
        for i in range(1,n): a[i,i-1] = 1
        power = a**m
        return sum(power[i,i] for i in range(n))

    def value(self,which,m,p):
        if which not in self.weights or m < 2 or not 0 <= p <= 1:
            raise ValueError('invalid probability, length or sector')
        return sum(c*self.trace_power(self.coefficients(name,p),m)
                   for name,c in self.weights[which].items())

    def cdf(self,m,p):
        return (1+self.value('M',m,p))/2

    def cylinder_root(self):
        if self.width == 2:
            return mp.findroot(lambda p:2*p**3+2*p**2-1,('.55','.58'))
        if self.width == 3:
            return mp.findroot(lambda p:p**6-3*p**5-5*p**4-4*p**3+p+1,('.58','.60'))
        t = mp.findroot(lambda t:p_eval(self.crossing_poly,t),('1.44','1.45'))
        return t/(1+t)

    def leading_value(self,which,p):
        name = self.leading[which]
        co = self.coefficients(name,p)
        if len(co)==2: return -co[1]
        if which == 0:
            return (-co[1]+mp.sqrt(co[1]**2-4*co[2]))/2
        # The isolated Perron branch near the certified crossing.
        return mp.findroot(lambda x:mp.polyval(co,x),mp.mpf('.9'))

    def branch_data(self):
        q = self.cylinder_root()
        l0,l2 = [self.leading_value(j,q) for j in (0,1)]
        if abs(l0-l2) > mp.mpf('1e-50'): raise AssertionError('leading branches do not cross')
        rates = [mp.diff(lambda p:mp.log(self.leading_value(j,p)),q) for j in (0,1)]
        if not rates[0] < 0 < rates[1] or not 0 < l0 < 1:
            raise AssertionError('central-layer hypotheses failed')
        return q,l0,rates[0],rates[1]

    def endpoint_certificates(self):
        out=[]
        for which,name in enumerate(self.leading):
            co=self.factors[name]
            if which: co=[one_minus(a) for a in co]
            defect=p_add(*co)  # f(1,p), or f(1,1-v)
            order=next(i for i,c in enumerate(defect) if c)
            derivative=sum((len(co)-1-i)*row[0] for i,row in enumerate(co))
            expected=1 if which==0 else sum(comb(self.width,2*k)*comb(2*k,k)
                                           for k in range(self.width//2+1))
            if order != self.width or defect[order] != expected or derivative != 1:
                raise AssertionError('exact endpoint eigenvalue expansion failed')
            out.append({'endpoint':'p=0' if which==0 else 'p=1',
                        'leading_factor':name,'defect_polynomial_ascending_small_parameter':defect,
                        'first_nonzero_degree':order,'first_coefficient':defect[order],
                        'partial_lambda_at_endpoint':derivative,
                        'conclusion':f'lambda=1-{expected}*epsilon^{self.width}+O(epsilon^{self.width+1})'})
        return out


def fmt(x,digits=28):
    return mp.nstr(x,digits)


def numerical_report(definition: dict, dps: int = 80) -> dict:
    started=time.perf_counter()
    if dps < 70: raise ValueError('at least 70 decimal digits required')
    result=[]
    with mp.workdps(dps):
        for width in (2,3,4):
            model=WidthSpectrum(width,definition)
            q,lstar,a0,a2=model.branch_data()
            tw=sum(comb(width,2*k)*comb(2*k,k) for k in range(width//2+1))
            row={'width':width,'cylinder_crossing':fmt(q,50),'lambda_star':fmt(lstar,50),
                 'log_lambda0_derivative':fmt(a0),'log_lambda2_derivative':fmt(a2),
                 'hprime':fmt(a2-a0),'endpoint_certificates':model.endpoint_certificates()}
            tails=[]
            for m in (100,1000,10000,1000000):
                p_lo=mp.mpf(m)**(-mp.mpf(1)/width)
                vacancy=mp.mpf(tw*m)**(-mp.mpf(1)/width)
                flo=model.cdf(m,p_lo);fhi=model.cdf(m,1-vacancy)
                target_lo=(1-mp.exp(-1))/2;target_hi=(1+mp.exp(-1))/2
                tails.append({'length':m,'lower_test_p':fmt(p_lo),'upper_test_p':fmt(1-vacancy),
                              'F_lower':fmt(flo),'F_upper':fmt(fhi),
                              'F_lower_limit':fmt(target_lo),'F_upper_limit':fmt(target_hi),
                              'lower_absolute_error':fmt(abs(flo-target_lo)),
                              'upper_absolute_error':fmt(abs(fhi-target_hi))})
            row['two_endpoint_diagnostics']=tails
            if width==4:
                quantiles=[]
                for m in (32,256,65536,1000000):
                    values=[]
                    for u in (mp.mpf('.25'),mp.mpf('.75')):
                        lo,hi=mp.mpf('0'),mp.mpf('1')
                        for _ in range(170):
                            mid=(lo+hi)/2
                            if model.cdf(m,mid)<u: lo=mid
                            else: hi=mid
                        values.append((lo+hi)/2)
                    quantiles.append({'length':m,'Q_025':fmt(values[0]),'Q_075':fmt(values[1]),
                                      'interquartile_width':fmt(values[1]-values[0]),
                                      'separate_median_limit_q':fmt(q),
                                      'lower_to_Weibull_prediction_ratio':fmt(values[0]/(mp.log(2)/m)**(mp.mpf(1)/4)),
                                      'upper_gap_to_Weibull_prediction_ratio':fmt((1-values[1])/(mp.log(2)/(19*m))**(mp.mpf(1)/4))})
                row['fixed_quartile_diagnostics']=quantiles
            center=[]
            for m in (16,64,256):
                for x in (mp.mpf('-.5'),mp.mpf('0'),mp.mpf('.5')):
                    actual=model.value('M',m,q+x/m)/lstar**m
                    activity=(model.value('P0',m,q+x/m)+model.value('P2',m,q+x/m))/lstar**m
                    target=mp.exp(a2*x)-mp.exp(a0*x)
                    center.append({'length':m,'p_offset_times_length':fmt(x),
                                   'M_divided_by_lambda_star_power':(fmt(actual) if abs(actual)>mp.mpf(10)**(-dps+15) else None),
                                   'M_resolved_at_working_precision':abs(actual)>mp.mpf(10)**(-dps+15),
                                   'limiting_central_profile':fmt(target),
                                   'absolute_error':fmt(abs(actual-target)),
                                   'conditional_rank2_given_nonzero_X':fmt((1+actual/activity)/2),
                                   'conditional_logistic_limit':fmt(1/(1+mp.exp(-(a2-a0)*x)))})
            row['central_probability_window_diagnostics']=center
            tilted=[]
            for m in (16,64,256):
                for source in (mp.mpf('-.5'),mp.mpf('.5')):
                    target=-2*source/(a2-a0)
                    def residual(x):
                        p=q+x/m
                        return mp.log(model.value('P2',m,p)/model.value('P0',m,p))+2*source
                    root=mp.findroot(residual,(target-mp.mpf('.05'),target+mp.mpf('.05')))
                    tilted.append({'length':m,'topological_source_s':fmt(source),
                                   'source_balance_root':fmt(q+root/m),
                                   'm_times_root_minus_q':fmt(root),'limit':fmt(target),
                                   'absolute_error':fmt(abs(root-target)),
                                   'log_odds_equation_residual':fmt(abs(residual(root)),12)})
            row['finite_source_balance_roots']=tilted
            row['asymptotic_m_times_d_balance_root_d_source']=fmt(-2/(a2-a0))
            snapshots=[]
            for m in (8,32,128,512):
                p0=model.value('P0',m,q);p2=model.value('P2',m,q)
                matching=model.value('M',m,q)
                activity=p0+p2
                slope=mp.diff(lambda p:model.value('M',m,p),q)
                variance=activity-matching**2
                if min(p0,p2,variance,slope) <= 0 or activity>=1:
                    raise AssertionError('invalid fixed-p snapshot statistics')
                n95=mp.log(mp.mpf('.05'))/mp.log1p(-activity)
                # Large m: subtracting 1-activity in float would lose the event.
                snapshots.append({'length':m,'at_cylinder_q_not_exact_finite_median':True,
                                  'P0':fmt(p0),'P2':fmt(p2),'activity_E_X2':fmt(activity),
                                  'matching_M':(fmt(matching) if abs(matching)>activity*mp.mpf(10)**(-dps+15) else None),
                                  'matching_M_resolved_relative_to_activity':abs(matching)>activity*mp.mpf(10)**(-dps+15),
                                  'matching_slope':fmt(slope),
                                  'n_95percent_nonzero_snapshot_diagnostic_before_ceiling':fmt(n95,16),
                                  'local_inverse_mean_SE_times_sqrt_samples':fmt(mp.sqrt(variance)/slope),
                                  'inverse_CDF_slope_at_F_q':fmt(2/slope),
                                  'log10_activity':fmt(mp.log10(activity))})
            row['independent_fixed_p_snapshot_diagnostics']=snapshots
            source=[]
            # Large-deviation source s=m*sigma, evaluated near the certified crossing.
            for m in (16,64,256):
                p0=model.value('P0',m,q);p2=model.value('P2',m,q);p1=1-p0-p2
                for multiplier in (mp.mpf('-1.5'),mp.mpf('0'),mp.mpf('1.5')):
                    sigma=multiplier*(-mp.log(lstar))
                    logs=[mp.log(p0)-m*sigma,mp.log(p1),mp.log(p2)+m*sigma]
                    top=max(logs)
                    free=(top+mp.log(sum(mp.exp(z-top) for z in logs)))/m
                    expected=max(0,mp.log(lstar)-sigma,mp.log(lstar)+sigma)
                    source.append({'length':m,'scaled_source_sigma':fmt(sigma),
                                   'log_Z_over_length':fmt(free),
                                   'limit':fmt(expected),'absolute_error':fmt(abs(free-expected))})
            row['topological_source_large_deviation_diagnostics']=source
            # Precision check probes cancellation-prone and endpoint values, not a new theory check.
            probes=[(32,q),(1000000,mp.mpf(1000000)**(-mp.mpf(1)/width)),
                    (1000000,1-mp.mpf(tw*1000000)**(-mp.mpf(1)/width))]
            errors=[]
            for m,p in probes:
                reference=model.value('M',m,p)
                with mp.workdps(dps+30):
                    repeat=model.value('M',m,p)
                    errors.append(abs(repeat-reference))
            row['same_algorithm_extra_precision_max_absolute_difference']=fmt(max(errors),12)
            result.append(row)
    return {'schema':'matching-one.thin-torus-spectral-controls.v1','decimal_digits':dps,
            'standing':'Deterministic high-precision diagnostics of proved finite-width formulas; not exact interval results.',
            'widths':result,'elapsed_seconds':time.perf_counter()-started,
            'limits':['No Monte Carlo or new evidence block.',
                      'The all-fixed-width probability limits are proved geometrically, not extrapolated from widths 2/3/4.',
                      'Snapshot costs concern independent Bernoulli rank observations at one p only.',
                      'No cost lower bound for histogram reconstruction, Rao-Blackwellization, importance sampling or exact transfer methods.',
                      'Finite m snapshot rows are at q_w, not exactly at p_(w,m).',
                      'There is no width-uniform CFT/root-shift rate or pc value claim.']}


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--definition',type=Path,default=DEFAULT_DEFINITION)
    parser.add_argument('--out',type=Path)
    parser.add_argument('--dps',type=int,default=80)
    args=parser.parse_args()
    if not args.definition.exists():
        parser.error('missing preceding parametric delivery definition; supply --definition PATH')
    raw=args.definition.read_bytes();definition=json.loads(raw)
    report=numerical_report(definition,args.dps)
    report['input']={'repository_path':'results/research-control-20260912/width4-parametric-definition.json',
                     'sha256':hashlib.sha256(raw).hexdigest(),
                     'standing':'Supplied preceding delivery, not asserted merged into main.'}
    text=json.dumps(report,indent=2,allow_nan=False)+'\n'
    if args.out:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        with args.out.open('x',encoding='utf-8') as f:f.write(text)
    else: print(text,end='')

if __name__=='__main__':
    main()

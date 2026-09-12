#!/usr/bin/env python3
"""Rational interval certificate for the width-four crossing and all-m root sign.

mpmath supplies diagnostic approximations and proposes isolating bands. Every
load-bearing band/sign/inequality is checked afterward using fractions.Fraction.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
from mpmath import mp
from width4_parametric_traces import polynomial_value, RESULTS

class Interval:
    def __init__(self,lo,hi=None):
        self.lo=F(lo);self.hi=F(hi if hi is not None else lo)
        if self.lo>self.hi:raise ValueError('reversed interval')
    @staticmethod
    def of(x):return x if isinstance(x,Interval) else Interval(x)
    def __add__(self,other):
        o=self.of(other);return Interval(self.lo+o.lo,self.hi+o.hi)
    __radd__=__add__
    def __neg__(self):return Interval(-self.hi,-self.lo)
    def __sub__(self,other):return self+-self.of(other)
    def __rsub__(self,other):return self.of(other)+-self
    def __mul__(self,other):
        o=self.of(other);v=[self.lo*o.lo,self.lo*o.hi,self.hi*o.lo,self.hi*o.hi]
        return Interval(min(v),max(v))
    __rmul__=__mul__
    def __truediv__(self,other):
        o=self.of(other)
        if o.lo<=0<=o.hi:raise ValueError('interval division through zero')
        return self*Interval(1/o.hi,1/o.lo)
    def __rtruediv__(self,other):return self.of(other)/self
    def __pow__(self,n):
        if not isinstance(n,int) or n<0:raise ValueError('nonnegative integer power required')
        out=Interval(1)
        for _ in range(n):out=out*self
        return out
    def positive(self):return self.lo>0
    def negative(self):return self.hi<0
    def json(self):return [str(self.lo),str(self.hi)]


def evaluate_x(coeff,t,x):
    out=0
    for c in coeff:out=out*x+polynomial_value(c,t)
    return out


def certificate():
    definition=json.loads((RESULTS/'width4-parametric-definition.json').read_text())
    result=json.loads((RESULTS/'width4-parametric-traces.json').read_text())
    R=definition['crossing_polynomial_ascending_t']
    signs=[1 if c>0 else -1 for c in reversed(R) if c]
    changes=sum(a!=b for a,b in zip(signs,signs[1:]))
    if changes!=1:raise AssertionError('Descartes count is not one')
    lo,hi=F(1),F(2)
    for _ in range(160):
        mid=(lo+hi)/2
        if polynomial_value(R,mid)<0:lo=mid
        else:hi=mid
    if not polynomial_value(R,lo)<0<polynomial_value(R,hi):raise AssertionError('root bracket failed')
    ti=Interval(lo,hi);pi=ti/(1+ti)
    factors=definition['factors'];bands={};diagnostics={}
    # Each real root band is verified uniformly for t in the isolating interval.
    starts={'closed_2':[30.315,.2396],'open_3a':[7.632,-2.983,.989],
            'open_3b':[4.660], 'open_5':[30.315,-2.874,-.548]}
    with mp.workdps(75):
        tn=mp.findroot(lambda t:polynomial_value(R,t),(mp.mpf('1.4'),mp.mpf('1.5')))
        pn=tn/(1+tn)
        for name,guesses in starts.items():
            coefficients=factors[name]['coefficients_desc_x_ascending_t'];bb=[];vals=[]
            for guess in guesses:
                rr=mp.findroot(lambda x:evaluate_x(coefficients,tn,x),guess)
                scale=10**12;a=F(int(mp.floor(rr*scale)),scale);b=a+F(1,scale)
                va=evaluate_x(coefficients,ti,Interval(a));vb=evaluate_x(coefficients,ti,Interval(b))
                if not ((va.positive() and vb.negative()) or (va.negative() and vb.positive())):
                    raise AssertionError(('unverified real root band',name))
                bb.append(Interval(a,b));vals.append(rr)
            bands[name]=bb;diagnostics[name]=[mp.nstr(v,40) for v in vals]
        # The unlisted roots of these factors form one nonreal conjugate pair.
        complex_controls={}
        for name in ('open_3b','open_5'):
            co=factors[name]['coefficients_desc_x_ascending_t'];degree=len(co)-1
            total=-polynomial_value(co[1],ti)
            product=(-1)**degree*polynomial_value(co[-1],ti)
            for band in bands[name]:total=total-band;product=product/band
            discr=total*total-4*product
            if not discr.negative() or not product.positive() or product.hi>=F('5.13')**2:
                raise AssertionError('complex-pair enclosure failed')
            complex_controls[name]={'sum':total.json(),'modulus_squared':product.json(),'discriminant':discr.json()}
        lam=bands['closed_2'][0];mu=bands['open_3a'][0]
        # Equality of the Perron roots follows from the unique resultant root,
        # positive equitable quotients and an endpoint sign change, not decimals.
        other=bands['closed_2'][1:]+bands['open_3a'][1:]+bands['open_3b']+bands['open_5'][1:]
        other += [ti*ti*(1+ti),ti*ti*(1-ti),-ti*ti]
        if max(max(abs(v.lo),abs(v.hi)) for v in other)>=F('5.13'):
            raise AssertionError('other-root modulus bound failed')
        theta=F(674,1000)
        if F('5.13')/mu.lo>=theta or 16*theta**6>=2:raise AssertionError('all-m positivity bound failed')
        def derivative_t(coeff):return [[j*c[j] for j in range(1,len(c))] or [0] for c in coeff]
        def derivative_x(coeff):return [[(len(coeff)-1-i)*v for v in c] for i,c in enumerate(coeff[:-1])]
        derivatives={}
        for name in ('closed_2','open_5'):
            co=factors[name]['coefficients_desc_x_ascending_t']
            derivatives[name]=-evaluate_x(derivative_t(co),ti,lam)/evaluate_x(derivative_x(co),ti,lam)
        hp=(derivatives['open_5']-derivatives['closed_2'])/lam/(1-pi)**2
        if not hp.positive():raise AssertionError('crossing slope not positive')
        # The shared sector's Perron root is rigorously smaller at the crossing.
        rows=result['blocks']['16']['perron_quotient']['rows'];n=len(rows)
        A=mp.zeros(n)
        for i,row in enumerate(rows):
            for k,j,c in row:A[i,j]+=c*tn**k
        eig,vec=mp.eig(A,left=False,right=True);idx=max(range(n),key=lambda i:mp.re(eig[i]));v=[mp.re(vec[i,idx]/vec[0,idx]) for i in range(n)]
        vf=[F(int(mp.nint(a*10**10)),10**10) for a in v]
        if not all(v>0 for v in vf):raise AssertionError('Collatz vector not positive')
        ratios=[sum(c*ti**k*vf[j] for k,j,c in row)/vf[i] for i,row in enumerate(rows)]
        lower=min(v.lo for v in ratios);upper=max(v.hi for v in ratios)
        if upper>=lam.lo:raise AssertionError('shared sector not strictly below leading')
        finite_checks=[]
        for row in result['small_integer_sector_polynomials']:
            m=row['length']
            if 2<=m<=5:
                coefficients=[b-a for a,b in zip(row['rank0_coefficients'],row['rank2_coefficients'])]
                value=polynomial_value(coefficients,ti)
                if not value.positive():raise AssertionError(('finite root sign failed',m))
                finite_checks.append({'m':m,'positive_unnormalized_M_interval':value.json()})
        # Diagnostic numbers use the same implicit derivatives with mp arithmetic.
        ln=mp.findroot(lambda x:evaluate_x(factors['closed_2']['coefficients_desc_x_ascending_t'],tn,x),30)
        mn=mp.findroot(lambda x:evaluate_x(factors['open_3a']['coefficients_desc_x_ascending_t'],tn,x),7.6)
        def dp(name):
            co=factors[name]['coefficients_desc_x_ascending_t']
            return -evaluate_x(derivative_t(co),tn,ln)/evaluate_x(derivative_x(co),tn,ln)
        hpn=(dp('open_5')-dp('closed_2'))/ln/(1-pn)**2
        pub=mp.mpf('0.5914171708531384817988341017359231779642')
        return {'schema':'matching-one.width4-cylinder-rational-certificate.v1',
          'root_t_interval':ti.json(),'root_p_interval':pi.json(),'one_positive_root_Descartes':True,
          'q4_diagnostic':mp.nstr(pn,60),'matches_Jacobsen_2015_table2_n4_within_1e_minus40':abs(pn-pub)<mp.mpf('1e-40'),
          'real_root_intervals':{k:[v.json() for v in vv] for k,vv in bands.items()},
          'real_root_diagnostics':diagnostics,'complex_pair_controls':complex_controls,
          'shared_Perron_Collatz':{'vector':[str(x) for x in vf],'lower':str(lower),'upper':str(upper),
                                 'ratio_to_leading_diagnostic':mp.nstr(mp.re(eig[idx])/ln,30)},
          'leading_prefactors':{'P0':1,'P2':1,'M_open':1,'M_closed':-1},
          'hprime_p_interval':hp.json(),'hprime_p_diagnostic':mp.nstr(hpn,45),
          'rho_diagnostic':mp.nstr(mn/ln,45),'relative_remainder_ratio_upper':str(theta),
          'all_m_root_sign':{'m2_to_5':finite_checks,'m_at_least6':'2 - 16*(674/1000)^m > 0',
             'conclusion':'M_(4,m)(q4)>0 and p_(4,m)<q4 for every integer m>=2'},
          'asymptotic':'p_(4,m)-q4 = -2*rho^m/(m*hprime_p)*(1+O(theta^m)+O(rho^m)) at fixed width four',
          'limits':['not a new infinite-square pc','width-specific, not a uniform-width estimate','numerical roots propose bands; rational inequalities certify all load-bearing signs']}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path);args=ap.parse_args()
    data=certificate();text=json.dumps(data,indent=2,allow_nan=False)+'\n'
    if args.out:
        with args.out.open('x',encoding='utf-8') as f:f.write(text)
    else:print(text,end='')

#!/usr/bin/env python3
"""Post-reveal original-U influence/umbrella diagnostic on pinned N25 m64 sums.

No configurations, RNG, Markov chains, or cloud jobs are used.
Original root brackets are reused without refinement; derivative QA solves
small artificial perturbed-law roots, not additional production targets.
Numerical results are high-precision diagnostics, not interval certificates.
"""
from __future__ import annotations
import argparse
import csv
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import mpmath as mp

ROOT = Path(__file__).resolve().parent
BASE = 'cae9c8997b5994c218bfe060f75656137f745755'
EXPECTED = {
    'axis': '2d23fecc98d276d9ad15ad1867199cd308f0570cb5040ef94eb6b923b4c53458',
    'tilted': '225031e612929ed922ba75c55e76703d59990f5283e7ac39b94f022841798da5',
}
BRACKETS = {
    'star': (
        '365464612308200343214368568562460173533457402519/365375409332725729550921208179070754913983135744',
        '730929224616400686428737137124920347066914805039/730750818665451459101842416358141509827966271488'),
    'drop': (
        '1048047337982522312534188669395795905951669035581/1461501637330902918203684832716283019655932542976',
        '524023668991261156267094334697897952975834517791/730750818665451459101842416358141509827966271488'),
}

def require(test: bool, message: str) -> None:
    if not test:
        raise ValueError(message)

def dec(x: str) -> mp.mpf:
    f = Fraction(x)
    return mp.mpf(f.numerator) / f.denominator

def load_tables() -> dict:
    tables = {}
    for name, expected in EXPECTED.items():
        data = (ROOT/'inputs'/f'{name}.csv').read_bytes()
        require(hashlib.sha256(data).hexdigest() == expected, f'hash mismatch: {name}')
        reader = csv.DictReader(data.decode('utf8').splitlines())
        require(reader.fieldnames == ['k','g','q','count'], f'wrong schema: {name}')
        rows = [tuple(int(row[key]) for key in reader.fieldnames) for row in reader]
        require(all(sum(c for k,_,_,c in rows if k == j) == math.comb(25,j)
                    for j in range(26)), f'count mismatch: {name}')
        tables[name] = rows
    return tables

def cell(rows, h, delta, amplitude=None, perturbation=None):
    weights = [c * h**k * mp.mpf(64)**(-cost+delta*(q+1))
               * (mp.exp(amplitude * perturbation(k,cost,q))
                  if perturbation is not None else 1)
               for k,cost,q,c in rows]
    total = mp.fsum(weights)
    prob = [w/total for w in weights]
    k = [mp.mpf(row[0]) for row in rows]
    q = [mp.mpf(row[2]) for row in rows]
    e = [x*x for x in q]
    def avg(v):
        return mp.fsum(p*x for p,x in zip(prob,v))
    mu, qb, eb = avg(k), avg(q), avg(e)
    ck = [x-mu for x in k]
    cq = avg([(x-qb)*y for x,y in zip(q,ck)])
    ce = avg([(x-eb)*y for x,y in zip(e,ck)])
    cqh = avg([(x-qb)*y*y for x,y in zip(q,ck)])
    ceh = avg([(x-eb)*y*y for x,y in zip(e,ck)])
    return dict(prob=prob,k=k,q=q,e=e,mu=mu,qb=qb,eb=eb,cq=cq,ce=ce,
                cqh=cqh,ceh=ceh,avg=avg,ck=ck)

def pair(tables, h, delta, amplitude=None, perturbations=(None,None)):
    cells = [cell(tables[name],h,delta,amplitude,f)
             for name,f in zip(('axis','tilted'),perturbations)]
    angular = mp.mpf(1152)/625
    ell = [1/angular,-1/angular]
    D = mp.fsum(c['cq'] for c in cells)/2
    B = mp.fsum(a*c['ce'] for a,c in zip(ell,cells))
    R = B/D
    Dh = mp.fsum(c['cqh'] for c in cells)/2
    Bh = mp.fsum(a*c['ceh'] for a,c in zip(ell,cells))
    Rh = (Bh-R*Dh)/D
    for a,c in zip(ell,cells):
        c['phi'] = [(a*((e-c['eb'])*x-c['ce'])
                     -R*((q-c['qb'])*x-c['cq'])/2
                     -Rh*(q-c['qb'])/2)/D
                    for e,q,x in zip(c['e'],c['q'],c['ck'])]
    return cells, R, Rh, D

def calculate(tables, law, dps=120):
    with mp.workdps(dps):
        lo,hi = map(dec,BRACKETS[law]); h=(lo+hi)/2
        delta = int(law=='drop')
        cells, R, Rh, D = pair(tables,h,delta)
        root_residual = mp.fsum(c['qb'] for c in cells)/2
        require(D>0,'nonpositive local thermal slope')
        require(abs(root_residual)<mp.mpf('1e-43'),'root midpoint residual too large')
        out={}
        for name,c in zip(('axis','tilted'),cells):
            probs, phi = c['prob'],c['phi']
            sectors=[]
            for j in range(3):
                selected=[i for i,q in enumerate(c['q']) if q+1==j]
                p=mp.fsum(probs[i] for i in selected)
                mass2=mp.fsum(probs[i]*phi[i]**2 for i in selected)
                sectors.append((p,mass2/p))
            def avar(bias):
                return mp.fsum(p*b for (p,_),b in zip(sectors,bias))*mp.fsum(
                    p*M/b for (p,M),b in zip(sectors,bias))
            opt=[mp.sqrt(M) for _,M in sectors]
            factor=mp.fsum(p*b for (p,_),b in zip(sectors,opt))
            minimum=factor**2
            require(abs(avar(opt)/minimum-1)<mp.mpf('1e-80'),'CS equality failed')
            constant_bias=[mp.mpf(1)]*3
            equal=[1/p for p,_ in sectors]
            # Fixed geometry/barrier-informed choice, not a fitted sampler:
            # axis minimum g=9, tilted minimum g=13.
            barrier=[mp.mpf(1),mp.mpf(64)**(9 if name=='axis' else 13),mp.mpf(1)]
            direct=c['avg']([x*x for x in phi])
            require(abs(avar(constant_bias)/direct-1)<mp.mpf('1e-80'),'direct check')
            # Verify SNIS reconstruction of both means and normalizer algebra.
            rbias=mp.fsum(p*b for (p,_),b in zip(sectors,barrier))
            nu=[p*barrier[int(q+1)]/rbias for p,q in zip(probs,c['q'])]
            iw=[1/barrier[int(q+1)] for q in c['q']]
            z=mp.fsum(p*w for p,w in zip(nu,iw))
            snis_errors=[abs(mp.fsum(p*w*x for p,w,x in zip(nu,iw,c[key]))/z-c['avg'](c[key]))
                         for key in ('k','q','e')]
            dyadic=[mp.mpf(2)**int(mp.floor(mp.log(b/opt[0],2)+mp.mpf('0.5'))) for b in opt]
            dyadic_var=avar(dyadic)
            require(dyadic_var/minimum <= mp.mpf(9)/8,'dyadic Kantorovich bound failed')
            mean_phi=c['avg'](phi)
            require(abs(mean_phi)<mp.mpf('1e-90'),'uncentered influence')
            require(max(snis_errors)<mp.mpf('1e-90'),'SNIS identity failed')
            def text(x): return mp.nstr(x,45)
            out[name]={'sector_P':[text(p) for p,_ in sectors],
                       'sector_conditional_phi2':[text(M) for _,M in sectors],
                       'optimal_sector_probabilities':[text(p*b/factor) for (p,_),b in zip(sectors,opt)],
                       'iid_original_variance':text(direct),
                       'iid_dyadic_bias_variance':text(dyadic_var),
                       'dyadic_bias_factors':[text(x) for x in dyadic],
                       'dyadic_over_sector_optimum':text(dyadic_var/minimum),
                       'iid_equal_sector_variance':text(avar(equal)),
                       'iid_fixed_barrier_bias_variance':text(avar(barrier)),
                       'iid_sector_optimum_variance':text(minimum),
                       'iid_unrestricted_Snis_infimum':text(c['avg']([abs(x) for x in phi])**2),
                       'fixed_barrier_bias':[text(x) for x in barrier],
                       'full_influence_mean':text(mean_phi),
                       'snis_reconstruction_max_error':text(max(snis_errors))}
        clock=mp.fsum(c['avg']([f*x for f,x in zip(c['phi'],c['ck'])]) for c in cells)
        scale=abs(R)+abs(Rh)+mp.mpf('1e-100')
        require(abs(clock)/scale<mp.mpf('1e-70'),'common thermal gauge not cancelled')
        # Ratios cancel A_N. These are asymptotic variance-equivalent iid counts,
        # not finite-sample SNR, CI, Markov-step budget or all-algorithm lower bound.
        budgets={}
        for field in ('iid_original_variance','iid_equal_sector_variance',
                      'iid_fixed_barrier_bias_variance','iid_sector_optimum_variance',
                      'iid_unrestricted_Snis_infimum','iid_dyadic_bias_variance'):
            value=9*mp.fsum(mp.mpf(v[field]) for v in out.values())/R**2
            budgets[field]=mp.nstr(value,45)
        return {'law':law,'input_root_bracket':BRACKETS[law],
                'root_midpoint':mp.nstr(h,55),'root_residual':mp.nstr(root_residual,15),
                'U_over_A':mp.nstr(R,45),'d_logh_U_over_A':mp.nstr(Rh,45),
                'pooled_cov_q_K':mp.nstr(D,45),'geometry':out,
                'common_clock_identity_residual':mp.nstr(clock,15),
                'asymptotic_variance_equivalent_n_per_geometry_for_3se':budgets}


def derivative_checks(tables, law):
    """Independent full recomputation at perturbed finite laws; no sampling."""
    with mp.workdps(130):
        h0=(dec(BRACKETS[law][0])+dec(BRACKETS[law][1]))/2
        delta=int(law=='drop')
        cells,R,Rh,D=pair(tables,h0,delta)
        f1=lambda k,g,q: mp.mpf(int(q==0))
        f2=lambda k,g,q: mp.mpf(k*k*q)/625
        directions=((f1,None),(f2,f2))
        output=[]
        step=mp.mpf('1e-26')
        for fs in directions:
            expected=mp.fsum(c['avg']([v*(f(*row[:3]) if f else 0)
                        for v,row in zip(c['phi'],tables[name])])
                        for name,c,f in zip(('axis','tilted'),cells,fs))
            values=[]
            for sign in (-1,1):
                x=mp.log(h0)
                for _ in range(15):
                    cs,r,_,d=pair(tables,mp.exp(x),delta,sign*step,fs)
                    root_residual=mp.fsum(c['qb'] for c in cs)/2
                    if abs(root_residual)<mp.mpf('1e-115'):
                        break
                    x-=root_residual/d
                else:
                    raise ValueError('perturbed-law root check did not converge')
                values.append(r)
            observed=(values[1]-values[0])/(2*step)
            relative=abs(observed-expected)/max(abs(expected),mp.mpf('1e-70'))
            require(relative<mp.mpf('1e-35'),'full recomputation derivative mismatch')
            output.append({'influence_derivative':mp.nstr(expected,45),
                           'full_recomputed_derivative':mp.nstr(observed,45),
                           'relative_error':mp.nstr(relative,12)})
        return output

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    require(not args.output.exists(),'output must not exist')
    tables=load_tables()
    results=[calculate(tables,law,120) for law in ('star','drop')]
    checks=[calculate(tables,law,170) for law in ('star','drop')]
    # Independent precision is a stability check, not an independent algorithm.
    for a,b in zip(results,checks):
        for key in a['asymptotic_variance_equivalent_n_per_geometry_for_3se']:
            require(a['asymptotic_variance_equivalent_n_per_geometry_for_3se'][key]
                    ==b['asymptotic_variance_equivalent_n_per_geometry_for_3se'][key],
                    '120/170-digit stability check failed')
    payload={'status':'POST_REVEAL_HIGH_PRECISION_IID_ENVELOPE_NOT_A_SAMPLER',
             'input_commit':BASE,'input_sha256':EXPECTED,
             'mpmath_version':mp.__version__,'precision_dps':[120,170],
             'geometry':'N25 axis (5,0) / tilted (4,3); same original DeltaCos4=1152/625',
             'claim_boundary':'Exact algebra; decimal diagnostics at existing root-bracket midpoints. '
                 'Not new simulation, not prospective evidence, not rational interval certification. '
                 'IID asymptotic variance only; root/mean/denominator influence INCLUDED. '
                 'Optimal biases use full old archive as an oracle. No sampling/mixing cost is established. '
                 'Zero-phi proposal infimum is restricted to regular single-proposal SNIS.',
             'new_samples':0,'new_enumerations':0,'cloud_jobs':0,'laws':results,
             'independent_numerical_derivative_checks':{law:derivative_checks(tables,law)
                                                       for law in ('star','drop')}}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(payload,ensure_ascii=False,indent=2)+'\n')
    for law in results:
        print(law['law'],law['U_over_A'])
        print(json.dumps(law['asymptotic_variance_equivalent_n_per_geometry_for_3se'],indent=2))
        for g,v in law['geometry'].items(): print(g,'optimal probabilities',v['optimal_sector_probabilities'])
if __name__=='__main__': main()

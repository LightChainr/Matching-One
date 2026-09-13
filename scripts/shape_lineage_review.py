#!/usr/bin/env python3
"""Read #706's stored Q/covariance: correct uncertainty, test width-scaled skew.

No histogram reconstruction or simulation. Original artifacts are immutable.
New nonlinear statistics use first-order propagation of the full Q covariance.
The two weightings are sensitivities of the same blocks, never independent votes.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from mpmath import mp

ROOT=Path(__file__).resolve().parents[1]
INPUT=ROOT/'results/probe-invariant-shape/quantile-shape-lineage-622.json'
SIZES=(145,290,725)
AI=np.array([0,2,3,4]); ZI=np.array([0,2,3,4,5,6,8])


def chart(q):
    q=np.asarray(q); w=q[7]-q[1]
    z=(q-q[1])/w
    return w,z,z+z[::-1]-1


def derivative(fun, x):
    """Complex-step Jacobian; functions here are rational analytic maps of Q."""
    x=np.asarray(x,dtype=float); cols=[]
    for i in range(x.size):
        shifted=x.astype(complex);shifted[i]+=1e-25j
        cols.append(np.imag(np.atleast_1d(fun(shifted)))/1e-25)
    return np.asarray(cols).T


def propagated(fun,x,cov):
    val=np.atleast_1d(fun(np.asarray(x,dtype=float)))
    jac=derivative(fun,x);c=jac@cov@jac.T
    return val,(c+c.T)/2,jac


def block_diagonal(blocks):
    d=sum(len(b) for b in blocks);out=np.zeros((d,d));at=0
    for block in blocks:
        n=len(block);out[at:at+n,at:at+n]=block;at+=n
    return out


def zero_test(value,cov):
    """60-dps solve; estimated covariance means only a nominal reference law."""
    value=np.asarray(value,dtype=float);cov=np.asarray(cov,dtype=float)
    if not np.all(np.isfinite(value)) or not np.all(np.isfinite(cov)):
        raise ValueError('nonfinite data')
    sd=np.sqrt(np.diag(cov));ev=np.linalg.eigvalsh(cov)
    out={'max_abs_marginal_z':float(np.max(np.abs(value)/sd)),
         'covariance_condition':float(ev[-1]/ev[0]) if ev[0]>0 else None,
         'dof':len(value),'reference':'nominal Gaussian; estimated covariance'}
    if ev[0]<=0:
        out.update(statistic=None,reason='nonpositive covariance; no pseudoinverse')
        return out
    with mp.workdps(60):
        c=mp.matrix([[mp.mpf(str(v)) for v in row] for row in cov])
        y=mp.matrix([mp.mpf(str(v)) for v in value]);sol=mp.lu_solve(c,y)
        d=(y.T*sol)[0]
        if d<0: raise ValueError('negative quadratic form')
        tail=mp.gammainc(mp.mpf(len(value))/2,d/2,mp.inf,regularized=True)
        out.update(statistic=float(d),log10_p_nominal=float(mp.log10(tail)))
    return out


def width_se(cov):
    cov=np.asarray(cov);var=cov[7,7]+cov[1,1]-2*cov[1,7]
    if var<=0:raise ValueError('nonpositive width variance')
    return float(np.sqrt(var))


def interval_norm_change(values,covs):
    """Difference of adjacent displacement norms, retaining the shared middle."""
    d1=values[1]-values[0];d2=values[2]-values[1]
    n1,n2=np.linalg.norm(d1),np.linalg.norm(d2)
    if min(n1,n2)<=0:raise ValueError('norm derivative undefined at zero')
    g1,g2=d1/n1,d2/n2
    independent=g1@(covs[0]+covs[1])@g1+g2@(covs[1]+covs[2])@g2
    shared=2*g1@covs[1]@g2
    variance=independent+shared
    if variance<=0:raise ValueError('nonpositive propagated variance')
    return {'value':float(n2-n1),'se':float(np.sqrt(variance)),
            'absolute_z_nominal':float(abs(n2-n1)/np.sqrt(variance)),
            'omitted_shared_middle_variance_term':float(shared),
            'old_independent_se':float(np.sqrt(independent))}


def normalized_skew(q):
    w,_,a=chart(q);return a[AI]/w


def normal_form_residual(q):
    """Necessary leading smooth-chart shape; median coefficient eliminated."""
    w,z,a=chart(q);x=(z-z[::-1])/2
    k=-4*a[4]/w
    return (a/w-k*(x*x-.25))[AI[:3]]


def beta_from_source(q):
    """One quadratic chart fixed by the N145 median only; not a target fit."""
    b=2*q[4]-q[1]-q[7]
    v=(q-.5)**2;c=2*v[4]-v[1]-v[7]
    return -b/c


def common_quadratic_residual(qall):
    beta=beta_from_source(qall[:9]);out=[]
    for i in range(3):
        q=qall[9*i:9*(i+1)];p=q+beta*(q-.5)**2
        _,_,a=chart(p)
        out.extend(a[AI[:3]] if i==0 else a[AI])
    return np.asarray(out)


def analyze(data):
    if data['deciles']!=[.1,.2,.3,.4,.5,.6,.7,.8,.9]:raise ValueError('unexpected grid')
    result={'schema':'matching-one.shape-lineage-reviewed.v1',
            'source_pr':706,'source_head':'b52e3a35d2c88dfc1bf0c17508d3dc4ce5d570f1',
            'source_artifact':str(INPUT.relative_to(ROOT)),
            'standing':'same-block C2; no new samples; nonlinear diagnostics use delta-method full-Q propagation',
            'weights':{}}
    for weighting in ('spin0','equal'):
        rows=[data['sizes'][str(n)]['analysis'][weighting] for n in SIZES]
        qs=[np.array(r['Q_pooled']) for r in rows]
        cs=[np.array(r['Q_covariance']) for r in rows]
        report={'sizes':{},'shared_middle_corrected':{}}
        skew,cskew=[],[]
        for n,r,q,c in zip(SIZES,rows,qs,cs):
            w,z,a=chart(q)
            if w<=0 or np.any(np.diff(q)<=0):raise ValueError('nonmonotone quantile')
            sk,sc,_=propagated(normalized_skew,q,c)
            nf,nc,_=propagated(normal_form_residual,q,c)
            az,acz,_=propagated(lambda x:chart(x)[2][AI],q,c)
            original=np.asarray(r['A_covariance'])[np.ix_(AI,AI)]
            diagdiff=float(np.max(abs(np.diag(acz)/np.diag(original)-1)))
            corrected=width_se(c)
            med,semed=sk[-1],np.sqrt(sc[-1,-1])
            skew.append(sk);cskew.append(sc)
            report['sizes'][str(n)]={
                'W':float(w),'W_se_corrected':corrected,'W_se_original':r['W_se'],
                'width_se_factor':corrected/r['W_se'],
                'A_over_W':sk.tolist(),'A_over_W_se':np.sqrt(np.diag(sc)).tolist(),
                'A_over_W_covariance':sc.tolist(),
                'A_norm_over_W':float(np.linalg.norm(a[AI])/w),
                'median_analytic_coefficient':float(-4*med),'median_coefficient_se':float(4*semed),
                'leading_normal_form_residual':nf.tolist(),
                'leading_normal_form_residual_over_W_squared':(nf/w**2).tolist(),
                'leading_normal_form_zero_diagnostic':zero_test(nf,nc),
                'A_zero_recheck_from_stored_A_covariance':zero_test(a[AI],original),
                'delta_vs_jackknife_A_covariance_max_diagonal_relative_difference':diagdiff}
        for name,indices in (('Z',ZI),('A',AI)):
            vals=[np.array(r[f'{name}_pooled'])[indices] for r in rows]
            covs=[np.array(r[f'{name}_covariance'])[np.ix_(indices,indices)] for r in rows]
            report['shared_middle_corrected'][name]=interval_norm_change(vals,covs)
        d=np.concatenate([skew[1]-skew[0],skew[2]-skew[1]])
        dc=np.block([[cskew[0]+cskew[1],-cskew[1]],[-cskew[1],cskew[1]+cskew[2]]])
        report['common_A_over_W_exact_equality_test']=zero_test(d,dc)
        ks=np.array([-4*v[-1] for v in skew]); kv=np.array([16*v[-1,-1] for v in cskew])
        report['median_coefficient_adjacent_changes']=[
            {'pair':f'{SIZES[i]}->{SIZES[i+1]}','difference':float(ks[i+1]-ks[i]),
             'se':float(np.sqrt(kv[i]+kv[i+1]))} for i in range(2)]
        qall=np.concatenate(qs);joint=block_diagonal(cs)
        beta,bc,_=propagated(lambda q:np.array([beta_from_source(q[:9])]),qall,joint)
        residual,rc,_=propagated(common_quadratic_residual,qall,joint)
        report['source_fixed_quadratic_chart']={
            'definition':'phi(p)=p+beta*(p-.5)^2; beta from N145 median only',
            'beta':float(beta[0]),'beta_se':float(np.sqrt(bc[0,0])),
            'minimum_derivative_on_unit_interval':float(1-abs(beta[0])),
            'residuals_source3_then_target4_then_target4':residual.tolist(),
            'held_target_joint_test_with_shared_source_uncertainty':zero_test(residual[3:],rc[3:,3:]),
            'source_other_coordinates_test':zero_test(residual[:3],rc[:3,:3]),
            'residual_norm_ratios_to_raw_A':[
                float(np.linalg.norm(residual[3:7])/np.linalg.norm(chart(qs[1])[2][AI])),
                float(np.linalg.norm(residual[7:])/np.linalg.norm(chart(qs[2])[2][AI]))],
            'interpretation':'test of ONE exact quadratic chart, not every smooth scaling field; retrospective, not preregistered'}
        result['weights'][weighting]=report
    result['limits']=[
        'A/W equality is a finite model, not an exponent or universality theorem.',
        'Smooth-coordinate normal form permits an O(W^2) residual in A/W; rejecting its exact zero does not refute a smooth chart.',
        'Quadratic-chart cancellation is a diagnostic, not permission to redefine the physical occupation parameter.',
        'Full chi-square smaller than a diagonal sum is not a robustness certificate.',
        'sqrt(N) is a Euclidean length, not an integer count of sites on a shortest winding path.']
    return result


def compact(result):
    out={k:v for k,v in result.items() if k!='weights'};out['weights']={}
    for name,r in result['weights'].items():
        rr={k:v for k,v in r.items() if k!='sizes'};rr['sizes']={}
        for n,row in r['sizes'].items():
            rr['sizes'][n]={k:v for k,v in row.items() if k not in ('A_over_W_covariance','A_over_W_se','A_zero_recheck_from_stored_A_covariance')}
        out['weights'][name]=rr
    return out


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,default=INPUT)
    ap.add_argument('--out',type=Path);ap.add_argument('--compact',action='store_true');args=ap.parse_args()
    result=analyze(json.loads(args.input.read_text()))
    if args.compact:result=compact(result)
    text=json.dumps(result,indent=2,allow_nan=False)+'\n'
    if args.out:
        args.out.parent.mkdir(parents=True,exist_ok=True)
        with args.out.open('x') as f:f.write(text)
    print(text,end='')

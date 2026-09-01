#!/usr/bin/env python3
"""Synthetic decision cases and old-archive derivative validation; no fresh RNG."""
import csv,json,math,subprocess,sys
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from archive_channel_split import ROOT,moments,kernels,channel_split,load
from score_production import decisions

def main():
    contract=json.loads((ROOT/'CONTRACT.json').read_text());checks={}
    cases={'weak':(0.,0.),'entry':(1.,0.),'completion':(0.,1.),'gap':(.45,.45),'cancel':(1.,-1.)}
    for name,(v1,v2) in cases.items():
        se=.01;est={}
        for n in (85,340):
            for key,value,error in (('entry',v1,se),('completion',v2,se),('total',v1+v2,math.sqrt(2)*se)):
                est[f'N{n}.{key}.v']={'value':value,'se':error}
        result=decisions(est,contract);checks[name]=result
    assert checks['weak']['remaining_templates']==['W_double_channel_weak']
    assert checks['entry']['remaining_templates']==['B_entry_selective']
    assert checks['completion']['remaining_templates']==['C_completion_selective']
    assert checks['gap']['all_three_numerical_templates_contradicted'] and checks['gap']['net_above_positive_target_at_both_N']
    assert checks['cancel']['all_three_numerical_templates_contradicted'] and checks['cancel']['stop_current_main_H4_priority_for_this_source']
    bases=np.load(ROOT/'inputs/old_profiles.npz');anchor=json.loads((ROOT/'inputs/anchors.json').read_text());fd=[];audit=[]
    for n in (85,340):
        base=bases[f'n{n}'].sum(axis=0).astype(float);events=load(n).sum(axis=0);m=100000 if n==85 else 1000000;p0=anchor['by_N'][str(n)]['p0'];t=kernels(base,events,m)
        profile=base.copy();profile[...,2]=0;profile[...,3]=t[...,0]+2*t[...,1]+t[...,2];profile[...,4]=-t[...,0]+t[...,2]
        analytic=channel_split(base,t,m,n,p0)
        def tilted(lam):
            def packet(p):
                out=[]
                for g in range(2):
                    v,d,_=moments(profile[g],m,p,n);out.append((v[:2]+lam*v[3:],d[:2]+lam*d[3:]))
                return out
            p=brentq(lambda p:sum(x[0][0] for x in packet(p))/2,.55,.65,xtol=5e-15,rtol=5e-15);x=packet(p)
            d=np.mean([r[1][0] for r in x]);delta=2304/1445;factor=n**(13/8)/2/d/delta
            entry=-factor*((x[0][1][0]-x[0][1][1])-(x[1][1][0]-x[1][1][1]))/2
            completion=factor*((x[0][1][0]+x[0][1][1])-(x[1][1][0]+x[1][1][1]))/2
            return np.array([entry,completion])
        numeric=(tilted(5e-6)-tilted(-5e-6))/1e-5;truth=np.array([analytic['entry.v'],analytic['completion.v']]);error=float(np.max(np.abs(numeric-truth)))
        assert error<2e-5,(n,error)
        fd.append({'N':n,'shared_root_channel_finite_difference_max_error':error})
        old_csv=ROOT/f'qa/old-batch0-n{n}.csv'
        if old_csv.exists():
            observed=np.zeros((2,n+1,11))
            with old_csv.open() as f:
                for r in csv.DictReader(f):observed[int(r['g']),int(r['k'])]=[int(r[k]) for k in ('q','e','s','qs','es','count01','count02','count12','s_previous01','s_previous02','s_previous12')]
            expected=np.concatenate((bases[f'n{n}'][0],load(n)[0]),axis=-1)
            assert np.array_equal(observed,expected),(n,np.max(np.abs(observed-expected)))
            audit.append({'N':n,'original_batch0_exact_all11_moments':True,'old_permutations':1000 if n==85 else 10000,'new_random_samples':0})
    process=subprocess.run([sys.executable,str(ROOT/'run_production.py'),'--n','85'],capture_output=True,text=True)
    assert process.returncode!=0 and 'NO FREEZE AUTHORIZATION' in process.stderr
    result={'status':'pre_freeze_QA_passed','synthetic_decisions':checks,'old_archive_derivative_checks':fd,'producer_old_batch_checks':audit,
            'unauthorized_production_refused_before_sampling':True,'prospective_samples':0,
            'producer_old_audit_scope':'Only old original batch0; missing old audit files means this part not yet run, not a passed producer check'}
    (ROOT/'PRE_FREEZE_QA.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'producer_old_batch_checks':audit,'fd':fd,'prospective_samples':0}))
if __name__=='__main__':main()

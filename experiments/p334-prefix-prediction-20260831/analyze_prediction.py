#!/usr/bin/env python3
"""Grouped cross-fitted prediction of local response from exact source Gram."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import time
import numpy as np

MODELS=['constant','BG']

def crossfit(data,omit=None,keep_arrays=False):
    take=np.ones(len(data['batch']),dtype=bool) if omit is None else data['batch']!=omit
    x=data['Jbar'][take];g=data['G'][take];u=data['energy_U'][take]
    batch=data['batch'][take];fold=data['fold'][take]
    preds=np.empty((len(x),2,x.shape[1],2,2))
    fits=[]
    for k in range(5):
        train=fold!=k;test=~train
        assert train.any() and test.any()
        normal=np.einsum('pij,pkj->ik',g[train],g[train])
        eigen=np.linalg.eigvalsh(normal)
        assert eigen[0]>0 and eigen[-1]/eigen[0]<1e8
        constant=x[train].mean(axis=0)
        rhs=np.einsum('poij,pkj->oik',x[train],g[train])
        coeff=np.stack([np.linalg.solve(normal,r.T).T for r in rhs])
        preds[test,0]=constant
        preds[test,1]=np.einsum('oij,pjk->poik',coeff,g[test])
        assert not set(batch[train]) & set(batch[test])
        fits.append({'fold':k,'training_prefixes':int(train.sum()),'test_prefixes':int(test.sum()),
                     'training_batches':np.unique(batch[train]).tolist(),'test_batches':np.unique(batch[test]).tolist(),
                     'constant_matrix':constant.tolist(),'B':coeff.tolist(),'normal_eigenvalues':eigen.tolist(),
                     'normal_condition_number':float(eigen[-1]/eigen[0])})
    dot=np.einsum('pmoij,poij->pmo',preds,x)
    prednorm=np.sum(preds**2,axis=(-2,-1))
    risk=u[:,None,:]-2*dot+prednorm
    naive=np.sum((preds-x[:,None])**2,axis=(-2,-1))
    stats={}
    n0=len(x)
    original_n=20000 if omit is None else 19000
    # Main conditional-cell risk pools held-out prefixes equally despite unequal cell counts per batch.
    for o,name in enumerate(data['output_names']):
        name=str(name)
        stats[name+'.energy_U']=float(u[:,o].mean())
        stats[name+'.MC_noise_in_naive_MSE']=float(np.mean(np.sum(x[:,o]**2,axis=(-2,-1))-u[:,o]))
        for m,model in enumerate(MODELS):
            stats[name+'.'+model+'.risk_U']=float(risk[:,m,o].mean())
            stats[name+'.'+model+'.risk_naive']=float(naive[:,m,o].mean())
            stats[name+'.'+model+'.gain_over_zero']=float((2*dot[:,m,o]-prednorm[:,m,o]).mean())
            stats[name+'.'+model+'.population_risk_contribution']=float(risk[:,m,o].sum()/original_n)
        gain=risk[:,0,o]-risk[:,1,o]
        # Energy is the same for both predictions and must cancel in their paired difference.
        cancellation=(prednorm[:,0,o]-prednorm[:,1,o])-2*(dot[:,0,o]-dot[:,1,o])
        assert np.allclose(gain,cancellation,rtol=1e-10,atol=1e-18)
        stats[name+'.BG_gain_over_constant']=float(gain.mean())
        stats[name+'.BG_gain_population_contribution']=float(gain.sum()/original_n)
        denom=stats[name+'.constant.risk_U']
        stats[name+'.relative_BG_gain']=float(gain.mean()/denom) if denom!=0 else float('nan')
    stats['cell00_prefixes']=float(n0)
    out={'statistics':stats,'fits':fits}
    if keep_arrays:
        out['arrays']={'counter':data['counter'][take],'batch':batch,'fold':fold,'prediction':preds,'risk_U':risk,'risk_naive':naive,'energy_U':u,
                       'Jbar':x,'G':g,'output_names':data['output_names'],'models':np.array(MODELS)}
    return out

def finite_risk_identity_check():
    # Exhaust all two-draw samples from a fixed finite matrix law; prediction is independent.
    law=np.array([[[1.,2.],[-1.,0.]],[[0.,-1.],[2.,3.]],[[2.,0.],[1.,-2.]]])
    prediction=np.array([[.3,-.2],[.1,.5]])
    expected=np.sum((law.mean(0)-prediction)**2)
    risks=[]
    for a in law:
        for b in law:
            risks.append(np.sum(a*b)-2*np.sum(prediction*(a+b)/2)+np.sum(prediction**2))
    error=float(abs(np.mean(risks)-expected))
    assert error<1e-12
    return error

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inputs',type=Path,default=Path(__file__).parent/'inputs')
    p.add_argument('--output',type=Path,default=Path(__file__).parent/'results')
    a=p.parse_args();a.output.mkdir(exist_ok=False);started=time.time()
    manifest=json.loads((a.inputs/'manifest.json').read_text())
    for f in manifest['derived_files']:
        assert hashlib.sha256((a.inputs/f['path']).read_bytes()).hexdigest()==f['sha256']
    out={'estimand':'Within-cell00 held-out squared Frobenius prediction risk for true local mean-response J(Z), N325/N425 fit separately.',
         'outputs':['A at p_ref=.59274605079','C/(N+1), derived exactly as -.5 times integral-A response'],
         'source_gram':'Exact G=E[HH^T|Z]=sum_a pi_a^3 Cov(L|a), from full label census; deterministic prefix feature.',
         'models':'constant K and BG, each4 parameters per output. No feature search, penalty, tuning, or new sampling.',
         'split':'5 grouped folds, batch mod5. All72 quartets of each prefix stay together. Fitting never uses held-out batches.',
         'risk':'U_energy-2<P,Jbar>+||P||². U_energy uses distinct quartets. Paired constant-minus-BG risk cancels common U_energy exactly.',
         'uncertainty':'Delete one original batch and refit all5 folds, then20-batch joint jackknife factor. Prefix-weighted within-cell means; also original-population contributions.',
         'interpretation':'Predictive risk reduction is not proof of mechanistic closure, causal uniqueness, or universality. No reduction limits only the fixed shared-B model.',
         'new_prefixes':0,'new_labels':0,'new_tails':0,'finite_risk_identity_error':finite_risk_identity_check(),'sizes':{}}
    for n in (325,425):
        data=dict(np.load(a.inputs/f'N{n}.npz'))
        assert len(data['counter'])=={325:1502,425:1551}[n]
        full=crossfit(data,keep_arrays=True)
        drops=[crossfit(data,omit=b) for b in range(20)]
        labels=list(full['statistics']);point=np.array(list(full['statistics'].values()))
        loo=np.array([list(z['statistics'].values()) for z in drops])
        ratio_gates={}
        for name in data['output_names']:
            lab=str(name)+'.relative_BG_gain';j=labels.index(lab)
            denom_j=labels.index(str(name)+'.constant.risk_U')
            valid=point[denom_j]>0 and np.all(loo[:,denom_j]>0) and np.all(np.isfinite(loo[:,j]))
            ratio_gates[str(name)]={'all_denominators_positive':bool(valid),'point_denominator':float(point[denom_j]),'LOO_min_denominator':float(loo[:,denom_j].min())}
            if not valid:
                # Ratios are omitted rather than clipping a noisy risk denominator.
                del labels[j];point=np.delete(point,j);loo=np.delete(loo,j,axis=1)
        factor=np.sqrt(19/20)*(loo-loo.mean(axis=0));se=np.linalg.norm(factor,axis=0)
        np.savez_compressed(a.output/f'heldout_predictions_N{n}.npz',**full['arrays'])
        out['sizes'][str(n)]={'labels':labels,'estimate':point.tolist(),'se':se.tolist(),'LOO':loo.tolist(),'factor':factor.tolist(),
           'omitted_batch_ids':list(range(20)),'full_crossfit_models':full['fits'],'delete_one_crossfit_models':[z['fits'] for z in drops],
           'relative_gain_gates':ratio_gates,'original_cell00_batch_counts':np.bincount(data['batch'],minlength=20).tolist()}
        for label,val,err in zip(labels,point,se):
            if any(k in label for k in ('risk_U','energy_U','BG_gain_over_constant','relative_BG_gain')):
                print(n,label,f'{val:.12g} +/- {err:.6g}',flush=True)
    (a.output/'score.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    receipt={'host':platform.node(),'python':platform.python_version(),'numpy':np.__version__,'started_unix':started,'finished_unix':time.time(),'elapsed_seconds':time.time()-started,
       'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'input_manifest_sha256':hashlib.sha256((a.inputs/'manifest.json').read_bytes()).hexdigest(),
       'plan_sha256':hashlib.sha256(Path(__file__).with_name('PLAN.md').read_bytes()).hexdigest(),'analysis_attempt':1,'new_sampling':0}
    (a.output/'run_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)

if __name__=='__main__':main()

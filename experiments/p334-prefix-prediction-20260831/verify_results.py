#!/usr/bin/env python3
"""Audit saved held-out predictions, leakage exclusion, and independent QR fits."""
import hashlib
import json
from pathlib import Path
import numpy as np

root=Path(__file__).parent
score=json.loads((root/'results/score.json').read_text())
manifest=json.loads((root/'inputs/manifest.json').read_text())
for f in manifest['derived_files']:
    assert hashlib.sha256((root/'inputs'/f['path']).read_bytes()).hexdigest()==f['sha256']
qa={'sizes':{},'finite_risk_identity_error':score['finite_risk_identity_error'],'new_sampling':0}
for n in (325,425):
    data=dict(np.load(root/'inputs'/f'N{n}.npz'))
    held=dict(np.load(root/'results'/f'heldout_predictions_N{n}.npz'))
    s=score['sizes'][str(n)]
    max_fit_error=0.;max_condition=0.
    coefficient_change=0.
    for omitted,models in [(None,s['full_crossfit_models'])]+list(enumerate(s['delete_one_crossfit_models'])):
        for fit in models:
            train=np.isin(data['batch'],fit['training_batches'])
            test=np.isin(data['batch'],fit['test_batches'])
            assert not set(fit['training_batches']) & set(fit['test_batches'])
            if omitted is not None:
                assert omitted not in fit['training_batches'] and omitted not in fit['test_batches']
            assert train.sum()==fit['training_prefixes'] and test.sum()==fit['test_prefixes']
            # Stack both source columns as observations; solve by QR/SVD least squares,
            # independently of the producer's2x2 normal-equation construction.
            design=data['G'][train].transpose(0,2,1).reshape(-1,2)
            target=data['Jbar'][train]
            coeff=np.array([[np.linalg.lstsq(design,target[:,o,i,:].reshape(-1),rcond=None)[0] for i in range(2)] for o in range(2)])
            max_fit_error=max(max_fit_error,float(np.max(abs(coeff-np.asarray(fit['B'])))))
            assert np.allclose(coeff,fit['B'],rtol=1e-10,atol=1e-12)
            assert np.allclose(target.mean(0),fit['constant_matrix'],rtol=1e-12,atol=1e-16)
            max_condition=max(max_condition,fit['normal_condition_number'])
            if omitted is not None:
                coefficient_change=max(coefficient_change,float(np.max(abs(coeff-np.asarray(s['full_crossfit_models'][fit['fold']]['B'])))))
    assert coefficient_change>0
    u=held['energy_U'];p=held['prediction'];x=held['Jbar']
    risk=u[:,None,:]-2*np.sum(p*x[:,None],axis=(-2,-1))+np.sum(p*p,axis=(-2,-1))
    assert np.allclose(risk,held['risk_U'],rtol=1e-10,atol=1e-18)
    assert np.array_equal(held['counter'],data['counter']) and np.array_equal(held['fold'],data['batch']%5)
    d={k:v for k,v in zip(s['labels'],s['estimate'])}
    for o,name in enumerate(data['output_names']):
        for m,model in enumerate(('constant','BG')):
            assert np.isclose(risk[:,m,o].mean(),d[str(name)+'.'+model+'.risk_U'],rtol=1e-12,atol=1e-20)
        assert np.isclose((risk[:,0,o]-risk[:,1,o]).mean(),d[str(name)+'.BG_gain_over_constant'],rtol=1e-12,atol=1e-20)
    loo=np.array(s['LOO']);factor=np.sqrt(19/20)*(loo-loo.mean(0))
    assert np.allclose(factor,s['factor'],rtol=1e-12,atol=1e-20)
    assert np.allclose(np.linalg.norm(factor,axis=0),s['se'],rtol=1e-12,atol=1e-20)
    qa['sizes'][str(n)]={'prefixes':len(x),'all_full_and_LOO_folds_no_leakage':True,'independent_lstsq_B_max_abs_difference':max_fit_error,
                        'maximum_training_normal_condition':max_condition,'maximum_LOO_coefficient_change':coefficient_change,
                        'heldout_risk_and_full_refit_factor_checks':'pass'}
(root/'QA.json').write_text(json.dumps(qa,indent=2)+'\n')
print(json.dumps(qa,indent=2))

#!/usr/bin/env python3
"""Old archived inputs and deterministic arithmetic QA only; no new sampling."""
import csv
from fractions import Fraction
import hashlib
import io
import json
from pathlib import Path
import subprocess

import numpy as np

from exact_policy import policy, verify_finite_contrast
from export_existing_model import export_model
from score_prospective import derive

ROOT=Path(__file__).resolve().parent


def main():
    receipt={'new_prefixes':0,'new_tails':0,'production_binary_executed':False,'new_fits':0}
    training=json.loads((ROOT/'inputs/existing_contact_score.json').read_text())
    model=export_model(training)
    saved=json.loads((ROOT/'existing_model.json').read_text())
    assert model['sizes']==saved['sizes']
    residual_errors=[];loading_errors=[]
    predictors=list(training['predictors'])
    for n,item in training['sizes'].items():
        d=dict(zip(item['labels'],item['estimate']))
        for oi,ori in enumerate(('first','second')):
            beta=np.array(model['sizes'][n]['point'][ori]['beta'])
            def k(i,j):
                i,j=sorted((i,j));return d[f'{ori}.K.{predictors[i]}|{predictors[j]}']
            captured=2*sum(k(4,i)*beta[i,2*oi] for i in range(4))-.5*sum(k(5,i)*beta[i,2*oi+1] for i in range(4))
            old_captured=d[f'{ori}.new64.contact.captured_signed_loading']
            residual=d[f'{ori}.new64.own_signed_loading']-captured
            old_residual=model['sizes'][n]['point'][ori]['R_old']
            loading_errors.append(abs(captured-old_captured));residual_errors.append(abs(residual-old_residual))
    assert max(loading_errors)<1e-22 and max(residual_errors)<1e-22
    receipt.update(max_old_captured_loading_reconstruction_error=max(loading_errors),max_old_residual_loading_reconstruction_error=max(residual_errors))
    # Two degree classes plus one unchanged outside-safe label; all rational.
    census=[(0,True,2,2,0,1),(1,True,2,2,1,0),(2,True,2,2,1,1),
            (3,True,3,3,0,2),(4,True,3,3,2,0),(5,False,4,4,0,0)]
    means={i:Fraction(i*i+2*i,7) for i in range(6)}
    checks=[]
    for source in (0,1):
        response=verify_finite_contrast(census,source,means)
        _,qp,qm,w=policy(census,source)
        checks.append({'source':source,'response':str(response),'min_probability':str(min(*qp.values(),*qm.values())),'contrast_mass_W':str(w)})
    receipt['exact_affine_and_RB_checks']=checks
    # Old-prefix checkpoints/censuses only, using original immutable counters.
    old_root=ROOT.parent/'p334-mechanism-response-20260831/inputs/prefix_archive'
    oldqa=ROOT/'qa';oldqa.mkdir(exist_ok=True)
    field_errors=[];rows=0;source_hashes={}
    for n in (325,425):
        original=old_root/f'N{n}.csv'
        raw=subprocess.check_output([str(ROOT/'qa_old'),str(n),str(original)])
        (oldqa/f'old_N{n}_first1024.csv').write_bytes(raw)
        blob=subprocess.check_output(['git','show',f'1cfa4ae892a2f7f4168e9a71690efd7a5560d4cd:results/p334-exact-prefix-structure/N{n}.npz'],cwd=ROOT)
        source_hashes[str(n)]=hashlib.sha256(blob).hexdigest()
        with np.load(io.BytesIO(blob),allow_pickle=False) as z:
            for r in csv.DictReader(io.StringIO(raw.decode())):
                i=int(r['index']);assert int(z['counter'][i])==int(r['counter'])
                actual=np.array([float(r[o+'_'+f]) for o in ('first','second') for f in ('mass','energy','degree','loop')])
                expected=z['features'][i,[0,6,2,4,0,8,3,5]]
                field_errors.append(float(np.max(np.abs(actual-expected))));rows+=1
    assert max(field_errors)<2e-14
    receipt.update(old_prefix_reconstruction_checks=2048,old_cell00_census_comparisons=rows,max_old_feature_absolute_error=max(field_errors),descriptor_npz_sha256=source_hashes)
    # Verify full-population covariance normalization on deterministic values.
    population=1000;m=37;x=np.arange(m,dtype=float)
    muC=.6+.02*np.sin(x);muW=.1+.01*np.cos(x)
    rc=1e-5*np.cos(.7*x);rw=2e-5*np.sin(.3*x)
    d={'cell00.mass':m/population}
    for ori in ('first','second'):
        for name,value in [('mu_C',muC),('mu_W',muW),('residual_C',rc),('residual_W',rw),('muC_residualC',muC*rc),('muW_residualW',muW*rw)]:
            d[ori+'.'+name]=float(value.sum()/population)
    got=derive(list(d.values()),list(d),population,model['sizes']['325']['point'])['R.primary_receiver_mean']
    expected=m/population*(2*np.cov(muC,rc,ddof=1)[0,1]-.5*np.cov(muW,rw,ddof=1)[0,1])
    assert abs(got-expected)<1e-21
    receipt['deterministic_covariance_normalization_error']=abs(got-expected)
    (ROOT/'QA.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))


if __name__=='__main__':
    main()

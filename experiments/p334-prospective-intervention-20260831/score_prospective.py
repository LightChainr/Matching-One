#!/usr/bin/env python3
"""The single prospective scorer, conditional on frozen archive coefficients."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.stats import binom, t

ROOT=Path(__file__).resolve().parent
OBS=("A_ref","E_ref","C","W","K1_over_Nplus1","K2_over_Nplus1","integral_A","integral_E")
ORIS=("first","second")
PREF=.59274605079


def outcomes(k,n):
    survival=binom.sf(np.arange(n+1)-1,n,PREF)
    f1,f2=survival[k[...,0]],survival[k[...,1]]
    x,y=k[...,0]/(n+1),k[...,1]/(n+1)
    c,w=(x+y)/2,y-x
    return np.stack((f1+f2-1,1-f1+f2,c,w,x,y,1-2*c,1-w),axis=-1)


def summarize_shard(stem,source):
    n=int(stem.split('.')[0][1:]);shard=int(stem.split('shard')[1])
    with gzip.open(source/(stem+'.prefix.csv.gz'),'rt') as f:
        all_rows=list(csv.DictReader(f))
    if len(all_rows)!=5000:
        raise ValueError('all-prefix denominator changed')
    selected=[r for r in all_rows if r['rank_first']=='0' and r['rank_second']=='0']
    ids=np.array([int(r['index']) for r in selected]);lookup={int(v):i for i,v in enumerate(ids)};m=len(ids)
    weights=np.array([[float(r['W_first']),float(r['W_second'])] for r in selected])
    if weights.shape!=(m,2) or np.any(weights<0):
        raise ValueError('invalid exact contrast weights')
    baseline=np.full((m,32,2,2),-1,dtype=np.int32)
    arms=np.full((m,2,64,2,2,2),-1,dtype=np.int32)
    seen_baseline=np.zeros((m,32),bool);seen_arms=np.zeros((m,2,64,2),bool)
    classes=np.full((m,2,64,2,2),-1,dtype=np.int32)
    with gzip.open(source/(stem+'.tails.csv.gz'),'rt') as f:
        for row in csv.DictReader(f):
            idx=lookup[int(row['index'])];draw=int(row['draw']);arm=int(row['arm'])
            if int(row['batch'])!=shard or int(row['counter'])!=int(selected[idx]['counter']):
                raise ValueError('tail/prefix identity mismatch')
            k=np.array([[int(row['first_k1']),int(row['first_k2'])],[int(row['second_k1']),int(row['second_k2'])]])
            if np.any(k[:,0]>k[:,1]) or np.any(k<1) or np.any(k>n):
                raise ValueError('invalid birth ranks')
            if row['group']=='B':
                if int(row['source'])!=-1 or arm!=0 or not 0<=draw<32 or seen_baseline[idx,draw]:
                    raise ValueError('baseline draw duplication/schema')
                baseline[idx,draw]=k;seen_baseline[idx,draw]=True
            elif row['group']=='I':
                s=int(row['source'])
                if s not in (0,1) or not 0<=draw<64 or arm not in (0,1) or seen_arms[idx,s,draw,arm]:
                    raise ValueError('contrast draw duplication/schema')
                arms[idx,s,draw,arm]=k;seen_arms[idx,s,draw,arm]=True
                classes[idx,s,draw,arm]=[int(row['e_first']),int(row['e_second'])]
                if np.any(k[:,0]<=int(selected[idx]['k0'])+1):
                    raise ValueError('rank was not preserved by intervention')
            else:
                raise ValueError('unknown group')
    if not np.all(seen_baseline):
        raise ValueError('incomplete independent baseline group')
    mu=outcomes(baseline,n).mean(axis=1)
    tau=np.zeros((m,2,2,len(OBS)))
    for s in (0,1):
        take=weights[:,s]>0
        if not np.all(seen_arms[take,s]) or np.any(seen_arms[~take,s]):
            raise ValueError('contrast budget/support mismatch')
        if np.any(classes[take,s,:,0]!=classes[take,s,:,1]):
            raise ValueError('paired degree classes differ')
        f=outcomes(arms[take,s],n)
        tau[take,s]=weights[take,s,None,None]*(f[:,:,0]-f[:,:,1]).mean(axis=1)
    with gzip.open(source/(stem+'.prediction.csv.gz'),'rt') as f:
        pred_rows={int(r['index']):r for r in csv.DictReader(f)}
    if len(pred_rows)!=5000:
        raise ValueError('prediction completeness changed')
    prediction=np.zeros((m,2,2,2))  # prefix,source,receiver,C/W
    for i,index in enumerate(ids):
        r=pred_rows[int(index)]
        if int(r['counter'])!=int(selected[i]['counter']):
            raise ValueError('prediction identity mismatch')
        for s,source_name in enumerate(ORIS):
            for o,ori in enumerate(ORIS):
                prediction[i,s,o]=[float(r[f'{ori}_{source_name}_{f}']) for f in ('C','W')]
    fields={'cell00.mass':m/5000}
    for o,ori in enumerate(ORIS):
        bc,bw=mu[:,o,2],mu[:,o,3]
        rc=tau[:,o,o,2]-prediction[:,o,o,0]
        rw=tau[:,o,o,3]-prediction[:,o,o,1]
        for name,value in [('mu_C',bc),('mu_W',bw),('residual_C',rc),('residual_W',rw),('muC_residualC',bc*rc),('muW_residualW',bw*rw)]:
            fields[f'{ori}.{name}']=float(value.sum()/5000)
        for s,source_name in enumerate(ORIS):
            for j,obs in enumerate(OBS):
                fields[f'response.{source_name}.{ori}.{obs}']=float(tau[:,s,o,j].sum()/5000)
            pc,pw=prediction[:,s,o,0],prediction[:,s,o,1]
            for obs,v in [('C',pc),('W',pw),('K1_over_Nplus1',pc-pw/2),('K2_over_Nplus1',pc+pw/2),('integral_A',-2*pc),('integral_E',-pw)]:
                fields[f'prediction.{source_name}.{ori}.{obs}']=float(v.sum()/5000)
    counts=[sum(int(r['rank_first'])*3+int(r['rank_second'])==j for r in all_rows) for j in range(9)]
    result={'N':n,'shard':shard,'all_prefixes':5000,'cell_counts':counts,'labels':list(fields),'batch_means':list(fields.values()),
            'baseline_paths':int(seen_baseline.sum()),'contrast_paths':int(seen_arms.sum()),'zero_source_counts':(weights==0).sum(axis=0).tolist(),
            'new_fit_calls':0,'conditioning':'archive beta, means and candidate numerical point values fixed'}
    target=source/(stem+'.sufficient.json')
    if target.exists():
        raise ValueError('refuse overwrite sufficient statistics')
    target.write_text(json.dumps(result,indent=2)+'\n')
    np.savez_compressed(source/(stem+'.prefix_scores.npz'),index=ids,mu=mu,tau=tau,prediction=prediction,weights=weights)
    return result


def derive(row,labels,population,old):
    d=dict(zip(labels,row));pi=d['cell00.mass'];count=int(round(pi*population))
    if count<2:
        raise ValueError('fewer than two00 prefixes; no repair')
    out={}
    for ori in ORIS:
        covc=d[f'{ori}.muC_residualC']/pi-d[f'{ori}.mu_C']*d[f'{ori}.residual_C']/pi**2
        covw=d[f'{ori}.muW_residualW']/pi-d[f'{ori}.mu_W']*d[f'{ori}.residual_W']/pi**2
        out[f'R.{ori}']=pi*count/(count-1)*(2*covc-.5*covw)
    out['R.primary_receiver_mean']=(out['R.first']+out['R.second'])/2
    old_mean=(old['first']['R_old']+old['second']['R_old'])/2
    out['R.primary_over_frozen_R_old']=out['R.primary_receiver_mean']/old_mean
    out['cell00.mass']=pi
    for key,value in d.items():
        if key.startswith(('response.','prediction.')):
            out[key]=value
    return out


def collect(source,output):
    model=json.loads((ROOT/'existing_model.json').read_text())
    freeze=json.loads((ROOT/'FREEZE.json').read_text())
    if hashlib.sha256((ROOT/'existing_model.json').read_bytes()).hexdigest()!=freeze['file_sha256']['existing_model.json']:
        raise ValueError('frozen model changed before validation')
    result={'status':'completed_prospective_fixed_prediction_test','model_sha256':hashlib.sha256((ROOT/'existing_model.json').read_bytes()).hexdigest(),
            'new_fit_calls':0,'primary_inference':'conditional on frozen training coefficients, means and R_old point forecasts; old observations do not enter validation score',
            'confidence':'97.5% marginal two-sided Student-t60batch intervals for the two prespecified size-specific primary scores; Bonferroni family coverage95%',
            'candidate_bands':{'C0_residual_projection_near_zero':[-.25,.25],'C1_old_point_residual_transfers':[.75,1.25]},'sizes':{}}
    for n in (325,425):
        paths=[source/f'N{n}.shard{s:03d}.sufficient.json' for s in range(60)]
        rows=[json.loads(p.read_text()) for p in paths]
        labels=rows[0]['labels']
        for s,r in enumerate(rows):
            if r['N']!=n or r['shard']!=s or r['all_prefixes']!=5000 or r['labels']!=labels:
                raise ValueError('missing/mismatched fixed shard')
        raw=np.array([r['batch_means'] for r in rows]);mean=raw.mean(axis=0)
        old=model['sizes'][str(n)]['point'];point=derive(mean,labels,300000,old)
        loo=np.array([list(derive((60*mean-r)/59,labels,295000,old).values()) for r in raw])
        factor=np.sqrt(59/60)*(loo-loo.mean(axis=0));se=np.linalg.norm(factor,axis=0)
        j=list(point).index('R.primary_over_frozen_R_old');estimate=list(point.values())[j]
        critical=float(t.ppf(.9875,59));ci=[float(estimate-critical*se[j]),float(estimate+critical*se[j])]
        decisions={}
        for name,(lo,hi) in result['candidate_bands'].items():
            decisions[name]='compatible_with_prespecified_tolerance' if lo<=ci[0] and ci[1]<=hi else ('excluded_by_primary_interval' if ci[1]<lo or ci[0]>hi else 'unresolved_at_fixed_budget')
        result['sizes'][str(n)]={'all_prefixes':300000,'cell_counts':np.array([r['cell_counts'] for r in rows]).sum(axis=0).tolist(),
                                'labels':list(point),'estimate':list(point.values()),'se':se.tolist(),'LOO':loo.tolist(),'factor':factor.tolist(),
                                'primary_ratio_interval':ci,'decisions':decisions,
                                'baseline_paths':sum(r['baseline_paths'] for r in rows),'contrast_paths':sum(r['contrast_paths'] for r in rows)}
    if output.exists():
        raise ValueError('refuse overwrite final analysis')
    output.mkdir(parents=True)
    (output/'latest.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--summarize-shard')
    p.add_argument('--input',type=Path,required=True)
    p.add_argument('--output',type=Path,default=ROOT/'results')
    args=p.parse_args()
    if args.summarize_shard:
        r=summarize_shard(args.summarize_shard,args.input)
        print(json.dumps({k:r[k] for k in ('N','shard','all_prefixes','cell_counts','baseline_paths','contrast_paths')}))
    else:
        r=collect(args.input,args.output)
        print(json.dumps({n:s['decisions'] for n,s in r['sizes'].items()},indent=2))

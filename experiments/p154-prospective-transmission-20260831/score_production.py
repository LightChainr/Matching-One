#!/usr/bin/env python3
"""Score complete fresh production once, with fixed model/stop rules."""
import argparse,csv,gzip,hashlib,json,math,time
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.stats import norm
from archive_channel_split import moments,kernels,channel_split
ROOT=Path(__file__).resolve().parent
NS=(85,340);BATCHES=200
PRIMARY=tuple(f'N{n}.{channel}.v' for n in NS for channel in ('entry','completion','total'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def decisions(estimates,contract):
    critical=float(norm.ppf(1-contract['inference']['alpha']/(2*len(PRIMARY))))
    intervals={k:[float(estimates[k]['value']-critical*estimates[k]['se']),float(estimates[k]['value']+critical*estimates[k]['se'])] for k in PRIMARY}
    models={}
    for name,model in contract['models'].items():
        contradictions=[];wholly_inside=True
        for n in NS:
            for channel in ('entry','completion'):
                lower,upper=model[channel+'_interval'];low,high=intervals[f'N{n}.{channel}.v']
                if lower is not None and high<lower or upper is not None and low>upper:contradictions.append(f'N{n}.{channel}.v')
                wholly_inside &= (lower is None or low>=lower) and (upper is None or high<=upper)
        models[name]={'contradicted':bool(contradictions),'contradicting_coordinates':contradictions,
                      'intervals_wholly_inside':bool(wholly_inside),'status':'numerical_prediction_contradicted' if contradictions else 'not_excluded'}
    eps=contract['net_transmission_primary']['equivalence_interval'][1]
    net_equiv=all(-eps<=intervals[f'N{n}.total.v'][0] and intervals[f'N{n}.total.v'][1]<=eps for n in NS)
    positive=all(intervals[f'N{n}.total.v'][0]>eps for n in NS)
    negative=all(intervals[f'N{n}.total.v'][1]<-eps for n in NS)
    return {'critical_value':critical,'simultaneous_intervals':intervals,'numerical_templates':models,
            'remaining_templates':[k for k,v in models.items() if not v['contradicted']],
            'all_three_numerical_templates_contradicted':all(v['contradicted'] for v in models.values()),
            'net_equivalent_within_declared_band_at_both_N':net_equiv,
            'stop_current_main_H4_priority_for_this_source':net_equiv,
            'net_above_positive_target_at_both_N':positive,'net_below_negative_target_at_both_N':negative,
            'interpretation':'fixed finite-size numerical restrictions; no forced winner, no physical-theory identification, no fourth template or sample extension'}
def load_fresh(directory,contract):
    arrays={n:np.zeros((BATCHES,2,n+1,11)) for n in NS};seen={n:set() for n in NS};receipts=[];freeze=None
    fields=('sum_q','sum_e','sum_s','sum_qs','sum_es','event_count01','event_count02','event_count12','sum_s_previous01','sum_s_previous02','sum_s_previous12')
    paths=sorted(directory.glob('n*-b*.csv.gz'))
    if len(paths)!=9:raise ValueError('Require exactly the frozen nine completed shards; no partial/interim scoring')
    for path in paths:
        receipt_path=path.with_name(path.name.replace('.csv.gz','.run.json'));r=json.loads(receipt_path.read_text())
        if r['status']!='completed' or r['gzip_sha256']!=sha(path):raise ValueError('Incomplete or altered fresh shard')
        if freeze is None:freeze=r['freeze_commit']
        if r['freeze_commit']!=freeze:raise ValueError('Mixed freeze commits')
        for name,digest in r['frozen_sha256'].items():
            if sha(ROOT/name)!=digest:raise ValueError('Frozen scorer/input changed: '+name)
        n=r['N'];per_batch=contract['samples_per_N'][str(n)]//BATCHES
        with gzip.open(path,'rt') as f:
            for row in csv.DictReader(f):
                b,g,k=int(row['batch']),('first','second').index(row['orientation']),int(row['k']);idx=(b,g,k)
                if idx in seen[n] or int(row['n'])!=n or not r['batch_begin']<=b<r['batch_end'] or not 0<=k<=n or int(row['samples'])!=per_batch:
                    raise ValueError('Duplicate, inconsistent or unplanned fresh sample row')
                seen[n].add(idx);arrays[n][idx]=[int(row[f]) for f in fields]
        receipts.append({'path':str(receipt_path),'sha256':sha(receipt_path),'samples':r['samples'],'N':n})
    for n in NS:
        if len(seen[n])!=BATCHES*2*(n+1):raise ValueError('Missing planned batch/K cells')
        per_batch=contract['samples_per_N'][str(n)]//BATCHES;c=arrays[n][...,5:8]
        q=-per_batch+np.cumsum(c[...,0]+2*c[...,1]+c[...,2],axis=-1)
        e=per_batch+np.cumsum(-c[...,0]+c[...,2],axis=-1)
        if not np.array_equal(q,arrays[n][...,0]) or not np.array_equal(e,arrays[n][...,1]):raise ValueError('First/completion event identity failed')
        if sum(r['samples'] for r in receipts if r['N']==n)!=contract['samples_per_N'][str(n)]:raise ValueError('Wrong frozen sample total')
    return arrays,receipts,freeze
def evaluate(data,samples,n,contract):
    base=data[...,:5];events=data[...,5:];bracket=contract['root_bracket']
    q_at=lambda p:sum(moments(base[g],samples,p,n)[0][0] for g in range(2))/2
    if q_at(bracket[0])*q_at(bracket[1])>0:raise ValueError('Frozen root bracket failed; do not retune after reveal')
    root=brentq(q_at,*bracket,xtol=5e-14,rtol=5e-14)
    return channel_split(base,kernels(base,events,samples),samples,n,root)
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--input-dir',type=Path,default=ROOT/'production');args=parser.parse_args()
    dest=ROOT/'PROSPECTIVE_RESULT.json'
    if dest.exists():raise RuntimeError('Final score already exists; no overwrite/rescore after changing rules')
    started=time.perf_counter();contract=json.loads((ROOT/'CONTRACT.json').read_text())
    raw,receipts,freeze=load_fresh(args.input_dir,contract);totals={n:raw[n].sum(axis=0) for n in NS}
    samples={n:contract['samples_per_N'][str(n)] for n in NS};points={n:evaluate(totals[n],samples[n],n,contract) for n in NS}
    vector=lambda rows:{f'N{n}.{k}':v for n in NS for k,v in rows[n].items()}
    central=vector(points);labels=list(central);cov=np.zeros((len(labels),len(labels)));groups={}
    for n in NS:
        vectors=[]
        for b in range(BATCHES):
            changed=dict(points);changed[n]=evaluate(totals[n]-raw[n][b],samples[n]*(BATCHES-1)/BATCHES,n,contract)
            vectors.append(list(vector(changed).values()))
        vectors=np.asarray(vectors);delta=vectors-vectors.mean(axis=0);factor=np.sqrt((BATCHES-1)/BATCHES)*delta;cov+=factor.T@factor
        groups[f'fresh_N{n}']={'Ns':[n],'delete_one_batch_ids':list(range(BATCHES)),'delete_one_vectors':vectors.tolist()}
    se=np.sqrt(np.maximum(0,cov.diagonal()));est={k:{'value':float(v),'se':float(s),'z':float(v/s) if s>0 else None} for k,v,s in zip(labels,central.values(),se)}
    selected=[labels.index(k) for k in PRIMARY];primary_cov=cov[np.ix_(selected,selected)]
    if not np.isfinite(primary_cov).all() or not np.isfinite(list(central.values())).all():raise ValueError('Nonfinite final score')
    variance_comparison={}
    for n in NS:
        plan=contract['power_plan']['old_variances'][str(n)]
        for channel in ('entry','completion','total'):
            expected=1.25*plan[channel]*plan['samples']/samples[n];actual=est[f'N{n}.{channel}.v']['se']**2
            variance_comparison[f'N{n}.{channel}.v']={'protected_planned_variance':expected,'observed_variance':actual,'ratio':actual/expected,'above_plan':actual>expected}
    result={'status':'FRESH_INDEPENDENT_FIXED_BUDGET_COMPLETED','freeze_commit':freeze,'contract':contract,'sample_totals':samples,
            'old_data_pooled':False,'labels':labels,'estimates':est,'covariance':cov.tolist(),'covariance_groups':groups,
            'primary_labels':list(PRIMARY),'primary_covariance':primary_cov.tolist(),'decision':decisions(est,contract),
            'variance_comparison':variance_comparison,'receipts':receipts,'elapsed_seconds':time.perf_counter()-started,
            'code_sha256':sha(Path(__file__)),'contract_sha256':sha(ROOT/'CONTRACT.json')}
    dest.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result['decision'],indent=2))
if __name__=='__main__':main()

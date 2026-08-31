#!/usr/bin/env python3
"""Additive soft-angular allocation of the existing norm-4 U/source response."""
from __future__ import annotations
import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import platform
import resource
import socket
import sys
import time
from pathlib import Path
from fractions import Fraction

import numpy as np
import scipy
from frozen_moments import binomial_moments

ROOT=Path(__file__).resolve().parent
NS=(65,85,130,170,260,340)
LINEAGES=((65,130,260),(85,170,340))
SOURCE_FIELDS=('sum_q','sum_e','sum_s','sum_qs','sum_es')
LINE_FIELDS=('sum_rank1','sum_rank1_s','sum_line4_re','sum_line4_im','sum_line4_s_re','sum_line4_s_im')


def payload(short):
    return gzip.decompress((ROOT/'inputs'/f'{short}.gz').read_bytes())


def read_profile(short,n,per_batch,run,fields):
    out=np.zeros((100,2,n+1,len(fields)))
    seen=set()
    for row in csv.DictReader(io.StringIO(payload(short).decode())):
        g=('first','second').index(row['orientation'])
        batch,k=int(row['batch']),int(row['k']);key=(batch,g,k)
        assert key not in seen and 0<=batch<100 and 0<=k<=n
        assert int(row['n'])==n and int(row['samples'])==per_batch
        assert [int(row['a']),int(row['b'])]==run[row['orientation']]
        seen.add(key);out[batch,g,k]=[float(row[f]) for f in fields]
    assert len(seen)==100*2*(n+1) and np.isfinite(out).all()
    return out


def evaluate(n,source,line,samples,run,old):
    p=old['p0'];delta=float(Fraction(run['delta_cos4']))
    packs=[binomial_moments(source[g],samples,p,n) for g in range(2)]
    means=np.array([x[0] for x in packs]);first=np.array([x[1] for x in packs]);second=np.array([x[2] for x in packs])
    q,e,s,qs,es=means.T;qp,ep,sp,qsp,esp=first.T
    jq=qs-q*s;jqp=qsp-qp*s-q*sp
    d=qp.mean();qpp=second[:,0].mean();pdot=old['rootdot_fugacity']
    assert abs(-jq.mean()/d-pdot)<1e-9
    a=n**(13/8)/2
    phases=np.array([complex(*run[g]) for g in ('first','second')]).conj()**4/n**2
    transformed=np.empty((2,n+1,4))
    within_source=np.zeros((2,n+1,2))
    for g in range(2):
        r,rs,ore,oim,osre,osim=line[g].T
        h=(phases[g]*(ore+1j*oim)).real
        hs=(phases[g]*(osre+1j*osim)).real
        transformed[g]=np.column_stack(((r+h)/2,(r-h)/2,(rs+hs)/2,(rs-hs)/2))
        product=np.zeros(n+1)
        np.divide(h*rs,r,out=product,where=r>0)
        angular_within=hs-product
        angular_within[r<=1]=0
        within_source[g]=np.column_stack((angular_within/2,-angular_within/2))
        assert np.min(transformed[g,:,:2])>-1e-7
    pp=[binomial_moments(transformed[g],samples,p,n) for g in range(2)]
    within_packets=[binomial_moments(within_source[g],samples,p,n) for g in range(2)]
    within_p=np.array([x[1] for x in within_packets])
    w=np.array([x[0] for x in pp]);wp=np.array([x[1] for x in pp]);wpp=np.array([x[2] for x in pp])
    result={'p0':p,'D':d,'rootdot':pdot,'U_old':old['U'],'v_old':old['Udot_fugacity']}
    for j,sign in enumerate(('plus','minus')):
        b=(wp[0,j]-wp[1,j])/delta
        bpp=(wpp[0,j]-wpp[1,j])/delta
        jw=w[:,j+2]-w[:,j]*s
        jwp=wp[:,j+2]-wp[:,j]*s-w[:,j]*sp
        bsource=(jwp[0]-jwp[1])/delta
        terms={'direct':-a*bsource/d,'rootmotion':-a*bpp*pdot/d,
               'slope_source':a*b*jqp.mean()/d**2,'slope_root':a*b*qpp*pdot/d**2}
        result[f'U_{sign}']=-a*b/d
        result[f'v_{sign}']=math.fsum(terms.values())
        result[f'v_{sign}_rank1_within_source']=-a*(within_p[0,j]-within_p[1,j])/(delta*d)
        result[f'v_{sign}_remaining_source']=result[f'v_{sign}']-result[f'v_{sign}_rank1_within_source']
        result.update({f'v_{sign}_{name}':value for name,value in terms.items()})
        for g,direction in enumerate(('first','second')):
            result[f'{direction}.W_{sign}']=w[g,j]
            result[f'{direction}.Wp_{sign}']=wp[g,j]
            result[f'{direction}.Wsource_{sign}']=jw[g]+pdot*wp[g,j]
    for base in ('U','v'):
        result[f'{base}_sum']=result[f'{base}_plus']+result[f'{base}_minus']
        result[f'{base}_angular_contrast']=result[f'{base}_plus']-result[f'{base}_minus']
        result[f'{base}_sum_residual']=result[f'{base}_sum']-result[f'{base}_old']
        tolerance=2e-7*max(1,abs(result[f'{base}_old']))
        assert abs(result[f'{base}_sum_residual'])<tolerance,(n,base,result[f'{base}_sum_residual'])
    result['v_rank1_within_sum']=result['v_plus_rank1_within_source']+result['v_minus_rank1_within_source']
    assert result['v_rank1_within_sum']==0
    return {key:float(value) for key,value in result.items()}


def vectorize(points):
    out={f'N{n}.{key}':v for n in NS for key,v in points[n].items()}
    for chain in LINEAGES:
        for model,coeff in (('q2',(1,-3,2)),('Jordan',(1,-2,1))):
            for quantity in ('U','v'):
                for part in ('plus','minus','sum','old','angular_contrast'):
                    out[f'{model}.{chain[0]}.{quantity}_{part}']=sum(c*points[n][f'{quantity}_{part}'] for c,n in zip(coeff,chain))
    return out


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,default=ROOT/'results');args=ap.parse_args()
    if (args.output/'latest.json').exists():raise RuntimeError('Existing result; refusing repeat or overwrite')
    started=time.perf_counter();manifest=json.loads((ROOT/'SOURCES.json').read_text())
    for item in manifest['inputs']:
        data=(ROOT/item['local']).read_bytes()
        assert hashlib.sha256(data).hexdigest()==item['compressed_sha256']
        assert hashlib.sha256(gzip.decompress(data)).hexdigest()==item['uncompressed_sha256']
    saved=json.loads(payload('source_result.json'));contract=json.loads(payload('line_contract.json'))
    runs={r['N']:r for r in json.loads(payload('geometry.json'))['runs']}
    profiles,lines,totals,line_totals,points,inputs_checks={},{},{},{},{},{}
    per_batch={n:contract['marked_permutations_by_N'][str(n)]//100 for n in NS}
    indices={label:i for i,label in enumerate(saved['labels'])}
    saved_fields=('p0','rootdot_fugacity','U','Udot_fugacity')
    for n in NS:
        profiles[n]=read_profile(f'source_n{n}.csv',n,1000,runs[n],SOURCE_FIELDS)
        if n in (260,340):profiles[n]+=read_profile(f'increment_n{n}.csv',n,9000,runs[n],SOURCE_FIELDS)
        lines[n]=read_profile(f'line_n{n}.csv',n,per_batch[n],runs[n],LINE_FIELDS)
        assert np.array_equal(lines[n][...,0],per_batch[n]-profiles[n][...,1])
        assert np.array_equal(lines[n][...,1],profiles[n][...,2]-profiles[n][...,4])
        totals[n]=profiles[n].sum(axis=0);line_totals[n]=lines[n].sum(axis=0)
        old={f:saved['estimates'][f'N{n}.{f}']['value'] for f in saved_fields}
        points[n]=evaluate(n,totals[n],line_totals[n],100*per_batch[n],runs[n],old)
        inputs_checks[n]={'R_equals_1_minus_E_every_batch_K':True,'Rs_equals_s_minus_Es_every_batch_K':True}
    central_map=vectorize(points);labels=list(central_map);central=np.array(list(central_map.values()))
    covariance=np.zeros((len(labels),len(labels)));contributions={};max_residual={'U':0.,'v':0.}
    for group in contract['dependency_groups']:
        gid,sizes=group['id'],group['Ns'];old_group=saved['covariance_contributions']['source:'+gid]
        assert old_group['Ns']==sizes and old_group['delete_one_batch_ids']==list(range(100))
        vectors=[]
        for batch in range(100):
            changed=dict(points)
            for n in sizes:
                assert old_group['batch_counts']==[per_batch[n]]*100
                old={f:old_group['delete_one_vectors'][batch][indices[f'N{n}.{f}']] for f in saved_fields}
                changed[n]=evaluate(n,totals[n]-profiles[n][batch],line_totals[n]-lines[n][batch],99*per_batch[n],runs[n],old)
                for quantity in ('U','v'):max_residual[quantity]=max(max_residual[quantity],abs(changed[n][f'{quantity}_sum_residual']))
            vectors.append(list(vectorize(changed).values()))
        vectors=np.array(vectors);centered=vectors-vectors.mean(axis=0);component=.99*centered.T@centered;covariance+=component
        contributions[gid]={'Ns':sizes,'batch_counts':old_group['batch_counts'],'delete_one_vectors':vectors.tolist(),'covariance':component.tolist()}
        print('completed aligned covariance',sizes,flush=True)
    errors=np.sqrt(np.maximum(0,np.diag(covariance)))
    estimates={key:{'value':float(v),'se':float(e),'z':float(v/e) if e else None} for key,v,e in zip(labels,central,errors)}
    allocations={}
    for n in NS:
        allocations[n]={}
        for quantity in ('U','v'):
            wanted=[f'N{n}.{quantity}_{part}' for part in ('plus','minus','sum')];ii=[labels.index(k) for k in wanted];cc=covariance[np.ix_(ii,ii)]
            allocations[n][quantity]={'labels':wanted,'estimates':[estimates[k] for k in wanted],'covariance':cc.tolist(),'plus_minus_correlation':float(cc[0,1]/math.sqrt(cc[0,0]*cc[1,1]))}
    result={'schema':'matching-one/p154-soft-angular-global-U/v1','status':'completed_existing_marks_no_replay',
            'source_commit':manifest['commit'],'labels':labels,'estimates':estimates,'covariance':covariance.tolist(),
            'covariance_contributions':contributions,'allocations':allocations,'by_N':points,'identity_checks':inputs_checks,
            'max_aligned_sum_residual':max_residual,'source_units':'bulk s=CB+CW, coupled by exp(t*s); no extra N multiplier',
            'uncertainty':'Same100 aligned batches, four cyclic N in one group and two separate endpoint groups. Each omission reuses the matching root/rootdot from the same source result and recomputes W, its jets, D and slope-source response. No independent complement anchor.',
            'scope':'Global-observable allocation on the marked100k/1M source subset; not full1.9B/1B precision. Soft nonnegative angular weights are not discrete winding classes; paired components and chain contrasts are dependent views, not new replications or an identified continuum field.',
            'source_result_sha256':hashlib.sha256(payload('source_result.json')).hexdigest(),
            'code_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),ROOT/'frozen_moments.py',ROOT/'SOURCES.json')},
            'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'machine':platform.machine(),'hostname':socket.gethostname()},
            'command':[sys.executable,*sys.argv],'elapsed_seconds':time.perf_counter()-started,'max_rss_native':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'configuration_replays':0,'new_random_counters':0,'root_solver_calls':0}
    args.output.mkdir(parents=True,exist_ok=True)
    with (args.output/'latest.json').open('x') as handle:json.dump(result,handle,indent=2,allow_nan=False);handle.write('\n')
    print(json.dumps({'elapsed_seconds':result['elapsed_seconds'],'max_sum_residual':max_residual,'allocations':allocations},indent=2))


if __name__=='__main__':main()

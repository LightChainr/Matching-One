#!/usr/bin/env python3
"""Fixed new64 and combined72 local-J readouts, using O(Q) fourth-order algebra."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
from itertools import combinations
import json
import math
from pathlib import Path
import platform
import time
import numpy as np
import scipy
from scipy.stats import binom
from analyze_prefix_local_rank import array_csv, determinant_polarization, P_REF

OBS=["p_ref.A","p_ref.E","p_integral.A","p_integral.E","score_G"]
NAMES=[obs+suffix for obs in OBS for suffix in [".E_det_JZ",".E_det_JZ_squared",".E_frobenius_JZ_squared"]+
       [f".mean_J[{o},{s}]" for o in ("first","second") for s in ("first","second")]]

def fast_moments(x):
    """All distinct-pair and disjoint-four-subset U statistics, without Q^4 work."""
    p,q,o=x.shape[:3]
    flat=x.reshape(p,q,o,4)
    sm=x.sum(axis=1)
    diagonal=np.linalg.det(x)
    edge_sum=(np.linalg.det(sm)-diagonal.sum(axis=1))/2
    row_sum=determinant_polarization(x,sm[:,None])-diagonal
    m2=np.einsum('pqoi,pqoj->poij',flat,flat)
    matrix=np.array([[0.,0.,0.,1.],[0.,0.,-1.,0.],[0.,-1.,0.,0.],[1.,0.,0.,0.]])
    mm=np.einsum('ij,pojk->poik',matrix,m2)
    all_ordered_squares=np.einsum('poij,poji->po',mm,mm)/4
    edge_square_sum=(all_ordered_squares-(diagonal**2).sum(axis=1))/2
    det=edge_sum/math.comb(q,2)
    det2=(edge_sum**2+edge_square_sum-(row_sum**2).sum(axis=1))/(6*math.comb(q,4))
    energy=((sm**2).sum(axis=(-2,-1))-(x**2).sum(axis=(1,3,4)))/(q*(q-1))
    mean=x.mean(axis=1)
    fields=[]
    for j in range(o):
        fields.extend([det[:,j],det2[:,j],energy[:,j]])
        fields.extend(mean[:,j,i,k] for i in range(2) for k in range(2))
    return np.column_stack(fields)

def check_fast_algebra():
    v=np.array([[[1.,2.],[3.,4.]],[[2.,-1.],[0.,3.]],[[0.,2.],[-1.,1.]],[[3.,0.],[2.,-2.]],
                [[1.,-2.],[2.,0.]],[[0.,1.],[3.,-1.]],[[2.,2.],[1.,4.]],[[-1.,0.],[2.,3.]]])
    x=v[None,:,None]
    pol={ij:determinant_polarization(x[:,ij[0]],x[:,ij[1]]) for ij in combinations(range(8),2)}
    direct=sum(pol.values())/28
    fourth=0.
    for a,b,c,d in combinations(range(8),4):
        fourth+=pol[(a,b)]*pol[(c,d)]+pol[(a,c)]*pol[(b,d)]+pol[(a,d)]*pol[(b,c)]
    fourth/=210
    fast=fast_moments(x)
    error=max(float(abs(fast[0,0]-direct[0,0])),float(abs(fast[0,1]-fourth[0,0])))
    assert error<1e-12,error
    return error

def source_classes(census):
    groups={}
    for row in census:
        n,b,ctr,k0,r0,r1,e0,e1,l0,l1,count=map(int,row)
        if r0 or r1: continue
        groups.setdefault(ctr,{}).setdefault((e0,e1),[]).append((l0,l1,count))
    result={}
    for ctr,classes in groups.items():
        result[ctr]={}
        for key,rows in classes.items():
            mass=sum(r[2] for r in rows)
            totals=np.array([sum(r[j]*r[2] for r in rows) for j in (0,1)],dtype=np.int64)
            assert np.all(sum(r[2]*(mass*np.array(r[:2])-totals) for r in rows)==0)
            result[ctr][key]=(mass,totals,{r[:2] for r in rows})
    return result

def matrices(raw,hi,contact,ci,classes,n,new=False):
    p,q=raw.shape[:2]
    counters=raw[:,0,0,0,hi['counter']]
    k0=int(raw[0,0,0,0,hi['k0']]);vacant=n-k0
    h=np.zeros((p,q,2,2))
    for i,ctr in enumerate(counters):
        for j in range(q):
            for g in (0,1):
                row=contact[i,j,g]
                suffix='next_rank' if new else 'rank_after'
                if any(row[ci[f'{o}_{suffix}']] != 0 for o in ('first','second')):continue
                degree=tuple(int(row[ci[f'{o}_e']]) for o in ('first','second'))
                marks=tuple(int(row[ci[f'{o}_e']]-row[ci[f'{o}_c']]) for o in ('first','second'))
                mass,totals,support=classes[int(ctr)][degree]
                assert marks in support
                h[i,j,g]=(mass*np.array(marks)-totals)/vacant
    dh=h[:,:,0]-h[:,:,1]
    tail=binom.sf(np.arange(n+1)-1,n,P_REF)
    yy=[]
    for o in ('first','second'):
        k1,k2=(raw[...,hi[f'{o}_{key}']] for key in ('k1','k2'))
        assert np.all((k0<k1)&(k1<=k2)&(k2<=n))
        yy.append(np.stack((tail[k1]+tail[k2]-1,1-tail[k1]+tail[k2],1-(k1+k2)/(n+1),1+(k1-k2)/(n+1)),axis=-1))
    y=np.stack(yy,axis=-1)
    dy=y[:,:,0].mean(axis=2)-y[:,:,1].mean(axis=2)
    j=dy[..., :, :,None]*dh[...,None,None,:]/2
    g=dh[..., :,None]*dh[...,None,:]/2
    return np.concatenate((j,g[:,:,None]),axis=2)

def batch(task):
    root,n,b,census=task;root=Path(root)
    hi,raw=array_csv(root/'inputs/forks'/f'N{n}.batch{b:02}.csv.gz')
    ci,contact=array_csv(root/'inputs/contact'/f'N{n}.batch{b:02}.csv.gz')
    raw=raw[np.lexsort(tuple(raw[:,hi[k]] for k in ('replica','group','quartet','counter')))].reshape(1000,8,2,2,-1)
    contact=contact[np.lexsort(tuple(contact[:,ci[k]] for k in ('group','quartet','counter')))].reshape(1000,8,2,-1)
    take=(raw[:,0,0,0,hi['first_rank']]==0)&(raw[:,0,0,0,hi['second_rank']]==0)
    raw,contact=raw[take],contact[take]
    for key in ('counter','quartet','group','next_label'):
        assert np.array_equal(raw[...,0,hi[key]],contact[...,ci[key]])
    nh,nraw=array_csv(root/'extension'/f'N{n}.batch{b:02}.csv.gz')
    nraw=nraw[np.lexsort(tuple(nraw[:,nh[k]] for k in ('replica','group','quartet','counter')))].reshape(len(raw),64,2,2,-1)
    assert np.array_equal(raw[:,0,0,0,hi['counter']],nraw[:,0,0,0,nh['counter']])
    assert np.array_equal(nraw[0,:,0,0,nh['quartet']],np.arange(8,72))
    for key in ('next_label','first_next_rank','second_next_rank','first_e','first_c','second_e','second_c'):
        assert np.array_equal(nraw[...,0,nh[key]],nraw[...,1,nh[key]])
    classes=source_classes(census)
    old=matrices(raw,hi,contact,ci,classes,n)
    new=matrices(nraw,nh,nraw[:,:,:,0],nh,classes,n,new=True)
    original=fast_moments(old)
    reference=np.load(root/'results-exact-score'/f'prefix_statistics_N{n}.npz')
    rt=(reference['batch']==b)&(reference['cell']==0)
    assert np.array_equal(reference['counter'][rt],raw[:,0,0,0,hi['counter']])
    cols=[list(reference['labels']).index(k) for k in NAMES]
    difference=original-reference['values'][rt][:,cols]
    assert np.allclose(original,reference['values'][rt][:,cols],rtol=1e-9,atol=1e-18),float(np.max(abs(difference)))
    return {'N':n,'batch':b,'counter':raw[:,0,0,0,hi['counter']],
            'old8':original,'new64':fast_moments(new),'combined72':fast_moments(np.concatenate((old,new),axis=1)),
            'old8_vs_enumerated_max_error':float(np.max(abs(difference)))}

def derive(v,labels):
    out=dict(zip(labels,v))
    for mode in ('old8','new64','combined72'):
        for group in ['all']+[f'cell{a}{b}' for a in range(3) for b in range(3)]:
            for obs in OBS:
                base=f'{mode}.{group}.{obs}'
                j=np.array([[out[base+f'.mean_J[{o},{s}]'] for s in ('first','second')] for o in ('first','second')])
                out[base+'.det_mean_JZ']=float(np.linalg.det(j))
    return out

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers',type=int,default=16)
    parser.add_argument('--output',type=Path,default=Path(__file__).parent/'results-extension')
    args=parser.parse_args();root=Path(__file__).parent;started=time.time()
    args.output.mkdir(exist_ok=False)
    ext_receipt=json.loads((root/'extension/run_receipt.json').read_text())
    assert ext_receipt['new_tail_paths']==781568 and ext_receipt['prefixes']==3053
    for f in ext_receipt['files']:
        assert hashlib.sha256((root/'extension'/f['name']).read_bytes()).hexdigest()==f['sha256']
    for where in ('inputs','census','inputs/prefix_archive'):
        m=json.loads((root/where/'manifest.json').read_text())
        for f in m['files']:
            assert hashlib.sha256((root/where/f['local_path']).read_bytes()).hexdigest()==f['sha256']
    algebra_error=check_fast_algebra();tasks=[]
    for n in (325,425):
        ch,census=array_csv(root/'census'/f'N{n}/census.csv.gz')
        for b in range(20):tasks.append((str(root),n,b,census[census[:,1]==b]))
    with ProcessPoolExecutor(max_workers=args.workers) as pool:blocks=list(pool.map(batch,tasks))
    old=json.loads((root/'results-exact-score/score.json').read_text())
    result={'target':'E_Z det J(Z) and E_Z[(det J(Z))^2] on original population; all non00 contributions are exact zero.',
       'fixed_extension':'64 new iid quartets only at3053 existing00 prefixes;781568 new suffix paths; no new independent prefixes. q8..71 bit31 stream domain.',
       'modes':'old8, new64, combined72. For all-population mean matrices/energies, non00 retains original8; local determinants outside00 are zero. Every cell uses original1000-per-batch denominator.',
       'dependence':'All modes share original prefixes/batches; combined72 includes old8 and new64. Joint20 delete-one factors retain their dependence. New64 is not an independent population replication.',
       'fast_identity':'D(q,r)=polarized determinant. S=sum_edges D, T=sum_edges D^2, R=sum_vertices(row_sum D)^2. U_det2=(S^2+T-R)/(6 choose(Q,4)). T computed from4x4 sample second-moment matrix; O(Q), not Q^4.',
       'algebra_vs_direct_four_subset_error':algebra_error,'sizes':{}}
    for n in (325,425):
        bb=[b for b in blocks if b['N']==n];s=old['sizes'][str(n)]
        bx=np.array(s['joint_20_batch_means']);bi={k:i for i,k in enumerate(s['base_labels'])}
        labels=[];matrix=[]
        for mode in ('old8','new64','combined72'):
            for group in ['all']+[f'cell{a}{b}' for a in range(3) for b in range(3)]:
                for k,key in enumerate(NAMES):
                    oldv=bx[:,bi[f'{group}.{key}']]
                    if mode=='old8' or group not in ('all','cell00'):val=oldv
                    else:
                        fresh=np.array([b[mode][:,k].sum()/1000 for b in bb])
                        val=fresh if group=='cell00' else oldv+fresh-bx[:,bi[f'cell00.{key}']]
                    labels.append(f'{mode}.{group}.{key}');matrix.append(val)
        x=np.column_stack(matrix);mean=x.mean(0);primary=derive(mean,labels)
        loo=np.array([list(derive((20*mean-row)/19,labels).values()) for row in x]);factor=np.sqrt(19/20)*(loo-loo.mean(0))
        comparisons={}
        for obs in OBS[:4]:
            for key in ('.E_det_JZ','.E_det_JZ_squared'):
                label=obs+key
                z=x[:,labels.index('new64.all.'+label)]-x[:,labels.index('old8.all.'+label)]
                comparisons[label]={'new64_minus_old8':float(z.mean()),'paired_se':float(z.std(ddof=1)/np.sqrt(20)),'batch_differences':z.tolist()}
        np.savez_compressed(args.output/f'prefix_statistics_N{n}.npz',counter=np.concatenate([b['counter'] for b in bb]),batch=np.concatenate([np.full(len(b['counter']),b['batch']) for b in bb]),
            labels=np.array(NAMES),old8=np.concatenate([b['old8'] for b in bb]),new64=np.concatenate([b['new64'] for b in bb]),combined72=np.concatenate([b['combined72'] for b in bb]))
        result['sizes'][str(n)]={'batch_ids':list(range(20)),'labels':list(primary),'estimate':list(primary.values()),'se':np.linalg.norm(factor,axis=0).tolist(),'LOO':loo.tolist(),'factor':factor.tolist(),
            'base_labels':labels,'joint_20_batch_means':x.tolist(),'cell_counts':s['cell_counts'],'new64_vs_old8':comparisons,'fast_vs_old8_enumeration_max_error':max(b['old8_vs_enumerated_max_error'] for b in bb)}
        for k,v,e in zip(primary,primary.values(),np.linalg.norm(factor,axis=0)):
            if k.startswith(('new64.all.','combined72.all.')) and ('E_det_JZ' in k or ('det_mean_JZ' in k and '.A.' in k)):
                print(n,k,f'{v:.12g} +/- {e:.6g}',flush=True)
    (args.output/'score.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    receipt={'host':platform.node(),'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'workers':args.workers,'started_unix':started,'finished_unix':time.time(),'elapsed_seconds':time.time()-started,
       'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'extension_receipt_sha256':hashlib.sha256((root/'extension/run_receipt.json').read_bytes()).hexdigest(),'new_sampling_in_this_analysis':0,'analysis_attempt':1}
    (args.output/'run_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)

if __name__=='__main__':main()

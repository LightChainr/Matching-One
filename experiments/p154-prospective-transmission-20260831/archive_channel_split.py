#!/usr/bin/env python3
"""Planning-only old-archive readout split; never generates random samples."""
import csv,gzip,hashlib,json,math
from fractions import Fraction
from pathlib import Path
import numpy as np
from scipy.special import gammaln

ROOT=Path(__file__).resolve().parent
NS=(65,85,130,170,260,340)
def moments(sums,samples,p,n):
    k=np.arange(n+1,dtype=float)
    w=np.exp(gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1)+k*math.log(p)+(n-k)*math.log1p(-p))
    score=k/p-(n-k)/(1-p);wp=w*score;wpp=w*(score**2-k/p**2-(n-k)/(1-p)**2)
    z,zp,zpp=w.sum(),wp.sum(),wpp.sum();m=w@sums/(samples*z)
    mp=(wp@sums/samples-m*zp)/z;mpp=(wpp@sums/samples-m*zpp-2*mp*zp)/z
    return m,mp,mpp
def kernels(base,events,samples):
    q,e,s,qs,es=np.moveaxis(base,-1,0)
    mu0=np.divide(es-qs,e-q,out=np.zeros_like(q),where=e>q)
    mu1=np.divide(s-es,samples-e,out=np.zeros_like(q),where=samples>e)
    shift=lambda x:np.pad(x[:,:-1],((0,0),(1,0)))
    t=events[...,3:].copy();t[...,0]-=shift(mu0)*events[...,0];t[...,1]-=shift(mu0)*events[...,1];t[...,2]-=shift(mu1)*events[...,2]
    return t
def channel_split(base,t,samples,n,p0):
    profile=base.copy();profile[...,2]=0;profile[...,3]=t[...,0]+2*t[...,1]+t[...,2];profile[...,4]=-t[...,0]+t[...,2]
    m,p,pp=(np.asarray(v) for v in zip(*(moments(profile[g],samples,p0,n) for g in range(2))))
    q,e,zero,jq,je=m.T;qp,ep,_,jqp,jep=p.T;qpp,epp=pp[:,:2].T
    delta=float(Fraction('1152/845' if n in (65,130,260) else '2304/1445'));p4=lambda v:(v[0]-v[1])/delta
    d=qp.mean();rootdot=-jq.mean()/d;ddot=jqp.mean()+rootdot*qpp.mean();a=n**(13/8)/2
    values={'p0':float(p0),'rootdot':float(rootdot),'rank1_rootdot':float(-je.mean()-rootdot*ep.mean())}
    for name,sign,f_p,f_pp,j_p in (
      ('entry',-1,(qp-ep)/2,(qpp-epp)/2,(jqp-jep)/2),
      ('completion',1,(qp+ep)/2,(qpp+epp)/2,(jqp+jep)/2)):
        terms=(sign*a*p4(j_p)/d,sign*a*rootdot*p4(f_pp)/d,-sign*a*p4(f_p)*ddot/d**2)
        values[name+'.U']=float(sign*a*p4(f_p)/d)
        values[name+'.v']=float(math.fsum(terms))
        for key,v in zip(('direct','root','slope'),terms):values[name+'.v_'+key]=float(v)
    values['total.U']=values['entry.U']+values['completion.U'];values['total.v']=values['entry.v']+values['completion.v']
    values['difference.v']=values['completion.v']-values['entry.v']
    return values
def load(n):
    out=np.zeros((100,2,n+1,6));fields=['event_count'+t for t in ('01','02','12')]+['sum_s_previous'+t for t in ('01','02','12')]
    seen=set()
    with gzip.open(ROOT/f'inputs/lag1-n{n}.csv.gz','rt') as f:
        for r in csv.DictReader(f):
            key=(int(r['batch']),('first','second').index(r['orientation']),int(r['k']))
            assert key not in seen and int(r['n'])==n;seen.add(key);out[key]=[int(r[k]) for k in fields]
    assert len(seen)==100*2*(n+1)
    return out
def main():
    anchors=json.loads((ROOT/'inputs/anchors.json').read_text());old=json.loads((ROOT/'inputs/lag1-latest.json').read_text())
    bases=np.load(ROOT/'inputs/old_profiles.npz');base={n:bases[f'n{n}'].astype(float) for n in NS};events={n:load(n) for n in NS}
    bs={n:base[n].sum(axis=0) for n in NS};ts={n:events[n].sum(axis=0) for n in NS};samples={n:1000000 if n in (260,340) else 100000 for n in NS}
    points={n:channel_split(bs[n],kernels(bs[n],ts[n],samples[n]),samples[n],n,anchors['by_N'][str(n)]['p0']) for n in NS}
    def vec(rows):return {f'N{n}.{k}':v for n in NS for k,v in rows[n].items()}
    central=vec(points);labels=list(central);cov=np.zeros((len(labels),len(labels)));groups={};err=0.
    for n in NS:err=max(err,abs(points[n]['total.v']-old['estimates'][f'N{n}.total.v']['value']))
    old_idx={l:i for i,l in enumerate(old['labels'])}
    for name,g in anchors['groups'].items():
        vectors=[]
        for b in range(100):
            changed=dict(points)
            for n in g['Ns']:
                bb=bs[n]-base[n][b];tt=ts[n]-events[n][b];m=samples[n]*.99;p0=g['by_N'][str(n)]['p0'][b]
                changed[n]=channel_split(bb,kernels(bb,tt,m),m,n,p0)
                saved=old['covariance_contributions'][name]['delete_one_vectors'][b][old_idx[f'N{n}.total.v']]
                err=max(err,abs(changed[n]['total.v']-saved))
            vectors.append(list(vec(changed).values()))
        vectors=np.asarray(vectors);d=vectors-vectors.mean(axis=0);cov+=.99*d.T@d
        groups[name]={'Ns':g['Ns'],'delete_one_vectors':vectors.tolist()}
    se=np.sqrt(np.maximum(0,cov.diagonal()));est={k:{'value':float(v),'se':float(s),'z':float(v/s) if s>0 else None} for k,v,s in zip(labels,central.values(),se)}
    assert err<1e-8,err
    result={'status':'OLD_ARCHIVE_PLANNING_ONLY_NO_NEW_RANDOM_SAMPLE','source_commit':'4daae57eef5c945aa050a95cd3d5d5d77582161b',
            'source':'one-activation previous-rank centered bulk cluster count','labels':labels,'estimates':est,'covariance':cov.tolist(),'groups':groups,
            'samples_per_N':samples,'max_total_v_reconstruction_error_central_and_LOO':err,'new_random_samples':0}
    (ROOT/'OLD_ARCHIVE_PLANNING.json').write_text(json.dumps(result,indent=2)+'\n')
    for n in NS:
        keys=[f'N{n}.{s}.v' for s in ('entry','completion','total','difference')]
        print(n,{k:est[k] for k in keys})
        ii=[labels.index(k) for k in keys[:2]];cc=cov[np.ix_(ii,ii)]
        print('channel covariance',cc.tolist(),'corr',cc[0,1]/math.sqrt(cc[0,0]*cc[1,1]))
    print('reconstruction error',err)
if __name__=='__main__':main()

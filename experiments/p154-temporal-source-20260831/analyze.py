#!/usr/bin/env python3
"""Original-U response to one explicitly lagged path source; no simulation."""
import csv,gzip,hashlib,json,math,platform,time
from fractions import Fraction
from pathlib import Path
import numpy as np
import scipy
from scipy.special import gammaln

ROOT=Path(__file__).resolve().parent
NS=(65,85,130,170,260,340)
COMPONENTS=('current','early_raw','early_centered','early_r0','early_r1','early_r2')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def moments(sums,samples,p,n):
    k=np.arange(n+1,dtype=float)
    w=np.exp(gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1)+k*math.log(p)+(n-k)*math.log1p(-p))
    score=k/p-(n-k)/(1-p);wp=w*score
    wpp=w*(score**2-k/p**2-(n-k)/(1-p)**2)
    z,zp,zpp=w.sum(),wp.sum(),wpp.sum()
    m=w@sums/(samples*z);mp=(wp@sums/samples-m*zp)/z
    mpp=(wpp@sums/samples-m*zpp-2*mp*zp)/z
    return m,mp,mpp
def response(profile,samples,n,p0):
    m,p,pp=(np.array(v) for v in zip(*(moments(profile[g],samples,p0,n) for g in range(2))))
    q,e,s,qs,es=m.T;qp,ep,sp,qsp,esp=p.T
    jq=qs-q*s;je=es-e*s;jqp=qsp-qp*s-q*sp;jep=esp-ep*s-e*sp
    delta=float(Fraction('1152/845' if n in (65,130,260) else '2304/1445'))
    d=qp.mean();b=(ep[0]-ep[1])/delta;h=(pp[0,1]-pp[1,1])/delta;t=pp[:,0].mean()
    root=-jq.mean()/d;scale=n**(13/8)/2
    pieces=(scale*(jep[0]-jep[1])/delta/d,scale*h*root/d,-scale*b*jqp.mean()/d**2,-scale*b*t*root/d**2)
    out={'U':scale*b/d,'v':math.fsum(pieces),'rootdot':root,'rank1_rootdot':-je.mean()-ep.mean()*root,
         'v_direct':pieces[0],'v_rootmotion':pieces[1],'v_slope_source':pieces[2],'v_slope_root':pieces[3]}
    for g,name in enumerate(('first','second')):
        out.update({f'{name}.Jq':jq[g],f'{name}.JE':je[g],f'{name}.entry_fixedp':(jq[g]-je[g])/2,
                    f'{name}.exit_fixedp':(jq[g]+je[g])/2,
                    f'{name}.rank1_fixedp':-je[g]})
    return {k:float(v) for k,v in out.items()}
def load(n):
    shape=(100,2,n+1,3,10);v=np.zeros(shape,dtype=np.float64);seen=set();lag=math.isqrt(n)+(math.isqrt(n)**2<n)
    keys=('count','sum_q','sum_e','sum_s_early','sum_qs_early','sum_es_early','sum_s2_early','sum_s_now','sum_qs_now','sum_es_now')
    path=ROOT/'results/raw'/f'n{n}.csv.gz'
    with gzip.open(path,'rt') as f:
        for row in csv.DictReader(f):
            b,g,k,r=int(row['batch']),('first','second').index(row['orientation']),int(row['k']),int(row['early_rank'])
            idx=(b,g,k,r)
            assert idx not in seen and int(row['n'])==n and int(row['lag'])==lag and int(row['source_k'])==max(0,k-lag)
            seen.add(idx);v[idx]=[int(row[x]) for x in keys]
    assert len(seen)==100*2*(n+1)*3
    assert np.all(v[...,0].sum(axis=-1)==(10000 if n in (260,340) else 1000))
    return v
def point(total,samples,n,p0):
    pooled=total.sum(axis=2);base=pooled[...,[1,2,7,8,9]]
    profiles={'current':base,'early_raw':pooled[...,[1,2,3,4,5]]}
    count=total[...,0];mu=np.divide(total[...,3],count,out=np.zeros_like(count),where=count>0)
    qs=total[...,4]-mu*total[...,1];es=total[...,5]-mu*total[...,2]
    # Mean source is zero within every early stratum; no cross-rank mean imputation.
    for r in range(3):
        z=base.copy();z[...,2]=0;z[...,3]=qs[...,r];z[...,4]=es[...,r]
        profiles[f'early_r{r}']=z
    z=base.copy();z[...,2]=0;z[...,3]=qs.sum(axis=-1);z[...,4]=es.sum(axis=-1)
    profiles['early_centered']=z
    values={'p0':float(p0)}
    for name in COMPONENTS:
        values.update({f'{name}.{k}':v for k,v in response(profiles[name],samples,n,p0).items()})
    # Exact constraints: r=2 is absorbing, and entry is already1 on r=1.
    check={'immediate_rank_source_sum':float(np.max(np.abs(total[...,3]-mu*count))),
           'absorbed_late_qs':float(np.max(np.abs(qs[...,2]))),
           'absorbed_late_es':float(np.max(np.abs(es[...,2]))),
           'r1_entry_source_sum':float(np.max(np.abs(qs[...,1]-es[...,1]))),
           'stratum_addback':max(abs(values['early_centered.'+f]-sum(values[f'early_r{r}.'+f] for r in range(3))) for f in ('v','rootdot','rank1_rootdot'))}
    # All sources share the same original unperturbed U; do not add it by strata.
    return values,check
def vector(points):
    vals={f'N{n}.{k}':v for n in NS for k,v in points[n].items()}
    for lineage in ((65,130,260),(85,170,340)):
        for model,c in (('q2',(1,-3,2)),('Jordan',(1,-2,1))):
            for source in ('early_raw','early_centered','early_r0','early_r1'):
                vals[f'{model}.{lineage[0]}.{source}.v']=math.fsum(x*points[n][f'{source}.v'] for n,x in zip(lineage,c))
    return vals
def main():
    t=time.perf_counter();dest=ROOT/'results/latest.json'
    if dest.exists():raise RuntimeError('Analysis output exists; refusing overwrite')
    anchor=json.loads((ROOT/'inputs/anchors.json').read_text())
    raw={n:load(n) for n in NS};total={n:v.sum(axis=0) for n,v in raw.items()}
    samples={n:1000000 if n in (260,340) else 100000 for n in NS}
    points={};identities={};saved_error={k:0. for k in ('U','v','rootdot','rank1_rootdot')}
    oldkeys={'U':'U','v':'Udot_fugacity','rootdot':'rootdot_fugacity','rank1_rootdot':'root_comoving_rank1_fugacity'}
    def evaluate(n,sums,size,p0,saved):
        value,checks=point(sums,size,n,p0)
        for k,x in checks.items():identities[k]=max(identities.get(k,0.),x)
        for k,old in oldkeys.items():saved_error[k]=max(saved_error[k],abs(value['current.'+k]-saved[old]))
        return value
    for n in NS:points[n]=evaluate(n,total[n],samples[n],anchor['by_N'][str(n)]['p0'],anchor['by_N'][str(n)])
    central=vector(points);labels=list(central);mean=np.array(list(central.values()));cov=np.zeros((len(labels),len(labels)))
    groups={}
    for name,group in anchor['groups'].items():
        assert group['delete_one_batch_ids']==list(range(100));vectors=[]
        for b in range(100):
            changed=dict(points)
            for n in group['Ns']:
                saved={k:v[b] for k,v in group['by_N'][str(n)].items()}
                changed[n]=evaluate(n,total[n]-raw[n][b],samples[n]*.99,saved['p0'],saved)
            vectors.append(list(vector(changed).values()))
        vectors=np.asarray(vectors);deviation=vectors-vectors.mean(axis=0);cov+=.99*deviation.T@deviation
        groups[name]={'Ns':group['Ns'],'delete_one_batch_ids':list(range(100)),'delete_one_vectors':vectors.tolist()}
    errors=np.sqrt(np.maximum(0,cov.diagonal()))
    estimates={key:{'value':float(x),'se':float(e),'z':float(x/e) if e>1e-13 else None} for key,x,e in zip(labels,mean,errors)}
    assert max(saved_error.values())<1e-6,saved_error
    assert max(identities.values())<1e-6,identities
    assert np.isfinite(mean).all() and np.isfinite(cov).all()
    result={'schema':'p154-temporal-source-result.v1','status':'completed','contract':json.loads((ROOT/'CONTRACT.json').read_text()),
            'labels':labels,'estimates':estimates,'covariance':cov.tolist(),'covariance_groups':groups,
            'identities':identities,'same_time_control_reconstruction_max_error':saved_error,
            'samples_per_N':samples,'elapsed_seconds':time.perf_counter()-t,
            'environment':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'hostname':platform.node()},
            'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in sorted((ROOT/'results/raw').glob('*.gz'))},
            'analysis_sha256':sha(Path(__file__))}
    dest.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'analysis_seconds':result['elapsed_seconds'],'dimensions':len(labels),'identities':identities,'control_errors':saved_error}),flush=True)
    for n in NS:print(n,{s:{f:estimates[f'N{n}.{s}.{f}'] for f in ('v','rootdot','rank1_rootdot')} for s in ('early_raw','early_centered','early_r0','early_r1')},flush=True)
if __name__=='__main__':main()

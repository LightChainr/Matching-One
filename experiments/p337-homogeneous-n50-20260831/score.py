#!/usr/bin/env python3
"""Fixed N50 homogeneous original-U and Sstar tangent, exact intervals."""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys

PACKAGE=Path(__file__).resolve().parent
BACKEND=PACKAGE.parent/'p337-finite-law-window-20260831/vendor/interval_backend.py'
BACKEND_SHA='001a4ec8d85934c11690c1948f47ea8bdc892ad854ee422882cbaa4053fd09db'
assert hashlib.sha256(BACKEND.read_bytes()).hexdigest()==BACKEND_SHA
spec=importlib.util.spec_from_file_location('fixed_interval',BACKEND)
backend=importlib.util.module_from_spec(spec);sys.modules[spec.name]=backend;spec.loader.exec_module(backend)
I=backend.Interval
N=50;DELTA=F(-1152,625)


def poly(c,x):
    value=0
    for a in reversed(c):value=value*x+a
    return value


def read(path,geometry):
    data=json.loads(path.read_text())
    if not data['complete'] or data['N']!=N or data['geometry']!=geometry:
        raise ValueError('complete exact fixed geometry required')
    arrays={s:[0]*(N+1) for s in ('one','q','e','s','qs','es')}
    seen=set()
    for row in data['histogram']:
        k,q,c,s=(row[x] for x in ('K','q','count','sum_S'))
        if any(type(v)!=int for v in (k,q,c,s)) or not(0<=k<=N and q in(-1,0,1) and c>0):
            raise ValueError('integer histogram required')
        if (k,q) in seen:raise ValueError('duplicate row')
        seen.add((k,q))
        for field,value in zip(arrays,(c,q*c,q*q*c,s,q*s,q*q*s)):
            arrays[field][k]+=value
    if arrays['one']!=[math.comb(N,k) for k in range(N+1)]:
        raise ValueError('incomplete binomial population')
    return arrays


def packet(arrays,h):
    norm=poly(arrays['one'],h)
    raw={name:[poly([k**j*c for k,c in enumerate(values)],h)/norm
               for j in range(3)] for name,values in arrays.items()}
    mu,mu2=raw['one'][1:];s,ks=raw['s'][:2]
    out={}
    for name in ('q','e'):
        o,ok,ok2=raw[name];os,oks=raw[name+'s'][:2]
        out[name]=dict(value=o,z=ok-o*mu,zz=ok2-2*mu*ok+(2*mu*mu-mu2)*o,
                       t=os-o*s,zt=oks-ok*s-os*mu-o*ks+2*o*mu*s)
    out['mean_s']=s
    return out


def evaluate(pair):
    coefficients=[a+b for a,b in zip(pair[0]['q'],pair[1]['q'])]
    signs=[1 if x>0 else -1 for x in coefficients if x]
    if sum(x!=y for x,y in zip(signs,signs[1:]))!=1:raise ValueError('root uniqueness gate')
    lo,hi=F(0),F(4)
    if not poly(coefficients,lo)<0<poly(coefficients,hi):raise ValueError('root bracket gate')
    for _ in range(160):
        mid=(lo+hi)/2;val=poly(coefficients,mid)
        if val==0:lo=hi=mid;break
        if val<0:lo=mid
        else:hi=mid
    h=I(lo,hi);ps=[packet(data,h) for data in pair]
    mean=lambda field:(ps[0]['q'][field]+ps[1]['q'][field])/2
    diff=lambda field:(ps[0]['e'][field]-ps[1]['e'][field])/DELTA
    d=mean('z')
    if d.lo<=0:raise ValueError('positive slope gate')
    b=diff('z');zt=-mean('t')/d
    terms=dict(direct=diff('zt')/d,root_motion=zt*diff('zz')/d,
               slope_source=-b*mean('zt')/(d*d),slope_root=-b*zt*mean('zz')/(d*d))
    v=sum(terms.values());u=b/d
    encode=backend.interval_json
    return dict(h_root=encode(h),p_root=encode(h/(1+h)),D_z=encode(d),
                reduced_U=encode(u),reduced_V=encode(v),root_z_t=encode(zt),
                reduced_terms={k:encode(x) for k,x in terms.items()},
                U_display=float(backend.middle(u))*N**(13/8)/2,
                V_display=float(backend.middle(v))*N**(13/8)/2,
                V_positive=v.lo>0,V_negative=v.hi<0,finite_zero_excluded=v.lo>0 or v.hi<0,
                endpoint_positive_sign_rejected=v.hi<0,
                geometry_mean_S=[encode(p['mean_s']) for p in ps])


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--freeze-commit',required=True)
    ap.add_argument('--first',type=Path,required=True);ap.add_argument('--second',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    repo=Path(subprocess.check_output(['git','-C',str(PACKAGE),'rev-parse','--show-toplevel'],text=True).strip())
    freeze=subprocess.check_output(['git','-C',str(repo),'rev-parse',args.freeze_commit+'^{commit}'],text=True).strip()
    if freeze!=args.freeze_commit:raise ValueError('full freeze SHA required')
    hashes={}
    for name in ('CONTRACT.md','score.py'):
        path=PACKAGE/name
        blob=subprocess.check_output(['git','-C',str(repo),'show',freeze+':'+str(path.relative_to(repo))])
        if blob!=path.read_bytes():raise ValueError('frozen score package changed')
        hashes[name]=hashlib.sha256(blob).hexdigest()
    if args.output.exists():raise ValueError('do not overwrite prior score')
    result=evaluate([read(args.first,[5,5]),read(args.second,[1,7])])
    result.update(freeze_commit=freeze,package_hashes=hashes,backend_sha256=BACKEND_SHA,
                  input_sha256=[hashlib.sha256(p.read_bytes()).hexdigest() for p in (args.first,args.second)],
                  dependency='single exact N50 finite population, not independent random evidence')
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ('U_display','V_display','endpoint_positive_sign_rejected')},indent=2))


if __name__=='__main__':main()

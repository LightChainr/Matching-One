#!/usr/bin/env python3
"""Frozen pooled-root score for the N65 kernel-changing contact x stage MC."""
import argparse,csv,hashlib,json,math
from collections import defaultdict
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.stats import binom

N=65; B=100; DELTA=1152/845; FIELDS=("count","sum_q0","sum_E0","sum_a16_0","sum_q0_a16_0","sum_E0_a16_0","sum_q1","sum_E1","sum_a16_1","sum_q1_a16_1","sum_E1_a16_1")
FROZEN_SAMPLES=20_000_000;FROZEN_SHARDS=4;FROZEN_SEED="20260901537";FROZEN_P=0.5927311266364432

def hist(path):
    x=np.zeros((B,2,2,N+1)); samples=np.zeros(B)
    with path.open(newline="") as f:
        for r in csv.DictReader(f):
            if int(r["n"])!=N: continue
            b=int(r["batch"]); g=("first","second").index(r["orientation"]); t=("minus","plus").index(r["kind"]); x[b,g,t,int(r["k"])]+=int(r["count"]); samples[b]=int(r["samples"])
    if np.any(samples==0): raise ValueError("incomplete N65 baseline")
    return x,samples

def baseline(raw,samples,p,omit=None):
    h=raw.sum(0)-(raw[omit] if omit is not None else 0); total=samples.sum()-(samples[omit] if omit is not None else 0)
    c=np.cumsum(h,axis=-1); q=(-total+c[:,0]+c[:,1])/total; e=(total-c[:,0]+c[:,1])/total
    k=np.arange(N+1); w=binom.pmf(k,N,p); s=k-N*p
    return {"q":q@w,"e":e@w,"qt":q@(w*s),"et":e@(w*s)}

def tables(paths):
    rows=defaultdict(lambda:np.zeros(len(FIELDS),dtype=np.int64)); metas=[]
    for path in paths:
        meta={}; lines=[]
        for line in path.read_text().splitlines():
            if line.startswith("# "): k,v=line[2:].split("=",1);meta[k]=v
            elif not line.startswith("#"): lines.append(line)
        metas.append(meta)
        for r in csv.DictReader(lines,delimiter="\t"):
            key=(int(r["batch"]),r["kind"],("axis","tilted").index(r["geometry"]),int(r["dx"]),int(r["dy"]),r["stage"],int(r["contact_mask"]),int(r["k"]))
            rows[key]+=np.array([int(r[f]) for f in FIELDS],dtype=np.int64)
    shards=int(metas[0]["shard_count"]); pstar=float(metas[0]["proposal_p"]); seed=metas[0]["seed"]; samples=int(metas[0]["samples"])
    if {int(m["shard_index"]) for m in metas}!=set(range(shards)) or any((int(m["samples"]),float(m["proposal_p"]),m["seed"])!=(samples,pstar,seed) for m in metas): raise ValueError("mixed/incomplete shards")
    for m in metas:
        i=int(m["shard_index"])
        if int(m.get("begin",-1))!=samples*i//shards or int(m.get("end",-1))!=samples*(i+1)//shards: raise ValueError("shard range does not partition the frozen counter domain")
    smoke=samples<=10_000
    if not smoke and (samples!=FROZEN_SAMPLES or shards!=FROZEN_SHARDS or seed!=FROZEN_SEED or pstar!=FROZEN_P): raise ValueError("non-smoke input does not match the frozen production contract")
    return rows,pstar,{"samples":samples,"shards":shards,"seed":seed,"proposal_p":pstar}

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--baseline",required=True,type=Path);ap.add_argument("--tables",nargs="+",required=True,type=Path);ap.add_argument("--output",required=True,type=Path);ap.add_argument("--audit-output",type=Path);a=ap.parse_args()
    if a.output.exists(): raise FileExistsError(a.output)
    if a.audit_output is not None and a.audit_output.exists(): raise FileExistsError(a.audit_output)
    H,HS=hist(a.baseline); rows,pstar,run=tables(a.tables); records=[(k,v) for k,v in rows.items()]
    disps=sorted({k[3:5] for k,v in records if k[1]=="global"}); probe=disps[0]
    batch_n=np.zeros(B)
    for (b,kind,g,dx,dy,stage,mask,k),v in records:
        if kind=="global" and g==0 and (dx,dy)==probe: batch_n[b]+=v[0]
    if np.any(batch_n==0): raise ValueError("every producer batch must be populated")

    def evaluate(base_omit=None,prod_omit=None,detail=False):
        root=brentq(lambda p:baseline(H,HS,p,base_omit)["q"].mean(),.58,.61); base=baseline(H,HS,root,base_omit)
        mt=base["qt"].mean(); yt=(base["et"][0]-base["et"][1])/DELTA; R=yt/mt; c=(1/DELTA,-1/DELTA); muh=[2*c[g]*base["e"][g]-R*base["q"][g] for g in range(2)]
        denom=batch_n.sum()-(batch_n[prod_omit] if prod_omit is not None else 0); packet=defaultdict(lambda:[0.,0.])
        def factor(k): return (1-root)*(root/pstar)**k*((1-root)/(1-pstar))**(N-2-k)/denom
        for (b,kind,g,dx,dy,stage,mask,k),v in records:
            if kind!="global" or b==prod_omit: continue
            f=factor(k); a0,a1=v[3]/(N*16),v[8]/(N*16);qa0,qa1=v[4]/(N*16),v[9]/(N*16)
            packet[(g,dx,dy)][0]+=f*((1-root)*a0+root*a1);packet[(g,dx,dy)][1]+=f*((1-root)*qa0+root*qa1)
        beta={}; pdata={}
        for d in disps:
            cov=[]
            for g in range(2):
                mu,mqa=packet[(g,*d)];cov.append(mqa-base["q"][g]*mu);pdata[(g,*d)]=(mu,cov[-1])
            beta[d]=sum(cov)/(2*mt)
        M=np.zeros((2,3));exposure=np.zeros((2,3))
        for (b,kind,g,dx,dy,stage,mask,k),v in records:
            if kind!="carrier" or b==prod_omit: continue
            f=factor(k);mu=pdata[(g,dx,dy)][0];bet=beta[(dx,dy)];st=("01","12").index(stage);sm=k-(N-1)*root
            for i in (0,1):
                q,e,a16,qa16,ea16=(v[j] for j in ((1,2,3,4,5) if i==0 else (6,7,8,9,10)));aa,qa,ea=a16/(N*16),qa16/(N*16),ea16/(N*16);w=(1-root,root)[i];u=i-root;bi=u*(sm+u)-root*(1-root)
                sh=2*c[g]*e-R*q-muh[g]*v[0];sha=2*c[g]*(ea-mu*e)-R*(qa-mu*q)-muh[g]*(aa-mu*v[0]);M[st,mask-1]+=2*f*w*(u*sha-bet*bi*sh);exposure[st,mask-1]+=2*f*w*v[0]/N
        extra={"p":root,"M_t":mt,"R":R}
        if detail: extra["per_displacement"]={f"{d[0]},{d[1]}":{"beta":beta[d],"axis_mu_a":pdata[(0,*d)][0],"tilted_mu_a":pdata[(1,*d)][0]} for d in disps}
        # Keep root and R at their historical indices so the frozen default
        # result remains byte-identical.  Exposure is an audit-only extension.
        return np.r_[M.ravel(),root,R,exposure.ravel()],extra

    central,detail=evaluate(detail=True);prod=np.array([evaluate(prod_omit=b)[0] for b in range(B)]);base=np.array([evaluate(base_omit=b)[0] for b in range(B)])
    fp=math.sqrt((B-1)/B)*(prod-prod.mean(0));fb=math.sqrt((B-1)/B)*(base-base.mean(0));se=np.sqrt((fp*fp).sum(0)+(fb*fb).sum(0))
    def primary(v):
        m=v[:6].reshape(2,3);c=np.c_[m[:,0]+m[:,1],m[:,2]];left,right=c[0,0]*c[1,1],c[0,1]*c[1,0];delta=left-right;theta=delta/(abs(left)+abs(right)) if left or right else 0.;return np.r_[c.ravel(),delta,theta]
    pc=primary(central);pp=np.array([primary(v) for v in prod]);pb=np.array([primary(v) for v in base]);fpp=math.sqrt((B-1)/B)*(pp-pp.mean(0));fpb=math.sqrt((B-1)/B)*(pb-pb.mean(0));pse=np.sqrt((fpp*fpp).sum(0)+(fpb*fpb).sum(0));ci=[pc[4]-1.96*pse[4],pc[4]+1.96*pse[4]]
    smoke=run["samples"]<=10000; signs=np.sign(pc[:4]).astype(int).tolist();decision="SMOKE_ONLY" if smoke else ("CONTACT_FUSION_COMPLETION_TRANSMITS" if ci[1]<0 and signs==[-1,-1,-1,1] else ("SIGN_ROTATION_REJECTED" if ci[0]>0 or signs!=[-1,-1,-1,1] else "UNRESOLVED_CONTACT_STAGE_GATE"))
    payload={"schema":"matching-one/p537-contact-stage-n65-score/v1","status":"SMOKE" if smoke else "COMPLETED","geometry_order":["axis(8,1)","tilted(7,4)"],"stage_order":["01","12"],"contact_mask_order":[1,2,3],"definition":"alternating rank-changing Bell edge with g16_before != g16_after; full pooled-root Schur allocation; beta fixed per common displacement before aggregation","global":detail,"matrix":{"estimate":central[:6].reshape(2,3).tolist(),"se":se[:6].reshape(2,3).tolist()},"primary":{"column_order":["single=mask1+mask2","double=mask3"],"matrix":pc[:4].reshape(2,2).tolist(),"matrix_se":pse[:4].reshape(2,2).tolist(),"Delta":pc[4],"Delta_se":pse[4],"Delta_95":ci,"theta":pc[5],"theta_se":pse[5],"theta_definition":"Delta/(abs(L01_single*L12_double)+abs(L01_double*L12_single))","decision":decision},"independent_covariance_groups":["P45_baseline_100_batches","new_MC_100_batches"],"run":run}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(payload,indent=2)+"\n")
    if a.audit_output is not None:
        # This completes the contract's retained covariance/exposure fields.
        # It does not alter the preregistered primary score above.
        fp6=fp[:,:6];fb6=fb[:,:6];covp=fp6.T@fp6;covb=fb6.T@fb6;cov=covp+covb
        T=np.array([[1,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,1,0],[0,0,0,0,0,1]],dtype=float)
        ccovp=T@covp@T.T;ccovb=T@covb@T.T;ccov=ccovp+ccovb
        exp=central[8:14];expp=prod[:,8:14];expb=base[:,8:14]
        ratio=central[:6]/exp;ratiop=prod[:,:6]/expp;ratiob=base[:,:6]/expb
        frp=math.sqrt((B-1)/B)*(ratiop-ratiop.mean(0));frb=math.sqrt((B-1)/B)*(ratiob-ratiob.mean(0));ratiose=np.sqrt((frp*frp).sum(0)+(frb*frb).sum(0))
        one=np.ones(6);total=float(one@central[:6]);vp=float(one@covp@one);vb=float(one@covb@one);vt=vp+vb;tse=math.sqrt(vt)
        z=central[:6]/np.sqrt(np.diag(cov));cellci=np.c_[central[:6]-1.96*np.sqrt(np.diag(cov)),central[:6]+1.96*np.sqrt(np.diag(cov))]
        sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
        audit={"schema":"matching-one/p537-contact-stage-n65-audit/v1","status":"COMPLETED_SAME_FROZEN_BLOCK","scope":{"confirmatory":"the preregistered four signs and Delta decision only","contract_completion":["positive exposure","full two-group delete-one covariance","selected-carrier total"],"exploratory":["conditional signed density","marginal cell z scores"],"nonclaims":["not six independent evidence votes","not a coordinate-invariant physical operator","not the complete full-T response","not a size exponent or CFT identification"]},"vector_order":["01_mask1","01_mask2","01_mask3","12_mask1","12_mask2","12_mask3"],"matrix":{"estimate":central[:6].reshape(2,3).tolist(),"se":np.sqrt(np.diag(cov)).reshape(2,3).tolist(),"cell_z_marginal":z.reshape(2,3).tolist(),"cell_ci95_marginal":cellci.reshape(2,3,2).tolist(),"covariance":{"production":covp.tolist(),"baseline":covb.tolist(),"total":cov.tolist()}},"positive_exposure":{"definition":"source-normalized importance-weighted selected-carrier background mass; not an event probability","estimate":exp.reshape(2,3).tolist(),"se":se[8:14].reshape(2,3).tolist()},"conditional_signed_density":{"status":"EXPLORATORY_NOT_PREREGISTERED","estimate":ratio.reshape(2,3).tolist(),"se_nonlinear_delete_one":ratiose.reshape(2,3).tolist()},"primary":{"collapse_transform":T.tolist(),"column_order":["01_single","01_double","12_single","12_double"],"estimate":pc[:4].tolist(),"covariance":{"production":ccovp.tolist(),"baseline":ccovb.tolist(),"total":ccov.tolist()},"Delta":pc[4],"Delta_se":pse[4],"Delta_95":ci,"decision":decision,"theta":{"value":pc[5],"se":pse[5],"saturated_by_sign_pattern":bool(signs==[-1,-1,-1,1]),"production_loo_minmax":[float(pp[:,5].min()),float(pp[:,5].max())],"baseline_loo_minmax":[float(pb[:,5].min()),float(pb[:,5].max())],"interpretation":"theta=-1 is algebraically forced throughout this open sign cone; zero SE is not infinite precision"}},"selected_carrier_total":{"estimate":total,"variance_production":vp,"variance_baseline":vb,"variance_total":vt,"se":tse,"ci95":[total-1.96*tse,total+1.96*tse],"z":total/tse},"jackknife":{"production_batches":B,"baseline_batches":B,"groups_independent":True,"factor":"99/100"},"inputs":{"baseline":{"path":str(a.baseline),"sha256":sha(a.baseline)},"tables":[{"path":str(p),"sha256":sha(p)} for p in a.tables]}}
        a.audit_output.parent.mkdir(parents=True,exist_ok=True);a.audit_output.write_text(json.dumps(audit,indent=2)+"\n")
    print(json.dumps({"status":payload["status"],"p":detail["p"],"matrix":payload["matrix"]}))
if __name__=="__main__": main()

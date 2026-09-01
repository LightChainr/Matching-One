#!/usr/bin/env python3
"""Frozen pooled-root score for the N65 kernel-changing contact x stage MC."""
import argparse,csv,json,math
from collections import defaultdict
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.stats import binom

N=65; B=100; DELTA=1152/845; FIELDS=("count","sum_q0","sum_E0","sum_a16_0","sum_q0_a16_0","sum_E0_a16_0","sum_q1","sum_E1","sum_a16_1","sum_q1_a16_1","sum_E1_a16_1")

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
    return rows,pstar,{"samples":samples,"shards":shards,"seed":seed,"proposal_p":pstar}

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--baseline",required=True,type=Path);ap.add_argument("--tables",nargs="+",required=True,type=Path);ap.add_argument("--output",required=True,type=Path);a=ap.parse_args()
    if a.output.exists(): raise FileExistsError(a.output)
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
        M=np.zeros((2,3))
        for (b,kind,g,dx,dy,stage,mask,k),v in records:
            if kind!="carrier" or b==prod_omit: continue
            f=factor(k);mu=pdata[(g,dx,dy)][0];bet=beta[(dx,dy)];st=("01","12").index(stage);sm=k-(N-1)*root
            for i in (0,1):
                q,e,a16,qa16,ea16=(v[j] for j in ((1,2,3,4,5) if i==0 else (6,7,8,9,10)));aa,qa,ea=a16/(N*16),qa16/(N*16),ea16/(N*16);w=(1-root,root)[i];u=i-root;bi=u*(sm+u)-root*(1-root)
                sh=2*c[g]*e-R*q-muh[g]*v[0];sha=2*c[g]*(ea-mu*e)-R*(qa-mu*q)-muh[g]*(aa-mu*v[0]);M[st,mask-1]+=2*f*w*(u*sha-bet*bi*sh)
        extra={"p":root,"M_t":mt,"R":R}
        if detail: extra["per_displacement"]={f"{d[0]},{d[1]}":{"beta":beta[d],"axis_mu_a":pdata[(0,*d)][0],"tilted_mu_a":pdata[(1,*d)][0]} for d in disps}
        return np.r_[M.ravel(),root,R],extra

    central,detail=evaluate(detail=True);prod=np.array([evaluate(prod_omit=b)[0] for b in range(B)]);base=np.array([evaluate(base_omit=b)[0] for b in range(B)])
    fp=math.sqrt((B-1)/B)*(prod-prod.mean(0));fb=math.sqrt((B-1)/B)*(base-base.mean(0));se=np.sqrt((fp*fp).sum(0)+(fb*fb).sum(0))
    def primary(v):
        m=v[:6].reshape(2,3);c=np.c_[m[:,0]+m[:,1],m[:,2]];left,right=c[0,0]*c[1,1],c[0,1]*c[1,0];delta=left-right;theta=delta/(abs(left)+abs(right)) if left or right else 0.;return np.r_[c.ravel(),delta,theta]
    pc=primary(central);pp=np.array([primary(v) for v in prod]);pb=np.array([primary(v) for v in base]);fpp=math.sqrt((B-1)/B)*(pp-pp.mean(0));fpb=math.sqrt((B-1)/B)*(pb-pb.mean(0));pse=np.sqrt((fpp*fpp).sum(0)+(fpb*fpb).sum(0));ci=[pc[4]-1.96*pse[4],pc[4]+1.96*pse[4]]
    cell95=np.c_[pc[:4]-1.96*pse[:4],pc[:4]+1.96*pse[:4]]
    expected=np.array([-1,-1,-1,1])
    observed=np.sign(pc[:4]).astype(int)
    opposite_excluded=np.array([
        lo>0 if sign<0 else hi<0
        for sign,(lo,hi) in zip(expected,cell95)
    ])
    smoke=run["samples"]<=10000
    if smoke:
        decision="SMOKE_ONLY"
    elif ci[1]<0 and np.array_equal(observed,expected):
        decision="CONTACT_FUSION_COMPLETION_TRANSMITS"
    elif ci[0]>0 or np.any(opposite_excluded):
        decision="SIGN_ROTATION_REJECTED"
    else:
        decision="UNRESOLVED_CONTACT_STAGE_GATE"
    matrix_cov_prod=fp[:,:6].T@fp[:,:6]
    matrix_cov_base=fb[:,:6].T@fb[:,:6]
    primary_cov_prod=fpp.T@fpp
    primary_cov_base=fpb.T@fpb
    payload={"schema":"matching-one/p537-contact-stage-n65-score/v1","status":"SMOKE" if smoke else "COMPLETED","geometry_order":["axis(8,1)","tilted(7,4)"],"stage_order":["01","12"],"contact_mask_order":[1,2,3],"definition":"alternating rank-changing Bell edge with g16_before != g16_after; full pooled-root Schur allocation; beta fixed per common displacement before aggregation","global":detail,"matrix":{"estimate":central[:6].reshape(2,3).tolist(),"se":se[:6].reshape(2,3).tolist(),"covariance_order":["01_mask1","01_mask2","01_mask3","12_mask1","12_mask2","12_mask3"],"covariance":{"new_MC_100_batches":matrix_cov_prod.tolist(),"P45_baseline_100_batches":matrix_cov_base.tolist(),"combined":(matrix_cov_prod+matrix_cov_base).tolist()}},"primary":{"column_order":["single=mask1+mask2","double=mask3"],"matrix":pc[:4].reshape(2,2).tolist(),"matrix_se":pse[:4].reshape(2,2).tolist(),"matrix_95":cell95.reshape(2,2,2).tolist(),"expected_cell_signs":[[-1,-1],[-1,1]],"observed_cell_signs":observed.reshape(2,2).tolist(),"opposite_sign_excluded_95":opposite_excluded.reshape(2,2).tolist(),"Delta":pc[4],"Delta_se":pse[4],"Delta_95":ci,"theta":pc[5],"theta_se":pse[5],"theta_definition":"Delta/(abs(L01_single*L12_double)+abs(L01_double*L12_single))","covariance_order":["01_single","01_double","12_single","12_double","Delta","theta"],"covariance":{"new_MC_100_batches":primary_cov_prod.tolist(),"P45_baseline_100_batches":primary_cov_base.tolist(),"combined":(primary_cov_prod+primary_cov_base).tolist()},"decision":decision},"independent_covariance_groups":["P45_baseline_100_batches","new_MC_100_batches"],"run":run}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(payload,indent=2)+"\n");print(json.dumps({"status":payload["status"],"p":detail["p"],"matrix":payload["matrix"]}))
if __name__=="__main__": main()

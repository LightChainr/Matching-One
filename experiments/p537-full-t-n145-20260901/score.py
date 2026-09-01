#!/usr/bin/env python3
"""Score the N145 complete canonical-pair thermal response."""
from __future__ import annotations
import argparse,csv,hashlib,json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
from scipy.stats import binom

N=145; B=100; DELTA=8064/4205
FIELDS=("count","a0","qa0","ea0","a1","qa1","ea1","nn0","qnn0","enn0","nn1","qnn1","enn1")
OUT=("p","T_t","T_t_over_M_t","J","M_t","R","R_t","jM","jM_t","jY","jY_t")
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def meta(p):
 d={}
 with p.open() as f:
  for line in f:
   if not line.startswith("# "): break
   k,v=line[2:].rstrip().split("=",1);d[k]=v
 return d
def baseline(p):
 raw=np.zeros((B,2,2,N+1)); samples=np.zeros(B)
 with p.open(newline="") as f:
  for r in csv.DictReader(f):
   if int(r["n"])!=N: continue
   b=int(r["batch"]);g=("first","second").index(r["orientation"]);t=("minus","plus").index(r["kind"]);raw[b,g,t,int(r["k"])]+=int(r["count"]);samples[b]=int(r["samples"])
 if np.any(samples<=0): raise ValueError("incomplete baseline")
 return raw,samples
def source(paths):
 ms=[meta(p) for p in paths]; shards=int(ms[0]["shard_count"]); samples=int(ms[0]["samples"]); seed=ms[0]["seed"]; proposal=float(ms[0]["proposal_p"])
 if len(paths)!=shards or {int(x["shard_index"]) for x in ms}!=set(range(shards)) or any((int(x["samples"]),x["seed"],float(x["proposal_p"]))!=(samples,seed,proposal) for x in ms): raise ValueError("mixed/incomplete shards")
 data=np.zeros((B,2,N,len(FIELDS)),dtype=np.int64)
 for p in paths:
  with p.open(newline="") as f:
   rows=(x for x in f if not x.startswith("#"))
   for r in csv.DictReader(rows,delimiter="\t"):
    data[int(r["batch"]),("axis","tilted").index(r["geometry"]),int(r["k"])]+=np.array([int(r[x]) for x in FIELDS],dtype=np.int64)
 counts=data[...,0].sum(axis=-1)
 if np.any(counts<=0) or np.any(counts!=counts[:,:1]) or int(counts[:,0].sum())!=samples: raise ValueError("source counts")
 return data,{"samples":samples,"shards":shards,"seed":seed,"proposal_p":proposal}
def main():
 ap=argparse.ArgumentParser();ap.add_argument("--baseline",required=True,type=Path);ap.add_argument("--tables",required=True,nargs="+",type=Path);ap.add_argument("--output",required=True,type=Path);a=ap.parse_args()
 if a.output.exists(): raise FileExistsError(a.output)
 bb,bs=baseline(a.baseline); sb,run=source(a.tables); bt=bb.sum(0);st=sb.sum(0);ks=np.arange(N+1);rest=np.arange(N)
 def evaluate(bo=None,so=None):
  h=bt-(bb[bo] if bo is not None else 0); ns=bs.sum()-(bs[bo] if bo is not None else 0);cum=np.cumsum(h,axis=-1);q=(-ns+cum[:,0]+cum[:,1])/ns;e=(ns-cum[:,0]+cum[:,1])/ns
  def base(p):
   w=binom.pmf(ks,N,p);s=ks-N*p;b=s*s-N*p*(1-p);return q@w,e@w,q@(w*s),e@(w*s),q@(w*b),e@(w*b)
  p=brentq(lambda x:base(x)[0].mean(),.58,.61);qv,ev,qt,et,qtt,ett=base(p);Mt=qt.mean();Mtt=qtt.mean();Yt=(et[0]-et[1])/DELTA;Ytt=(ett[0]-ett[1])/DELTA;R=Yt/Mt;Rt=(Ytt-R*Mtt)/Mt
  s=st-(sb[so] if so is not None else 0);den=s[0,:,0].sum();pr=run["proposal_p"];ratio=(p/pr)**rest*((1-p)/(1-pr))**(N-2-rest);pack=np.zeros((2,6))
  for z,(cols,ncols) in enumerate((((1,2,3),(7,8,9)),((4,5,6),(10,11,12)))):
   score=rest+z-N*p;w=(1-p)*(p if z else 1-p)*ratio/den;vals=(s[...,list(cols)]+s[...,list(ncols)]/3)/(16*N);pack[:,:3]+=np.einsum("gkc,k->gc",vals,w);pack[:,3:]+=np.einsum("gkc,k->gc",vals,w*score)
  av,qav,eav,at,qat,eat=pack.T;jmg=qav-qv*av;jmtg=qat-qt*av-qv*at;jyg=eav-ev*av;jytg=eat-et*av-ev*at;jm=jmg.mean();jmt=jmtg.mean();jy=(jyg[0]-jyg[1])/DELTA;jyt=(jytg[0]-jytg[1])/DELTA;T=jyt-R*jmt-Rt*jm;J=(N**(13/8)/2)*T/Mt
  return np.array([p,T,T/Mt,J,Mt,R,Rt,jm,jmt,jy,jyt])
 center=evaluate();sj=np.array([evaluate(so=i) for i in range(B)]);bj=np.array([evaluate(bo=i) for i in range(B)]);f=(B-1)/B;se=np.sqrt(f*np.square(sj-sj.mean(0)).sum(0)+f*np.square(bj-bj.mean(0)).sum(0));values={k:{"value":float(v),"se":float(z),"ci95":[float(v-1.96*z),float(v+1.96*z)]} for k,v,z in zip(OUT,center,se)}
 payload={"schema":"matching-one/p537-full-t-n145/v1","status":"SMOKE" if run["samples"]<1000000 else "COMPLETED","N":N,"geometries":[[12,1],[9,8]],"delta_cos4":"8064/4205","definition":"complete canonical b16 sum with C4 reconstruction of omitted +e1 NN column","values":values,"jackknife":{"source_batches":B,"baseline_batches":B,"groups_independent":True},"run":run,"inputs":{"baseline":{"path":str(a.baseline),"sha256":sha(a.baseline)},"tables":[{"path":str(p),"sha256":sha(p)} for p in a.tables]}}
 a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(payload,indent=2)+"\n");print(json.dumps({"status":payload["status"],"T_t":values["T_t"],"J":values["J"]},indent=2))
if __name__=="__main__": main()

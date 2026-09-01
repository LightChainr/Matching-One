#!/usr/bin/env python3
"""Score the frozen held-out square-L6 root-conditioned G4 production."""
from __future__ import annotations
import argparse,csv,json,math
from collections import defaultdict
from pathlib import Path
import mpmath as mp


LOCAL_G6_BAND=(mp.mpf("0.0022905"),mp.mpf("0.0027995"))


def read(path):
    meta={}; lines=[]
    for line in path.read_text().splitlines():
        if line.startswith("# ") and "=" in line:
            k,v=line[2:].split("=",1);meta[k]=v
        elif not line.startswith("#"):lines.append(line)
    return meta,list(csv.DictReader(lines,delimiter="\t"))


def load(paths):
    metas=[]; g=defaultdict(lambda:[0,0,0,0,0]);a=defaultdict(lambda:[0,0,0])
    for path in paths:
        meta,rows=read(path);metas.append(meta)
        for row in rows:
            b,k=int(row["batch"]),int(row["k"])
            if row["kind"]=="global":
                z=g[(b,k)]
                for i,name in enumerate(("count","sum_q0","sum_q1","sum_source16_0","sum_source16_1")):
                    z[i]+=int(row[name])
            else:
                z=a[(b,row["transition"],k)]
                for i,name in enumerate(("signed_count","signed_source_mid16","unsigned_count")):
                    z[i]+=int(row[name])
    def unique(key,cast):
        values={cast(m[key]) for m in metas}
        if len(values)!=1:raise ValueError(f"mixed {key}")
        return values.pop()
    L=unique("L",int);N=unique("N",int);batches=unique("batches",int)
    proposal=unique("proposal_p",mp.mpf);samples=unique("samples",int)
    shards=unique("shard_count",int)
    if {int(m["shard_index"]) for m in metas}!=set(range(shards)):raise ValueError("incomplete shards")
    if sum(g[(b,k)][0] for b in range(batches) for k in range(N))!=samples:raise ValueError("sample mass mismatch")
    return L,N,batches,proposal,samples,g,a


def main():
    ap=argparse.ArgumentParser();ap.add_argument("tables",nargs="+",type=Path)
    ap.add_argument("--output",required=True,type=Path);ap.add_argument("--dps",type=int,default=60)
    ap.add_argument("--combined-output",type=Path)
    args=ap.parse_args()
    if args.output.exists():raise SystemExit(f"refusing to overwrite {args.output}")
    mp.mp.dps=args.dps
    L,N,B,pstar,samples,G,A=load(args.tables)
    if args.combined_output:
        if args.combined_output.exists():raise SystemExit(f"refusing to overwrite {args.combined_output}")
        with args.combined_output.open("w") as out:
            out.write("# schema=matching-one/p537-aggregate-wedge-mc-combined/v1\n")
            out.write(f"# L={L}\n# N={N}\n# samples={samples}\n# shard_index=0\n# shard_count=1\n# batches={B}\n# proposal_p={pstar}\n")
            out.write("batch\tkind\ttransition\tk\tcount\tsum_q0\tsum_q1\tsum_source16_0\tsum_source16_1\tsigned_count\tsigned_source_mid16\tunsigned_count\n")
            for b in range(B):
                for k in range(N):
                    g=G[(b,k)]
                    out.write(f"{b}\tglobal\t-\t{k}\t{g[0]}\t{g[1]}\t{g[2]}\t{g[3]}\t{g[4]}\t0\t0\t0\n")
                    for tr in ("01","12"):
                        a=A[(b,tr,k)]
                        out.write(f"{b}\tlanding\t{tr}\t{k}\t0\t0\t0\t0\t0\t{a[0]}\t{a[1]}\t{a[2]}\n")

    def one(exclude=None):
        kept=[b for b in range(B) if b!=exclude]
        n=mp.mpf(sum(G[(b,k)][0] for b in kept for k in range(N)))
        def sums_global(k):return [sum(G[(b,k)][j] for b in kept) for j in range(5)]
        def sums_landing(tr,k):return [sum(A[(b,tr,k)][j] for b in kept) for j in range(3)]
        def ratio(k,p):return (p/pstar)**k*((1-p)/(1-pstar))**(N-1-k)
        def matching(p):
            return mp.fsum(((1-p)*sums_global(k)[1]+p*sums_global(k)[2])*ratio(k,p) for k in range(N))/n
        root=mp.findroot(matching,(mp.mpf("0.58"),mp.mpf("0.61")))
        def mean_source(p):
            return mp.fsum(((1-p)*sums_global(k)[3]+p*sums_global(k)[4])*ratio(k,p) for k in range(N))/(n*16*N*N)
        def cell(tr,p):
            mu=mean_source(p)
            sc=[sums_landing(tr,k)[0] for k in range(N)]
            sm=[sums_landing(tr,k)[1] for k in range(N)]
            b=mp.fsum(sc[k]*ratio(k,p) for k in range(N))/n
            t=mp.fsum(sc[k]*(mp.mpf(k)+mp.mpf("0.5")-N*p)*ratio(k,p) for k in range(N))/n
            raw=mp.fsum(sm[k]*ratio(k,p) for k in range(N))/(n*32*N*N)
            return t,raw-mu*b,b
        def psi_sum(p):
            t1,a1,_=cell("01",p);t2,a2,_=cell("12",p)
            return t1*a2-t2*a1,t1+t2
        t1,a1,b1=cell("01",root);t2,a2,b2=cell("12",root)
        psi=t1*a2-t2*a1;s=t1+t2
        dp=mp.diff(lambda p:psi_sum(p)[0],root);ds=mp.diff(lambda p:psi_sum(p)[1],root)
        wr=dp*s-psi*ds;c4=2*psi/s;g4=2*root*(1-root)*wr/s**3
        norm2=t1*t1+t2*t2
        return dict(root=root,mean_source=mean_source(root),T01=t1,T12=t2,A01=a1,A12=a2,
                    B01=b1,B12=b2,Psi4=psi,S=s,C4=c4,G4=g4,scaled=L**4*g4,
                    thermal_norm2=norm2,chi_perp=psi/norm2,L4_chi_perp=L**4*psi/norm2,n=int(n))

    full=one(); loo_scores=[one(b) for b in range(B)]
    loo=[x["scaled"] for x in loo_scores]
    center=mp.fsum(loo)/B
    se=mp.sqrt(mp.mpf(B-1)/B*mp.fsum((x-center)**2 for x in loo))
    loo_chi=[x["L4_chi_perp"] for x in loo_scores]
    center_chi=mp.fsum(loo_chi)/B
    se_chi=mp.sqrt(mp.mpf(B-1)/B*mp.fsum((x-center_chi)**2 for x in loo_chi))
    lo=full["scaled"]-mp.mpf("1.96")*se;hi=full["scaled"]+mp.mpf("1.96")*se
    g_lo=full["G4"]-mp.mpf("1.96")*(se/L**4)
    g_hi=full["G4"]+mp.mpf("1.96")*(se/L**4)
    if g_lo>=LOCAL_G6_BAND[0] and g_hi<=LOCAL_G6_BAND[1]:decision="LOCAL_L_MINUS_4_CONTINUATION"
    elif g_hi<LOCAL_G6_BAND[0] or g_lo>LOCAL_G6_BAND[1]:decision="LOCAL_L_MINUS_4_REJECTED"
    else:decision="UNRESOLVED_OVERLAP"
    g4_L4=mp.mpf("0.014008416865993059839566875264032706578302403993935")
    g4_L5=mp.mpf("0.0052774014587435467246546864566856462110310512070205")
    se_g=se/L**4
    predictions={str(power):g4_L5*(mp.mpf(5)/6)**mp.mpf(power) for power in (4,mp.mpf("4.5"),5)}
    z_power={power:(full["G4"]-value)/se_g for power,value in predictions.items()}
    effective_56=-mp.log(full["G4"]/g4_L5)/mp.log(mp.mpf(6)/5)
    effective_46=-mp.log(full["G4"]/g4_L4)/mp.log(mp.mpf(6)/4)
    s=lambda x:mp.nstr(x,40)
    payload={"schema":"matching-one/p537-aggregate-wedge-l6-score/v1","L":L,"N":N,
      "samples":samples,"batches":B,"proposal_p":s(pstar),
      "primary":"L^4*G4","decision":decision,
      "frozen_prediction":{"G4_point":"0.0025450","G4_acceptance_band":[s(x) for x in LOCAL_G6_BAND],"basis":"(5/6)^4 continuation with fixed 10 percent tolerance"},
      "score":{k:(v if isinstance(v,int) else s(v)) for k,v in full.items()},
      "jackknife":{"L4G4_estimate":s(full["scaled"]),"L4G4_se":s(se),"L4G4_normal_95":[s(lo),s(hi)],
                    "G4_normal_95":[s(g_lo),s(g_hi)],"L4_chi_perp_se":s(se_chi)},
      "exploratory_fixed_power_comparison":{
        "predicted_G4_L6_from_L5":{power:s(value) for power,value in predictions.items()},
        "z_from_prediction":{power:s(value) for power,value in z_power.items()},
        "effective_power_L5_to_L6":s(effective_56),"effective_power_L4_to_L6":s(effective_46),
        "L_power_4p5_G4":s(L**mp.mpf("4.5")*full["G4"])},
      "boundary":"single held-out L6 square-torus Bernoulli production; no radius/source/minor scan"}
    args.output.write_text(json.dumps(payload,indent=2)+"\n");print(json.dumps(payload,indent=2))


if __name__=="__main__":main()

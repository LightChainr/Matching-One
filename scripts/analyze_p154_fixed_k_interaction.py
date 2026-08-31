#!/usr/bin/env python3
"""New sector/interaction response on deterministic reobservations of old blocks."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import platform
import subprocess
from pathlib import Path

FIELDS = ["A", "E", "C", "W", "Q", "R", "H", "J_Q", "J_R", "J_H", "J_R0", "J_R2"]
ROOT = Path(__file__).resolve().parents[1]
P = 0.59274605079
SOURCE = ROOT / "results/p154-phase-e-mixed-plane-pilot/raw"

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def read(path):
    with path.open(newline="") as f:
        rows = list(csv.DictReader(f))
    return {(int(r["batch"]), r["orientation"]): {k: (v if k == "orientation" else int(v)) for k,v in r.items()} for r in rows}

def cos4(a,b):
    return (a**4-6*a*a*b*b+b**4)/(a*a+b*b)**2

def total(rows, direction, omit=None):
    selected = [r for (b,d),r in rows.items() if d == direction and b != omit]
    return {k:sum(r[k] for r in selected) for k in selected[0] if k.startswith("sum_") or k == "samples"}

def point(s,n):
    count=s["samples"]
    m=lambda name:s[name]/count
    z0,z2=m("sum_i0"),m("sum_i2")
    k,kk,t=m("sum_k"),m("sum_kk"),m("sum_edges")
    q=(t-4*P*k+2*n*P*P)/(2*n)
    r=(t-2*kk/(n-1))/(2*n)
    h=q-r
    eq=((m("sum_i0edges")+m("sum_i2edges"))-4*P*(m("sum_i0k")+m("sum_i2k"))+2*n*P*P*(z0+z2))/(2*n)
    r0=(m("sum_i0edges")-2*m("sum_i0kk")/(n-1))/(2*n)
    r2=(m("sum_i2edges")-2*m("sum_i2kk")/(n-1))/(2*n)
    unbiased=count/(count-1)
    j0=unbiased*(r0-z0*r)
    j2=unbiased*(r2-z2*r)
    jq=unbiased*(eq-(z0+z2)*q)
    jr=j0+j2
    values=[z2-z0,z2+z0,(m("sum_k1")+m("sum_k2"))/(2*(n+1))-.5,
            (m("sum_k2")-m("sum_k1"))/(n+1),q,r,h,jq,jr,jq-jr,j0,j2]
    return values

def project(rows,n,omit=None):
    a=rows[(0,"first")];b=rows[(0,"second")]
    delta=cos4(a["a"],a["b"])-cos4(b["a"],b["b"])
    x=point(total(rows,"first",omit),n);y=point(total(rows,"second",omit),n)
    return [(u-v)/delta for u,v in zip(x,y)]

def analyze(raw):
    output={}; inputs=[]; verified=0
    for n in (65,130):
        path=raw/f"n{n}.csv"; parent=SOURCE/f"n{n}_mixed.batches.csv"
        rows,old=read(path),read(parent)
        expected={(b,d) for b in range(100) for d in ("first","second")}
        if set(rows)!=expected or set(old)!=expected: raise ValueError("incomplete aligned block")
        identity_fields=["n","a","b","batch","samples","sum_k1","sum_k2","sum_i0","sum_i1","sum_i2"]
        for key in expected:
            if any(rows[key][k]!=old[key][k] for k in identity_fields):
                raise ValueError(f"{n}/{key}: replay is not the archived sample block")
            verified+=1
        for b in range(100):
            for k in ("sum_k","sum_kk"):
                if rows[(b,"first")][k]!=rows[(b,"second")][k]: raise ValueError("unpaired occupation counts")
        central=project(rows,n)
        loo=[project(rows,n,b) for b in range(100)]
        avg=[sum(x[j] for x in loo)/100 for j in range(len(FIELDS))]
        cov=[[99/100*sum((x[i]-avg[i])*(x[j]-avg[j]) for x in loo) for j in range(len(FIELDS))] for i in range(len(FIELDS))]
        est={f:{"value":central[i],"se":math.sqrt(max(0,cov[i][i])),
                "z":central[i]/math.sqrt(cov[i][i]) if cov[i][i]>0 else None} for i,f in enumerate(FIELDS)}
        output[str(n)]={"estimates":est,"covariance":cov,"delete_one_vectors":loo,
                        "orientation_values":{d:dict(zip(FIELDS,point(total(rows,d),n))) for d in ("first","second")}}
        inputs.extend({"path":str(p.relative_to(ROOT)),"sha256":sha(p)} for p in (path,parent,SOURCE/f"n{n}_mixed.metadata.json"))
    chi=sum(output[str(n)]["estimates"]["J_R"]["z"]**2 for n in (65,130))
    return {"schema":"matching-one.p154-fixed-k-interaction-result.v1",
            "analysis_type":"retrospective deterministic replay, new mixed observable, zero new samples",
            "fields":FIELDS,"by_N":output,"inputs":inputs,"verified_archived_batch_orientation_rows":verified,
            "J_R_joint_zero":{"chi_square":chi,"df":2,"nominal_p_value":math.exp(-chi/2),"basis":"independent N blocks, estimated delete-one covariance, not an exact confidence certificate"},
            "dependency_groups":[{"N":n,"same_as":f"P154 Phase-E mixed n{n}: original 20000 replicas / 100 batches"} for n in (65,130)],
            "identities":{"J_Q":"J_R+J_H","R_conditional_mean":"E[R|K]=0 exactly under the uniform fixed-K measure",
                          "H_response":"Cov(O,H)=p^2(1-p)^2/(N(N-1))*d^2 E_p[O]/dp^2 for any p-independent O"},
            "limitations":["fixed-K orthogonality is not RG thermal orthogonality","R has a global K counterterm; Q is the local canonical potential",
                           "no new independent replication","no exponent or continuum field assignment","no arbitrary M1/M2d/M2j certificate"],
            "environment":{"python":platform.python_version(),"machine":platform.machine(),"platform":platform.platform()},
            "code":[{"path":str(p.relative_to(ROOT)),"sha256":sha(p)} for p in (Path(__file__),ROOT/"src/p154_fixed_k_interaction_replay.cpp",ROOT/"src/threshold_rank_integer_period_mc.cpp",ROOT/"analysis/p154_fixed_k_interaction_replay.json")]}

def report(x):
    e=lambda n,f:x["by_N"][str(n)]["estimates"][f]
    lines=["# Norm-4 topology under a fixed-occupation local interaction","",
           "## The new interaction response remains unresolved at 20k","",
           "Neither the full local interaction response nor its fixed-K part is resolved on these inherited blocks. The fixed-K z-scores are 0.152 and 0.891; this supplies no positive evidence for a second physical direction, and does not show that the response is explained by the thermal-count term.","",
           "This reanalysis recovers the old Phase-E configurations and adds the missing `E_top × edge interaction` moment. It does not rerun the former connectivity-B scorer or claim a new energy field.","",
           "For each of the two archived sizes, the reported direction difference is divided only by the exact cosine-four difference. Values below are the response per mean NN edge; uncertainty is aligned delete-one-batch standard error.",""]
    for n in (65,130):
        lines.append(f"### N={n}")
        lines.append("")
        for f in ("J_Q","J_R","J_H","J_R0","J_R2"):
            r=e(n,f);lines.append(f"- `{f}` = {r['value']:.9g} ± {r['se']:.3g} (z={r['z']:.3f}).")
        lines.append("")
    j=x["J_R_joint_zero"]
    lines += [f"The two-size fixed-K mixed-response joint zero statistic is `{j['chi_square']:.6g}/2`, nominal p=`{j['nominal_p_value']:.6g}`. This uses an estimated covariance and is not an exact finite-sample model certificate.","",
              "## What distinguishes the three sources","",
              "Write K for occupied sites, T for occupied NN edges and p=.59274605079. There are 2N simple NN edges in each archived torus.","",
              "```text","Q = (T − 4pK + 2Np²)/(2N)","R = (T − 2K(K−1)/(N−1))/(2N)","H = Q−R", "J_X = Cov(E_top,X), E_top=I0+I2", "J_Q = J_R + J_H", "```","",
              "Q is a local pair interaction with positive finite-volume measure proportional to Bernoulli(p) × exp(λQ). Its λ derivative of E_top is J_Q. R removes the entire occupation-count conditional mean: E[R|K]=0. Conditional on K, Q and R generate exactly the same interaction. R includes a global K counterterm and is not itself a strictly local canonical field.","",
              "For any p-independent observable O, `Cov(O,H)=p²(1−p)²/[N(N−1)] × d²E_p[O]/dp²`. Thus H is an explicit second-thermal-score projection, whereas J_R measures how the sector distinguishes local edge arrangements at the same occupied count. This finite-product-measure separation does not prove vanishing RG thermal overlap.","",
              "J_R0 and J_R2 retain the rank-zero and rank-two contributions separately. Their sum is J_R. Neither a cancellation nor an unresolved sum means the underlying two mixed responses vanish.","",
              "## What this adds to the existing pilot","",
              "The old B asks whether nearby endpoints are connected through the complete torus; it is a legitimate global-connectivity readout, not a finite local occupation-cylinder function. Its unresolved mixed H4 response did not prove algebraic dependence. P40 stores the q × fixed-K motif second moment, where q=I2−I0 is configurationwise. E_top=q² requires a third moment absent from its Gram archive; it is not the square of the expected matching function. The present E_top × R measurement is therefore a different scientific question, not a renamed variance-reduction score.","",
              "Matching parity is stated for the fixed Euclidean stencil under the full graph-pair/complement transformation `(G_black,G_white,η,p)→(G_white,G_black,1−η,1−p)`. It does not identify the NN edge set with that of whichever graph becomes black.","",
              "## Sampling and reproducibility","",
              "The original 0578105 backend blob is unchanged. The original N65/N130 seeds and 20000 counter ranges are specified in the manifest. All 400 batch/orientation rows reproduce the original sample counts, K1/K2 sums and I0/I1/I2 sums exactly. Only local edge sufficient statistics were added. Parent and child results share their random streams; they are not independent evidence.","",
              "Pooled unbiased sample covariances are recomputed after removing each aligned batch from both orientations. The JSON retains all leave-one-out vectors and the full covariance, including the exact J_Q/J_R/J_H linear dependence. No inverse of the redundant full matrix is used. Sizes have different original seeds and remain separate covariance blocks.","",
              "## Interpretation and next physical question","",
              "Use the measured split to decide whether an observed local interaction response is explained by the occupation-count second derivative or requires fixed-K geometry. An unresolved result is a completed measurement at these inherited 20k blocks, not a prohibition on other singlets or scales. A resolved fixed-K row would still need transport across geometry/scale and overlap with the original norm-4 residual before a field could be named. This analysis does not fit an exponent, select a continuum operator or start additional production.","",
              "Reproduce by compiling `src/p154_fixed_k_interaction_replay.cpp`, replaying N65 and N130 to `raw/n65.csv` and `raw/n130.csv`, then running this script. The runner refuses to overwrite outputs. Input/output and code hashes are in `latest.json`.",""]
    return "\n".join(lines)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output-dir",type=Path,default=ROOT/"results/p154-fixed-k-interaction")
    args=parser.parse_args();out=args.output_dir
    for name in ("latest.json","REPORT.md"):
        if (out/name).exists(): raise ValueError("refusing to overwrite existing analysis")
    expected="22058703c12b168e844088277c9b61d64b9c1d2c"
    blob=subprocess.check_output(["git","hash-object","src/threshold_rank_integer_period_mc.cpp"],cwd=ROOT,text=True).strip()
    if blob!=expected: raise ValueError("backend differs from original replay source")
    result=analyze(out/"raw")
    (out/"REPORT.md").write_text(report(result))
    result["report_sha256"]=sha(out/"REPORT.md")
    (out/"latest.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"J_R_joint_zero":result["J_R_joint_zero"],"by_N":{n:v["estimates"]["J_R"] for n,v in result["by_N"].items()}},indent=2))

if __name__=="__main__": main()

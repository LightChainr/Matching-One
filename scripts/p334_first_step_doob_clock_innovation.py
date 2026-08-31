#!/usr/bin/env python3
"""Read saved full clocks; quantify the direct-event Doob variance, with no DP."""
from __future__ import annotations

import os
for _key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_key] = "1"

import argparse
import csv
import gzip
import hashlib
import json
from math import comb
from pathlib import Path

import numpy as np
from scipy.stats import binom

ROOT=Path(__file__).resolve().parents[1]
CONTRACT=ROOT/"analysis/p334_first_step_doob_clock_innovation.json"
OUT=ROOT/"results/p334-first-step-doob-clock-innovation"
FIELDS=["risk_count","solved_count","partial_pair_solved_count","h_zero_count","h_all_count",
        "sum_hazard","sum_mu","sum_second_wait_moment","sum_wait_variance",
        "integrated_B","integrated_V","canonical_B","canonical_V",
        "h_zero_integrated_V","h_zero_canonical_V"]


def summarize(sums):
    values=[]
    for orientation in range(2):
        row=dict(zip(FIELDS,sums[orientation]))
        for readout in ("integrated","canonical"):
            b,v=row[readout+"_B"],row[readout+"_V"]
            values.extend((b/row["solved_count"],v/row["solved_count"],b/v))
    return np.asarray(values)


def score(contract):
    result={"schema":contract["schema"],"parent":contract["parent"],"sizes":{},"input_sha256":{}}
    batch_rows=[]
    largest_roundoff_bound_excess=0.
    for n in contract["sizes"]:
        batches=[]
        for batch in range(contract["batches_per_size"]):
            path=ROOT/contract["input_directory"]/f"N{n}.batch{batch:02d}.json.gz"
            result["input_sha256"][str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
            with gzip.open(path,"rt") as stream:
                records=json.load(stream)["records"]
            sums=np.zeros((2,len(FIELDS)))
            for record in records:
                for o,(source,clock) in enumerate(zip(record["source_rows"],record["clocks"])):
                    if source is None:
                        continue
                    sums[o,0]+=1
                    if clock is None or "safe_coefficients" not in clock:
                        continue
                    d=n-source["k0"]
                    h=source["H2"]
                    coeff=clock["safe_coefficients"]
                    survival=np.array([value/comb(d,j) for j,value in enumerate(coeff)])
                    pmf=survival[:-1]-survival[1:]
                    t=np.arange(1,d+1,dtype=float)
                    mu=float(survival[:-1].sum())
                    second=float((survival[:-1]*(2*np.arange(d)+1)).sum())
                    wait_variance=float(pmf@((t-mu)**2))
                    g=binom.sf(source["k0"]+t-1,n,contract["p_ref"])
                    mean=float(pmf@g)
                    v_can=float(pmf@((g-mean)**2))
                    v_int=wait_variance/(n+1)**2
                    if h==0 or h==d:
                        b_int=b_can=0.
                    else:
                        b_int=h/(d-h)*(mu-1)**2/(n+1)**2
                        b_can=h/(d-h)*(g[0]-mean)**2
                    largest_roundoff_bound_excess=max(largest_roundoff_bound_excess,b_int-v_int,b_can-v_can)
                    row=[1,record["status"]=="whole_pair_fallback",h==0,h==d,h/d,mu,second,wait_variance,
                         b_int,v_int,b_can,v_can,v_int if h==0 else 0.,v_can if h==0 else 0.]
                    sums[o,1:]+=row
            batches.append(sums)
            for o,name in enumerate(("first","second")):
                batch_rows.append({"N":n,"batch":batch,"orientation":name,
                                   "counter_first":records[0]["counter"],"counter_last":records[-1]["counter"],
                                   **dict(zip(FIELDS,sums[o].tolist()))})
        batches=np.asarray(batches)
        total=batches.sum(axis=0)
        estimate=summarize(total)
        leave=np.array([summarize(total-b) for b in batches])
        centered=leave-leave.mean(axis=0)
        covariance=(len(batches)-1)/len(batches)*centered.T@centered
        labels=[f"{o}_{r}_{stat}" for o in ("first","second") for r in ("integrated","canonical")
                for stat in ("mean_bound","mean_conditional_variance","fraction")]
        size={"statistic_labels":labels,"estimates":estimate.tolist(),"joint_delete_one_batch_covariance":covariance.tolist(),
              "delete_one_batch_estimates":leave.tolist(),"orientations":{}}
        for o,name in enumerate(("first","second")):
            totals=dict(zip(FIELDS,total[o].tolist()))
            row={"counts":{key:int(totals[key]) for key in FIELDS if key.endswith("count")},
                 "solved_fraction_of_rank_one_rows":totals["solved_count"]/totals["risk_count"],
                 "mean_direct_next_probability":totals["sum_hazard"]/totals["solved_count"],
                 "mean_wait":totals["sum_mu"]/totals["solved_count"],
                 "mean_conditional_wait_variance":totals["sum_wait_variance"]/totals["solved_count"]}
            for r,readout in enumerate(("integrated","canonical")):
                j=6*o+3*r
                row[readout]={"mean_bound":estimate[j],"mean_conditional_suffix_variance":estimate[j+1],
                              "direct_event_fraction":estimate[j+2],"fraction_batch_se":float(np.sqrt(covariance[j+2,j+2])),
                              "remaining_safe_event_fraction":1-estimate[j+2],
                              "h_zero_share_of_total_suffix_variance":totals["h_zero_"+readout+"_V"]/totals[readout+"_V"]}
            size["orientations"][name]=row
        result["sizes"][str(n)]=size
        print(n,json.dumps(size["orientations"],indent=2),flush=True)
    result["largest_numeric_bound_minus_variance_if_positive"]=largest_roundoff_bound_excess
    result["scope"]="Solved single-orientation clock pool only; binary direct-event innovation, not paired H4 or variance necessarily surviving the full next label."
    return result,batch_rows


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,default=OUT)
    args=parser.parse_args()
    contract=json.loads(CONTRACT.read_text())
    result,batches=score(contract)
    result["input_sha256"].update({str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (CONTRACT,Path(__file__))})
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/"score.json").write_text(json.dumps(result,indent=2)+"\n")
    with (args.output/"batches.csv").open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(batches[0]))
        writer.writeheader();writer.writerows(batches)


if __name__=="__main__":
    main()

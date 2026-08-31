#!/usr/bin/env python3
"""Evaluate existing exact real-prefix F and all D_v at equal survival1/2.

No reliability solve or new source data: only high-precision polynomial
evaluation of the two previously completed1c06230b marked-site archives.
"""
import json
from pathlib import Path

import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/"results/p334-real-prefix-median-noise-energy"
COUNTERS=(43042514269,43042505280)
mp.mp.dps=80


def record(x):
    return {"decimal":mp.nstr(x,50),"value":float(x)}


def bernoulli_count_polynomial(coefficients,u,dimension):
    z=u/(1-u)
    value=mp.mpf(0)
    for coefficient in reversed(coefficients):
        value=value*z+coefficient
    return (1-u)**dimension*value


def main():
    clocks=json.loads((ROOT/"results/p334-contracted-full-clock/full_physical_birth_clock.json").read_text())["records"]
    clocks={r["counter"]:r for r in clocks}
    rows=[]
    for counter in COUNTERS:
        path=ROOT/f"results/p334-exact-marked-birth/marked_birth_{counter}.json"
        marked=json.loads(path.read_text())
        counts=clocks[counter]["true_safe_counts"]
        d=len(counts)-1
        left,right=mp.mpf(0),mp.mpf(1)
        for _ in range(180):
            middle=(left+right)/2
            if bernoulli_count_polynomial(counts,middle,d)>mp.mpf("0.5"):
                left=middle
            else:
                right=middle
        u=(left+right)/2
        survival=bernoulli_count_polynomial(counts,u,d)
        derivative_counts=[(d-k)*counts[k]-(k+1)*counts[k+1] for k in range(d)]
        slope=bernoulli_count_polynomial(derivative_counts,u,d-1)
        sites=[]
        influences=[]
        for site in marked["site_records"]:
            influence=bernoulli_count_polynomial(site["pivotal_count_by_prior_size"],u,d-1)
            influences.append(influence)
            sites.append({"site":site["site"],"type":site["type"],"pivotal_probability":record(influence),
                          "singleton_energy":record(u*(1-u)*influence**2)})
        square_sum=sum(i*i for i in influences)
        concentration=square_sum/slope**2
        clock_prefactor=u*(1-u)*slope**2
        E1=u*(1-u)*square_sum
        total=survival*(1-survival)
        symmetric=clock_prefactor/d
        row={"counter":counter,"source_marked_artifact":str(path.relative_to(ROOT)),
             "N":marked["N"],"k0":marked["k0"],"remaining_sites":d,
             "median_u":record(u),"root_bracket":[record(left),record(right)],
             "safe_probability":record(survival),"clock_slope_magnitude":record(slope),
             "sum_site_pivotal_probabilities":record(sum(influences)),
             "sum_squared_site_pivotal_probabilities":record(square_sum),
             "clock_slope_squared_prefactor_u1mu_slope2":record(clock_prefactor),
             "pivotal_concentration":record(concentration),
             "inverse_pivotal_concentration_effective_sites":record(1/concentration),
             "E1":record(E1),"total_indicator_variance":record(total),
             "degree2plus_energy":record(total-E1),"E1_fraction_total_variance":record(E1/total),
             "clock_determined_uniform_singleton_projection_energy":record(symmetric),
             "spatial_singleton_excess_energy":record(E1-symmetric),
             "spatial_excess_fraction_E1":record((E1-symmetric)/E1),
             "E1_to_clock_uniform_projection_ratio":record(E1/symmetric),
             "positive_pivotal_support":sum(any(s["pivotal_count_by_prior_size"]) for s in marked["site_records"]),
             "sites":sites}
        rows.append(row)
        print(counter, {k:row[k]["decimal"] for k in ("median_u","clock_slope_magnitude","pivotal_concentration","E1","degree2plus_energy","clock_determined_uniform_singleton_projection_energy","spatial_singleton_excess_energy")})
    A,B=rows
    ratio_keys=("clock_slope_squared_prefactor_u1mu_slope2","pivotal_concentration","E1")
    result={"parent_commit":"614eedb2429d74d6b4de7ebf15d6c8f918b54e3c",
            "marked_source_commit":"1c06230b8f7e13be98f128361ad72b23c0c425ae",
            "full_clock_source_commit":"6358ba49ef390c10a3f501b589ba7ba1d4e05b09",
            "mpmath_decimal_precision":80,"bisection_steps":180,"new_samples":0,"new_network_DP_solves":0,
            "selection":"Only the original two solved real N425 witnesses, not the later15-direct-gate pair or a147-row estimate",
            "estimand":"At each fixed prefix choose u with S_prefix(u)=1/2 among its d173 remaining sites; E1=u(1-u)sum I_v(u)^2. This u is not the original full-N canonical occupation p.",
            "records":rows,
            "B_over_A_multiplicative_decomposition":{k:record(mp.mpf(B[k]["decimal"])/mp.mpf(A[k]["decimal"])) for k in ratio_keys},
            "scope":"Deterministic polynomial evaluations on two posthoc real witnesses, same existing source block. Not paired/stratum-weighted cohort inference; no standard errors or population claim."}
    OUTPUT.mkdir(parents=True,exist_ok=True)
    (OUTPUT/"iso_survival_median_energy.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("ratios",result["B_over_A_multiplicative_decomposition"])


if __name__=="__main__":
    main()

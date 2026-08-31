#!/usr/bin/env python3
"""Eight prescribed middle-site blockades on two archived physical prefixes.

Only affected factor DPs are evaluated. Blocked sites remain inert dummies
in the same d=173 uniform insertion permutation.
"""
from fractions import Fraction
import json
from math import comb
from pathlib import Path
import time

from p334_full_birth_reliability import safety_polynomial
from p334_pair_only_survival import contiguous, frac, multiply

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"results/p334-all147-prefix-clocks/prefixes"
OUTPUT=ROOT/"results/p334-middle-bridge-interventions"
COUNTERS=(43042508631,43042514803)
BLOCKS={43042508631:((198,),),
        43042514803:((24,),(25,),(184,),(340,),(94,),(361,),(24,25,94,184,340,361))}
D=173


def summary(coefficients):
    survival=[Fraction(c,comb(D,k)) for k,c in enumerate(coefficients)]
    return {"safe_coefficients":coefficients,"survival":[frac(s) for s in survival],
            "mean_birth_step":frac(sum(survival[:-1])),
            "tail40":frac(survival[40]),"tail20":frac(survival[20]),
            "tail65":frac(survival[65]),
            "direct_share":frac(15*sum(Fraction(coefficients[k],D*comb(D-1,k)) for k in range(D))),
            "maximum_safe_k":max(k for k,c in enumerate(coefficients) if c)}


def main():
    marked=json.loads((ROOT/"results/p334-marked-triple-bridges/middle_site_bridges.json").read_text())
    marked={r["counter"]:r for r in marked["records"]}
    previous_path=OUTPUT/"physical_clock_blockades.json"
    previous=json.loads(previous_path.read_text()) if previous_path.exists() else {}
    factor_cache={(r["counter"],r["factor"],tuple(r["blocked_sites"])):r
                  for r in previous.get("affected_factor_polynomials",[])}
    new_factor_solves=0
    records=[]
    started=time.monotonic()
    for counter in COUNTERS:
        source=json.loads((SOURCE/f"{counter}.json").read_text())
        clock=source["clock"]
        baseline=summary(clock["true_safe_counts"])
        networks=[c["two_terminal_network"] for c in source["mapping"]["port_components"] if "two_terminal_network" in c]
        all_middle={m["middle_site_id"] for m in marked[counter]["positive_middle_sites"]}
        weights={m["middle_site_id"]:m["minimal_triple_contribution"] for m in marked[counter]["positive_middle_sites"]}
        base_g3=marked[counter]["genuine_minimal_triples"]
        rows=[]
        for prescribed in BLOCKS[counter]:
            blocked=set(prescribed)
            coefficients=[1]
            factors=[]
            for fi,(network,frozen_factor) in enumerate(zip(networks,source["factors"])):
                local=blocked & set(network["vacant_sites"])
                if not local:
                    factor_coefficients=frozen_factor["safe_coefficients"]
                    factors.append({"factor":fi,"reused_baseline":True})
                else:
                    key=(counter,fi,tuple(sorted(local)))
                    if key not in factor_cache:
                        reduced={**network,
                                 "vertices":[v for v in network["vertices"] if v not in local],
                                 "edges":[e for e in network["edges"] if not (set(e)&local)],
                                 "vacant_sites":[v for v in network["vacant_sites"] if v not in local]}
                        active,stats=safety_polynomial(reduced,set(reduced["vacant_sites"]))
                        new_factor_solves+=1
                        factor_cache[key]={"counter":counter,"factor":fi,"blocked_sites":sorted(local),
                                           "active_polynomial":active,
                                           "dummy_polynomial":[comb(len(local),k) for k in range(len(local)+1)],
                                           "safe_polynomial_on_original_sites":multiply(active,[comb(len(local),k) for k in range(len(local)+1)]),
                                           "dp":{k:v for k,v in stats.items() if k!="bag_state_counts"}}
                    data=factor_cache[key]
                    factor_coefficients=data["safe_polynomial_on_original_sites"]
                    factors.append({"factor":fi,"reused_baseline":False,"cache_key":list(key[:2])+[list(key[2])]})
                coefficients=multiply(coefficients,factor_coefficients)
            free=clock["structure"]["off_core_random_sites"]
            coefficients=multiply(coefficients,[comb(free,k) for k in range(free+1)])
            readout=summary(coefficients)
            delta=[off-on for off,on in zip(coefficients,clock["true_safe_counts"])]
            removed=sum(weights[y] for y in blocked)
            readout.update(blocked_sites=sorted(blocked),all_middle_blocked=(blocked==all_middle),
                           aliases=["single","all_middle"] if len(all_middle)==1 else ["all_middle" if blocked==all_middle else "single"],
                           unchanged_original_d=D,g3_remaining=base_g3-removed,g3_removed=removed,
                           f3_increase=delta[3],
                           mean_increase=frac(Fraction(readout["mean_birth_step"]["exact"])-Fraction(baseline["mean_birth_step"]["exact"])),
                           tail40_increase=frac(Fraction(readout["tail40"]["exact"])-Fraction(baseline["tail40"]["exact"])),
                           tail40_ratio=frac(Fraction(readout["tail40"]["exact"])/Fraction(baseline["tail40"]["exact"])),
                           survival_strictly_above_baseline=contiguous([k for k,v in enumerate(delta) if v>0]),
                           first_changed_coefficient=next((k for k,v in enumerate(delta) if v),None),
                           factors=factors)
            if len(blocked)==1:
                final_law=[Fraction(0)]+[Fraction(delta[k],k*comb(D,k)) for k in range(1,D+1)]
                pi=sum(final_law)
                readout["baseline_final_site_readout_from_knockout"]={
                    "site":next(iter(blocked)),"probability_final_birth_site":frac(pi),
                    "mean_birth_step_conditional_on_this_site_final":frac(sum(k*p for k,p in enumerate(final_law))/pi),
                    "joint_T_final_site_probability":[frac(p) for p in final_law],
                    "identity":"Delta S(k)=k P_baseline(T=k,Vfinal=blocked_site); Delta E[T]=pi_v E[T|Vfinal=v]"}
            rows.append(readout)
            print(counter,sorted(blocked),"g3",readout["g3_remaining"],"mean",readout["mean_birth_step"]["value"],
                  "delta",readout["mean_increase"]["value"],"S40",readout["tail40"]["value"],"ratio",readout["tail40_ratio"]["value"],flush=True)
        joint=next(r for r in rows if r["all_middle_blocked"])
        single=[r for r in rows if len(r["blocked_sites"])==1]
        nonadd=[joint["safe_coefficients"][k]-baseline["safe_coefficients"][k]
                -sum(r["safe_coefficients"][k]-baseline["safe_coefficients"][k] for r in single)
                for k in range(D+1)]
        baseline["g3"]=base_g3
        records.append({"counter":counter,"source_artifact":str((SOURCE/f"{counter}.json").relative_to(ROOT)),
                        "baseline_reused":baseline,"interventions":rows,
                        "joint_mean_increase_minus_sum_single_increases":frac(Fraction(joint["mean_increase"]["exact"])-sum(Fraction(r["mean_increase"]["exact"]) for r in single)),
                        "joint_tail40_increase_minus_sum_single_increases":frac(Fraction(joint["tail40_increase"]["exact"])-sum(Fraction(r["tail40_increase"]["exact"]) for r in single)),
                        "joint_minus_sum_single_polynomial":nonadd,
                        "nonadditivity_first_k":next((k for k,v in enumerate(nonadd) if v),None),
                        "nonadditivity_positive_k":contiguous([k for k,v in enumerate(nonadd) if v>0]),
                        "nonadditivity_negative_k":contiguous([k for k,v in enumerate(nonadd) if v<0])})
    comparisons={}
    for label in ("baseline","all_middle_blocked"):
        a,b=[r["baseline_reused"] if label=="baseline" else next(i for i in r["interventions"] if i["all_middle_blocked"]) for r in records]
        diff=[y-x for x,y in zip(a["safe_coefficients"],b["safe_coefficients"])]
        comparisons[label]={"mean_B_minus_A":frac(Fraction(b["mean_birth_step"]["exact"])-Fraction(a["mean_birth_step"]["exact"])),
                            "tail40_B_minus_A":frac(Fraction(b["tail40"]["exact"])-Fraction(a["tail40"]["exact"])),
                            "B_survival_above_A":contiguous([k for k,v in enumerate(diff) if v>0]),
                            "B_survival_below_A":contiguous([k for k,v in enumerate(diff) if v<0]),
                            "equal_survival":contiguous([k for k,v in enumerate(diff) if not v]),
                            "safe_coefficients_first6_A":a["safe_coefficients"][:6],
                            "safe_coefficients_first6_B":b["safe_coefficients"][:6]}
    # The saved pair-edge lists identify the same two small tree components,
    # P3 and a five-vertex double star. Their independence polynomials are
    # closed-form counts, not another network DP or hyperedge enumeration.
    pair_core=multiply([1,3,1],[1,5,6,2])
    pair_only=multiply([comb(150,k) for k in range(151)],pair_core)+[0]*18
    pair_reduction={"safe_non_direct_sites":158,"pair_incident_sites":8,"safe_pair_irrelevant_sites":150,
                    "component_independence_polynomials":[[1,3,1],[1,5,6,2]],
                    "component_types":["P3","five-vertex double star with degrees 3,2,1,1,1"],
                    "pair_edges_by_counter":{str(c):[e for f in marked[c]["factors"] for e in f["minimal_pair_edges"]] for c in COUNTERS},
                    "pair_graph_isomorphism_A_to_B":{129:204,130:202,286:360,40:27,39:180,196:338,354:71,355:72},
                    "common_pair_only_polynomial":pair_only,
                    "common_formula":"(1+z)^150 (1+3z+z^2) (1+5z+6z^2+2z^3)",
                    "equivalent_algebraic_factorization":"(1+z)^151 (1+7z+15z^2+10z^3+2z^4); the extra (1+z) is algebraic, not another Boolean-irrelevant site",
                    "joint_polynomial_equals_pair_only_by_counter":{str(r["counter"]):next(i for i in r["interventions"] if i["all_middle_blocked"])["safe_coefficients"]==pair_only for r in records},
                    "event_consequence":"Every pair-trigger is a physical trigger; all joint-off physical safe coefficients equal the pair-only coefficients. Inclusion plus equal counts at every cardinality makes the finite events identical. Thus each original minimal trigger of order>=3 intersects the prescribed middle-site set, for these two prefixes only."}
    OUTPUT.mkdir(parents=True,exist_ok=True)
    output={"parent_commit":"fd96cd95f50a2c2ed7020d180acd0b6f6e0c9747",
            "new_samples":0,"baseline_network_solves":0,"prescribed_interventions":8,
            "factor_dp_solves":len(factor_cache),"new_factor_dp_solves_this_invocation":new_factor_solves,
            "elapsed_seconds":previous.get("elapsed_seconds",time.monotonic()-started),
            "last_summary_invocation_seconds":time.monotonic()-started,
            "clock_semantics":"Original d173 uniform insertion permutation. Every blocked site is a permanently vacant inert dummy retaining its insertion label. Each active factor polynomial is multiplied by (1+z)^number_of_blocked_factor_sites.",
            "records":records,"two_prefix_comparisons":comparisons,
            "common_pair_core_reduction":pair_reduction,
            "affected_factor_polynomials":list(factor_cache.values())}
    (OUTPUT/"physical_clock_blockades.json").write_text(json.dumps(output,indent=2,sort_keys=True)+"\n")
    print("comparison",comparisons,"factor solves",len(factor_cache),"seconds",output["elapsed_seconds"])


if __name__=="__main__":
    main()

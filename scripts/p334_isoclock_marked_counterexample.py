#!/usr/bin/env python3
"""Two prescribed five-vertex pair-trigger graphs, 32 subsets each.

Recover exact (birth step, final vertex) counts among all120 insertion orders
by weighting pivotal preceding subsets, not by sampling or graph-family search.
"""
from fractions import Fraction
import json
from math import comb, factorial
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/"results/p334-isoclock-marked-counterexample"
GRAPHS={"double_star":{"edges":[[0,1],[0,2],[0,3],[1,4]],
                       "roles":["degree3 hub","degree2 hub","degree3-side leaf","degree3-side leaf","degree2-side leaf"]},
        "C4_plus_isolate":{"edges":[[0,1],[1,2],[2,3],[3,0]],
                           "roles":["cycle vertex","cycle vertex","cycle vertex","cycle vertex","inert isolate"]}}
D=5


def frac(x):
    x=Fraction(x)
    return {"exact":str(x),"value":float(x)}


def score(name,graph):
    safe={mask:not any(mask&(1<<u) and mask&(1<<v) for u,v in graph["edges"])
          for mask in range(1<<D)}
    independent=[sum(ok and mask.bit_count()==k for mask,ok in safe.items()) for k in range(D+1)]
    survival=[Fraction(independent[k],comb(D,k)) for k in range(D+1)]
    mean=sum(survival[:-1])
    joint=[[0]*(D+1) for _ in range(D)]
    pivotal=[[0]*D for _ in range(D)]
    for mask,ok in safe.items():
        if not ok:continue
        k=mask.bit_count()+1
        for v in range(D):
            if not mask&(1<<v) and not safe[mask|(1<<v)]:
                pivotal[v][k-1]+=1
                joint[v][k]+=factorial(k-1)*factorial(D-k)
    vertices=[]
    for v in range(D):
        p=[Fraction(count,factorial(D)) for count in joint[v]]
        pi=sum(p)
        knockout=sum(k*p[k] for k in range(D+1))
        vertices.append({"vertex":v,"role":graph["roles"][v],
                         "degree":sum(v in e for e in graph["edges"]),
                         "pivotal_preceding_subset_counts":pivotal[v],
                         "joint_T_final_vertex_sequence_counts_out_of120":joint[v],
                         "joint_T_final_vertex_probabilities":[frac(x) for x in p],
                         "pi":frac(pi),"inert_knockout_mean_increase":frac(knockout),
                         "normalized_mean_knockout_weight":frac(knockout/mean),
                         "mean_T_conditional_final_vertex":frac(knockout/pi) if pi else None})
    totals=[sum(joint[v][k] for v in range(D)) for k in range(D+1)]
    return {"name":name,"edges":graph["edges"],"vertices":vertices,
            "subsets_evaluated":1<<D,"uniform_sequence_denominator":factorial(D),
            "independence_polynomial":independent,"survival":[frac(s) for s in survival],
            "birth_sequence_counts_by_T":totals,"birth_probabilities_by_T":[frac(Fraction(n,factorial(D))) for n in totals],
            "mean_birth_step":frac(mean),"positive_pivotal_support":sum(Fraction(v["pi"]["exact"])>0 for v in vertices),
            "sum_pi_squared":frac(sum(Fraction(v["pi"]["exact"])**2 for v in vertices)),
            "two_independent_continuations_same_final_vertex_probability":frac(sum(Fraction(v["pi"]["exact"])**2 for v in vertices)),
            "time_conditional_final_vertex_collision":{str(k):frac(sum(Fraction(joint[v][k],totals[k])**2 for v in range(D))) for k in range(D+1) if totals[k]},
            "normalized_mean_knockout_concentration":frac(sum(Fraction(v["normalized_mean_knockout_weight"]["exact"])**2 for v in vertices)),
            "pi_multiset":sorted([v["pi"]["exact"] for v in vertices],key=Fraction),
            "safe_subset_masks":[mask for mask,ok in safe.items() if ok]}


def main():
    records=[score(name,graph) for name,graph in GRAPHS.items()]
    OUTPUT.mkdir(parents=True,exist_ok=True)
    result={"parent_commit":"0143632db59d867cfb658a6ad4465e5036684fff",
            "origin":"Exact independence-polynomial factor discovered in the actual P334 residual pair core; the alternate C4+isolate graph is a constructed finite mechanism counterexample, not another N425 checkpoint.",
            "graph_selection":"Exactly the prescribed double star and C4+one inert isolate; no family search",
            "new_N425_samples":0,"new_Monte_Carlo_samples":0,"records":records,
            "polynomial_identity":"I_double-star=1+5z+6z^2+2z^3=(1+z)(1+4z+2z^2)=I_C4+isolated-site",
            "common_birth_probabilities":{"2":"2/5","3":"2/5","4":"1/5"},
            "mark_interpretation":"V is the actual final insertion site at the first pair-trigger birth; all five sites retain their insertion labels, including the inert isolate."}
    (OUTPUT/"exact_marked_clock.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    for r in records:
        print(r["name"],"I",r["independence_polynomial"],"T",r["birth_sequence_counts_by_T"],
              "support",r["positive_pivotal_support"],"pi2",r["sum_pi_squared"],
              "knockout concentration",r["normalized_mean_knockout_concentration"])
        for v in r["vertices"]:
            print(v["vertex"],v["role"],v["joint_T_final_vertex_sequence_counts_out_of120"],
                  "pi",v["pi"]["exact"],"mean knockout",v["inert_knockout_mean_increase"]["exact"])


if __name__=="__main__":
    main()

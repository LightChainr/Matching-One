#!/usr/bin/env python3
"""All fixed-row clock distributions, state widths and exact crossing census."""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations
import json
from math import comb
from pathlib import Path

from p334_pair_only_survival import frac,contiguous

ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/"results/p334-all147-prefix-clocks"


def distribution(values):
    ordered=sorted(values)
    out={"n":len(ordered),"min":frac(ordered[0]),"max":frac(ordered[-1]),"mean":frac(sum(ordered)/len(ordered))}
    for q in (Fraction(1,4),Fraction(1,2),Fraction(3,4)):
        at=q*(len(ordered)-1); low=at.numerator//at.denominator; part=at-low
        out[str(q)]=frac(ordered[low]*(1-part)+ordered[min(low+1,len(ordered)-1)]*part)
    return out


def group_width(rows):
    low=min(rows,key=lambda r:Fraction(r["mean_true_birth_step"]["exact"]))
    high=max(rows,key=lambda r:Fraction(r["mean_true_birth_step"]["exact"]))
    tails=[Fraction(r["tail_after_40"]["exact"]) for r in rows]
    direct=[Fraction(r["direct_share"]["exact"]) for r in rows]
    return {"count":len(rows),"counters":[r["counter"] for r in rows],
            "distinct_full_polynomials":len({tuple(r["true_safe_counts"]) for r in rows}),
            "mean_min":low["mean_true_birth_step"],"mean_max":high["mean_true_birth_step"],
            "mean_range":frac(Fraction(high["mean_true_birth_step"]["exact"])-Fraction(low["mean_true_birth_step"]["exact"])),
            "mean_extreme_counters":[low["counter"],high["counter"]],
            "tail40_min":frac(min(tails)),"tail40_max":frac(max(tails)),
            "direct_share_min":frac(min(direct)),"direct_share_max":frac(max(direct))}


def main():
    batch=json.loads((OUTPUT/"full_clocks.json").read_text())
    rows=[r for r in batch["records"] if r["status"]=="solved_full_physical"]
    by_h2,by_pair,by_branch=defaultdict(list),defaultdict(list),defaultdict(list)
    for r in rows:
        h,b2=int(r["original_row"]["H2"]),int(r["original_row"]["checkpoint_b2_safe_pairs"])
        s_mean=sum(Fraction(r["true_safe_counts"][k],173*comb(172,k)) for k in range(173))
        r["mean_log_uniform_birth_time"]=frac(s_mean)
        r["direct_share"]=frac(h*s_mean)
        r["collective_share"]=frac(1-h*s_mean)
        by_h2[h].append(r);by_pair[h,b2].append(r)
        by_branch[h,b2,int(r["original_row"]["checkpoint_sum_child_b1_sq"])].append(r)
    crossings=[];dominance_counts=defaultdict(int)
    for A,B in combinations(rows,2):
        differences=[b-a for a,b in zip(A["true_safe_counts"],B["true_safe_counts"])]
        pos=[k for k,v in enumerate(differences) if v>0];neg=[k for k,v in enumerate(differences) if v<0]
        if pos and neg:
            previous=None;switches=[]
            for k,value in enumerate(differences):
                if not value:continue
                sign=1 if value>0 else -1
                if previous is not None and sign!=previous:
                    switches.append(k)
                previous=sign
            crossings.append({"A":A["counter"],"B":B["counter"],"B_above":contiguous(pos),"B_below":contiguous(neg),
                              "sign_switch_steps":switches,"first_switch_step":switches[0],
                              "first_switch_difference_count":differences[switches[0]],
                              "first_switch_survival_A":A["true_survival"][switches[0]],
                              "first_switch_survival_B":B["true_survival"][switches[0]]})
        else:
            dominance_counts["B_dominates" if pos else "A_dominates" if neg else "identical"]+=1
    all_collision_pairs=[]
    for key,group in by_pair.items():
        for A,B in combinations(group,2):
            if A["true_safe_counts"]!=B["true_safe_counts"]:
                all_collision_pairs.append((A,B))
    both_new=[pair for pair in all_collision_pairs if all(not r["reused_original_twelve"] for r in pair)]
    candidates=both_new or all_collision_pairs
    witness=None
    if candidates:
        A,B=min(candidates,key=lambda pair:(pair[0]["counter"],pair[1]["counter"]))
        k=next(k for k,(a,b) in enumerate(zip(A["true_safe_counts"],B["true_safe_counts"])) if a!=b)
        witness={"selection":"lexicographically first differing pair among newly evaluated135; fallback all147 only if no such new/new pair",
                 "both_new135":bool(both_new),"counters":[A["counter"],B["counter"]],
                 "H2":int(A["original_row"]["H2"]),"b2":int(A["original_row"]["checkpoint_b2_safe_pairs"]),
                 "first_differing_k":k,"first_differing_safe_counts":[A["true_safe_counts"][k],B["true_safe_counts"][k]],
                 "source_sum_child_b1_sq":[int(r["original_row"]["checkpoint_sum_child_b1_sq"]) for r in (A,B)],
                 "means":[r["mean_true_birth_step"] for r in (A,B)],
                 "mean_B_minus_A":frac(Fraction(B["mean_true_birth_step"]["exact"])-Fraction(A["mean_true_birth_step"]["exact"])),
                 "tail40":[r["tail_after_40"] for r in (A,B)],"direct_share":[r["direct_share"] for r in (A,B)],
                 "prefix_artifacts":[r.get("prefix_artifact",r.get("mapping_artifact")) for r in (A,B)]}
    branch_collisions=[(A,B) for group in by_branch.values() for A,B in combinations(group,2)
                       if A["true_safe_counts"]!=B["true_safe_counts"]]
    branch_witness=None
    if branch_collisions:
        A,B=min(branch_collisions,key=lambda pair:(pair[0]["counter"],pair[1]["counter"]))
        h=int(A["original_row"]["H2"]); a=173-h
        b2=int(A["original_row"]["checkpoint_b2_safe_pairs"])
        squares=int(A["original_row"]["checkpoint_sum_child_b1_sq"])
        m=comb(a,2)-b2
        trigger_squares=squares-a*(a-1)**2+4*(a-1)*m
        wedges=(trigger_squares-2*m)//2
        # The existing NN-torus rank-one two-cycle lemma makes the minimal
        # pair-trigger graph bipartite. Therefore its triangle term is zero.
        pair_safe3=comb(a,3)-m*(a-2)+wedges
        branch_witness={"selection":"lexicographically first exact (H2,b2,sum_child_b1_sq) collision with different full clocks; only one such pair in fixed147",
                        "counters":[A["counter"],B["counter"]],"H2":h,"safe_singletons":a,"b2":b2,
                        "source_sum_child_b1_sq":squares,"minimal_pair_edges":m,
                        "trigger_degree_square_sum":trigger_squares,"trigger_wedges":wedges,
                        "trigger_triangles":0,"pair_only_safe3":pair_safe3,
                        "true_safe3":[r["true_safe_counts"][3] for r in (A,B)],
                        "genuine_minimal_triples":[pair_safe3-r["true_safe_counts"][3] for r in (A,B)],
                        "safe3_difference_A_minus_B":A["true_safe_counts"][3]-B["true_safe_counts"][3],
                        "survival3_difference_A_minus_B":frac(Fraction(A["true_safe_counts"][3]-B["true_safe_counts"][3],comb(173,3))),
                        "safe_coefficients_first6":[r["true_safe_counts"][:6] for r in (A,B)],
                        "means":[r["mean_true_birth_step"] for r in (A,B)],
                        "tail40":[r["tail_after_40"] for r in (A,B)],
                        "direct_share":[r["direct_share"] for r in (A,B)],
                        "prefix_artifacts":[r["prefix_artifact"] for r in (A,B)],
                        "method":"Exact coefficient and archived-degree arithmetic; no triple enumeration or new samples. Source squares are the sum of squared safe successor counts (safe-pair graph degrees). Triangle=0 uses the existing rank-one black-NN torus bipartite-trigger lemma."}
    table=[{"counter":r["counter"],"H2":int(r["original_row"]["H2"]),"b2":int(r["original_row"]["checkpoint_b2_safe_pairs"]),
            "sum_child_b1_sq":int(r["original_row"]["checkpoint_sum_child_b1_sq"]),
            "mean":r["mean_true_birth_step"],"tail40":r["tail_after_40"],"direct_share":r["direct_share"],"collective_share":r["collective_share"],
            "mean_log_uniform_birth_time":r["mean_log_uniform_birth_time"],
            "parallel_factors":r["structure"]["parallel_two_port_factors"],"core_sites":r["structure"]["core_random_sites"],
            "treewidth":r["structure"]["treewidth_upper_bound"],"new_seconds":r["new_evaluation_seconds"]} for r in rows]
    summary={"manifest_commit":"8d7ac0e9","selected":batch["selected_rows"],"completed":len(rows),"reused":batch["reused_rows"],
             "new_evaluated":batch["new_rows"],"new_samples":0,"batch_wall_seconds":batch["batch_wall_seconds"],
             "incomplete_rows":[r for r in batch["records"] if r["status"]!="solved_full_physical"],
             "mean_distribution":distribution([Fraction(r["mean_true_birth_step"]["exact"]) for r in rows]),
             "tail40_distribution":distribution([Fraction(r["tail_after_40"]["exact"]) for r in rows]),
             "direct_share_distribution":distribution([Fraction(r["direct_share"]["exact"]) for r in rows]),
             "within_H2":{str(h):group_width(group) for h,group in sorted(by_h2.items())},
             "repeated_H2_b2_groups":[{"H2":h,"b2":b2,**group_width(group)} for (h,b2),group in sorted(by_pair.items()) if len(group)>1],
             "same_H2_b2_different_clock_pairs":len(all_collision_pairs),
             "same_H2_b2_different_clock_new_new_pairs":len(both_new),"new_collision_witness":witness,
             "same_H2_b2_degree_square_different_clock_pairs":len(branch_collisions),
             "branching_scalar_collision_witness":branch_witness,
             "total_completed_pairs":len(rows)*(len(rows)-1)//2,"crossing_pairs":len(crossings),
             "crossings_first_switch_le40":sum(r["first_switch_step"]<=40 for r in crossings),
             "noncrossing_pair_counts":dict(dominance_counts),"crossing_certificates":crossings,"table":table}
    (OUTPUT/"scientific_summary.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")
    print("completed",len(rows),"crossings",len(crossings),"early",summary["crossings_first_switch_le40"],"same-state pairs",len(all_collision_pairs),"new/new",len(both_new))
    print("mean range",summary["mean_distribution"]["min"]["value"],summary["mean_distribution"]["max"]["value"],"median",summary["mean_distribution"]["1/2"]["value"])
    print("tail40 range",summary["tail40_distribution"]["min"]["value"],summary["tail40_distribution"]["max"]["value"])
    print("witness",witness)
    for h,g in summary["within_H2"].items():
        if g["count"]>1: print("H2",h,"count",g["count"],"mean",g["mean_min"]["value"],g["mean_max"]["value"],"tail",g["tail40_min"]["value"],g["tail40_max"]["value"])


if __name__=="__main__":
    main()

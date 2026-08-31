#!/usr/bin/env python3
"""Zero-new-work scientific comparisons of the frozen twelve full clocks."""
from fractions import Fraction
from itertools import combinations
import json
from math import comb
from pathlib import Path

from p334_pair_only_survival import contiguous, frac

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results/p334-twelve-prefix-clocks"


def compare(A, B):
    difference = [b-a for a,b in zip(A["true_safe_counts"], B["true_safe_counts"])]
    cross = {k: A["true_safe_counts"][k]*B["true_safe_counts"][k-1] - B["true_safe_counts"][k]*A["true_safe_counts"][k-1]
             for k in range(1,174) if A["true_safe_counts"][k-1] and B["true_safe_counts"][k-1]}
    return {"A_counter": A["counter"], "B_counter": B["counter"],
            "survival_difference_count_B_minus_A": difference,
            "survival_B_above": contiguous([k for k,v in enumerate(difference) if v>0]),
            "survival_B_below": contiguous([k for k,v in enumerate(difference) if v<0]),
            "hazard_difference_crossproducts_B_minus_A": cross,
            "hazard_B_above": contiguous([k for k,v in cross.items() if v>0]),
            "hazard_B_below": contiguous([k for k,v in cross.items() if v<0]),
            "mean_B_minus_A": frac(Fraction(B["mean_true_birth_step"]["exact"])-Fraction(A["mean_true_birth_step"]["exact"]))}


def main():
    records = json.loads((OUTPUT / "full_clocks.json").read_text())["records"]
    solved = [r for r in records if r["status"] == "solved_full_physical"]
    crossings = []
    for A,B in combinations(solved,2):
        comparison = compare(A,B)
        if comparison["survival_B_above"] and comparison["survival_B_below"]:
            crossings.append(comparison)
    by_counter = {r["counter"]: r for r in solved}
    representative = compare(by_counter[43042501006],by_counter[43042500083])
    A, B = by_counter[43042501006], by_counter[43042500083]
    aA, aB = 173-int(A["original_row"]["H2"]), 173-int(B["original_row"]["H2"])
    decomposition = {"A_counter": A["counter"], "B_counter": B["counter"],
                     "single_safe_site_counts_A_B": [aA,aB],
                     "minimal_pair_triggers_among_safe_sites_A_B": [comb(aA,2)-int(A["original_row"]["checkpoint_b2_safe_pairs"]),
                                                                     comb(aB,2)-int(B["original_row"]["checkpoint_b2_safe_pairs"])],
                     "rows": []}
    for k in (1,5,10,11,20,40):
        gate_A,gate_B = Fraction(comb(aA,k),comb(173,k)),Fraction(comb(aB,k),comb(173,k))
        coll_A,coll_B = Fraction(A["true_safe_counts"][k],comb(aA,k)),Fraction(B["true_safe_counts"][k],comb(aB,k))
        decomposition["rows"].append({"k": k, "direct_avoidance_A": frac(gate_A), "direct_avoidance_B": frac(gate_B),
                                       "collective_survival_A": frac(coll_A), "collective_survival_B": frac(coll_B),
                                       "direct_ratio_B_over_A": frac(gate_B/gate_A),
                                       "collective_ratio_B_over_A": frac(coll_B/coll_A),
                                       "total_survival_ratio_B_over_A": frac((gate_B*coll_B)/(gate_A*coll_A))})
    table = [{"counter": r["counter"], "H2": int(r["original_row"]["H2"]),
              "b2": int(r["original_row"]["checkpoint_b2_safe_pairs"]),
              "parallel_factors": r["structure"]["parallel_two_port_factors"],
              "core_random_sites": r["structure"]["core_random_sites"],
              "treewidth_upper_bound": r["structure"]["treewidth_upper_bound"],
              "mean": r["mean_true_birth_step"], "S40": r["tail_after_40"],
              "median_step": r["birth_quantiles"]["1/2"],
              "max_safe_k": r["maximum_true_safe_k"], "seconds": r["single_process_seconds"]}
             for r in solved]
    summary = {"selection_manifest_commit": "b9cbe13e", "selected_count": len(records), "solved_count": len(solved),
               "crossing_pairs": len(crossings), "total_pairs": len(solved)*(len(solved)-1)//2,
               "table": table, "all_crossing_certificates": crossings,
               "representative_posthoc_clock_inversion": representative,
               "exact_direct_collective_decomposition": decomposition,
               "representative_choice_status": "descriptive mechanism example selected after all frozen clocks; not a preregistered pair hypothesis",
               "new_samples": 0, "marked_birth_probabilities_evaluated": False}
    (OUTPUT / "scientific_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True)+"\n")
    print("solved",len(solved),"crossings",len(crossings),"of",summary["total_pairs"])
    print({k:v for k,v in representative.items() if "count_B_minus" not in k and "crossproducts" not in k})
    print("safe-site minimal pairs A/B",decomposition["minimal_pair_triggers_among_safe_sites_A_B"])
    for row in decomposition["rows"]:
        print(row["k"],*[row[k]["value"] for k in ("direct_ratio_B_over_A","collective_ratio_B_over_A","total_survival_ratio_B_over_A")])


if __name__ == "__main__":
    main()

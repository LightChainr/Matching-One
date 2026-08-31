#!/usr/bin/env python3
"""Exact full-polynomial grouping of the fixed147 already-solved prefixes."""
from collections import defaultdict
import json
from math import comb
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"results/p334-all147-prefix-clocks/full_clocks.json"
OUTPUT=ROOT/"results/p334-archive-isoclock-groups"


def main():
    source=json.loads(SOURCE.read_text())
    groups=defaultdict(list)
    for row in source["records"]:
        if row["status"]=="solved_full_physical":
            groups[tuple(row["true_safe_counts"])].append(row["counter"])
    classes=[{"representative_counter":min(counters),"members":sorted(counters),
              "size":len(counters),"coefficient_count":len(coefficients)}
             for coefficients,counters in groups.items()]
    classes.sort(key=lambda group:group["representative_counter"])
    nontrivial=[g for g in classes if g["size"]>1]
    result={"parent_commit":"250c589958fc09c52380feb4c99276c8e9c9455b",
            "full_clock_source_commit":"9cca7bc6",
            "source_artifact":str(SOURCE.relative_to(ROOT)),
            "selection":"Exactly the existing fixed147 cohort, same N425 second/k0=252/age10/ell(12,-19); no new rows",
            "equality":"Equality of the entire ordered tuple of all174 exact integer safe coefficients, including trailing zeros; no float curves, approximate tolerance or digest-only comparison",
            "source_rows":len(source["records"]),"completed_rows_grouped":sum(g["size"] for g in classes),
            "exact_isoclock_group_count":len(classes),"singleton_group_count":sum(g["size"]==1 for g in classes),
            "nontrivial_group_count":len(nontrivial),"nontrivial_groups":nontrivial,
            "equal_clock_pairs":sum(comb(g["size"],2) for g in classes),
            "all_groups":classes,
            "new_samples":0,"network_DP_solves":0,"marked_site_DP_solves":0,
            "structural_comparison_status":"not_applicable_no_nontrivial_exact_isoclock_group" if not nontrivial else "candidate_group_available",
            "interpretation":"No real isoclock witness in this fixed147 cohort. This does not make the unmarked clock a generally identifiable graph invariant; the five-site constructed counterexample remains separate from production-prefix evidence."}
    OUTPUT.mkdir(parents=True,exist_ok=True)
    (OUTPUT/"exact_groups.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("rows",result["completed_rows_grouped"],"groups",len(classes),"nontrivial",len(nontrivial),"equal pairs",result["equal_clock_pairs"])


if __name__=="__main__":
    main()

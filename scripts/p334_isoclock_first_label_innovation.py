#!/usr/bin/env python3
"""First-label Doob information in two fixed five-site isoclock mechanisms."""
from fractions import Fraction
import hashlib
from itertools import permutations
import json
from pathlib import Path

import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
CONTRACT=ROOT/"analysis/p334_isoclock_first_label_innovation.json"
OUT=ROOT/"results/p334-isoclock-first-label-innovation"


def value(x):
    return {"exact":str(x),"value":float(x)}


def score(name,graph,vertices):
    d=len(vertices)
    times={v:[] for v in vertices}
    for order in permutations(vertices):
        occupied=set()
        for k,v in enumerate(order,1):
            occupied.add(v)
            if any(a in occupied and b in occupied for a,b in graph["edges"]):
                times[order[0]].append(k)
                break
        else:
            raise ValueError("The prescribed graph unexpectedly has no absorbing edge")
    all_times=[t for row in times.values() for t in row]
    denominator=len(all_times)
    mean=Fraction(sum(all_times),denominator)
    second=Fraction(sum(t*t for t in all_times),denominator)
    variance=second-mean*mean
    conditional_means=[]
    survival_rows=[]
    rows=[]
    for v in vertices:
        ts=times[v]
        m=Fraction(sum(ts),len(ts))
        survival=[Fraction(sum(t>k for t in ts),len(ts)) for k in range(d+1)]
        local_var=Fraction(sum(t*t for t in ts),len(ts))-m*m
        conditional_means.append(m);survival_rows.append(survival)
        rows.append({"first_label":v,"role":graph["roles"][v],
                     "conditional_order_count":len(ts),
                     "conditional_counts_by_T":[sum(t==k for t in ts) for k in range(d+1)],
                     "conditional_mean_T":value(m),"conditional_mean_child_wait":value(m-1),
                     "conditional_variance_T":value(local_var),
                     "conditional_survival_vector":[str(s) for s in survival]})
    mean_survival=[sum(row[k] for row in survival_rows)/d for k in range(d+1)]
    gram=[[sum((row[k]-mean_survival[k])*(row[l]-mean_survival[l]) for row in survival_rows)/d
           for l in range(d+1)] for k in range(d+1)]
    active=[k for k in range(d+1) if gram[k][k]]
    matrix=sp.Matrix([[sp.Rational(gram[k][l].numerator,gram[k][l].denominator) for l in active] for k in active])
    eigenvalues=matrix.eigenvals()
    first_variance=sum((m-mean)**2 for m in conditional_means)/d
    return {"name":name,"edges":graph["edges"],"orders_enumerated":denominator,
            "first_labels":rows,"mean_T":value(mean),"variance_T":value(variance),
            "first_label_innovation_variance":value(first_variance),
            "first_label_fraction_of_clock_variance":value(first_variance/variance),
            "mean_remaining_variance_after_first_label":value(sum(Fraction(row["conditional_variance_T"]["exact"]) for row in rows)/d),
            "direct_absorbing_count_at_empty_prefix":0,"binary_direct_safe_Doob_bound":0,
            "mean_survival_vector":[str(s) for s in mean_survival],
            "first_label_survival_covariance_Gram":[[str(x) for x in row] for row in gram],
            "nonconstant_survival_cuts":active,
            "active_Gram":[[str(matrix[i,j]) for j in range(len(active))] for i in range(len(active))],
            "active_Gram_rank":matrix.rank(),"active_Gram_determinant":str(matrix.det()),
            "active_Gram_characteristic_polynomial":str(matrix.charpoly().as_expr()),
            "nonzero_Gram_eigenvalues":[{"exact":str(x),"value":float(x),"multiplicity":multiplicity}
                                         for x,multiplicity in sorted(eigenvalues.items(),key=lambda pair:float(pair[0])) if x!=0],
            "full_Gram_zero_eigenvalue_multiplicity":d+1-matrix.rank(),
            "variance_from_sum_of_Gram_entries":value(sum(sum(row) for row in gram))}


def main():
    contract=json.loads(CONTRACT.read_text())
    records=[score(name,graph,contract["vertices"]) for name,graph in contract["graphs"].items()]
    result={"schema":contract["schema"],"graph_source_commit":contract["graph_source_commit"],
            "records":records,
            "exact_difference_in_first_label_variance_fraction":value(Fraction(records[0]["first_label_fraction_of_clock_variance"]["exact"])-Fraction(records[1]["first_label_fraction_of_clock_variance"]["exact"])),
            "interpretation":"Constructed exact finite pair-trigger example, not a new N425 prefix or a proof of arbitrary graph identification",
            "input_sha256":{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (CONTRACT,Path(__file__))}}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/"exact_first_label.json").write_text(json.dumps(result,indent=2)+"\n")
    for row in records:
        print(row["name"],"means",[v["conditional_mean_T"]["exact"] for v in row["first_labels"]],
              "innovation",row["first_label_innovation_variance"],"fraction",row["first_label_fraction_of_clock_variance"])
        print("Gram",row["active_Gram"],"rank",row["active_Gram_rank"],"eigenvalues",row["nonzero_Gram_eigenvalues"])


if __name__=="__main__":
    main()

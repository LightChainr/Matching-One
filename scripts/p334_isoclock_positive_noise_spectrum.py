#!/usr/bin/env python3
"""Exact Walsh degree energies of the two specified five-site safe indicators.

The already-saved32-entry truth table is transformed once per graph. No new
graph/source samples, network solves, or correlated-pair enumeration.
"""
from fractions import Fraction
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"results/p334-isoclock-marked-counterexample/exact_marked_clock.json"
OUTPUT=ROOT/"results/p334-isoclock-positive-noise-spectrum"


def frac(value):
    value=Fraction(value)
    return {"exact":str(value),"value":float(value)}


def transform(values):
    output=list(values)
    width=1
    while width<len(output):
        for base in range(0,len(output),2*width):
            for j in range(width):
                a,b=output[base+j],output[base+j+width]
                output[base+j],output[base+j+width]=a+b,a-b
        width*=2
    return output


def main():
    source=json.loads(SOURCE.read_text())
    rho,p=Fraction(1,2),Fraction(1,2)
    u=p+rho*(1-p)
    a=p/u
    records=[]
    for graph in source["records"]:
        safe=set(graph["safe_subset_masks"])
        values=[int(mask in safe) for mask in range(32)]
        numerators=transform(values)
        energy=[sum(Fraction(n*n,1024) for mask,n in enumerate(numerators) if mask.bit_count()==k) for k in range(6)]
        mean=Fraction(sum(values),32)
        variance=mean*(1-mean)
        mask_variance=sum(rho**k*energy[k] for k in range(1,6))
        row={"name":graph["name"],"edges":graph["edges"],
             "truth_table_size":32,"safe_mean":frac(mean),"safe_variance":frac(variance),
             "walsh_basis":"chi_S(x)=(-1)^sum_{i in S}x_i; p=1/2 orthonormal product basis",
             "walsh_numerators_over32":numerators,
             "singleton_walsh_coefficients":[frac(Fraction(numerators[1<<v],32)) for v in range(5)],
             "pivotal_probabilities_at_p_half":[frac(Fraction(numerators[1<<v],16)) for v in range(5)],
             "degree_energies_including_E0":[frac(e) for e in energy],
             "degree_energy_numerators_over256":[int(256*e) for e in energy],
             "nonconstant_total_energy":frac(sum(energy[1:])),
             "same_mask_conditional_response_variance_rho_half":frac(mask_variance),
             "same_mask_two_replica_safe_product_expectation":frac(mean*mean+mask_variance)}
        records.append(row)
        print(row["name"],"energies over256",row["degree_energy_numerators_over256"],
              "totalvar",variance,"variance rhohalf",mask_variance)
    difference=[Fraction(x["exact"])-Fraction(y["exact"]) for x,y in zip(records[0]["degree_energies_including_E0"],records[1]["degree_energies_including_E0"])]
    result={"parent_commit":"a11d649991340bb77af850bb59713a86ffe90c9b",
            "truth_table_source_commit":"250c589958fc09c52380feb4c99276c8e9c9455b",
            "source_artifact":str(SOURCE.relative_to(ROOT)),
            "new_samples":0,"new_N425_samples":0,"network_DP_solves":0,
            "selection":"Only the prescribed double-star and C4+isolate safe indicators; no graph-family scan",
            "fixed_p":frac(p),"selected_rho":frac(rho),"mask_active_probability_a":frac(a),"replica_proposal_probability_u":frac(u),
            "records":records,"degree_energy_difference_double_star_minus_cycle":[frac(e) for e in difference],
            "variance_difference_at_rho_half":frac(sum(rho**k*difference[k] for k in range(1,6))),
            "exact_variance_difference_for_fixed_p_half":"rho*(1-rho)^2/64",
            "strict_variance_order":"double_star > C4_plus_isolate for every 0<rho<1; equal at rho0 and rho1",
            "scope":"Exact constructed finite pair-trigger example. It is not an observed N425 spectrum or the P437 N112 C3 observer; only the positive noise-degree identity is shared."}
    OUTPUT.mkdir(parents=True,exist_ok=True)
    (OUTPUT/"exact_positive_spectra.json").write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    print("energy difference",[str(e) for e in difference],"variance gap",result["variance_difference_at_rho_half"])


if __name__=="__main__":
    main()

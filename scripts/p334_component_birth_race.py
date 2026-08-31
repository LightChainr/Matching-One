#!/usr/bin/env python3
"""Export small factor polynomials and resolve exact component birth races."""
from fractions import Fraction
import json
from math import comb
from pathlib import Path
import time

from p334_full_birth_reliability import safety_polynomial
from p334_pair_only_survival import frac, multiply

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results/p334-twelve-prefix-clocks"
OUTPUT = ROOT / "results/p334-component-birth-race"


def divide_one_plus_z(coefficients, times):
    p = list(coefficients)
    for _ in range(times):
        q = [0] * (len(p)-1)
        q[-1] = p[-1]
        for j in range(len(q)-1,0,-1):
            q[j-1] = p[j]-q[j]
        if p[0] != q[0]:
            raise ValueError("free-site factor is not exact")
        p = q
    return p


def main():
    clocks = json.loads((SOURCE / "full_clocks.json").read_text())["records"]
    records = []
    OUTPUT.mkdir(parents=True,exist_ok=True)
    for clock in clocks:
        counter = clock["counter"]
        mapping = json.loads((SOURCE / "maps" / f"{counter}.json").read_text())
        components = [c for c in mapping["port_components"] if "two_terminal_network" in c]
        started = time.monotonic()
        previous_path = OUTPUT / f"{counter}.json"
        previous = json.loads(previous_path.read_text()) if previous_path.exists() else None
        factors = []
        for index,c in enumerate(components):
            network = c["two_terminal_network"]
            sites = network["vacant_sites"]
            if previous is not None:
                coefficients = previous["factors"][index]["safe_coefficients"]
                method = previous["factors"][index]["coefficient_source"]
            elif len(components) == 1:
                coefficients = divide_one_plus_z(clock["true_safe_counts"],173-len(sites))
                method = "exact division of saved whole polynomial by free-site factor"
            else:
                coefficients,_ = safety_polynomial(network,sites)
                method = "new export of previously computed but unretained small-factor polynomial"
            n = len(sites)
            boundary = [(n-k)*coefficients[k]-(k+1)*coefficients[k+1] for k in range(n)]
            if any(x<0 for x in boundary):
                raise ValueError("negative local monotone boundary coefficient")
            factors.append({"component": index, "site_labels": sites, "n_sites": n,
                            "addresses": c["addresses"], "H2": n-coefficients[1],
                            "safe_coefficients": coefficients, "boundary_coefficients": boundary,
                            "coefficient_source": method})
        prefix = [[1]]
        for f in factors:
            prefix.append(multiply(prefix[-1],f["safe_coefficients"]))
        suffix = [[1] for _ in range(len(factors)+1)]
        for i in range(len(factors)-1,-1,-1):
            suffix[i] = multiply(factors[i]["safe_coefficients"],suffix[i+1])
        free = 173-sum(f["n_sites"] for f in factors)
        free_poly = [comb(free,k) for k in range(free+1)]
        total_boundary = [0]*173
        direct_unit_joint = [Fraction(clock["true_safe_counts"][k],173*comb(172,k)) for k in range(173)]
        direct_unit = sum(direct_unit_joint)
        direct_unit_tail = sum(direct_unit_joint[40:])
        tail40 = Fraction(clock["tail_after_40"]["exact"])
        for i,f in enumerate(factors):
            other = multiply(multiply(prefix[i],suffix[i+1]),free_poly)
            counts = multiply(f["boundary_coefficients"],other)
            counts += [0]*(173-len(counts))
            total_boundary = [a+b for a,b in zip(total_boundary,counts)]
            joint = [Fraction(v,173*comb(172,k)) for k,v in enumerate(counts)]
            probability = sum(joint)
            f.update({"full_boundary_coefficients": counts, "winning_probability": frac(probability),
                      "mean_step_given_winning": frac(sum((k+1)*v for k,v in enumerate(joint))/probability) if probability else None,
                      "winning_and_T_gt40": frac(sum(joint[40:])),
                      "winning_share_given_T_gt40": frac(sum(joint[40:])/Fraction(clock["tail_after_40"]["exact"])),
                      "winning_share_given_T_gt20": frac(sum(joint[20:])/Fraction(clock["tail_after_20"]["exact"])),
                      "direct_singleton_winning_probability": frac(f["H2"]*direct_unit),
                      "collective_winning_probability": frac(probability-f["H2"]*direct_unit),
                      "direct_singleton_share_given_T_gt40": frac(f["H2"]*direct_unit_tail/tail40),
                      "collective_share_given_T_gt40": frac((sum(joint[40:])-f["H2"]*direct_unit_tail)/tail40)})
        global_boundary = [(173-k)*clock["true_safe_counts"][k]-(k+1)*clock["true_safe_counts"][k+1] for k in range(173)]
        if total_boundary != global_boundary:
            raise ValueError("component boundary sum does not equal full physical boundary")
        row = {"counter": counter, "original_H2": int(clock["original_row"]["H2"]),
               "mean_true_birth_step": clock["mean_true_birth_step"], "factors": factors,
               "free_sites": free, "component_winning_sum": frac(sum(Fraction(f["winning_probability"]["exact"]) for f in factors)),
               "component_sum_equals_full_clock_each_step": True,
               "original_direct_gate_winning_probability": frac(int(clock["original_row"]["H2"])*direct_unit),
               "original_direct_gate_share_given_T_gt40": frac(int(clock["original_row"]["H2"])*direct_unit_tail/tail40),
               "original_factor_export_seconds": previous["factor_export_and_score_seconds"] if previous else None,
               "factor_coefficients_reused_this_pass": previous is not None,
               "factor_export_and_score_seconds": time.monotonic()-started}
        records.append(row)
        (OUTPUT / f"{counter}.json").write_text(json.dumps(row,indent=2,sort_keys=True)+"\n")
        print(counter,[(f["component"],f["n_sites"],f["H2"],round(f["winning_probability"]["value"],6),round(f["winning_share_given_T_gt40"]["value"],6)) for f in factors],flush=True)
    summary = {"parent_commit":"bd95f2a048d5780568b689bd42e0a684daf74315", "new_samples":0,
               "source":"same frozen twelve real prefixes; small factor coefficients newly exported, not claimed as zero-DP archival readout",
               "records":[{**{k:v for k,v in r.items() if k!='factors'},"factors":[{k:v for k,v in f.items() if k not in ('safe_coefficients','boundary_coefficients','full_boundary_coefficients','site_labels')} for f in r['factors']]} for r in records]}
    (OUTPUT / "summary.json").write_text(json.dumps(summary,indent=2,sort_keys=True)+"\n")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Exact final-birth site probabilities from the two physical reliability maps."""
from fractions import Fraction
import json
from math import comb
from pathlib import Path
import time

from p334_full_birth_reliability import safety_polynomial
from p334_pair_only_survival import frac

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results/p334-contracted-full-clock"
OUTPUT = ROOT / "results/p334-exact-marked-birth"


def main():
    mapped = json.loads((SOURCE / "whole_event_networks.json").read_text())["records"]
    clocks = {r["counter"]: r for r in json.loads((SOURCE / "full_physical_birth_clock.json").read_text())["records"]}
    records = []
    for mapping in mapped:
        counter, sites = mapping["counter"], mapping["vacant_sites"]
        component = next(c for c in mapping["port_components"] if "two_terminal_network" in c)
        network = component["two_terminal_network"]
        core = set(network["vacant_sites"])
        port_class = {v: f"port_{component['addresses'].index(a)}" for v, a in component["ports"]}
        d = len(sites)
        denominators = [d * comb(d - 1, k) for k in range(d)]
        started = time.monotonic()
        site_records = []
        total_by_step = [Fraction(0) for _ in range(d)]
        for index, site in enumerate(sites):
            if site in core:
                others = set(sites) - {site}
                off, _ = safety_polynomial(network, others, forced_site=site, forced_value=False)
                on, _ = safety_polynomial(network, others, forced_site=site, forced_value=True)
                pivotal = [a - b for a,b in zip(off, on)]
                if any(value < 0 for value in pivotal):
                    raise ValueError("nonmonotone physical pivotal coefficient")
            else:
                pivotal = [0] * d
            joint = [Fraction(value, denominators[k]) for k,value in enumerate(pivotal)]
            total_by_step = [a + b for a,b in zip(total_by_step, joint)]
            probability = sum(joint)
            time_mass = sum((k + 1) * value for k,value in enumerate(joint))
            site_records.append({"site": site, "type": port_class.get(site, "interior" if site in core else "outside_core"),
                                 "pivotal_count_by_prior_size": pivotal,
                                 "birth_probability": frac(probability),
                                 "mean_birth_step_conditional_on_site": frac(time_mass / probability) if probability else None,
                                 "birth_after_40_joint_probability": frac(sum(joint[40:])),
                                 "birth_after_65_joint_probability": frac(sum(joint[65:]))})
            if index % 30 == 0:
                print(counter, "sites", index + 1, "/", d, "seconds", time.monotonic() - started, flush=True)
        clock = clocks[counter]
        target = [Fraction(clock["true_survival"][k]["exact"]) - Fraction(clock["true_survival"][k + 1]["exact"]) for k in range(d)]
        if total_by_step != target:
            raise ValueError("site-resolved joint law does not sum to the physical clock")
        probability_by_type = {}
        for key in ("port_0", "port_1", "interior", "outside_core"):
            group = [s for s in site_records if s["type"] == key]
            probability_by_type[key] = {"sites": [s["site"] for s in group],
                                        "birth_probability": frac(sum(Fraction(s["birth_probability"]["exact"]) for s in group)),
                                        "after_40_share": frac(sum(Fraction(s["birth_after_40_joint_probability"]["exact"]) for s in group) / Fraction(clock["true_survival"][40]["exact"]))}
        top = sorted(site_records, key=lambda s: Fraction(s["birth_probability"]["exact"]), reverse=True)
        late = sorted(site_records, key=lambda s: Fraction(s["birth_after_40_joint_probability"]["exact"]), reverse=True)
        summary = {"positive_birth_sites": sum(Fraction(s["birth_probability"]["exact"]) > 0 for s in site_records),
                   "top5_sites": [s["site"] for s in top[:5]],
                   "top5_probability": frac(sum(Fraction(s["birth_probability"]["exact"]) for s in top[:5])),
                   "late_top5_sites": [s["site"] for s in late[:5]],
                   "late_top5_share": frac(sum(Fraction(s["birth_after_40_joint_probability"]["exact"]) for s in late[:5]) / Fraction(clock["true_survival"][40]["exact"])),
                   "inverse_simpson_effective_sites": frac(1 / sum(Fraction(s["birth_probability"]["exact"])**2 for s in site_records)),
                   "probability_by_type": probability_by_type,
                   "exact_total_birth_probability": frac(sum(Fraction(s["birth_probability"]["exact"]) for s in site_records)),
                   "exact_site_sum_matches_clock_every_step": True}
        row = {"counter": counter, "seed": mapping["seed"], "N": mapping["N"], "k0": mapping["k0"],
               "site_records": site_records, "summary": summary, "single_process_seconds": time.monotonic() - started}
        records.append(row)
        OUTPUT.mkdir(parents=True, exist_ok=True)
        (OUTPUT / f"marked_birth_{counter}.json").write_text(json.dumps(row, indent=2, sort_keys=True) + "\n")
        print("DONE", counter, {k:v for k,v in summary.items() if k != "probability_by_type"}, flush=True)
    (OUTPUT / "summary.json").write_text(json.dumps({"parent_commit": "6358ba49ef390c10a3f501b589ba7ba1d4e05b09", "new_samples": 0,
                                                    "records": [{k:v for k,v in r.items() if k != "site_records"} for r in records]}, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()

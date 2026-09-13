#!/usr/bin/env python3
"""Exact inert-site knockout effects read from already saved marked birth laws."""
from fractions import Fraction
from math import comb
from pathlib import Path
import hashlib
import json
import subprocess

from scipy.stats import binom

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/p334-clock-knockout-from-marked-law"
SOURCE = "1c06230b8f7e13be98f128361ad72b23c0c425ae"
COUNTERS = (43042514269, 43042505280)
P_REF = .59274605079


def value(x):
    return {"exact": str(x), "value": float(x)}


def main():
    result = {"source_commit": SOURCE, "new_MC": 0, "new_network_solves": 0,
              "clock": "Retain all d original insertion positions; a blocked site is an inert dummy.",
              "identity": "Delta S_v(k)=k P(T=k,V_final=v); Delta E[T]=E[T 1(V_final=v)] if eventual birth remains possible.",
              "sources": {}}
    report = ["# The complete single-site knockout landscape was already in the marked law", ""]
    for counter in COUNTERS:
        path = f"results/p334-exact-marked-birth/marked_birth_{counter}.json"
        raw = subprocess.check_output(["git", "show", SOURCE+":"+path], cwd=ROOT)
        source = json.loads(raw)
        d = source["N"]-source["k0"]
        thermal_kernel = binom.pmf([source["k0"]+k for k in range(d+1)], source["N"], P_REF)
        rows, exact_losses = [], []
        for row in source["site_records"]:
            counts = row["pivotal_count_by_prior_size"]
            increments = [Fraction(0)] + [Fraction(counts[k-1], comb(d, k)) for k in range(1, d+1)]
            loss = sum(increments, Fraction(0))
            pi = Fraction(row["birth_probability"]["exact"])
            late = sum(increments[41:], Fraction(0))
            rows.append({"site": row["site"], "type": row["type"],
                         "birth_probability": value(pi), "mean_wait_increase": value(loss),
                         "survival40_increase": value(increments[40]),
                         "survival65_increase": value(increments[65]),
                         "canonical_F2_change_at_pref": -float(sum(float(x)*b for x, b in zip(increments, thermal_kernel))),
                         "integrated_F2_change": value(-loss/(source["N"]+1)),
                         "mean_loss_from_birth_after40": value(late),
                         "permanent_nonarrival_probability": value(increments[-1])})
            exact_losses.append(loss)
        total = sum(exact_losses, Fraction(0))
        pi_order = sorted(rows, key=lambda r: (-Fraction(r["birth_probability"]["exact"]), r["site"]))
        loss_order = sorted(rows, key=lambda r: (-Fraction(r["mean_wait_increase"]["exact"]), r["site"]))
        tail_order = sorted(rows, key=lambda r: (-Fraction(r["survival40_increase"]["exact"]), r["site"]))
        thermal_order = sorted(rows, key=lambda r: (r["canonical_F2_change_at_pref"], r["site"]))
        thermal_total = sum(r["canonical_F2_change_at_pref"] for r in rows)
        for row in rows:
            row["fraction_of_total_single_site_mean_impact"] = value(Fraction(row["mean_wait_increase"]["exact"])/total)
        by_type = {}
        for typ in sorted({r["type"] for r in rows}):
            selected = [r for r in rows if r["type"] == typ]
            by_type[typ] = {
                "birth_probability_share": value(sum((Fraction(r["birth_probability"]["exact"]) for r in selected), Fraction(0))),
                "mean_impact_share": value(sum((Fraction(r["mean_wait_increase"]["exact"]) for r in selected), Fraction(0))/total),
                "canonical_F2_loss_share": sum(r["canonical_F2_change_at_pref"] for r in selected)/thermal_total}
        top5pi = [r["site"] for r in pi_order[:5]]
        top5loss = [r["site"] for r in loss_order[:5]]
        summary = {"counter": counter, "d": d, "source_path": path,
                   "source_sha256": hashlib.sha256(raw).hexdigest(), "site_records": rows,
                   "mean_wait_equal_total_individual_knockout_impact": value(total),
                   "all_individual_knockouts_still_eventually_birth": all(r["permanent_nonarrival_probability"]["exact"] == "0" for r in rows),
                   "top5_birth_probability_sites": top5pi, "top5_mean_impact_sites": top5loss,
                   "top5_tail40_impact_sites": [r["site"] for r in tail_order[:5]],
                   "p_ref": P_REF, "top5_canonical_F2_impact_sites": [r["site"] for r in thermal_order[:5]],
                   "sum_of_individual_canonical_F2_changes": thermal_total,
                   "top5_mean_impact_share": value(sum((Fraction(r["mean_wait_increase"]["exact"]) for r in loss_order[:5]), Fraction(0))/total),
                   "inverse_simpson_impact_sites": value(total*total/sum((x*x for x in exact_losses), Fraction(0))),
                   "shares_by_type": by_type}
        result["sources"][str(counter)] = summary
        report += [f"## Prefix {counter}", "",
                   f"Total individual mean extension = original mean wait = {float(total):.10g}.",
                   f"Top five by birth frequency: {top5pi}; by mean impact: {top5loss}; by step40 tail impact: {summary['top5_tail40_impact_sites']}.",
                   f"Effective mean-impact sites: {summary['inverse_simpson_impact_sites']['value']:.9g}.", "",
                   "| site | type | birth probability | mean-wait extension | step40 tail extension | canonical F2 change |", "|---|---|---:|---:|---:|---:|"]
        for row in loss_order[:10]:
            report.append(f"| {row['site']} | {row['type']} | {row['birth_probability']['value']:.8g} | {row['mean_wait_increase']['value']:.8g} | {row['survival40_increase']['value']:.8g} | {row['canonical_F2_change_at_pref']:.8g} |")
        report += ["", f"Top five canonical impact sites: {summary['top5_canonical_F2_impact_sites']}.", "",
                   "Type shares of birth probability / mean impact / canonical F2 loss:", ""]
        for typ, share in by_type.items():
            report.append(f"- {typ}: {share['birth_probability_share']['value']:.8g} / {share['mean_impact_share']['value']:.8g} / {share['canonical_F2_loss_share']:.8g}")
        report.append("")
    OUT.mkdir(exist_ok=True)
    (OUT/"score.json").write_text(json.dumps(result, indent=2)+"\n")
    (OUT/"REPORT.md").write_text("\n".join(report))
    print("\n".join(report))


if __name__ == "__main__":
    main()

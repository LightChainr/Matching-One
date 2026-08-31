#!/usr/bin/env python3
"""Resolve support versus side-capacity effects on already saved trigger graphs."""
from collections import defaultdict
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
TERMS = ("localization", "bipartition_constraint", "residual_organization")


def frac(x):
    return {"exact": str(x), "numerator": x.numerator, "denominator": x.denominator, "decimal": float(x)}


def simple_stars(s, m):
    return F(2*m*(m-1), s+1)


def analyze(record):
    structure, source = record["structure"], record["source_row"]
    if not structure["bipartite"]:
        raise ValueError("This fixed saved census contains a nonbipartite graph")
    a, m = structure["safe_vertices"], structure["trigger_edges"]
    global_null = simple_stars(a, m)
    local_null = bip_null = F(0)
    components = []
    for component in structure["components_with_edges"]:
        L, R = map(len, component["sides"])
        edges, s = component["edges"], L+R
        simple = simple_stars(s, edges)
        bip = F(edges*(edges-1)*(s-2), 2*(L*R-1)) if L*R > 1 else F(0)
        side_formula = F(edges*(edges-1)*((L-R)**2-s+2), 2*(L*R-1)*(s+1)) if L*R > 1 else F(0)
        if bip-simple != side_formula:
            raise ValueError("Closed-form bipartition increment differs")
        local_null += simple
        bip_null += bip
        components.append({"L": L, "R": R, "m": edges, "s": s,
                           "imbalance_squared_minus_s_plus_2": (L-R)**2-s+2,
                           "simple_expected_stars": frac(simple), "bipartite_expected_stars": frac(bip),
                           "side_constraint_increment": frac(side_formula)})
    observed = F(structure["trigger_wedges"])
    parts = (local_null-global_null, bip_null-local_null, observed-bip_null)
    if sum(parts) != observed-global_null:
        raise ValueError("Exact additive allocation failed")
    d = int(source["n"])-int(source["k0"])
    factor = F(2, d*(d-1)**2)
    return {"selection": record["selection"], "N": int(source["n"]), "orientation": source["orientation"],
            "replica": int(source["replica"]), "a": a, "d": d, "m": m,
            "support_vertices": a-structure["isolated_vertices"], "components": components,
            "observed_W2": frac(observed), "global_null_W2": frac(global_null),
            "block_simple_null_W2": frac(local_null), "block_bipartite_null_W2": frac(bip_null),
            "excess_over_global_W2": frac(observed-global_null),
            "allocation_W2": {k: frac(v) for k, v in zip(TERMS, parts)},
            "allocation_cooperative_excess": {k: frac(factor*v) for k, v in zip(TERMS, parts)},
            "W2_to_cooperative_factor": frac(factor), "exact_addition_residual": "0"}


def value(row, key):
    return F(row[key]["exact"])


def main():
    manifest = json.loads((ROOT/"analysis/p334_trigger_capacity_allocation_manifest.json").read_text())
    blob = subprocess.check_output(["git", "cat-file", "blob", f"{manifest['source_commit']}:{manifest['source_path']}"], cwd=ROOT)
    if hashlib.sha256(blob).hexdigest() != manifest["source_sha256"]:
        raise ValueError("Saved graph census differs from pinned input")
    source = json.loads(blob)
    rows = [analyze(r) for r in source["records"]]
    witnesses = {r["selection"]: r for r in rows if r["selection"].startswith("existing_")}
    A, B = witnesses["existing_N425_A"], witnesses["existing_N425_B"]
    pair_parts = {k: F(B["allocation_W2"][k]["exact"])-F(A["allocation_W2"][k]["exact"]) for k in TERMS}
    observed_difference = value(B, "observed_W2")-value(A, "observed_W2")
    if sum(pair_parts.values()) != observed_difference:
        raise ValueError("Saved-pair contrast does not close")
    groups = defaultdict(list)
    for row in rows:
        if row["selection"].startswith("first5_"):
            groups[f"N{row['N']}_{row['orientation']}"].append(row)
    summaries = {}
    for name, selected in groups.items():
        total = sum(value(r, "excess_over_global_W2") for r in selected)
        sums = {k: sum(F(r["allocation_W2"][k]["exact"]) for r in selected) for k in TERMS}
        summaries[name] = {"graphs": len(selected), "excess_over_global_W2_sum": frac(total),
                           "allocation_W2_sums": {k: frac(v) for k, v in sums.items()},
                           "positive_bipartition_increment_graphs": sum(F(r["allocation_W2"]["bipartition_constraint"]["exact"]) > 0 for r in selected),
                           "positive_residual_graphs": sum(F(r["allocation_W2"]["residual_organization"]["exact"]) > 0 for r in selected)}
    result = {"schema": "matching-one.p334-trigger-capacity-allocation.v1", "provenance": manifest,
              "environment": {"python": sys.version, "machine": platform.machine()}, "new_samples": 0,
              "exact_identities": {"simple_W2": "2m(m-1)/(s+1)", "bipartite_W2": "m(m-1)(s-2)/(2(LR-1))", "bipartite_minus_simple": "m(m-1)((L-R)^2-s+2)/(2(LR-1)(s+1))", "cooperative_increment": "2 DeltaW2/[d(d-1)^2]"},
              "rows": rows, "first_five_fixed_census": summaries,
              "saved_N425_B_minus_A": {"observed_W2_difference": frac(observed_difference),
                                        "allocation_W2": {k: frac(v) for k, v in pair_parts.items()},
                                        "allocation_cooperative_difference": {k: frac(F(A["W2_to_cooperative_factor"]["exact"])*v) for k, v in pair_parts.items()},
                                        "exact_addition_residual": "0"},
              "uncertainty": "Exact arithmetic on fixed selected graphs, not estimated population effects. No new sampling covariance, p-values or prevalence estimates are asserted."}
    lines = ["# P334: side imbalance, not bipartiteness alone, explains the selected overlap contrast", "",
             "**The saved N425 pair separates two effects that the previous 84% capacity figure combined.** Expanding the active support from 26 to 34 vertices lowers the fixed-edge simple-graph expected two-star count; the much more unequal 5×29 bipartition raises it strongly enough to reverse that change. Both real graphs retain overlap beyond their capacity benchmarks.", "",
             "This is exact zero-new-sample arithmetic on the 22 already saved trigger graphs. It does not repeat graph replay, bipartite/Ferrers screening, the full-production overlap scorer, or the minimal-triple census.", "",
             "## Three additive structural comparisons", "",
             "Edges are minimal rank-two-triggering pairs on individually safe sites; W2 counts shared-endpoint edge pairs. Keep a safe sites and m trigger edges. Compare: (i) G(a,m); (ii) uniform simple edges independently within each observed nonisolated component's vertex block with its observed edge count; (iii) uniform cross edges within that same block's observed L/R sides; (iv) the observed graph.", "",
             "The last two benchmarks may disconnect a block or introduce new isolates. They do not preserve connectedness or the complete component decomposition. The first increment is support/block/edge-allocation localization, not purely support size or a causal physical intervention.", "",
             "For one block with s=L+R, exact expectations are `E_simple W2=2m(m-1)/(s+1)` and `E_bip W2=m(m-1)(s-2)/[2(LR-1)]`. Hence", "",
             "```text", "E_bip W2 - E_simple W2", "  = m(m-1) [(L-R)^2-s+2] / [2(LR-1)(s+1)].", "```", "",
             "For m≥2, bipartition constraints increase expected overlap only when `(L-R)^2 > s-2`. Near-balanced bipartite capacities can reduce overlap. The exceptional one-slot case has both expectations zero.", "",
             "## Matched real N425 witnesses", "",
             "A/B remain the previously selected counters 43042514269 / 43042505280, with identical a=d=173, m=108, k0=252, age=10 and ell=(12,-19). Their single nonisolated blocks are 14×12 and 5×29.", "",
             "| W2 quantity | A | B | B−A |", "|---|---:|---:|---:|"]
    for key in ("global_null_W2", "block_simple_null_W2", "block_bipartite_null_W2", "observed_W2"):
        va, vb = value(A,key), value(B,key)
        lines.append(f"| {key} | {float(va):.6f} | {float(vb):.6f} | {float(vb-va):+.6f} |")
    lines += ["", "| Contribution to observed B−A = 540 | exact | decimal |", "|---|---:|---:|"]
    for name, v in pair_parts.items():
        lines.append(f"| {name} | {v} | {float(v):+.6f} |")
    lines += ["", "Thus the already reported +453.628743 capacity contribution is **−195.657143 localization +649.285885 side constraint**. Residual organization contributes +86.371257. This refines, rather than repeats, the earlier 84.0053% arithmetic: two-sided imbalance overcompensates an opposing support effect. A's side constraint lowers its expected W2 by 25.628743; B's raises it by 623.657143.", "",
              "Multiplying each contribution by `2/[173·172²]` gives an exact decomposition of the saved double-clone probability difference `135/639754`; the JSON stores these rational contributions. There is no new independent branching evidence.", "",
              "## Existing first-five extension: fixed-set descriptive check", "",
              "The two specially chosen witnesses are excluded from this table. Each row uses the already fixed five lowest eligible counters in that environment. Entries are sums of W2 increments across those five saved graphs, not population means, prevalence estimates or a new sampling test.", "",
              "| fixed environment | support/block localization | side constraint | residual organization | positive side increments / 5 | positive residuals / 5 |", "|---|---:|---:|---:|---:|---:|"]
    for name, summary in summaries.items():
        vals = [summary["allocation_W2_sums"][k]["decimal"] for k in TERMS]
        lines.append(f"| {name} | {vals[0]:+.4f} | {vals[1]:+.4f} | {vals[2]:+.4f} | {summary['positive_bipartition_increment_graphs']}/5 | {summary['positive_residual_graphs']}/5 |")
    lines += ["", f"Only {sum(s['positive_residual_graphs'] for s in summaries.values())}/20 of these fixed extension graphs have strictly positive residual W2 beyond the bipartite-capacity expectation. The N325-second residual sum is negative, whereas the other three are positive. This does not estimate population prevalence; it shows why the full-production excess above G(a,m) cannot be interpreted as excess beyond every geometry-aware benchmark. The selected A/B positive residuals are not representative evidence for that stronger statement."]
    lines += ["", "## Scientific consequence and next output", "",
              "The microscopic target is now sharper than 'a bipartite graph' or 'more hidden memory': explain the physical two-sided capacity and the residual nonexchangeable degree organization. Ferrers nesting is already false; none of this proves universal bipartiteness. The separately measured genuine triple-trigger counts 583/509 remain needed for three-step survival and are not supplied by this pair-graph benchmark.", "",
              "**Next actual output:** use the saved graph/site labels and period matrix to identify the two sides with explicit topological cut or boundary-landing classes on the existing matched configurations, and predict their capacities before consulting degrees. That would turn the presently observed L/R partition into a physical explanation. Another graph census, generic third-clone moment, or repeated support-baseline score is not required to reach this question.", "",
              f"Input: `{manifest['source_commit']}` / `{manifest['source_path']}`, SHA256 `{manifest['source_sha256']}`. Background full-production overlap source: `{manifest['background_overlap_commit']}`. All rows remain in the original P334 N325/N425 dependency groups.", "",
              "All additions and benchmark identities use exact rational arithmetic. No graph simulations, topology runner, parent scorer, or test suite were run. One script produces both outputs:", "", "```bash", "/Users/lc/python-envs/research-py311/bin/python scripts/p334_trigger_capacity_allocation.py", "```", ""]
    for kind, text in (("json", json.dumps(result, indent=2)+"\n"), ("markdown", "\n".join(lines))):
        path = ROOT/manifest["outputs"][kind]
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
    print(json.dumps({"saved_pair": result["saved_N425_B_minus_A"], "first_five": summaries}, indent=2))


if __name__ == "__main__":
    main()

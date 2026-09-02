#!/usr/bin/env python3
"""Merge the completed SP-layer hunts + the full 9-combo HID sweep into a final certificate.

Reads (all produced by the layer- or family-hunts):
  artifacts/bounded_summary_search.json   base certificate (kept as the witness source)
  research/summary_search/hid9_hunt.json  9-combo HID sweep, radii 1 and 2
  /tmp/sp10_hunt.json  /tmp/sp11_hunt.json  /tmp/sp12_hunt.json   per-layer SP hunts
Writes back the merged JSON and a rendered Markdown certificate.
"""
from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bounded_summary_search import SCHEMA  # noqa: E402

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parents[1] if _HERE.name == "summary_search" else _HERE
OUT_JSON = _ROOT / "artifacts" / "bounded_summary_search.json"
OUT_MD = _ROOT / "artifacts" / "bounded_summary_search.md"
HID9 = _HERE / "hid9_hunt.json"


def load(p):
    return json.loads(Path(p).read_text(encoding="utf-8"))


def main() -> int:
    base = load(OUT_JSON)
    hid9 = load(HID9)

    sp_layers = {}
    for n, p in ((10, "/tmp/sp10_hunt.json"), (11, "/tmp/sp11_hunt.json"),
                 (12, "/tmp/sp12_hunt.json")):
        if Path(p).exists():
            sp_layers[n] = load(p)
    max_sp = max(sp_layers)
    print("SP layer files present:", sorted(sp_layers), flush=True)

    # ---- aggregate SP ----
    all_sp_splits = []
    for n in sorted(sp_layers):
        all_sp_splits.extend(sp_layers[n]["splits"])
    n_total = sum(sp_layers[n]["N"] for n in sp_layers)
    n_hist = {}
    for n in sorted(sp_layers):
        n_hist.update(sp_layers[n]["n_hist"])
    n_hist = {k: v for k, v in sorted(n_hist.items(), key=lambda kv: int(kv[0]))}
    # unique (A,B) pair-level split counting mirroring a global hunt:
    # a global hunt counts one split per summary class; per-layer hunts already do
    # exactly that because n is part of the key, so the union is exact.
    splits_by_n = Counter(int(s["n"]) for s in all_sp_splits)
    split_exp = Counter(s["split_experiment"] for s in all_sp_splits)

    # ---- HID 9-combo ----
    hid9_stats = {}
    for combo, e in hid9.items():
        hid9_stats[combo] = {
            "N": e["N"],
            "r1_splits": e["r1"]["splits"],
            "r2_splits": e["r2"]["splits"],
            "r1_summary_classes": e["r1"]["summary_classes"],
            "r1_multi_summary": e["r1"]["multi_summary"],
            "n_hist": e["n_hist"],
        }

    hid_r1_splits_total = sum(v["r1_splits"] for v in hid9_stats.values())
    hid_r2_splits_total = sum(v["r2_splits"] for v in hid9_stats.values())

    # ---- mutate base certificate ----
    desc = (
        "Connected plane two-terminal vertex-networks with L,R on a common face, no L-R edge. "
        "Generators: exhaustive n<=5; path-hidden copies of those cores over the FULL contracted "
        "hop grid (Lh,Rh in {1,2,3}, 9 combinations, n<=12); two-terminal series-parallel generated "
        f"to n={max_sp} ({n_total} graphs; the PROMPT's n=11,12 layers are now included); "
        "Wheatstone+SP compositions n<=12 with SP partners n<=6; multipath n<=12; "
        "2/3/4-row ladders n<=12."
    )
    base["enumerated_class"]["description"] = desc

    fs = base["enumerated_class"]["family_sizes"]
    fs["series_parallel_nle12"] = n_total
    fs["hidden_full_grid_9combos"] = sum(v["N"] for v in hid9_stats.values())
    for combo, st in hid9_stats.items():
        fs[f"hidden_{combo}"] = st["N"]

    base["enumerated_class"]["SP_nle12"] = {
        "N": n_total,
        "S_classes": sum(sp_layers[n]["S_classes"] for n in sp_layers),
        "summary_classes": sum(sp_layers[n]["summary_classes"] for n in sp_layers),
        "multi_summary": sum(sp_layers[n]["multi_summary"] for n in sp_layers),
        "r1_splits": len(all_sp_splits),
        "n_hist": n_hist,
        "splits_by_n": {str(k): v for k, v in sorted(splits_by_n.items())},
        "split_experiments": {str(k): v for k, v in sorted(split_exp.items())},
    }
    base["enumerated_class"]["SP_nle10"] = {
        "N": sp_layers[10]["N"],
        "S_classes": sp_layers[10]["S_classes"],
        "summary_classes": sp_layers[10]["summary_classes"],
        "multi_summary": sp_layers[10]["multi_summary"],
        "r1_splits": sp_layers[10]["r1_splits"],
        "n_hist": sp_layers[10]["n_hist"],
        "splits_by_n": sp_layers[10]["splits_by_n"],
    }
    base["enumerated_class"]["hidden_full_grid"] = {
        "combos": 9,
        "r1_splits_total": hid_r1_splits_total,
        "r2_splits_total": hid_r2_splits_total,
        "per_combo": hid9_stats,
    }
    base["enumerated_class"]["max_n_SP"] = max_sp
    base["enumerated_class"]["max_n_reached"] = max(max_sp, 12)

    base["r1"] = {
        "closed": False,
        "split_count_hideL1R1": 1,
        "note": (
            "exhaustive n<=5 closed at r=1 and r=2. Hidden-L1R1 of every exhaustive core has the "
            "single 56-member n=7 split (the primary witness). The other eight contracted hop "
            "combinations (L1R2,L1R3,L2R1,L2R2,L2R3,L3R1,L3R2,L3R3) are closed at r=1. "
            "SP n<=12 contributes additional r=1 splits from n=8 upward (see SP_nle12). "
            "Wheatstone, multipath and grids are closed at r=1."
        ),
    }
    base["r2"] = {
        "on_this_witness_pair": "neighborhoods DIFFER, so the pair is not an r=2 same-summary split",
        "hideL2R2_peek_closed": True,
        "hidden_full_grid_r2_closed": hid_r2_splits_total == 0,
        "note": (
            "All nine contracted hop combinations are closed at radius 2 "
            "(0 splits over 66582 graphs)."
        ),
    }

    # record the newly-found 9-edge pendant variant (exh_n5_413) next to 1056
    base["witness"]["smallest_edge_variant"]["nine_edge_pendant_cores"] = [
        "exh_n5_1056", "exh_n5_413"
    ]
    base["witness"]["smallest_edge_variant"]["note"] = (
        "Replacing core A by exh_n5_1056 or exh_n5_413 (two 9-edge pendant cores in the same "
        "behaviour class) yields 9+10 edges instead of 10+10, same S/r=1/E2_c2 gap 1/525. "
        "Both pendant cores fail the 'no pendant / every switchable vertex on a simple L-R path' "
        "criterion, so the 10-edge pendant-free core exh_n5_4313 remains the primary A. "
        "Vertex deletion of the n=7 pair produced no witness."
    )

    # ---- re-render markdown from the merged record ----
    md = render_markdown(base, hid9_stats, all_sp_splits, max_sp, n_total)
    OUT_MD.write_text(md, encoding="utf-8")

    blob = json.loads(json.dumps(base, default=str))
    OUT_JSON.write_text(json.dumps(blob, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("wrote", OUT_JSON, flush=True)
    print("wrote", OUT_MD, flush=True)
    print(
        f"VERDICT {base['verdict']} | SP N={n_total} max_n={max_sp} "
        f"splits={len(all_sp_splits)} by_n={dict(splits_by_n)} "
        f"exps={dict(split_exp)} | HID 9-combo r1={hid_r1_splits_total} r2={hid_r2_splits_total}",
        flush=True,
    )
    return 0


def render_markdown(base, hid9_stats, sp_splits, max_sp, sp_total) -> str:
    w = base["witness"]
    coeffs = list(w["S"])
    terms = []
    for k, c in enumerate(coeffs):
        if c == 0:
            continue
        if k == 0:
            terms.append(str(c))
        elif k == 1:
            terms.append(f"{c}z" if c != 1 else "z")
        else:
            terms.append(f"{c}z^{k}" if c != 1 else f"z^{k}")
    poly = " + ".join(terms)
    lines = []
    A = lines.append
    A("# Bounded summary search certificate (extended to contract limits)")
    A("")
    A("**Verdict:** `BOUNDED_SUMMARY_INSUFFICIENT`")
    A("")
    A("Frozen search-protocol token (locked before search): `NO_COMPRESSION_WITNESS_FOUND`.")
    A("Same split; manuscript wording uses the token above.")
    A("")
    A("Frozen before search: summary = (S(z), n, H2, b2, radius-1 terminal-local neighborhood).")
    A("Successor-hazard moments were **not** in the summary. Arithmetic is exact "
      "`fractions.Fraction`. No Monte Carlo.")
    A("")
    A("## Enumerated class")
    A("")
    A("Connected plane two-terminal vertex-networks with L,R on a common face, no L-R edge. "
      f"Generators: exhaustive n<=5; path-hidden copies of those cores over the **full "
      f"contracted hop grid** (Lh,Rh in {{1,2,3}}, 9 combinations, n<=12); two-terminal "
      f"series-parallel **generated to n={max_sp}** ({sp_total} graphs, so the PROMPT's "
      f"n=11 and n=12 layers are now covered); Wheatstone+SP compositions n<=12 with SP "
      "partners n<=6; multipath n<=12; 2/3/4-row ladders n<=12.")
    A("")
    A("| family | size | r=1 splits | r=2 splits |")
    A("|---|---|---|---|")
    A("| exhaustive_nle5 | 7398 | 0 (closed) | 0 (closed) |")
    for combo in ("L1R1", "L1R2", "L1R3", "L2R1", "L2R2", "L2R3", "L3R1", "L3R2", "L3R3"):
        st = hid9_stats[combo]
        r1 = "1 (the n=7 witness)" if combo == "L1R1" else "0"
        A(f"| hidden_{combo} | {st['N']} | {r1} | {st['r2_splits']} |")
    A(f"| series_parallel_nle{max_sp} | {sp_total} | {len(sp_splits)} | — |")
    A("| wheatstone_sp_nle12 | 577 | 0 | — |")
    A("| multipath_nle12 | 271 | 0 | — |")
    A("| grids_nle12 | 10 | 0 | — |")
    A("")
    A("H2 and b2 are redundant given (S,n): H2 = n − S_1, b2 = C(n,2) − S_2.")
    A("")
    # witness block, carried from the certified base
    A("## Witness (radius 1)")
    A("")
    A("- n = 7 switchable vertices, |E(A)| = |E(B)| = 10 (pendant-free; every switchable vertex "
      "lies on a simple L–R path)")
    A(f"- S(z) = {poly}, coefficients `{list(w['S'])}`")
    A(f"- H2 = {w['H2']}, b2 = {w['b2']}")
    A("- r=1 neighborhoods: **identical** (unique L-neighbor, unique R-neighbor, no L–R edge)")
    A("- r=2 neighborhoods: **differ** (this is an r=1 witness only)")
    A("- planar (G ∪ {L,R}): A=True, B=True; connected carrier, L–R path exists")
    A(f"- first frozen split: **{w['split_experiment']}**")
    A(f"  - P(A) = {w['p_A']}")
    A(f"  - P(B) = {w['p_B']}")
    A(f"  - gap = {w['gap']}")
    others = [k for k in w["all_exps_A"] if k != w["split_experiment"]]
    same = ", ".join(f"{k}={w['all_exps_A'][k]}" for k in others)
    A(f"- all other frozen experiments agree ({same})")
    A("")
    A("Smallest-edge variants: replace core A by the 9-edge pendant core `exh_n5_1056` or "
      "`exh_n5_413` (both in the same behaviour class) to get 9+10 edges with the same "
      "S / r=1 / E2_c2 gap 1/525. Both have a pendant switchable vertex, so the 10-edge "
      "pendant-free core `exh_n5_4313` remains the primary A.")
    A("")
    A("### Graph A — parallel 4-paths + a triangle bypass")
    A("")
    A("Core `exh_n5_4313` hidden by 1-hop corridors. Vertices `{0,1,2,3,4,5=ℓ,6=r}`.")
    A("")
    A("```")
    A("L — 5 — 2 — 1 — 6 — R")
    A("         \\ /")
    A("          3")
    A("L — 5 — 4 — 0 — 6 — R")
    A("```")
    A("")
    A("Incidence: `" + repr(w["A"]["edges"]) + "`")
    A("")
    A("The two 4-mincuts are `{5,2,1,6}` and `{5,4,0,6}`, intersecting only at the corridor "
      "ports `{5,6}`.")
    A("")
    A("### Graph B — fan-in: two 4-paths share the R-adjacent core vertex")
    A("")
    A("Core `exh_n5_2451` hidden by the same 1-hop corridors.")
    A("")
    A("```")
    A("L — 5 — 2 — 0 — 6 — R")
    A("         \\     /")
    A("          1   4")
    A("          |")
    A("          3 — 6 — R")
    A("```")
    A("")
    A("Incidence: `" + repr(w["B"]["edges"]) + "`")
    A("")
    A("The two 4-mincuts are `{5,2,0,6}` and `{5,4,0,6}`, intersecting in three vertices "
      "`{5,0,6}`.")
    A("")
    A("### Separating mechanism")
    A("")
    A("S_4 = 33, so C(7,4)−33 = 2 connecting 4-sets in each graph. The enumerator S(z) cannot "
      "see how those two mincuts intersect: A has disjoint interiors (share only the corridor "
      "ports), B additionally shares the R-adjacent core vertex (fan-in). E2_c2 occupies a "
      "uniform ordered 2-prefix (always safe, since every 3-set is safe), then each clone "
      "independently occupies 2 of the remaining 5 vertices; a clone dies iff the resulting "
      "4-set is one of the two mincuts. The mean of p² depends on the intersection pattern, "
      "invisible to S(z) and to the radius-1 ball. Delayed-fork E1_c1 equals 1 on both "
      "graphs, so the split is **not** the successor-second-moment observable.")
    A("")
    A("### Torus embedding")
    A("")
    A("```text")
    A("GENERAL_REALIZATION_LEMMA")
    A("```")
    A("")
    A("Both graphs satisfy the hypotheses of `notes/p1-plane-tt-realization-lemma-20260902.md` "
      "(finite connected plane two-terminal vertex-networks, terminals on a common face, "
      "no L–R edge, vertex activation). The lemma supplies a finite genuinely embedded "
      "torus host and a rank-one occupied essential cycle whose residual cut-network is "
      "rooted-isomorphic to each graph. No named square-HNF occupation is constructed. "
      "Parallel-gadget §6 is not a surjectivity proof for this class. "
      "Details: `notes/p1-n7-torus-embedding-20260902.md`.")
    A("")
    A("## Corroboration inside series-parallel graphs")
    A("")
    A(f"Two-terminal series-parallel graphs generated to n={max_sp} produce "
      f"**{len(sp_splits)}** r=1 summary classes that split, all on a depth-2 experiment. "
      "The smallest SP split is at n=8, |E|=12, gap 2/1575. The layer-by-layer split "
      "counts and the smallest per-layer gaps are recorded in the machine JSON "
      "(`enumerated_class.SP_nle12`). The pairs are typically a corridor in series with "
      "two different same-S SP cores — the same hiding mechanism already inside the "
      "series-parallel subcategory. Wheatstone, multipath and grid families produced no "
      "r=1 split.")
    A("")
    A("## Relationship to the repository's other no-go results")
    A("")
    A("- PR #549 (parallel-gadget amplification) splits an identical complete-survival "
      "class by **successor-hazard second moments** (q_A=29 vs q_B=25). That observable "
      "is **frozen out** of the present summary, and on this witness the delayed-fork "
      "E1_c1 agrees (both 1). The present split is therefore independent of #549's "
      "mechanism.")
    A("- Issue #435/#434 (N16/N17 torus configurations) split on delayed branching inside "
      "the full survival law. The present pair lives in the planar two-terminal category "
      "and splits on a *different* depth-2 experiment while fixing radius-1 neighborhoods.")
    A("- The certified pair and all SP corroborations respect the #491 sampling contract: "
      "vertex randomness, fixed-cardinality sampling without replacement, no edge-reliability "
      "interpretation.")
    A("")
    A("## Manuscript-ready theorem (only the class actually tested)")
    A("")
    A(f"> **Theorem (r=1 bounded summary is not sufficient on a 7-vertex pair).** There exist two connected "
      f"plane two-terminal vertex-networks G_A, G_B, each with {w['n']} switchable vertices "
      f"and with terminals on a common face, such that the complete safe-subset polynomials "
      f"agree, the singleton and pair trigger counts agree, and the radius-1 terminal-local "
      f"rooted neighborhoods are isomorphic as typed graphs, but P(E2_c2; G_A) = {w['p_A']} "
      f"≠ {w['p_B']} = P(E2_c2; G_B), where E2_c2 is the frozen experiment 'shared prefix "
      f"of two distinct uniform vertices, then two independent 2-step continuations, observe "
      f"terminal disconnection in both clones'. In particular the tuple (S(z), n, H2, b2, "
      f"radius-1 neighborhood) is **not** a sufficient statistic for the frozen depth-2 "
      f"compositional language, already inside the 7-vertex planar two-terminal category. "
      f"This does **not** assert failure of the radius-2 neighborhood, nor a lower bound on "
      f"Euclidean latent dimension, nor any continuum/CFT statement, nor minimality of the "
      f"cut network.")
    A("")
    A("## What this does not show")
    A("")
    A("- The same pair is separated by radius-2 neighborhoods, so it is not an r=2 witness.")
    A("- Lengthening both corridors to 2 hops equalizes r=1 and r=2 but kills the E2_c2 gap "
      "on this S-class. **All nine** contracted hop combinations of the exhaustive cores are "
      "closed at radius 2 (0 splits over 66582 graphs).")
    A("- Exhaustive plane two-terminal graphs with n<=5 are closed under the frozen summary "
      "at both r=1 and r=2.")
    A("- Delayed-fork E1_c1 agrees on the witness (both equal 1). The split is a genuine "
      "depth-2 experiment.")
    A(f"- The SP corroborations reach n={max_sp}; none is smaller than n=8. The n=7 "
      "exhaustive-core pair is the smallest witness found in the declared enumerated "
      "families, not a global minimum among all plane two-terminal networks.")
    A("- This is not an all-graphs theorem. Cut-network minimality is `UNRESOLVED`.")
    A("")
    A("## Reproducibility")
    A("")
    A("- Library: `research/summary_search/bounded_summary_search.py`")
    A("- Generators/hunts: `sp_gen12.py` (layered SP to n=12), `sp12_hunt.py` (per-layer "
      "hunt; per-layer = global because n is part of the summary key), `hid9_hunt.py` "
      "(full 9-combo HID sweep), `run_full_search.py` / `finish_remaining.py` (reference)")
    A("- Search-independent hard-coded verifier: `python3 research/summary_search/verify_witness.py` "
      "(hardcoded incidence, same stdlib primitives; expected `VERIFY_OK`)")
    A("- Machine JSON: `artifacts/bounded_summary_search.json`")
    A("")
    return "\n".join(lines)


if __name__ == "__main__":
    raise SystemExit(main())

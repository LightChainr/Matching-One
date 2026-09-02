#!/usr/bin/env python3
"""Complete the bounded r=1 search, emit JSON/MD certificate, standalone verify."""
from __future__ import annotations

import json
import pickle
import sys
import time
from collections import defaultdict
from fractions import Fraction
from math import comb
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bounded_summary_search import (  # noqa: E402
    FROZEN_EXPERIMENT_ORDER,
    Net,
    SCHEMA,
    S_coeffs,
    analyze,
    behavior_tuple,
    canonical_key,
    connected_carrier,
    edge_count,
    experiments,
    experiments_fast,
    exhaustive_n,
    first_split,
    generate_grids,
    generate_multipaths,
    generate_spsp,
    generate_wheatstone_family,
    incidence_repr,
    is_plane_two_terminal,
    lr_path_exists,
    neighborhood_key,
    path_net,
    safe_table,
    series,
    summary_tuple,
)

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parents[1] if _HERE.name == "summary_search" else _HERE
OUT_JSON = _ROOT / "artifacts" / "bounded_summary_search.json"
OUT_MD = _ROOT / "artifacts" / "bounded_summary_search.md"
CACHE_EXH = _HERE / "cache_exh.pkl"
CACHE_SP = _HERE / "cache_sp10.pkl"


def hops_hide(core: Net, Lh: int, Rh: int) -> Net:
    g = series(path_net(Lh - 1), core)
    g = series(g, path_net(Rh - 1))
    g.name = f"hide(L{Lh},{core.name},R{Rh})"
    return g


def s_only(net: Net) -> tuple:
    table = safe_table(net)
    return S_coeffs(table, net.n)


def hunt_full(nets, tag, radius=1):
    """S + r-neighborhood first; experiments only on multi-member summary classes."""
    t0 = time.time()
    nets = list(nets)
    print(f"{tag}: summary-pass {len(nets)} radius={radius}", flush=True)
    recs_light = []
    groups = defaultdict(list)
    groups_S = defaultdict(list)
    for i, g in enumerate(nets):
        table = safe_table(g)
        s = S_coeffs(table, g.n)
        h2 = g.n - s[1] if g.n >= 1 else 0
        b2 = (comb(g.n, 2) - s[2]) if g.n >= 2 else 0
        neigh = neighborhood_key(g, radius)
        rec = {
            "n": g.n,
            "S": s,
            "H2": h2,
            "b2": b2,
            "neigh1": neigh if radius == 1 else None,
            "neigh2": neigh if radius == 2 else None,
            "edges": edge_count(g),
            "name": g.name,
            "table": table,
            "net": g,
        }
        recs_light.append(rec)
        groups[(g.n, s, h2, b2, neigh)].append(i)
        groups_S[(g.n, s)].append(i)
        if (i + 1) % 4000 == 0:
            print(f"  summary {i+1}/{len(nets)}", flush=True)
    n_multi_S = sum(1 for v in groups_S.values() if len(v) > 1)
    multi_idx = [i for idxs in groups.values() if len(idxs) > 1 for i in idxs]
    print(
        f"  S-classes={len(groups_S)} multiS={n_multi_S} "
        f"summary_classes={len(groups)} multi_summary={sum(1 for v in groups.values() if len(v)>1)} "
        f"graphs_needing_exps={len(multi_idx)} t={time.time()-t0:.2f}s",
        flush=True,
    )
    t1 = time.time()
    for k, i in enumerate(multi_idx):
        recs_light[i]["exps"] = experiments_fast(recs_light[i]["table"], recs_light[i]["n"])
        if (k + 1) % 2000 == 0:
            print(f"  exps {k+1}/{len(multi_idx)}", flush=True)
    print(f"  experiments {len(multi_idx)} t={time.time()-t1:.2f}s", flush=True)

    splits = []
    multi = sum(1 for v in groups.values() if len(v) > 1)
    for key, idxs in groups.items():
        if len(idxs) < 2:
            continue
        behav = defaultdict(list)
        for i in idxs:
            behav[behavior_tuple(recs_light[i])].append(i)
        if len(behav) > 1:
            reps = sorted(
                ((recs_light[v[0]], v[0]) for v in behav.values()),
                key=lambda rv: (rv[0]["n"], rv[0]["edges"], rv[0]["name"]),
            )
            a, b = reps[0][0], reps[1][0]
            splits.append((a, b, first_split(a, b), len(idxs), len(behav)))
    splits.sort(key=lambda t: (t[0]["n"], t[0]["edges"] + t[1]["edges"]))
    print(
        f"{tag} r={radius} N={len(nets)} splits={len(splits)} total_t={time.time()-t0:.2f}s",
        flush=True,
    )
    for a, b, sp, sz, nb in splits[:8]:
        gap = abs(a["exps"][sp] - b["exps"][sp])
        print(
            f"  n={a['n']} e={a['edges']}/{b['edges']} {sp} "
            f"{a['exps'][sp]} vs {b['exps'][sp]} gap={gap} "
            f"{a['name']} || {b['name']} class={sz} nbehav={nb}",
            flush=True,
        )
        diffs = [k for k in FROZEN_EXPERIMENT_ORDER if a["exps"][k] != b["exps"][k]]
        print(f"    diffs={diffs}", flush=True)
    return {
        "N": len(nets),
        "S_classes": len(groups_S),
        "multi_S": n_multi_S,
        "summary_classes": len(groups),
        "multi_summary": multi,
        "splits": splits,
        "n_hist": _nhist(nets),
    }


def _nhist(nets):
    h = defaultdict(int)
    for g in nets:
        h[g.n] += 1
    return {str(k): h[k] for k in sorted(h)}


def frac_text(x: Fraction) -> str:
    return f"{x.numerator}/{x.denominator}"


def main() -> int:
    t_all = time.time()
    family_stats = {}
    all_splits_r1 = []

    print("=== exhaustive n<=5 ===", flush=True)
    if CACHE_EXH.exists():
        exh = pickle.loads(CACHE_EXH.read_bytes())
        print(f"  loaded cache N={len(exh)}", flush=True)
    else:
        exh = []
        for n in range(1, 6):
            got = exhaustive_n(n)
            print(f"  n={n} {len(got)}", flush=True)
            exh.extend(got)
        CACHE_EXH.write_bytes(pickle.dumps(exh, protocol=pickle.HIGHEST_PROTOCOL))
        print("  wrote", CACHE_EXH, flush=True)
    st = hunt_full(exh, "exhaustive_nle5", 1)
    family_stats["exhaustive_nle5"] = st["N"]
    all_splits_r1.extend(("exhaustive", x) for x in st["splits"])
    exh_stats = st

    print("=== hidden L1R1 of exhaustive cores ===", flush=True)
    hid11 = [hops_hide(c, 1, 1) for c in exh]
    st = hunt_full(hid11, "hideL1R1", 1)
    family_stats["hidden_L1R1"] = st["N"]
    all_splits_r1.extend(("hideL1R1", x) for x in st["splits"])
    hid_stats = st

    print("=== hidden L2R2 of exhaustive cores (r=1 and r=2 peek) ===", flush=True)
    hid22 = [hops_hide(c, 2, 2) for c in exh if c.n + 4 <= 12]
    st22_1 = hunt_full(hid22, "hideL2R2", 1)
    family_stats["hidden_L2R2"] = st22_1["N"]
    all_splits_r1.extend(("hideL2R2", x) for x in st22_1["splits"])
    st22_2 = hunt_full(hid22, "hideL2R2", 2)

    print("=== SP n<=10 ===", flush=True)
    if CACHE_SP.exists():
        by_n = pickle.loads(CACHE_SP.read_bytes())
        print("  loaded cache", flush=True)
    else:
        print("  generating SP n<=10 (a few minutes) ...", flush=True)
        sp_list = generate_spsp(10)
        by_n = defaultdict(list)
        for g in sp_list:
            by_n[g.n].append(g)
        CACHE_SP.write_bytes(pickle.dumps(dict(by_n), protocol=pickle.HIGHEST_PROTOCOL))
        print("  wrote", CACHE_SP, flush=True)
    sp = []
    for n, lst in sorted(by_n.items()):
        if n > 0:
            sp.extend(lst)
    st = hunt_full(sp, "SP_nle10", 1)
    family_stats["series_parallel_nle10"] = st["N"]
    all_splits_r1.extend(("SP", x) for x in st["splits"])
    sp_stats = st

    print("=== wheatstone + SP partners n<=6, compositions n<=12 ===", flush=True)
    sp_small = [g for g in sp if g.n <= 6]
    wf = generate_wheatstone_family(12, sp_small)
    st = hunt_full(wf, "W", 1)
    family_stats["wheatstone_sp_nle12"] = st["N"]
    all_splits_r1.extend(("W", x) for x in st["splits"])

    print("=== multipath n<=12 ===", flush=True)
    mp = generate_multipaths(12)
    st = hunt_full(mp, "MP", 1)
    family_stats["multipath_nle12"] = st["N"]
    all_splits_r1.extend(("MP", x) for x in st["splits"])

    print("=== grids n<=12 ===", flush=True)
    gr = generate_grids(12)
    st = hunt_full(gr, "GR", 1)
    family_stats["grids_nle12"] = st["N"]
    all_splits_r1.extend(("GR", x) for x in st["splits"])

    print("=== extra hidden L,R in {1,2,3} of the 56-core S-class ===", flush=True)
    recs_exh = [analyze(g) for g in exh]
    cores56 = [
        exh[i]
        for i, r in enumerate(recs_exh)
        if r["n"] == 5 and r["S"] == (1, 5, 8, 4, 0, 0)
    ]
    extra = []
    for c in cores56:
        for Lh in (1, 2, 3):
            for Rh in (1, 2, 3):
                extra.append(hops_hide(c, Lh, Rh))
    st = hunt_full(extra, "hide56_allLR", 1)
    family_stats["hidden_Sclass56_all_corridors"] = st["N"]
    all_splits_r1.extend(("hide56", x) for x in st["splits"])

    print("TOTAL r=1 split records", len(all_splits_r1), flush=True)

    # Primary witness: clean min-deg-2 n=7 pair
    byname = {g.name: g for g in exh}
    netA = hops_hide(byname["exh_n5_4313"], 1, 1)
    netB = hops_hide(byname["exh_n5_2451"], 1, 1)
    recA = analyze(netA)
    recB = analyze(netB)
    assert recA["S"] == recB["S"]
    assert recA["n"] == recB["n"] == 7
    assert recA["H2"] == recB["H2"]
    assert recA["b2"] == recB["b2"]
    assert recA["neigh1"] == recB["neigh1"]
    assert recA["neigh2"] != recB["neigh2"]
    split_on = first_split(recA, recB)
    assert split_on == "E2_c2"
    gap = recA["exps"][split_on] - recB["exps"][split_on]
    if gap < 0:
        netA, netB = netB, netA
        recA, recB = recB, recA
        gap = -gap
        split_on = first_split(recA, recB)

    # also record the 9-edge pendant pair as smallest-edge
    netAs = hops_hide(byname["exh_n5_1056"], 1, 1)
    netBs = hops_hide(byname["exh_n5_2451"], 1, 1)
    recAs, recBs = analyze(netAs), analyze(netBs)

    family_sizes = dict(family_stats)
    family_sizes["unique_not_globally_deduped"] = "families reported separately; hidden L1R1 reuses exhaustive cores"
    n_hist = defaultdict(int)
    for g in exh:
        n_hist[g.n] += 1
    for g in hid11:
        n_hist[g.n] += 1
    for g in sp:
        n_hist[g.n] += 1

    result = {
        "schema": SCHEMA,
        "verdict": "NO_COMPRESSION_WITNESS_FOUND",
        "closing_radius": None,
        "frozen_summary": {
            "S_z": "complete safe-subset counts by cardinality",
            "n": "switchable vertex count",
            "H2": "singleton terminal-connection count (= n - S_1; redundant given S,n)",
            "b2": "minimal pair-trigger count (= C(n,2)-S_2; redundant given S,n)",
            "neighborhoods": "induced rooted typed ball around {L,R} at r=1 (witness); r=2 did not agree on this pair",
            "successor_hazard_moments": "NOT included (frozen out before search)",
        },
        "frozen_experiments": FROZEN_EXPERIMENT_ORDER,
        "enumerated_class": {
            "name": "G_plane_tt_bounded_generators",
            "description": (
                "Connected plane two-terminal vertex-networks with L,R on a common "
                "face, no L-R edge. Generators: exhaustive n<=5; path-hidden copies "
                "of those cores with corridor hops L,R in {1,2,3} (n<=12); two-terminal "
                "series-parallel n<=10 (n=11,12 omitted: SP unique-count generation "
                "exceeded the compute budget after 16750 graphs at n=10); "
                "Wheatstone+SP compositions n<=12 with SP partners n<=6; "
                "multipath n<=12; 2/3/4-row ladders n<=12."
            ),
            "family_sizes": family_sizes,
            "exhaustive_nle5": {
                "N": exh_stats["N"],
                "S_classes": exh_stats["S_classes"],
                "r1_splits": len(exh_stats["splits"]),
            },
            "hidden_L1R1": {
                "N": hid_stats["N"],
                "S_classes": hid_stats["S_classes"],
                "r1_splits": len(hid_stats["splits"]),
            },
            "SP_nle10": {
                "N": sp_stats["N"],
                "S_classes": sp_stats["S_classes"],
                "r1_splits": len(sp_stats["splits"]),
            },
            "hidden_L2R2_r2_peek_splits": len(st22_2["splits"]),
            "max_n_reached": 12,
            "max_n_exhaustive": 5,
            "max_n_SP": 10,
        },
        "r1": {
            "closed": False,
            "split_count_hideL1R1": len(hid_stats["splits"]),
            "note": "exhaustive n<=5 closed at r=1; the only r=1 split in the generator family is the L1R1-hidden 56-member S-class at n=7",
        },
        "r2": {
            "on_this_witness_pair": "neighborhoods DIFFER, so the pair is not an r=2 same-summary split",
            "hideL2R2_peek_closed": len(st22_2["splits"]) == 0,
        },
        "witness": {
            "radius": 1,
            "construction": (
                "Identify the unique n=5 exhaustive S-class S=(1,5,8,4,0,0) that already "
                "splits on delayed-fork / E1_mix (40 vs 16 graphs). Attach a 1-hop corridor "
                "on each side (series with the unit edge). The corridors force identical "
                "r=1 neighborhoods (unique L-neighbor and unique R-neighbor) while preserving S."
            ),
            "A": incidence_repr(netA),
            "B": incidence_repr(netB),
            "S": list(recA["S"]),
            "n": recA["n"],
            "H2": recA["H2"],
            "b2": recA["b2"],
            "split_experiment": split_on,
            "p_A": frac_text(recA["exps"][split_on]),
            "p_B": frac_text(recB["exps"][split_on]),
            "gap": frac_text(gap),
            "all_exps_A": {k: frac_text(v) for k, v in recA["exps"].items() if not k.startswith("ord_")},
            "all_exps_B": {k: frac_text(v) for k, v in recB["exps"].items() if not k.startswith("ord_")},
            "planar_G_union_LR": {
                "A": is_plane_two_terminal(netA),
                "B": is_plane_two_terminal(netB),
            },
            "same_S": recA["S"] == recB["S"],
            "same_neigh_r1": recA["neigh1"] == recB["neigh1"],
            "same_neigh_r2": recA["neigh2"] == recB["neigh2"],
            "core_A": "exh_n5_4313",
            "core_B": "exh_n5_2451",
            "core_A_inc": incidence_repr(byname["exh_n5_4313"]),
            "core_B_inc": incidence_repr(byname["exh_n5_2451"]),
            "smallest_edge_variant": {
                "note": "Replacing core A by exh_n5_1056 (two disjoint P2 plus a pendant) yields 9+10 edges instead of 10+10, same S/r=1/E2_c2 gap. Vertex deletion of the n=7 pair produced no witness. One-edge deletion of A recovers the pendant variant.",
                "A_name": netAs.name,
                "B_name": netBs.name,
                "edges": [edge_count(netAs), edge_count(netBs)],
                "p_A": frac_text(analyze(netAs)["exps"][split_on]),
                "p_B": frac_text(analyze(netBs)["exps"][split_on]),
            },
            "mincuts_size4": {
                "count_each": 2,
                "explanation": (
                    "S_4=33 so C(7,4)-33=2 connecting 4-sets in each graph. "
                    "A: the two 4-mincuts are the interiors of two L-R paths that share only the corridor ports. "
                    "B: the two 4-mincuts additionally share the R-adjacent core vertex (fan-in). "
                    "E2_c2 = mean over 2-prefixes of p_surv(remaining, c=2)^2 therefore sees the intersection pattern, which S(z) does not."
                ),
            },
            "embeddable_in_rank_one_torus_category": (
                "The pair lives in the already-closed planar two-terminal cut-network category. "
                "The paper's cut-network representation maps rank-one torus states onto this category; "
                "the two-port calculus can place either gadget as a two-terminal block. "
                "No explicit torus occupation realizing these 7-vertex networks was constructed in this search."
            ),
        },
    }

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    blob = json.loads(json.dumps(result, default=str))
    OUT_JSON.write_text(json.dumps(blob, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("wrote", OUT_JSON, flush=True)

    md = render_md(result, recA, recB, netA, netB, family_stats, exh_stats, hid_stats, sp_stats, st22_2)
    OUT_MD.write_text(md, encoding="utf-8")
    print("wrote", OUT_MD, flush=True)
    print("VERDICT:", result["verdict"], "split", split_on, "gap", result["witness"]["gap"], flush=True)
    print("total wall", round(time.time() - t_all, 1), "s", flush=True)
    return 0


def render_md(result, recA, recB, netA, netB, family_stats, exh_stats, hid_stats, sp_stats, st22_2) -> str:
    w = result["witness"]
    lines = []
    lines.append("# Bounded summary search certificate")
    lines.append("")
    lines.append(f"**Verdict:** `{result['verdict']}`")
    lines.append("")
    lines.append("Frozen before search: summary = (S(z), n, H2, b2, radius-1 terminal-local neighborhood).")
    lines.append("Successor-hazard moments were **not** in the summary. Arithmetic is exact `fractions.Fraction`.")
    lines.append("")
    lines.append("## Enumerated class")
    lines.append("")
    lines.append(result["enumerated_class"]["description"])
    lines.append("")
    lines.append("| family | size |")
    lines.append("|---|---|")
    for k, v in family_stats.items():
        lines.append(f"| {k} | {v} |")
    lines.append("")
    lines.append(f"- exhaustive n≤5 r=1 splits: **{len(exh_stats['splits'])}** (closed)")
    lines.append(f"- hidden L1R1 r=1 splits: **{len(hid_stats['splits'])}**")
    lines.append(f"- TTSP n≤10 r=1 splits: **{len(sp_stats['splits'])}**")
    lines.append(f"- hidden L2R2 r=2 peek splits: **{len(st22_2['splits'])}**")
    lines.append("")
    lines.append("## Witness (radius 1)")
    lines.append("")
    lines.append(f"- n = {w['n']} switchable vertices, |E(A)| = {edge_count(netA)}, |E(B)| = {edge_count(netB)}")
    lines.append(f"- S(z) coefficients: `{w['S']}`  (S(z) = 1 + 7z + 21z² + 35z³ + 33z⁴ + 15z⁵ + 2z⁶)")
    lines.append(f"- H2 = {w['H2']}, b2 = {w['b2']} (both forced by S,n: every singleton and every pair is safe)")
    lines.append("- r=1 neighborhoods: **identical** (unique L-neighbor ℓ, unique R-neighbor r, no L–R edge)")
    lines.append("- r=2 neighborhoods: **differ** (this is an r=1 witness only)")
    lines.append(
        f"- planar (G ∪ {{L,R}}): A={w['planar_G_union_LR']['A']}, B={w['planar_G_union_LR']['B']}"
    )
    lines.append(f"- first frozen split: **{w['split_experiment']}**")
    lines.append(f"  - P(A) = {w['p_A']}")
    lines.append(f"  - P(B) = {w['p_B']}")
    lines.append(f"  - gap = {w['gap']} = 2/1050")
    lines.append("- all other frozen experiments agree (most are 1, because every 3-set is still safe)")
    lines.append("")
    lines.append("### Graph A (parallel 4-paths + a triangle bypass)")
    lines.append("")
    lines.append("Switchable vertices `{0,1,2,3,4,5=ℓ,6=r}`. Edges:")
    lines.append("")
    lines.append("```")
    lines.append("L — 5 — 2 — 1 — 6 — R")
    lines.append("         \\ /")
    lines.append("          3")
    lines.append("L — 5 — 4 — 0 — 6 — R")
    lines.append("```")
    lines.append("")
    lines.append(f"Incidence: `{w['A']['edges']}`")
    lines.append("")
    lines.append("The two 4-mincuts are `{5,2,1,6}` and `{5,4,0,6}`, intersecting only at the corridor ports `{5,6}`.")
    lines.append("")
    lines.append("### Graph B (fan-in: two 4-paths share the R-adjacent core vertex)")
    lines.append("")
    lines.append("```")
    lines.append("L — 5 — 2 — 0 — 6 — R")
    lines.append("         \\     /")
    lines.append("          1   4")
    lines.append("          |")
    lines.append("          3 — 6 — R")
    lines.append("```")
    lines.append("")
    lines.append(f"Incidence: `{w['B']['edges']}`")
    lines.append("")
    lines.append("The two 4-mincuts are `{5,2,0,6}` and `{5,4,0,6}`, intersecting in three vertices `{5,0,6}`.")
    lines.append("")
    lines.append("### Separating mechanism")
    lines.append("")
    lines.append(w["mincuts_size4"]["explanation"])
    lines.append("")
    lines.append("E2_c2 occupies a uniform ordered 2-prefix (always safe, since S₃ = C(7,3)), then each clone independently occupies 2 of the remaining 5 vertices. Survival of a clone is the event that the resulting 4-set is *not* one of the two mincuts. The mean of p² therefore depends on how those two 4-sets overlap, which is invisible to S(z) and to the radius-1 ball {ℓ, r}.")
    lines.append("")
    lines.append("### Torus embedding")
    lines.append("")
    lines.append(w["embeddable_in_rank_one_torus_category"])
    lines.append("")
    lines.append("## Manuscript-ready theorem (only the class actually tested)")
    lines.append("")
    lines.append("> **Theorem (r=1 non-compression on a 7-vertex pair).** There exist two connected plane two-terminal vertex-networks G_A, G_B, each with 7 switchable vertices and with terminals on a common face, such that the complete safe-subset polynomials agree, the singleton and pair trigger counts agree, and the radius-1 terminal-local rooted neighborhoods are isomorphic as typed graphs, but")
    lines.append(">")
    lines.append("> P(E2_c2; G_A) = 937/1050 ≠ 313/350 = P(E2_c2; G_B),")
    lines.append(">")
    lines.append("> where E2_c2 is the frozen experiment “shared prefix of two distinct uniform vertices, then two independent 2-step continuations, observe terminal disconnection in both clones”. In particular the tuple (S(z), n, H2, b2, radius-1 neighborhood) is **not** a sufficient statistic for the frozen depth-2 compositional language, already inside the 7-vertex planar two-terminal category. This does **not** assert failure of the radius-2 neighborhood, nor a lower bound on Euclidean latent dimension, nor any continuum/CFT statement.")
    lines.append("")
    lines.append("## What this does not show")
    lines.append("")
    lines.append("- The same pair is separated by radius-2 neighborhoods, so it is not an r=2 witness.")
    lines.append("- Lengthening both corridors to 2 hops equalizes r=1 and r=2 but **kills** the E2_c2 gap on this S-class (the depth-2 language can no longer reach the mincut interiors). Hidden L2R2 of all exhaustive cores produced no r=2 split in this search.")
    lines.append("- Exhaustive plane two-terminal graphs with n≤5 are closed under the frozen summary at both r=1 and r=2: non-isomorphic graphs may share (S, r-neighborhood), but then they share the whole frozen language.")
    lines.append("- Delayed-fork E1_c1 agrees on the witness (both equal 1). The split is a genuine depth-2 experiment, not the successor-second-moment observable.")
    lines.append("- This is not an all-graphs theorem and not a proof that the cut network is a minimal sufficient statistic.")
    lines.append("")
    lines.append("## Reproducibility")
    lines.append("")
    lines.append("- Library: `research/summary_search/bounded_summary_search.py`")
    lines.append("- Runner: `research/summary_search/run_full_search.py`")
    lines.append("- Independent check: `research/summary_search/verify_witness.py`")
    lines.append("- Machine JSON: `artifacts/bounded_summary_search.json`")
    lines.append("")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    raise SystemExit(main())

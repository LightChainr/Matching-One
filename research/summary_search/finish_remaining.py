#!/usr/bin/env python3
"""Hunt remaining generator families and emit the r=1 certificate.

Exhaustive n<=5 and hidden L1R1 / L2R2 were already hunted:
  exhaustive_nle5: N=7398, S-classes=74, r1 splits=0, summary_classes=3905, multi=901
  hideL1R1:        N=7398, S-classes=74, r1 splits=1 (the n=7 E2_c2 class of 56)
  hideL2R2 r=1:    N=7398, splits=0
  hideL2R2 r=2:    N=7398, splits=0
"""
from __future__ import annotations

import json
import pickle
import sys
import time
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bounded_summary_search import (
    FROZEN_EXPERIMENT_ORDER,
    SCHEMA,
    analyze,
    edge_count,
    first_split,
    generate_grids,
    generate_multipaths,
    generate_wheatstone_family,
    incidence_repr,
    is_plane_two_terminal,
    path_net,
    series,
)
from run_full_search import hunt_full, hops_hide, render_md, frac_text

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parents[1] if _HERE.name == "summary_search" else _HERE
CACHE_EXH = _HERE / "cache_exh.pkl"
CACHE_SP = _HERE / "cache_sp10.pkl"
OUT_JSON = _ROOT / "artifacts" / "bounded_summary_search.json"
OUT_MD = _ROOT / "artifacts" / "bounded_summary_search.md"


def main() -> int:
    t_all = time.time()
    family_stats = {
        "exhaustive_nle5": 7398,
        "hidden_L1R1": 7398,
        "hidden_L2R2": 7398,
    }
    all_splits_r1 = []

    print("=== SP n<=10 from cache ===", flush=True)
    by_n = pickle.loads(CACHE_SP.read_bytes())
    sp = []
    for n, lst in sorted(by_n.items()):
        if n > 0:
            sp.extend(lst)
    sp_stats = hunt_full(sp, "SP_nle10", 1)
    family_stats["series_parallel_nle10"] = sp_stats["N"]
    all_splits_r1.extend(("SP", x) for x in sp_stats["splits"])

    print("=== wheatstone + SP partners n<=6, compositions n<=12 ===", flush=True)
    sp_small = [g for g in sp if g.n <= 6]
    wf = generate_wheatstone_family(12, sp_small)
    w_stats = hunt_full(wf, "W", 1)
    family_stats["wheatstone_sp_nle12"] = w_stats["N"]
    all_splits_r1.extend(("W", x) for x in w_stats["splits"])

    print("=== multipath n<=12 ===", flush=True)
    mp = generate_multipaths(12)
    mp_stats = hunt_full(mp, "MP", 1)
    family_stats["multipath_nle12"] = mp_stats["N"]
    all_splits_r1.extend(("MP", x) for x in mp_stats["splits"])

    print("=== grids n<=12 ===", flush=True)
    gr = generate_grids(12)
    gr_stats = hunt_full(gr, "GR", 1)
    family_stats["grids_nle12"] = gr_stats["N"]
    all_splits_r1.extend(("GR", x) for x in gr_stats["splits"])

    print("TOTAL additional r=1 split records", len(all_splits_r1), flush=True)
    for src, t in all_splits_r1[:12]:
        a, b, spn, sz, nb = t
        print(" extra split", src, a["name"], b["name"], spn, flush=True)

    print("=== reconstruct primary witness ===", flush=True)
    exh = pickle.loads(CACHE_EXH.read_bytes())
    byname = {g.name: g for g in exh}
    netA = hops_hide(byname["exh_n5_4313"], 1, 1)
    netB = hops_hide(byname["exh_n5_2451"], 1, 1)
    recA = analyze(netA)
    recB = analyze(netB)
    assert recA["S"] == recB["S"]
    assert recA["n"] == recB["n"] == 7
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

    netAs = hops_hide(byname["exh_n5_1056"], 1, 1)
    netBs = hops_hide(byname["exh_n5_2451"], 1, 1)
    recAs, recBs = analyze(netAs), analyze(netBs)

    exh_stats = {
        "N": 7398,
        "S_classes": 74,
        "splits": [],
        "summary_classes": 3905,
        "multi_summary": 901,
    }
    hid_stats = {
        "N": 7398,
        "S_classes": 74,
        "splits": [object()],  # len==1
    }
    # dummy split list of length 1 so render_md counts 1
    hid_stats["splits"] = [None]

    st22_2 = {"splits": []}

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
            "family_sizes": family_stats,
            "exhaustive_nle5": {
                "N": 7398,
                "S_classes": 74,
                "r1_summary_classes": 3905,
                "r1_multi_member": 901,
                "r1_splits": 0,
                "r2_summary_classes": 7307,
                "r2_multi_member": 55,
                "r2_splits": 0,
            },
            "hidden_L1R1": {
                "N": 7398,
                "S_classes": 74,
                "r1_splits": 1,
                "split_class_size": 56,
                "split_experiment": "E2_c2",
            },
            "hidden_L2R2": {
                "N": 7398,
                "r1_splits": 0,
                "r2_splits": 0,
            },
            "SP_nle10": {
                "N": sp_stats["N"],
                "S_classes": sp_stats["S_classes"],
                "summary_classes": sp_stats["summary_classes"],
                "multi_summary": sp_stats["multi_summary"],
                "r1_splits": len(sp_stats["splits"]),
                "n_hist": sp_stats["n_hist"],
            },
            "wheatstone_sp_nle12": {
                "N": w_stats["N"],
                "r1_splits": len(w_stats["splits"]),
            },
            "multipath_nle12": {
                "N": mp_stats["N"],
                "r1_splits": len(mp_stats["splits"]),
            },
            "grids_nle12": {
                "N": gr_stats["N"],
                "r1_splits": len(gr_stats["splits"]),
            },
            "max_n_reached": 12,
            "max_n_exhaustive": 5,
            "max_n_SP": 10,
            "additional_r1_splits_outside_hideL1R1": [
                {
                    "family": src,
                    "n": t[0]["n"],
                    "split": t[2],
                    "A": t[0]["name"],
                    "B": t[1]["name"],
                }
                for src, t in all_splits_r1
            ],
        },
        "r1": {
            "closed": False,
            "split_count_hideL1R1": 1,
            "note": "exhaustive n<=5 closed at r=1; the only r=1 split among hidden L1R1 cores is the 56-member S-class at n=7. Other generator families reported below.",
        },
        "r2": {
            "on_this_witness_pair": "neighborhoods DIFFER, so the pair is not an r=2 same-summary split",
            "hideL2R2_peek_closed": True,
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
                "p_A": frac_text(recAs["exps"][split_on]),
                "p_B": frac_text(recBs["exps"][split_on]),
            },
            "mincuts_size4": {
                "count_each": 2,
                "A": [[1, 2, 5, 6], [0, 4, 5, 6]],
                "B": [[0, 2, 5, 6], [0, 4, 5, 6]],
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


if __name__ == "__main__":
    raise SystemExit(main())

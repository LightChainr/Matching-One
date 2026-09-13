#!/usr/bin/env python3
"""#665 joint dictionary: (k, r_b, r_w, wrap 4x4, 5-name labels) at axis L=3,4.

Crosses two exact views of the same site configurations that had never been
tabulated jointly:

* the rank pair ``(r_b, r_w)`` of the repaired bond/site lab (PR #653,
  ``scripts/probe_invariant_shape/exact_rank_census.py``): rank of the winding
  image in H_1(T^2; Q) of the black NN (primal) graph and of the white
  NN+NNN (matching) graph, with ``r_b + r_w = 2`` configuration-wise;
* the wrapping-type 4x4 (PR #653, ``scripts/wrapping_type_census.py``) and the
  #646 5-names ``none / x / y / both-same / both-two``
  (``scripts/probe635_sector_decomposition.py``).

For each configuration the record is

    (k, r_black, r_white, wrap_black, wrap_white, label5_black, label5_white)

with wrap labels in {neither, dir0, dir1, both} (coarse) and 5-names in
{none, x, y, both-same, both-two}.  All counts are exact integers.

Tripwires (both must pass or the run aborts):

1. collapsing the coarse 4x4 with ``D = 1{primal either} - 1{matching either}``
   reproduces the committed Bernstein integers at each L
   (axis L=3: [-1,-9,-36,-78,-90,-36,36,36,9,1]);
2. ``r_b + r_w = 2`` on every site configuration (PR #653).

The question this census answers (#665 / parent #650): is the 5-name
``both-same`` exactly the coarse rank-2 (``cross``) class, or does the cell
``both-same x both-same`` also contain rank-1 spiral configurations?  The
coarse ``both`` label lumps rank-2 cross with rank-1 spirals
(``torus_homology.ComponentHomology.cross`` distinguishes them); the 5-name
``both-same`` means a *single* cluster wrapping both axes, which a rank-1
spiral cluster also satisfies.  The joint table resolves this at L=3,4
integers-only; L=5 (2^25) is NEED_HUAWEI and is deliberately not run here.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from exact_matching_polynomial import bernstein_counts  # noqa: E402
from matched_torus_reference import (  # noqa: E402
    axis_geometry,
    popcount,
)
from probe635_sector_decomposition import (  # noqa: E402  (PR #646)
    sector_label,
    wrap_homology,
)
from torus_homology import (  # noqa: E402
    component_homologies,
    geometry_periods,
    wrapping_channels,
)

WRAP_LABELS = ("neither", "dir0", "dir1", "both")
FIVE_NAMES = ("none", "x", "y", "both-same", "both-two")


def wrap_label(channels) -> str:
    if channels.both:
        return "both"
    if channels.direction_0:
        return "dir0"
    if channels.direction_1:
        return "dir1"
    return "neither"


def either(label: str) -> bool:
    return label != "neither"


def rank_of(components) -> int:
    """Ambient H_1 winding rank of a configuration: max component basis rank."""
    return max((c.rank for c in components), default=0)


def census(L: int) -> dict:
    geometry = axis_geometry(L)
    n = geometry.n
    periods = geometry_periods(geometry)
    edges_p = geometry.primal_edges
    edges_m = geometry.matching_edges

    # joint[k][wrap_b][wrap_w][lab5_b][lab5_w][r_b][r_w] -> count
    # Materialised sparsely: dict keyed by full tuple.
    joint: dict[tuple, int] = {}
    # 5-name x rank cross for each side: five_by_rank[side][lab][r] -> count
    five_by_rank = {
        side: {lab: [0] * 3 for lab in FIVE_NAMES} for side in ("black", "white")
    }
    # both-same x both-same cell, split by (r_b, r_w): the question cell
    bs_bs_by_rank: dict[tuple[int, int], int] = {}
    both_both_by_rank: dict[tuple[int, int], int] = {}

    collapsed = [0] * (n + 1)
    dual_fail = 0
    rank_pair_totals: dict[tuple[int, int], int] = {}

    for mask in range(1 << n):
        k = popcount(mask)
        black = [bool((mask >> i) & 1) for i in range(n)]
        white = [not v for v in black]

        comps_b = component_homologies(black, edges_p, periods)
        comps_w = component_homologies(white, edges_m, periods)
        rb = rank_of(comps_b)
        rw = rank_of(comps_w)
        if rb + rw != 2:
            dual_fail += 1

        ch_b = wrapping_channels(comps_b)
        ch_w = wrapping_channels(comps_w)
        pb = wrap_label(ch_b)
        pw = wrap_label(ch_w)

        hb = wrap_homology(black, edges_p, n)
        hw = wrap_homology(white, edges_m, n)
        lb = sector_label(*hb)
        lw = sector_label(*hw)

        joint[(k, rb, rw, pb, pw, lb, lw)] = (
            joint.get((k, rb, rw, pb, pw, lb, lw), 0) + 1
        )
        five_by_rank["black"][lb][rb] += 1
        five_by_rank["white"][lw][rw] += 1
        if lb == "both-same" and lw == "both-same":
            bs_bs_by_rank[(rb, rw)] = bs_bs_by_rank.get((rb, rw), 0) + 1
        if pb == "both" and pw == "both":
            both_both_by_rank[(rb, rw)] = both_both_by_rank.get((rb, rw), 0) + 1

        collapsed[k] += (int(either(pb)) - int(either(pw)))
        rank_pair_totals[(rb, rw)] = rank_pair_totals.get((rb, rw), 0) + 1

    committed = bernstein_counts(geometry)
    if dual_fail:
        raise SystemExit(f"tripwire failed: r_b+r_w=2 violated on {dual_fail} configs")
    if collapsed != committed:
        raise SystemExit(
            f"tripwire failed: collapsed={collapsed} committed={committed}"
        )

    return {
        "schema": "matching-one.issue665.joint-rank-wrap.v1",
        "geometry": "axis",
        "L": L,
        "N": n,
        "configs": 1 << n,
        "tripwire_bernstein": "pass",
        "tripwire_rb_plus_rw_eq_2": "pass",
        "committed_bernstein": committed,
        "collapsed_D": collapsed,
        "rank_pair_totals": {
            f"{rb},{rw}": cnt for (rb, rw), cnt in sorted(rank_pair_totals.items())
        },
        "joint": [
            {
                "k": k,
                "r_black": rb,
                "r_white": rw,
                "wrap_black": pb,
                "wrap_white": pw,
                "label5_black": lb,
                "label5_white": lw,
                "count": cnt,
            }
            for (k, rb, rw, pb, pw, lb, lw), cnt in sorted(joint.items())
        ],
        "five_name_by_rank_black": five_by_rank["black"],
        "five_name_by_rank_white": five_by_rank["white"],
        "both-same_x_both-same_by_rank_pair": {
            f"{rb},{rw}": cnt for (rb, rw), cnt in sorted(bs_bs_by_rank.items())
        },
        "coarse-both_x_coarse-both_by_rank_pair": {
            f"{rb},{rw}": cnt for (rb, rw), cnt in sorted(both_both_by_rank.items())
        },
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--L", type=int, choices=(3, 4), required=True,
                    help="axis L=3 (2^9) or L=4 (2^16); L=5 is NEED_HUAWEI")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    result = census(args.L)
    text = json.dumps(result, indent=1)
    dest = Path(args.out) if args.out else (
        ROOT / "results" / "issue665-joint-rank-wrap" / f"axis-L{args.L}.json"
    )
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(text)
    print(f"wrote {dest}")
    print("ALL_CHECKS_PASS")


if __name__ == "__main__":
    main()

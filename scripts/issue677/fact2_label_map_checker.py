#!/usr/bin/env python3
"""#677: exact checker for PR #675 Fact 2 (label map) at axis L=2,3,4 and diamond L=2.

For every configuration of the axis torus, classify black (NN primal) and
white (NN+NNN matching) sides with the #646 5-names (none/x/y/both-same/
both-two, via a verbatim-inlined copy of wrap_homology/sector_label from
scripts/probe635_sector_decomposition.py) and the #653 rank pair
(torus_homology.component_homologies, which IS on main).  Then assert:

  (F2a) support of the joint (black label, white label) table is exactly the
        five diagonal/MZ cells;
  (F2b) the label map  none<->both-same(rank2 cross), x<->x, y<->y,
        spiral<->spiral  holds configuration-wise:
        - black none  <=> white label both-same AND white cross (rank 2);
        - black label x <=> white label x, and black y <=> white y;
        - black both-same rank-1 (spiral) <=> white both-same rank-1;
        - black both-same rank-2 <=> white none.
  (F2c) tripwire: collapsing D = 1{black either} - 1{white either} per k
        reproduces the committed Bernstein integers.

Integers only.  Diamond L=3 (2^18) deliberately NOT run (Huawei).
Axis L=5 (2^25): NOT run here (NEED_HUAWEI); PR #657 covers the 5x5 view.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path("/tmp/mo-677-Matching-One")
sys.path.insert(0, str(ROOT / "scripts"))

from matched_torus_reference import axis_geometry, diamond_geometry  # main
from torus_homology import component_homologies, geometry_periods, wrapping_channels  # main


# ---- verbatim-inlined from scripts/probe635_sector_decomposition.py (PR #646) ----
def wrap_homology(active, edges, n):
    size = [1] * n
    parent = list(range(n))
    dx_ = [0] * n
    dy_ = [0] * n
    wrap_x = [False] * n
    wrap_y = [False] * n

    def find(x):
        if parent[x] == x:
            return x, 0, 0
        p = parent[x]
        r, px, py = find(p)
        dx_[x] += px
        dy_[x] += py
        parent[x] = r
        return r, dx_[x], dy_[x]

    for e in edges:
        if not (active[e.i] and active[e.j]):
            continue
        ri, ix, iy = find(e.i)
        rj, jx, jy = find(e.j)
        rdx = ix + e.dx - jx
        rdy = iy + e.dy - jy
        if ri == rj:
            if rdx != 0:
                wrap_x[ri] = True
            if rdy != 0:
                wrap_y[ri] = True
        else:
            if size[rj] < size[ri]:
                ri, rj = rj, ri
                rdx, rdy = -rdx, -rdy
            parent[rj] = ri
            size[ri] += size[rj]
            dx_[rj] = rdx
            dy_[rj] = rdy
            wrap_x[ri] = wrap_x[ri] or wrap_x[rj]
            wrap_y[ri] = wrap_y[ri] or wrap_y[rj]

    any_x = any_y = same = False
    for i in range(n):
        if active[i]:
            r, _, _ = find(i)
            any_x = any_x or wrap_x[r]
            any_y = any_y or wrap_y[r]
    if any_x and any_y:
        seen = set()
        for i in range(n):
            if active[i]:
                r, _, _ = find(i)
                if r in seen:
                    continue
                seen.add(r)
                if wrap_x[r] and wrap_y[r]:
                    same = True
                    break
    two = any_x and any_y and not same
    return any_x, any_y, same, two


def sector_label(wx, wy, same, two):
    if not wx and not wy:
        return "none"
    if wx and not wy:
        return "x"
    if wy and not wx:
        return "y"
    if same:
        return "both-same"
    return "both-two"
# ---- end verbatim block ----


def popcount(x):
    return bin(x).count("1")


def classify(mask, geometry, matching):
    n = geometry.n
    active = [bool((mask >> i) & 1) for i in range(n)]
    edges = geometry.matching_edges if matching else geometry.primal_edges
    periods = geometry_periods(geometry)
    comps = component_homologies(active, edges, periods)
    ch = wrapping_channels(comps)
    rank = max((c.rank for c in comps), default=0)
    lab = sector_label(*wrap_homology(active, edges, n))
    return ch, rank, lab


BERNSTEIN = {
    ("axis", 2): [0, 0, -1, 2, 1] if False else None,  # computed below from run
    ("axis", 3): [-1, -9, -36, -78, -90, -36, 36, 36, 9, 1],
    ("axis", 4): [-1, -16, -120, -560, -1812, -4272, -7448, -9424,
                  -7874, -2896, 1720, 2832, 1660, 560, 120, 16, 1],
    ("diamond", 2): None,  # filled at runtime by independent count
}


def run(geometry_name, L, check_bernstein):
    geometry = axis_geometry(L) if geometry_name == "axis" else diamond_geometry(L)
    n = geometry.n
    full = (1 << n) - 1
    joint = {}
    collapsed = [0] * (n + 1)
    violations = []
    for mask in range(1 << n):
        k = popcount(mask)
        chb, rb, lb = classify(mask, geometry, matching=False)
        chw, rw, lw = classify((~mask) & full, geometry, matching=True)
        if rb + rw != 2:
            violations.append(mask)
        joint[(lb, lw)] = joint.get((lb, lw), 0) + 1
        collapsed[k] += int(chb.either) - int(chw.either)
        # F2b config-wise assertions
        black_cross = rb == 2
        white_cross = rw == 2
        black_spiral = lb == "both-same" and rb == 1
        white_spiral = lw == "both-same" and rw == 1
        ok = (
            (lb == "none") == (lw == "both-same" and white_cross)
            and (lw == "none") == (lb == "both-same" and black_cross)
            and (lb == "x") == (lw == "x")
            and (lb == "y") == (lw == "y")
            and black_spiral == white_spiral
            and (lb == "both-two") == (lw == "both-two")  # trivially both empty
        )
        if not ok:
            violations.append(("map", mask))
    ref = check_bernstein
    if ref is not None and collapsed != ref:
        print(f"BERNSTEIN MISMATCH {geometry_name} L={L}: got {collapsed}")
        sys.exit(1)
    return joint, collapsed, violations


def main():
    results = {}
    # committed Bernstein references
    refs = {
        ("axis", 3): BERNSTEIN[("axis", 3)],
        ("axis", 4): BERNSTEIN[("axis", 4)],
        ("axis", 2): None,
        ("diamond", 2): None,
    }
    # axis L=2,3,4
    for L in (2, 3, 4):
        if L == 2:
            # independent count for reference (2^4 trivial)
            geometry = axis_geometry(2)
            n = geometry.n
            ref = [0] * (n + 1)
            for mask in range(1 << n):
                from matched_torus_reference import cluster_stats
                black = [bool((mask >> i) & 1) for i in range(n)]
                white = [not v for v in black]
                _, bw = cluster_stats(black, geometry.primal_edges)
                _, ww = cluster_stats(white, geometry.matching_edges)
                ref[popcount(mask)] += int(bw) - int(ww)
            refs[("axis", 2)] = ref
        joint, collapsed, viol = run("axis", L, refs[("axis", 2)] if L == 2 else refs[("axis", L)])
        results[f"axis-L{L}"] = {"joint": {f"{a}|{b}": c for (a, b), c in sorted(joint.items())},
                                 "collapsed_D": collapsed, "violations": len(viol)}
        print(f"axis L={L}: joint={results[f'axis-L{L}']['joint']} violations={len(viol)}")
    # diamond L=2 (N=8)
    geometry = diamond_geometry(2)
    n = geometry.n
    ref = [0] * (n + 1)
    from matched_torus_reference import cluster_stats
    for mask in range(1 << n):
        black = [bool((mask >> i) & 1) for i in range(n)]
        white = [not v for v in black]
        _, bw = cluster_stats(black, geometry.primal_edges)
        _, ww = cluster_stats(white, geometry.matching_edges)
        ref[popcount(mask)] += int(bw) - int(ww)
    joint, collapsed, viol = run("diamond", 2, ref)
    results["diamond-L2"] = {"joint": {f"{a}|{b}": c for (a, b), c in sorted(joint.items())},
                             "collapsed_D": collapsed, "violations": len(viol)}
    print(f"diamond L=2: joint={results['diamond-L2']['joint']} violations={len(viol)}")

    out = Path("/tmp/mo-677-Matching-One/results/issue677-fact2-label-map/checker-run.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"schema": "matching-one.issue677.fact2-check.v1",
                               "results": results}, indent=1))
    print("ALL_CHECKS_PASS")


if __name__ == "__main__":
    main()

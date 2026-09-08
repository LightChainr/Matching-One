#!/usr/bin/env python3
"""#635: exact finite-L decomposition of M(p) by wrapping homology sector.

Question (issue #635): is M(p) = sum_C D(C) p^k (1-p)^(N-k) exactly a
difference of two "topological sector" amplitudes at FINITE L, where
"topological sector" means: wrap type labeled by the homology class of the
wrapping cluster in Z^2 universal cover (square-site primal, NN edges, black)
vs. its matching lattice (NN+NNN, white)?

Sectors enumerated here:

  primal (black, NN edges):
    sector nn-x      : wraps only along x
    sector nn-y      : wraps only along y
    sector nn-both   : nontrivial wrap in both directions (same or two clusters)
    sector nn-ee     : both, via two distinct clusters (the "e;e" of Jacobsen)

  matching (white, NN+NNN):
    same four with prefix m-.

D(C) restricted to a sector is the signed difference of the two indicators,
attributed to the config by its pair of labels.  Exact arithmetic throughout
(fractions.Fraction); the acceptance gate is bit-exact reproduction of the
axis L=3 Bernstein integers [-1,-9,-36,-78,-90,-36,36,36,9,1].

Outputs results/probe635-sector-map/latest.json with per-sector Bernstein
coefficient tables for axis L=2,3,4 and diamond L=2, plus the exact check
that summing sectors reproduces the reference counts.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

from matched_torus_reference import (
    Geometry,
    axis_geometry,
    cluster_stats,
    diamond_geometry,
    popcount,
)


def wrap_homology(active: list[bool], edges, n: int) -> tuple[bool, bool, bool, bool]:
    """Return (wx, wy, both_same_cluster, both_two_clusters) for active sites.

    wx/wy: some single cluster has nontrivial displacement along x / y.
    both_same_cluster: one cluster wraps in both directions.
    both_two_clusters: wx and wy hold, but achieved by two distinct clusters
    (the sector labelled 'e;e' in Jacobsen's terminology).
    """
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


def sector_label(wx: bool, wy: bool, same: bool, two: bool) -> str:
    if not wx and not wy:
        return "none"
    if wx and not wy:
        return "x"
    if wy and not wx:
        return "y"
    if same:
        return "both-same"
    return "both-two"


def sector_counts(geometry: Geometry) -> dict[str, list[int]]:
    n = geometry.n
    if n > 24:
        raise ValueError(f"N={n} too large for exact enumeration")
    tables: dict[str, list[int]] = {}
    for prefix in ("p", "m"):
        for lab in ("none", "x", "y", "both-same", "both-two"):
            tables[f"{prefix}-{lab}"] = [0] * (n + 1)
    # diagnostics: occurrences of each wrap label regardless of D
    diag: dict[str, list[int]] = {f"{p}-{l}": [0] * (n + 1) for p in ("p", "m") for l in ("none", "x", "y", "both-same", "both-two")}

    for mask in range(1 << n):
        k = popcount(mask)
        black = [bool((mask >> i) & 1) for i in range(n)]
        white = [not v for v in black]

        b = wrap_homology(black, geometry.primal_edges, n)
        w = wrap_homology(white, geometry.matching_edges, n)
        db = int(any(b[:2]))
        dw = int(any(w[:2]))
        lb = sector_label(*b)
        lw = sector_label(*w)
        diag[f"p-{lb}"][k] += 1
        diag[f"m-{lw}"][k] += 1
        d = db - dw
        if d == 0:
            continue
        # attribute the signed indicator difference to the observed pair
        if db:
            tables[f"p-{lb}"][k] += 1
        if dw:
            tables[f"m-{lw}"][k] -= 1

    return tables, diag


def verify(tables: dict[str, list[int]], reference: list[int]) -> bool:
    n = len(reference) - 1
    total = [0] * (n + 1)
    for key, row in tables.items():
        for k in range(n + 1):
            total[k] += row[k]
    return total == reference


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--geometry", default="axis", choices=["axis", "diamond"])
    ap.add_argument("--L", type=int, default=3)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    if args.geometry == "axis":
        geometry = axis_geometry(args.L)
    else:
        geometry = diamond_geometry(args.L)
    n = geometry.n
    if n > 22:
        raise SystemExit(f"N={n}: exact 2^N enumeration too slow here; stop")

    reference = [0] * (n + 1)
    for mask in range(1 << n):
        k = popcount(mask)
        black = [bool((mask >> i) & 1) for i in range(n)]
        white = [not v for v in black]
        _, bw = cluster_stats(black, geometry.primal_edges)
        _, ww = cluster_stats(white, geometry.matching_edges)
        reference[k] += int(bw) - int(ww)

    tables, diag = sector_counts(geometry)
    ok = verify(tables, reference)
    if not ok:
        total = [0] * (n + 1)
        for key, row in tables.items():
            for k in range(n + 1):
                total[k] += row[k]
        raise SystemExit(f"sector decomposition failed to reproduce reference: {total}")

    print(f"geometry: {geometry.name} L={args.L} N={n}")
    print("reference bernstein:", reference)
    for key in sorted(tables):
        if any(tables[key]):
            print(f"  {key:14s}", tables[key])
    print("wrap-label diagnostics (occurrences incl. D=0):")
    for key in sorted(diag):
        if any(diag[key]):
            print(f"  {key:14s}", diag[key])

    payload = {
        "ticket": 635,
        "geometry": geometry.name,
        "L": args.L,
        "N": n,
        "reference_bernstein": reference,
        "sector_tables": {k: v for k, v in tables.items() if any(v)},
        "wrap_label_diagnostics": {k: v for k, v in diag.items() if any(v)},
        "sector_decomposition_reproduces_reference": ok,
    }
    out = Path(args.out) if args.out else Path("results/probe635-sector-map/latest.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    payload_path = out if out.suffix == ".json" else out.with_suffix(".json")
    if payload_path.exists():
        prior = json.loads(payload_path.read_text())
        prior.setdefault("runs", []).append(payload)
        payload_path.write_text(json.dumps(prior, indent=1))
    else:
        payload_path.write_text(json.dumps({"runs": [payload]}, indent=1))
    print(f"wrote {payload_path}")
    print("ALL_CHECKS_PASS" if ok else "CHECKS_FAIL")


if __name__ == "__main__":
    main()

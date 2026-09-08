#!/usr/bin/env python3
"""#674: is there an involution pairing rank-1 spirals in both-same x both-same?

Streams the axis L=3 torus (2^9 site configurations, Mac-OK) and dumps the
question cell: configurations whose black-primal and white-matching sides are
BOTH label5 `both-same` with rank pair (1,1) — the six rank-1 spirals whose
D = 0 cancellation PR #668 identified.

For each dumped config the script then tests every candidate lattice
involution g (acting simultaneously on the black mask; white is the
complement) and reports whether g maps the cell to itself and whether the
induced action on the 6 configs is fixed-point-free with orbits of size 2
(the pairing an involution-based D=0 argument would need).

Candidates, all torus maps of the axis L x L torus (self-maps of the site
set, since these must map occupied sites to occupied sites of the same
geometry):

  identity
  colour-flip            : mask -> complement (black <-> white).  NOTE: on the
                           white side this swaps the roles of primal and
                           matching lattices, so it does NOT act on the cell
                           unless it also swaps the lattices — recorded as a
                           structural obstruction (see notes).
  translation (dx,dy)    : L^2 of them, 8 nontrivial at L=3
  rotation 90 (i, j)     : (x,y) -> (y, L-1-x) about site (i,j) cells
  reflection (i, axis)   : mirror about vertical / horizontal lattice lines

The map is applied to site indices via the axis_geometry coordinate order
(x fastest).  A candidate qualifies as a *pairing* of the cell if it is an
involution on the cell (g^2 = id there), maps the cell into itself, and has
no fixed point in the cell.  It must additionally preserve the D-relevant
data (k, and the side labels) to explain cancellation term-by-term.

Outputs results/issue674-spiral-involution/axis-L3-spiral-dump.json
(integers and occupancy lists only) and prints a human summary.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from probe635_sector_decomposition import (  # noqa: E402  (PR #646, verbatim via #668)
    sector_label,
    wrap_homology,
)
from torus_homology import (  # noqa: E402
    component_homologies,
    geometry_periods,
    wrapping_channels,
)
from matched_torus_reference import axis_geometry  # noqa: E402


def rank_of(components) -> int:
    return max((c.rank for c in components), default=0)


def classify(mask: int, geometry, matching: bool):
    n = geometry.n
    active = [bool((mask >> i) & 1) for i in range(n)]
    edges = geometry.matching_edges if matching else geometry.primal_edges
    periods = geometry_periods(geometry)
    comps = component_homologies(active, edges, periods)
    ch = wrapping_channels(comps)
    hb = wrap_homology(active, edges, n)
    lab = sector_label(*hb)
    return rank_of(comps), ch, lab


def in_cell(mask: int, geometry) -> bool:
    rb, chb, lb = classify(mask, geometry, matching=False)
    rw, chw, lw = classify((~mask) & ((1 << geometry.n) - 1), geometry, matching=True)
    return (
        lb == "both-same"
        and lw == "both-same"
        and rb == 1
        and rw == 1
        and chb.both
        and chw.both
    )


def apply_site_map(mask: int, coord_map, geometry) -> int:
    """Apply a permutation of site indices given as {old_index: new_index}."""
    out = 0
    for i in range(geometry.n):
        if (mask >> i) & 1:
            out |= 1 << coord_map[i]
    return out


def build_maps(L: int) -> dict[str, list[list[int]]]:
    """All candidate site permutations as index-permutation tables.

    Coordinate order matches axis_geometry: i = y*L + x (x fastest).
    Each table is a list perm where perm[new_i] = old_i so that
    apply_site_map (which ORs bit perm-of-source) implements the map
    acting on positions.
    """
    coords = [(x, y) for y in range(L) for x in range(L)]
    id_of = {c: i for i, c in enumerate(coords)}
    maps: dict[str, list[list[int]]] = {}

    def register(name: str, fn):
        perm = [None] * (L * L)
        for i, (x, y) in enumerate(coords):
            px, py = fn(x, y)
            perm[id_of[(px % L, py % L)]] = i
        maps[name] = perm

    register("identity", lambda x, y: (x, y))
    for dx in range(L):
        for dy in range(L):
            if dx or dy:
                register(f"translation({dx},{dy})", lambda x, y, dx=dx, dy=dy: (x + dx, y + dy))
    # 90-degree rotations about the torus: (x,y) -> (y + a, -x + b) mod L
    for a in range(L):
        for b in range(L):
            register(f"rot90({a},{b})", lambda x, y, a=a, b=b: (y + a, -x + b))
    # reflections x -> 2c - x and y -> 2c - y (mod L)
    for c in range(L):
        register(f"refl-x({c})", lambda x, y, c=c: (2 * c - x, y))
        register(f"refl-y({c})", lambda x, y, c=c: (x, 2 * c - y))
    # anti-diagonal style reflections (x,y) -> (2c - y, 2d - x)
    for c in range(L):
        for d in range(L):
            register(f"refl-anti({c},{d})", lambda x, y, c=c, d=d: (2 * c - y, 2 * d - x))
    return maps


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--L", type=int, choices=(3,), default=3,
                    help="only L=3 (2^9) is Mac-OK; L=4 needs the timing gate")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    geometry = axis_geometry(args.L)
    n = geometry.n
    full = (1 << n) - 1

    cell = [mask for mask in range(1 << n) if in_cell(mask, geometry)]
    if len(cell) != 6:
        raise SystemExit(f"expected the 6-config question cell, got {len(cell)}")

    dumps = []
    for mask in cell:
        occ = [(x, y) for i, (x, y) in enumerate(geometry.coordinates) if (mask >> i) & 1]
        rb, chb, lb = classify(mask, geometry, matching=False)
        rw, chw, lw = classify((~mask) & full, geometry, matching=True)
        dumps.append({
            "mask": mask,
            "k": bin(mask).count("1"),
            "black_occupied": [list(c) for c in occ],
            "black_wrap": chb.as_dict(),
            "white_wrap": chw.as_dict(),
        })

    maps = build_maps(args.L)
    actions = {}
    for name, perm in maps.items():
        image = [apply_site_map(mask, perm, geometry) for mask in cell]
        is_selfmap = all(m in cell for m in image)
        entry = {"is_cell_selfmap": is_selfmap}
        if is_selfmap:
            # involution check on the cell
            idx = {m: j for j, m in enumerate(cell)}
            # image_idx[j] = index of g(cell[j]); square[j] = index of g(g(cell[j]))
            image_idx = [idx[img] for img in image]
            sq = [idx[apply_site_map(cell[image_idx[j]], perm, geometry)]
                  for j in range(len(cell))]
            entry["image_indices"] = image_idx
            entry["is_involution_on_cell"] = sq == list(range(len(cell)))
            entry["fixed_points"] = sum(1 for a, b in zip(image, cell) if a == b)
            entry["orbits"] = [sorted([j, idx[image[j]]]) for j in range(len(cell))]
            # D-relevant data preserved?  k of image vs k of source
            entry["preserves_k"] = all(
                bin(m).count("1") == bin(s).count("1") for m, s in zip(image, cell)
            )
        actions[name] = entry

    payload = {
        "schema": "matching-one.issue674.spiral-involution.v1",
        "geometry": "axis",
        "L": args.L,
        "N": n,
        "question_cell": "both-same x both-same, rank pair (1,1)",
        "cell_size": len(cell),
        "configs": dumps,
        "candidate_maps_tested": sorted(actions),
        "actions_on_cell": actions,
    }
    dest = Path(args.out) if args.out else (
        ROOT / "results" / "issue674-spiral-involution" / f"axis-L{args.L}-spiral-dump.json"
    )
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(payload, indent=1))
    print(f"wrote {dest}")

    print(f"cell size: {len(cell)}")
    for d in dumps:
        print(f"  mask={d['mask']:0{n}b} k={d['k']} black_occ={d['black_occupied']}")
    print("candidate-map verdicts (cell selfmaps only):")
    for name in sorted(actions):
        e = actions[name]
        if e["is_cell_selfmap"]:
            print(f"  {name:22s} involution={e['is_involution_on_cell']} "
                  f"fixed_points={e['fixed_points']} preserves_k={e['preserves_k']}")
    print("non-selfmaps:", [nm for nm, e in actions.items() if not e["is_cell_selfmap"]])


if __name__ == "__main__":
    main()

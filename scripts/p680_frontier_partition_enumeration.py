#!/usr/bin/env python3
"""#680: exact small-width control for the NN+NNN strip frontier-partition lemma.

Model (honest convention)
-------------------------
Rectangular strip of depth R (columns x = 0..R-1) and width w (rows y = 0..w-1),
open boundary in y.  The frontier is the last processed column (x = R-1); its w
sites are the frontier sites.  The frontier occupancy partition of an occupancy
configuration is the set partition of the OCCUPIED frontier sites induced by
occupied-site connectivity through the processed rectangle, under adjacency:
  NN      = von-Neumann (4-neighbour),
  NN+NNN  = king (8-neighbour: adds the two diagonals of each unit face).
A partition class is recorded as (occupied-position tuple, canonical
restricted-growth string) over positions in increasing row order.

Enumeration is an exact BFS over transfer states (column-by-column, union-find
transitions).  A brute-force per-configuration union-find enumeration validates
the BFS wherever w*R <= 16.

Parts
-----
1. main table: for w = 2..6, R = 1..4, both edge sets: number of reachable
   partition classes, whether the NN and NN+NNN class SETS coincide, whether all
   classes are noncrossing (linear frontier), brute-force cross-check.
2. convention reconciliation: the delayed-closing / fully-occupied-frontier
   variant (the convention under which #675 quotes 28 states at w=5 over 4 rows
   and 66 at w=6 over 3 rows).
3. cluster-footprint witness: the recorded converse-false fact (diagonal relay
   footprints that NN cannot make with the same occupied set).

No floating point; exact integers throughout.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

# ---------------------------------------------------------------- utilities


def canonical_rgs(labels) -> tuple[int, ...]:
    mapping: dict[int, int] = {}
    out: list[int] = []
    for lab in labels:
        if lab not in mapping:
            mapping[lab] = len(mapping)
        out.append(mapping[lab])
    return tuple(out)


def canonical_rgs_keep_vacant(labels) -> list[int]:
    mapping: dict[int, int] = {}
    out: list[int] = []
    for lab in labels:
        if lab < 0:
            out.append(-1)
        else:
            if lab not in mapping:
                mapping[lab] = len(mapping)
            out.append(mapping[lab])
    return out


def is_noncrossing(rgs: tuple[int, ...]) -> bool:
    """Linear-frontier test: False iff blocks alternate a<b<c<d."""
    w = len(rgs)
    for a in range(w):
        for b in range(a + 1, w):
            if rgs[a] == rgs[b]:
                continue
            for c in range(b + 1, w):
                if rgs[c] != rgs[a]:
                    continue
                for d in range(c + 1, w):
                    if rgs[d] == rgs[b]:
                        return False
    return True


# ------------------------------------------------------- honest-convention BFS


def column_edges(w: int, king: bool):
    """Within-frontier-column edges.  King adjacency is l_infinity distance 1,
    so within one column only consecutive sites are adjacent for BOTH edge
    sets.  (The distance-2 bonds along the frontier that diagonals create go
    THROUGH the previous column and are handled by `between_edges` + the
    previous state, exactly as in the row-by-row build of PR #645 Q2.)"""
    return [(y, y + 1) for y in range(w - 1)]


def between_offsets(king: bool):
    return [0, 1, -1] if king else [0]


def step_state(w: int, state, mask: int, king: bool):
    """Append one column with occupancy `mask` to state (labels, -1 = vacant).

    The new state records the partition of the NEW column's occupied sites,
    including the new column's own within-column bonds (honest convention)."""
    parent: dict[tuple[int, int], tuple[int, int]] = {}

    def find(v):
        while parent[v] != v:
            parent[v] = parent[parent[v]]
            v = parent[v]
        return v

    def union(u, v):
        ru, rv = find(u), find(v)
        if ru != rv:
            parent[ru] = rv

    for y in range(w):
        if state[y] >= 0:
            parent[(0, y)] = (0, y)
        if (mask >> y) & 1:
            parent[(1, y)] = (1, y)
    by_label: dict[int, tuple[int, int]] = {}
    for y in range(w):
        lab = state[y]
        if lab >= 0:
            if lab in by_label:
                union((0, y), by_label[lab])
            else:
                by_label[lab] = (0, y)
    for y, y2 in column_edges(w, king):
        if (mask >> y) & 1 and (mask >> y2) & 1:
            union((1, y), (1, y2))
    for dy in between_offsets(king):
        for y in range(w):
            y2 = y + dy
            if 0 <= y2 < w and state[y] >= 0 and (mask >> y2) & 1:
                union((0, y), (1, y2))
    old_root_to_lab: dict[tuple[int, int], int] = {}
    for y in range(w):
        if state[y] >= 0:
            old_root_to_lab[find((0, y))] = state[y]
    new_labels = []
    for y in range(w):
        if (mask >> y) & 1:
            r = find((1, y))
            new_labels.append(old_root_to_lab.get(r, ("new", r)))
    rgs = canonical_rgs(new_labels)
    labels = [-1] * w
    occ = [y for y in range(w) if (mask >> y) & 1]
    for pos, lab in zip(occ, rgs):
        labels[pos] = lab
    return tuple(labels)


def reachable_states(w: int, R: int, king: bool) -> set[tuple[int, ...]]:
    states: set[tuple[int, ...]] = {(-1,) * w}
    for _ in range(R):
        nxt: set[tuple[int, ...]] = set()
        for st in states:
            for mask in range(1 << w):
                nxt.add(step_state(w, st, mask, king))
        states = nxt
    return states


def state_to_class(st: tuple[int, ...]):
    occ = tuple(y for y in range(len(st)) if st[y] >= 0)
    rgs = canonical_rgs(st[y] for y in occ)
    return occ, rgs


# --------------------------------------------------- brute-force validation


def brute_force_classes(w: int, R: int, king: bool):
    """Direct per-configuration union-find over the whole rectangle."""
    if w * R > 16:
        return None
    offs = [(0, 1), (1, 0)] + ([(1, 1), (1, -1)] if king else [])
    out = set()
    for cfg in range(1 << (w * R)):
        occ = lambda x, y: (cfg >> (x * w + y)) & 1
        parent = {}

        def find(v):
            while parent[v] != v:
                parent[v] = parent[parent[v]]
                v = parent[v]
            return v

        def union(u, v):
            ru, rv = find(u), find(v)
            if ru != rv:
                parent[ru] = rv

        for x in range(R):
            for y in range(w):
                if occ(x, y):
                    parent[(x, y)] = (x, y)
        for x in range(R):
            for y in range(w):
                if not occ(x, y):
                    continue
                for dx, dy in offs:
                    x2, y2 = x + dx, y + dy
                    if 0 <= x2 < R and 0 <= y2 < w and occ(x2, y2):
                        union((x, y), (x2, y2))
        root_label: dict = {}
        labels = []
        for y in range(w):
            if occ(R - 1, y):
                r = find((R - 1, y))
                if r not in root_label:
                    root_label[r] = len(root_label)
                labels.append(root_label[r])
        out.add(
            (tuple(y for y in range(w) if occ(R - 1, y)), canonical_rgs(labels))
        )
    return out


# ------------------------------------- delayed-closing / full-frontier variant


def reachable_delay_full(w: int, R: int, king: bool):
    """Delayed-closing convention (the one under which #675's 28/66 arise).

    Row-by-row build; when a row is added, the PREVIOUS row is closed (its
    within-row bonds applied) first, then the vertical (dy=0) and, for king,
    diagonal (dy=+-1) bonds between prev and new are applied; the NEW row's own
    within-row bonds stay pending.  Only configurations whose FINAL row is
    fully occupied are accepted, so classes are partitions of the full set [w]."""
    full = (1 << w) - 1
    states: set[tuple[int, ...]] = {(-1,) * w}
    for step in range(R):
        nxt: set[tuple[int, ...]] = set()
        for st in states:
            for mask in range(1 << w):
                if step == R - 1 and mask != full:
                    continue
                # close previous row
                parent = {y: y for y in range(w) if st[y] >= 0}

                def find(v):
                    while parent[v] != v:
                        parent[v] = parent[parent[v]]
                        v = parent[v]
                    return v

                for y in range(w - 1):
                    if st[y] >= 0 and st[y + 1] >= 0:
                        ry, ry1 = find(y), find(y + 1)
                        if ry != ry1:
                            parent[ry] = ry1
                root_lab = {}
                for y in range(w):
                    if st[y] >= 0:
                        root_lab[find(y)] = st[y]
                prev = canonical_rgs_keep_vacant(
                    [-1 if st[y] < 0 else root_lab[find(y)] for y in range(w)]
                )
                # add new row with vertical + diagonal bonds, no own bonds
                parent2: dict[tuple[int, int], tuple[int, int]] = {}

                def find2(v):
                    while parent2[v] != v:
                        parent2[v] = parent2[parent2[v]]
                        v = parent2[v]
                    return v

                def union2(u, v):
                    ru, rv = find2(u), find2(v)
                    if ru != rv:
                        parent2[ru] = rv

                for y in range(w):
                    if prev[y] >= 0:
                        parent2[(0, y)] = (0, y)
                    if (mask >> y) & 1:
                        parent2[(1, y)] = (1, y)
                by_label = {}
                for y in range(w):
                    lab = prev[y]
                    if lab >= 0:
                        if lab in by_label:
                            union2((0, y), by_label[lab])
                        else:
                            by_label[lab] = (0, y)
                for dy in between_offsets(king):
                    for y in range(w):
                        y2 = y + dy
                        if 0 <= y2 < w and prev[y] >= 0 and (mask >> y2) & 1:
                            union2((0, y), (1, y2))
                old_root = {}
                for y in range(w):
                    if prev[y] >= 0:
                        old_root[find2((0, y))] = prev[y]
                new_labels = []
                for y in range(w):
                    if (mask >> y) & 1:
                        r = find2((1, y))
                        new_labels.append(old_root.get(r, ("new", r)))
                rgs = canonical_rgs(new_labels)
                labels = [-1] * w
                occ = [y for y in range(w) if (mask >> y) & 1]
                for pos, lab in zip(occ, rgs):
                    labels[pos] = lab
                nxt.add(tuple(labels))
        states = nxt
    return {state_to_class(st)[1] for st in states}


# ------------------------------------------------------ witness (converse)


def cluster_footprint_witness() -> dict:
    """Converse-false at the CLUSTER level: on a 2x3 rectangle (depth 2,
    width 3), the occupied set {(1,0), (0,1), (1,2)} makes one king cluster
    touching frontier positions {0,2} with position 1 vacant; with the same
    three occupied sites the NN graph gives three singletons, so NN cannot
    produce that footprint.  (The PARTITION {0,2} is still NN-reachable via a
    longer path, which is why the class sets coincide.)"""
    occ = {(1, 0), (0, 1), (1, 2)}
    king_offs = [(1, 0), (0, 1), (-1, 0), (0, -1), (1, 1), (1, -1), (-1, 1), (-1, -1)]
    nn_offs = [(1, 0), (0, 1), (-1, 0), (0, -1)]

    def components(offs):
        parent = {p: p for p in occ}

        def find(v):
            while parent[v] != v:
                parent[v] = parent[parent[v]]
                v = parent[v]
            return v

        for (x, y) in occ:
            for dx, dy in offs:
                p2 = (x + dx, y + dy)
                if p2 in occ:
                    r1, r2 = find((x, y)), find(p2)
                    if r1 != r2:
                        parent[r1] = r2
        comps: dict[tuple, set] = {}
        for p in occ:
            comps.setdefault(find(p), set()).add(p)
        return list(comps.values())

    king_comps = components(king_offs)
    nn_comps = components(nn_offs)
    king_footprints = [
        tuple(y for (x, y) in sorted(c) if x == 1) for c in king_comps
    ]
    return {
        "occupied": sorted(occ),
        "king_components": [sorted(c) for c in king_comps],
        "king_frontier_footprints": sorted(f for f in king_footprints if f),
        "nn_components": [sorted(c) for c in nn_comps],
        "nn_frontier_footprints": sorted(
            f
            for f in (
                tuple(y for (x, y) in sorted(c) if x == 1) for c in nn_comps
            )
            if f
        ),
    }


# ---------------------------------------------------------------- main


def main() -> None:
    global W
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    rows = []
    for w in range(2, 7):
        for R in range(1, 5):
            cl_nnn = {state_to_class(s) for s in reachable_states(w, R, True)}
            cl_nn = {state_to_class(s) for s in reachable_states(w, R, False)}
            bf_nnn = brute_force_classes(w, R, True)
            bf_nn = brute_force_classes(w, R, False)
            rows.append(
                dict(
                    w=w,
                    R=R,
                    n_classes_nn=len(cl_nn),
                    n_classes_nnn=len(cl_nnn),
                    sets_equal=bool(cl_nnn == cl_nn),
                    all_noncrossing_nn=all(is_noncrossing(r) for _, r in cl_nn),
                    all_noncrossing_nnn=all(is_noncrossing(r) for _, r in cl_nnn),
                    brute_force_match_nn=None if bf_nn is None else bool(bf_nn == cl_nn),
                    brute_force_match_nnn=None if bf_nnn is None else bool(bf_nnn == cl_nnn),
                )
            )
            print(
                f"w={w} R={R}: NN={len(cl_nn):4d} NNN={len(cl_nnn):4d} "
                f"sets_equal={cl_nnn == cl_nn} "
                f"noncrossing={rows[-1]['all_noncrossing_nn']}/"
                f"{rows[-1]['all_noncrossing_nnn']} "
                f"bf={rows[-1]['brute_force_match_nn']}/{rows[-1]['brute_force_match_nnn']}"
            )
            if args.out is not None:
                d = args.out / f"w{w}_R{R}"
                d.mkdir(parents=True, exist_ok=True)
                with open(d / "nn.txt", "w") as f:
                    for occ, rgs in sorted(cl_nn):
                        f.write(f"{occ} {rgs}\n")
                with open(d / "nnn.txt", "w") as f:
                    for occ, rgs in sorted(cl_nnn):
                        f.write(f"{occ} {rgs}\n")

    print("\nClass-count sequences (w=2..6), honest convention:")
    for R in range(1, 5):
        seq_nn = [r["n_classes_nn"] for r in rows if r["R"] == R]
        seq_nnn = [r["n_classes_nnn"] for r in rows if r["R"] == R]
        print(f"  R={R}: NN={seq_nn}  NN+NNN={seq_nnn}")

    print("\nDelayed-closing / full-frontier variant (class counts, partitions of [w]):")
    delay_rows = []
    for w in range(2, 7):
        for R in range(1, 5):
            f_nn = reachable_delay_full(w, R, False)
            f_nnn = reachable_delay_full(w, R, True)
            delay_rows.append(
                dict(w=w, R=R, full_nn=len(f_nn), full_nnn=len(f_nnn),
                     sets_equal=bool(f_nn == f_nnn),
                     all_noncrossing_nn=all(is_noncrossing(r) for r in f_nn),
                     all_noncrossing_nnn=all(is_noncrossing(r) for r in f_nnn))
            )
            print(
                f"  w={w} R={R}: NN={len(f_nn):3d} NNN={len(f_nnn):3d} "
                f"eq={f_nn == f_nnn} nc={delay_rows[-1]['all_noncrossing_nn']}/"
                f"{delay_rows[-1]['all_noncrossing_nnn']}"
            )

    wit = cluster_footprint_witness()
    print("\nCluster-footprint witness:", json.dumps(wit))

    if args.out is not None:
        with open(args.out / "summary.json", "w") as f:
            json.dump(
                dict(main_table=rows, delay_full_table=delay_rows, witness=wit),
                f,
                indent=1,
            )


if __name__ == "__main__":
    main()

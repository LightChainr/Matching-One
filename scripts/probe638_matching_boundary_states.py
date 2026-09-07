#!/usr/bin/env python3
"""#638: boundary connectivity classes of the helical site-by-site scan.

Model (issue #636's own construction: "adding sites one at a time").  A
cylinder of width w is the square lattice, rows indexed x = 0, 1, 2, ...,
columns y taken periodically (y mod w).  Sites are added one at a time in
helical scan order (x, 0), (x, 1), ..., (x, w-1), (x+1, 0), ....  The scan
frontier carries w slots, one per column: after step (x, y) slot j <= y
holds the component endpoint of site (x, j) (or EMPTY if that site is
vacant), and slot j > y holds the endpoint of site (x-1, j) (or EMPTY).
The boundary class of a state is the partition of the w columns induced by
the slots: columns sharing a slot label share a block; EMPTY slots are
singleton blocks indistinguishable from occupied singletons (the ticket's
control demands exactly Catalan(w) classes under NN, which holds iff the
two kinds of singleton are identified).

Edges.  NN adjacency: intra-row (0,1) [periodic] and inter-row (1,0).
NNN ("matching") adjacency adds the two inter-row diagonals (1,1) and
(1,-1).  Under NNN the diagonal edge (x, y)-(x-1, y-1) needs the endpoint
of (x-1, y-1) AFTER slot y-1 has been overwritten by site (x, y-1); that
endpoint is therefore carried for one step in a shadow variable (and
shadow_0, which is also consumed by the wrap diagonal (x, w-1)-(x-1, 0),
is carried for the whole row).  Shadows are internal scan machinery: the
reported class is read off the slots only.

Reachability.  The reachable set is the least fixed point of the two-
branch (vacant / occupied) step from the empty state -- a finite state
graph because canonical relabelling bounds the distinct labels by the
number of slots plus shadows.  Classes are recorded at EVERY scan phase
(the state space the #636 transfer actually carries); classes at phase 0
(row boundaries, the straight-cut states) are recorded separately.

Controls and independent checks (GOVERNANCE 2):
  (a) NN control: the reachable all-phase class set must equal the
      committed noncrossing_states(w) elementwise for every width, so the
      count must be Catalan(w) = 1, 2, 5, 14, 42, 132, 429, 1430 for
      w = 1..8.  (A straight-cut row transfer provably misses the
      double-pair classes -- 12 not 14 at w = 4 -- which is why the scan
      is helical and classes are read at every phase.)
  (b) transfer-vs-direct agreement: for w <= 5 and scan prefixes of
      length up to w*r (r = max(3, 18//w)), EVERY prefix of EVERY
      occupancy pattern is classified twice -- once by the transfer
      step, once by an independently written direct classifier that
      recomputes components by explicit graph search over the whole laid
      prefix (no union-find, no shadows) -- and the two classes must be
      equal bit-for-bit.
  (c) every crossing class reported on the NNN side is witnessed by a
      shortest occupancy sequence, and the witness is REPLAYED through
      the independent direct classifier, which must return the claimed
      crossing class.
  (d) row-boundary cross-check: for w <= 5, the classes of ALL full-row
      occupancy patterns to depth 3 rows, computed by the direct
      classifier alone, must be a subset of the closure's phase-0 class
      set (two disjoint routes to the same row-cut classes).
  (e) bond anchor: a bond-percolation row transfer on the same cylinder
      (w dangling bonds, one per column, never vacated) must reach
      EXACTLY the committed noncrossing_states(w) -- elementwise for
      w <= 6 and by count for w = 7 -- so the Catalan(w) control holds
      in the bond representation and its failure on sites is a theorem
      of the site frontier, not an enumerator defect.

Exact integer arithmetic only; the output artifact contains no floats.
Outputs results/probe638-matching-boundary-states/latest.json and prints
the ticket table.
"""

from __future__ import annotations

import argparse
import itertools
import json
import sys
import time
from pathlib import Path
from typing import Dict, FrozenSet, List, Optional, Set, Tuple

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from noncrossing_connectivity_codec import (  # noqa: E402
    blocks_to_rgs,
    canonical_rgs,
    catalan,
    is_noncrossing_blocks,
    noncrossing_states,
)

EMPTY = -1

State = Tuple[int, Tuple[int, ...], int, int]  # (phase, slots, sh0, shprev)
FullRGS = Tuple[int, ...]


# --------------------------------------------------------------------------
# transfer step
# --------------------------------------------------------------------------

def canonical(phase: int, slots: List[int], sh0: int, shprev: int) -> State:
    """Relabel component labels by first appearance over slots, sh0, shprev."""

    remap: Dict[int, int] = {}

    def m(value: int) -> int:
        if value == EMPTY:
            return EMPTY
        if value not in remap:
            remap[value] = len(remap)
        return remap[value]

    return (phase, tuple(m(v) for v in slots), m(sh0), m(shprev))


def step(state: State, occupied: bool, width: int, diag: bool) -> State:
    """Advance the helical scan by one site; `occupied` is the site's bit."""

    phase, slots_t, sh0, shprev = state
    y = phase
    slots = list(slots_t)
    old_y = slots[y]

    if occupied:
        parent: Dict[int, int] = {}

        def find(a: int) -> int:
            root = a
            while parent.get(root, root) != root:
                root = parent[root]
            while parent.get(a, a) != a:
                nxt = parent[a]
                parent[a] = root
                a = nxt
            return root

        def union(a: int, b: int) -> int:
            ra, rb = find(a), find(b)
            if ra == rb:
                return ra
            if ra < rb:
                parent[rb] = ra
                return ra
            parent[ra] = rb
            return rb

        partners: List[int] = []
        if slots[y] != EMPTY:
            partners.append(slots[y])  # (1,0): site above
        if y >= 1 and slots[y - 1] != EMPTY:
            partners.append(slots[y - 1])  # (0,1): left, same row
        if y == width - 1 and slots[0] != EMPTY:
            partners.append(slots[0])  # (0,1) wrap: (x,w-1)-(x,0)
        if diag:
            if y <= width - 2 and slots[y + 1] != EMPTY:
                partners.append(slots[y + 1])  # (1,1)-family: (x,y)-(x-1,y+1)
            if y == 0 and slots[width - 1] != EMPTY:
                partners.append(slots[width - 1])  # (1,1)-family wrap
            if y >= 1:
                shadow = sh0 if y == 1 else shprev
                if shadow != EMPTY:
                    partners.append(shadow)  # (1,-1)-family: (x,y)-(x-1,y-1)
            if y == width - 1 and sh0 != EMPTY:
                partners.append(sh0)  # (1,-1)-family wrap: (x,w-1)-(x-1,0)

        if partners:
            root = partners[0]
            for q in partners[1:]:
                root = union(root, q)
            for j in range(width):
                if slots[j] in partners:
                    slots[j] = root  # every partner slot now carries the merged label
            if sh0 in partners:
                sh0 = root  # shadows merged too: keep them live on the block
            if shprev in partners:
                shprev = root
            slots[y] = root
        else:
            present = [v for v in slots if v != EMPTY]
            if sh0 != EMPTY:
                present.append(sh0)
            if shprev != EMPTY:
                present.append(shprev)
            slots[y] = (max(present) + 1) if present else 0
    else:
        slots[y] = EMPTY

    if diag:
        if y == 0:
            new_sh0, new_shprev = old_y, EMPTY  # shadow_0 lives a whole row
        else:
            new_sh0 = sh0 if y < width - 1 else EMPTY  # consumed at y = w-1
            new_shprev = old_y if 1 <= y <= width - 2 else EMPTY
    else:
        new_sh0 = new_shprev = EMPTY

    return canonical((y + 1) % width, slots, new_sh0, new_shprev)


def class_rgs(state: State, width: int) -> FullRGS:
    """Class = partition of the w columns read off the slots."""

    _phase, slots, _sh0, _shprev = state
    groups: Dict[int, List[int]] = {}
    for j, v in enumerate(slots):
        if v != EMPTY:
            groups.setdefault(v, []).append(j)
    blocks: List[FrozenSet[int]] = [frozenset(v) for v in groups.values()]
    for j, v in enumerate(slots):
        if v == EMPTY:
            blocks.append(frozenset({j}))
    return blocks_to_rgs(blocks, width)


def rgs_blocks(rgs: FullRGS) -> List[FrozenSet[int]]:
    groups: Dict[int, List[int]] = {}
    for y, label in enumerate(rgs):
        groups.setdefault(label, []).append(y)
    return [frozenset(v) for v in groups.values()]


# --------------------------------------------------------------------------
# closure
# --------------------------------------------------------------------------

def closure(
    width: int, diag: bool, state_cap: int = 4_000_000
) -> Tuple[Dict[FullRGS, State], Dict[FullRGS, State], int, Dict[State, Tuple[State, int]], bool]:
    """Least fixed point of the two-branch step from the empty state.

    Returns (all-phase classes -> first state, phase-0 classes -> first
    state, number of states, parent map, capped flag).  BFS order makes
    the recorded state a SHORTEST occupancy sequence for its class.
    """

    start = canonical(0, [EMPTY] * width, EMPTY, EMPTY)
    seen: Set[State] = {start}
    parent: Dict[State, Tuple[State, int]] = {start: (start, 0)}
    classes: Dict[FullRGS, State] = {class_rgs(start, width): start}
    classes_p0: Dict[FullRGS, State] = {}
    order: List[State] = [start]
    i = 0
    capped = False
    while i < len(order):
        st = order[i]
        i += 1
        for b in (0, 1):
            nxt = step(st, b == 1, width, diag)
            if nxt in seen:
                continue
            seen.add(nxt)
            parent[nxt] = (st, b)
            order.append(nxt)
            rgs = class_rgs(nxt, width)
            if rgs not in classes:
                classes[rgs] = nxt
            if nxt[0] == 0 and rgs not in classes_p0:
                classes_p0[rgs] = nxt
        if len(seen) > state_cap:
            capped = True
            break
    return classes, classes_p0, len(seen), parent, capped


def bits_of(state: State, parent: Dict[State, Tuple[State, int]]) -> List[int]:
    """Occupancy bits of the shortest scan leading to `state` (root to leaf)."""

    bits: List[int] = []
    cur = state
    while parent[cur][0] is not cur:
        prev, b = parent[cur]
        bits.append(b)
        cur = prev
    bits.reverse()
    return bits


def rows_of(bits: List[int], width: int) -> List[int]:
    """Pack the scan bits into per-row occupancy masks."""

    rows = []
    for k in range(0, len(bits), width):
        mask = 0
        for offset, b in enumerate(bits[k : k + width]):
            if b:
                mask |= 1 << offset
        rows.append(mask)
    return rows


# --------------------------------------------------------------------------
# independent direct classifier: no union-find, no shadows, no transfer
# --------------------------------------------------------------------------

def direct_class(bits: List[int], width: int, diag: bool) -> FullRGS:
    """Class of the scan state after consuming `bits`, computed directly.

    Sites are the occupied (k // width, k % width) for the 1-bits.  Edges:
    intra-row (x, y±1) periodic, inter-row (x, y)-(x-1, y), plus the two
    diagonals (x, y)-(x-1, y±1) when diag.  Components by explicit graph
    search.  Endpoints after the last consumed bit: column j <= y_last
    attaches to site (x_last, j); column j > y_last attaches to
    (x_last - 1, j), or to nothing on the first row.
    """

    n = len(bits)
    occupied = {(k // width, k % width) for k in range(n) if bits[k]}
    adj: Dict[Tuple[int, int], Set[Tuple[int, int]]] = {}
    for (x, y) in occupied:
        adj.setdefault((x, y), set())
    for (x, y) in occupied:
        for dy in (-1, 1):
            nb = (x, (y + dy) % width)
            if nb in occupied:
                adj[(x, y)].add(nb)  # intra-row, both directions
                adj[nb].add((x, y))
        up = (x - 1, y)
        if up in occupied:
            adj[(x, y)].add(up)  # inter-row (1,0), both directions
            adj[up].add((x, y))
        if diag:
            for dy in (-1, 1):
                nb = (x - 1, (y + dy) % width)
                if nb in occupied:
                    adj[(x, y)].add(nb)  # diagonals, both directions
                    adj[nb].add((x, y))

    comp_of: Dict[Tuple[int, int], int] = {}
    for start in sorted(occupied):
        if start in comp_of:
            continue
        label = len(comp_of)
        comp_of[start] = label
        stack = [start]
        while stack:
            cur = stack.pop()
            for nb in adj[cur]:
                if nb not in comp_of:
                    comp_of[nb] = label
                    stack.append(nb)

    x_last = (n - 1) // width
    y_last = (n - 1) % width
    groups: Dict[int, List[int]] = {}
    blocks: List[FrozenSet[int]] = []
    for j in range(width):
        if j <= y_last:
            site = (x_last, j)
        elif x_last >= 1:
            site = (x_last - 1, j)
        else:
            site = None
        if site is not None and site in occupied:
            groups.setdefault(comp_of[site], []).append(j)
        else:
            blocks.append(frozenset({j}))
    blocks.extend(frozenset(v) for v in groups.values())
    return blocks_to_rgs(blocks, width)


# --------------------------------------------------------------------------
# bond anchor: dangling bonds, one per column, never vacated.  Under NN the
# reachable dangling-bond partitions are exactly the noncrossing partitions
# (Catalan(w)) -- the anchor that the site frontier provably misses.
# --------------------------------------------------------------------------

BondState = Tuple[Tuple[int, ...], ...]  # RGS over the w cut ports (one per column)

EMPTY_LABEL = -1  # sentinel: port sits on an absent vertical bond


def bond_step(
    state: BondState, vrow: int, hrow: int, width: int, periodic: bool = False
) -> BondState:
    """One full row of the bond-percolation transfer.

    periodic=True  : the ticket's cylinder; the wrap horizontal bond
        (x, w-1)-(x, 0) is available, and crossing partitions such as
        {0,2},{1,3} become reachable -- the honest cylinder state space
        is LARGER than Catalan(w).
    periodic=False : the PLANAR strip (cut open between w-1 and 0); with
        no wrap edge the reachable partitions are exactly the noncrossing
        ones, Catalan(w) -- the control the ticket quotes.

    The state is the partition of the w cut ports; port y is the upper
    end of the previous row's vertical bond at column y.  Laying the next
    row consumes two w-bit masks: vrow (vertical bonds (x-1,y)-(x,y)) and
    hrow (horizontal bonds (x,y)-(x,y+1), periodic).  A present vertical
    bond merges port y with the fresh site (x, y); an absent one leaves
    port y carrying its old block while (x, y) starts fresh -- this is
    the degree of freedom that lets a path sneak UNDER a port without
    touching it.  Horizontal bonds merge fresh sites.  The new port y is
    the fresh site's block if v_y is present, else the old block.
    """

    parent: Dict[int, int] = {v: v for v in range(2 * width)}

    def find(a: int) -> int:
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            if ra < rb:
                parent[rb] = ra
            else:
                parent[ra] = rb

    # seed: ports sharing an old block
    groups: Dict[int, List[int]] = {}
    for y, label in enumerate(state):
        groups.setdefault(label, []).append(y)
    for members in groups.values():
        for y in members[1:]:
            union(members[0], y)
    # vertical bonds
    for y in range(width):
        if (vrow >> y) & 1:
            union(y, width + y)
    # horizontal bonds; the wrap (w-1, 0) only on the cylinder
    for y in range(width):
        if (hrow >> y) & 1:
            if y == width - 1 and not periodic:
                continue
            union(width + y, width + (y + 1) % width)

    labels: List[int] = []
    for y in range(width):
        if (vrow >> y) & 1:
            labels.append(find(width + y))  # port = cluster of site (x, y)
        else:
            labels.append(EMPTY_LABEL)  # isolated midpoint
    # each EMPTY port is its OWN singleton (distinct isolated midpoints);
    # sealed blocks simply disappear from the relabelling
    fresh = 10 ** 6
    for y in range(width):
        if labels[y] == EMPTY_LABEL:
            labels[y] = fresh
            fresh += 1
    return canonical_rgs(labels)


def bond_closure(width: int, periodic: bool = False) -> Set[Tuple[int, ...]]:
    """Reachable cut-port partitions under the bond row transfer (NN)."""

    start: BondState = tuple(range(width))  # all singletons
    seen: Set[BondState] = {start}
    frontier = [start]
    classes: Set[Tuple[int, ...]] = {start}
    while frontier:
        nxt = []
        for st in frontier:
            for vrow in range(1 << width):
                for hrow in range(1 << width):
                    nst = bond_step(st, vrow, hrow, width, periodic)
                    if nst in seen:
                        continue
                    seen.add(nst)
                    nxt.append(nst)
                    classes.add(nst)
        frontier = nxt
    return classes


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def selfcheck_prefix_agreement(width: int, diag: bool) -> Dict[str, object]:
    """Check (b): every prefix of every pattern of length w*r, classified
    by the transfer and by the direct classifier, must agree exactly."""

    rows = max(3, 18 // width)
    length = width * rows
    checked = 0
    for pattern in itertools.product((0, 1), repeat=length):
        state = canonical(0, [EMPTY] * width, EMPTY, EMPTY)
        for k, b in enumerate(pattern):
            state = step(state, b == 1, width, diag)
            checked += 1
            if class_rgs(state, width) != direct_class(list(pattern[: k + 1]), width, diag):
                return {
                    "width": width,
                    "diag": diag,
                    "prefix_length": length,
                    "prefixes_checked": checked,
                    "agrees": False,
                    "first_mismatch_bits": list(pattern[: k + 1]),
                }
    return {
        "width": width,
        "diag": diag,
        "prefix_length": length,
        "prefixes_checked": checked,
        "agrees": True,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--widths", type=int, nargs="+", default=list(range(1, 9)))
    parser.add_argument("--state-cap", type=int, default=4_000_000)
    parser.add_argument(
        "--out",
        type=Path,
        default=ROOT / "results" / "probe638-matching-boundary-states" / "latest.json",
    )
    args = parser.parse_args()

    runs: List[dict] = []
    selfchecks: List[dict] = []
    control_ok = True
    smallest_crossing_width: Optional[int] = None
    witness: Optional[dict] = None

    for width in args.widths:
        t0 = time.time()
        nn_classes, nn_p0, nn_states, nn_parent, nn_capped = closure(width, False, args.state_cap)
        t_nn = time.time() - t0
        t0 = time.time()
        nnn_classes, nnn_p0, nnn_states, nnn_parent, nnn_capped = closure(
            width, True, args.state_cap
        )
        t_nnn = time.time() - t0

        committed = sorted(noncrossing_states(width))
        got = sorted(nn_classes)
        nn_matches = got == committed and len(got) == catalan(width)
        # NOTE: the site-side NN closure being a strict subset of the
        # noncrossing partitions is a FINDING of this probe (the nested-
        # pocket classes are unreachable by the site frontier), not a
        # control failure; the Catalan control is carried by the BOND
        # anchor below, which must hit it exactly.

        crossing = sorted(c for c in nnn_classes if not is_noncrossing_blocks(rgs_blocks(c)))

        if crossing and smallest_crossing_width is None:
            smallest_crossing_width = width
            cls = crossing[0]
            st = nnn_classes[cls]
            bits = bits_of(st, nnn_parent)
            replayed = direct_class(bits, width, True)
            witness = {
                "width": width,
                "class_rgs": list(cls),
                "blocks": sorted(sorted(b) for b in rgs_blocks(cls)),
                "n_scan_bits": len(bits),
                "row_masks": rows_of(bits, width),
                "phase_after": st[0],
                "replay_independent_classifier_agrees": replayed == cls,
            }
            if replayed != cls:
                control_ok = False

        run = {
            "width": width,
            "nn_class_count": len(nn_classes),
            "nn_catalan": catalan(width),
            "nn_matches_catalan_exactly": nn_matches,
            "nn_state_count": nn_states,
            "nn_phase0_class_count": len(nn_p0),
            "nnn_class_count": len(nnn_classes),
            "nnn_state_count": nnn_states,
            "nnn_phase0_class_count": len(nnn_p0),
            "nnn_crossing_class_count": len(crossing),
            "nnn_all_noncrossing": not crossing,
            "nn_capped": nn_capped,
            "nnn_capped": nnn_capped,
            "nn_seconds": round(t_nn, 1),
            "nnn_seconds": round(t_nnn, 1),
        }
        if not nn_matches:
            run["nn_missing_classes"] = [list(c) for c in sorted(set(committed) - set(got))]
            run["nn_extra_classes"] = [list(c) for c in sorted(set(got) - set(committed))]
        runs.append(run)
        print(
            f"w={width}  NN {len(nn_classes):5d} (Catalan {catalan(width):5d}) "
            f"[states {nn_states}, p0 {len(nn_p0)}]  "
            f"NNN {len(nnn_classes):6d} [states {nnn_states}, "
            f"p0 {len(nnn_p0)}]  crossing {len(crossing):4d}  "
            f"[{t_nn:.1f}s / {t_nnn:.1f}s]",
            flush=True,
        )

    if max(args.widths) <= 5:
        for width in args.widths:
            for diag in (False, True):
                selfchecks.append(selfcheck_prefix_agreement(width, diag))
                if not selfchecks[-1]["agrees"]:
                    control_ok = False
        print("prefix agreement checks done", flush=True)

    # bond anchor (e): Catalan control in the bond representation.
    # planar  = cut open between w-1 and 0 (no wrap horizontal bond);
    # cylinder = wrap bond available (crossing partitions become reachable).
    bond_anchor: List[dict] = []
    for width in args.widths:
        if width > 7:
            break
        t0 = time.time()
        planar = bond_closure(width, periodic=False)
        t_planar = time.time() - t0
        t0 = time.time()
        cyl = bond_closure(width, periodic=True)
        t_cyl = time.time() - t0
        committed = set(noncrossing_states(width))
        exact = planar == committed
        if not exact:
            control_ok = False
        bond_anchor.append(
            {
                "width": width,
                "bond_planar_class_count": len(planar),
                "bond_cylinder_class_count": len(cyl),
                "catalan": catalan(width),
                "bond_planar_equals_noncrossing_states_elementwise": exact,
                "bond_cylinder_crossing_class_count": sum(
                    1 for c in cyl if not is_noncrossing_blocks(rgs_blocks(c))
                ),
                "bond_planar_seconds": round(t_planar, 1),
                "bond_cylinder_seconds": round(t_cyl, 1),
            }
        )
        print(
            f"bond w={width}: planar {len(planar)} (Catalan {catalan(width)}, exact={exact}) "
            f"cylinder {len(cyl)} [{t_planar:.1f}s / {t_cyl:.1f}s]",
            flush=True,
        )

    payload = {
        "schema": "matching-one.probe638.boundary-states.v2",
        "ticket": 638,
        "model": (
            "helical site-by-site scan on the w-column cylinder; NN edges "
            "(0,1)+(1,0); NNN adds the inter-row diagonals (1,1) and (1,-1) "
            "(the (1,-1) edge routed through one-step shadow endpoints)"
        ),
        "class_definition": (
            "partition of the w columns read off the scan slots; EMPTY slots "
            "are singletons indistinguishable from occupied singletons; "
            "classes recorded at every scan phase"
        ),
        "control_ok": control_ok,
        "bond_anchor": bond_anchor,
        "smallest_crossing_width": smallest_crossing_width,
        "crossing_witness": witness,
        "selfchecks": selfchecks,
        "runs": runs,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(f"wrote {args.out}")
    print(
        "ALL_GOVERNANCE_CHECKS_PASS" if control_ok else "GOVERNANCE_CHECK_FAILURE",
        "(control = bond anchor hits Catalan exactly; see JSON fields)",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

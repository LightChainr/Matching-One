#!/usr/bin/env python3
"""Independent audit of the bond ambient rank used by PR #628 (issue #635).

Suspicion under audit: ``bond_ambient_rank`` in
``scripts/probe_invariant_shape/exact_rank_census.py`` (PR #628) computes the
winding of each fundamental cycle as ``tree_lift(u, v) + (dx, dy)``, i.e. it
*adds* the stored edge displacement to the tree-path lift displacement.  The
correct winding of the fundamental cycle (tree path u->v, then the stored edge
back v->u, whose lift displacement is -(dx, dy)) is

    (ax, ay) - (dx, dy).

Because the stored displacement is not negated, the computed "winding" depends
on the chosen spanning tree, so the function is not a function of the edge set
at all: the same subgraph can score rank 1 or 2 depending on edge order (first
observed on the L=3 mask 63 dual graph: ascending order -> 1, e628 iteration
order -> 2).  Every (r_b, r_w) number PR #628 published for the bond census
(``dual_fail = 118133`` of 262144 and its rank-pair table) was computed with
this function, so the audit must decide which of those numbers survive.

Three rank implementations:

1. ``bond_ambient_rank_e628`` -- the PR #628 function, copied VERBATIM.  The
   census it produces must reproduce the published 118133 / pair table
   bit-for-bit; that closes the causal chain bug -> published numbers.
2. ``bond_ambient_rank_fixed`` -- same algorithm with the sign corrected
   (``ax - dx``).  This is the candidate truth used for the full census.
3. ``winding_image_rank_linalg`` -- an independent algorithm that never builds
   a spanning tree: for an edge subset S the winding image is
   ``w(ker B_S)`` with B_S the incidence matrix and w the Z^2 winding map;
   its rank is ``rank([B_S; W]) - rank(B_S)``, both computed by exact rational
   Gaussian elimination.  Certified on a 20000-mask random sample plus
   structured samples; agreement there certifies implementation 2, which then
   carries the full 2^18 census.

Exactness: every rank here is an exact integer/Fraction computation; the only
floats in this file are none.  Published #628 bond numbers reproduced for
comparison:  dual_fail = 118133; pair_counts (r_b, r_w) ->
(0,2):51704, (1,2):49263, (0,1):4397, (1,1):41839, (1,0):5633,
(2,2):12183, (2,1):46657, (2,0):50468.
"""

from __future__ import annotations

import json
import multiprocessing as mp
import random
import sys
from collections import deque
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from square_bond_kappa3 import square_bond_pairs  # noqa: E402

NB = 18  # 2 L^2, L = 3
NCONF = 1 << NB
PUBLISHED_E628 = {
    "dual_fail": 118133,
    "pair_counts": {"0,2": 51704, "1,2": 49263, "0,1": 4397, "1,1": 41839,
                    "1,0": 5633, "2,2": 12183, "2,1": 46657, "2,0": 50468},
}

# ---------------------------------------------------------------- bond table

_PAIRS = square_bond_pairs(3)
BOND_UV = [(_PAIRS[i].primal[0], _PAIRS[i].primal[1]) for i in range(NB)]
BOND_STEP = [(_PAIRS[i].primal[2], _PAIRS[i].primal[3]) for i in range(NB)]
NV = max(max(u, v) for (u, v) in BOND_UV) + 1
# Crossing-dual map, exactly as PR #628 resolves it: the dual edge of pair i
# coincides (under the crossing-dual identification) with a unique primal
# edge, whose index is dual_to_primal[i].
_dual_index = {}
for _i in range(NB):
    _du, _dv = _PAIRS[_i].dual[0], _PAIRS[_i].dual[1]
    _key = (min(_du, _dv), max(_du, _dv))
    _dual_index[_i] = _key
_primal_index = {}
for _i in range(NB):
    _u, _v = BOND_UV[_i]
    _primal_index[(min(_u, _v), max(_u, _v))] = _i
DUAL_TO_PRIMAL = tuple(_primal_index[_dual_index[_i]] for _i in range(NB))
assert len(set(DUAL_TO_PRIMAL)) == NB  # derangement of bond indices


# ------------------------------------------- 1. PR #628 rank, VERBATIM copy

def bond_ambient_rank_e628(edges: list[tuple[int, int, int, int]]) -> int:
    """Verbatim copy of exact_rank_census.bond_ambient_rank (PR #628).

    Kept byte-for-byte (modulo the name) so the census below reproduces the
    published numbers with the published algorithm, bug included.
    """
    if not edges:
        return 0
    n_v = max(max(u, v) for (u, v, _, _) in edges) + 1
    parent = list(range(n_v))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    tree: list[tuple[int, int, int, int]] = []
    nontree: list[tuple[int, int, int, int]] = []
    for (u, v, dx, dy) in edges:
        ru, rv = find(u), find(v)
        if ru != rv:
            parent[ru] = rv
            tree.append((u, v, dx, dy))
        else:
            nontree.append((u, v, dx, dy))

    # Tree adjacency with lift steps.
    from collections import defaultdict
    tadj: dict[int, list[tuple[int, tuple[int, int]]]] = defaultdict(list)
    for (u, v, dx, dy) in tree:
        tadj[u].append((v, (dx, dy)))
        tadj[v].append((u, (-dx, -dy)))

    def tree_lift(src: int, dst: int) -> tuple[int, int]:
        q = deque([src])
        prev: dict[int, tuple[int, tuple[int, int] | None]] = {src: (src, None)}
        while q:
            x = q.popleft()
            if x == dst:
                break
            for (y, step) in tadj[x]:
                if y not in prev:
                    prev[y] = (x, step)
                    q.append(y)
        sx = sy = 0
        x = dst
        while prev[x][0] != x:
            p, step = prev[x]
            sx += step[0]
            sy += step[1]
            x = p
        return sx, sy

    windings: list[tuple[int, int]] = []
    for (u, v, dx, dy) in nontree:
        (ax, ay) = tree_lift(u, v)
        # fundamental cycle: lift path u->v in tree, then the edge v->u with
        # displacement -(dx,dy) as stored directed u->v.
        sx, sy = ax + dx, ay + dy
        if sx or sy:
            windings.append((sx, sy))
    indep: list[tuple[int, int]] = []
    for (a, b) in windings:
        if len(indep) == 0:
            indep.append((a, b))
        elif len(indep) == 1:
            (a0, b0) = indep[0]
            if a0 * b - b0 * a != 0:
                indep.append((a, b))
        else:
            break
    return len(indep)


# ---------------------------------------- 2. sign-corrected rank, same scheme

def bond_ambient_rank_fixed(edges: list[tuple[int, int, int, int]]) -> int:
    """Same algorithm as above with the winding sign corrected.

    The fundamental cycle is tree path u->v followed by the stored edge
    v->u, so its winding is (ax, ay) - (dx, dy); the published function adds
    (dx, dy) instead, which makes the result depend on the spanning tree.
    """
    if not edges:
        return 0
    n_v = max(max(u, v) for (u, v, _, _) in edges) + 1
    parent = list(range(n_v))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    tree: list[tuple[int, int, int, int]] = []
    nontree: list[tuple[int, int, int, int]] = []
    for (u, v, dx, dy) in edges:
        ru, rv = find(u), find(v)
        if ru != rv:
            parent[ru] = rv
            tree.append((u, v, dx, dy))
        else:
            nontree.append((u, v, dx, dy))

    from collections import defaultdict
    tadj: dict[int, list[tuple[int, tuple[int, int]]]] = defaultdict(list)
    for (u, v, dx, dy) in tree:
        tadj[u].append((v, (dx, dy)))
        tadj[v].append((u, (-dx, -dy)))

    def tree_lift(src: int, dst: int) -> tuple[int, int]:
        q = deque([src])
        prev: dict[int, tuple[int, tuple[int, int] | None]] = {src: (src, None)}
        while q:
            x = q.popleft()
            if x == dst:
                break
            for (y, step) in tadj[x]:
                if y not in prev:
                    prev[y] = (x, step)
                    q.append(y)
        sx = sy = 0
        x = dst
        while prev[x][0] != x:
            p, step = prev[x]
            sx += step[0]
            sy += step[1]
            x = p
        return sx, sy

    windings: list[tuple[int, int]] = []
    for (u, v, dx, dy) in nontree:
        (ax, ay) = tree_lift(u, v)
        # CORRECTED: the cycle closes along the stored edge backwards, so the
        # edge contributes -(dx, dy), not +(dx, dy).
        sx, sy = ax - dx, ay - dy
        if sx or sy:
            windings.append((sx, sy))
    indep: list[tuple[int, int]] = []
    for (a, b) in windings:
        if len(indep) == 0:
            indep.append((a, b))
        elif len(indep) == 1:
            (a0, b0) = indep[0]
            if a0 * b - b0 * a != 0:
                indep.append((a, b))
        else:
            break
    return len(indep)


# --------------------------------------- 3. independent exact-linear-algebra

def _gauss_rank(rows: list[list[Fraction]]) -> int:
    """Exact rank by rational Gaussian elimination."""
    m = [row[:] for row in rows]
    nr, nc = len(m), len(m[0]) if m else 0
    r = 0
    for c in range(nc):
        piv = -1
        for i in range(r, nr):
            if m[i][c] != 0:
                piv = i
                break
        if piv < 0:
            continue
        m[r], m[piv] = m[piv], m[r]
        pr = m[r]
        prc = pr[c]
        for i in range(r + 1, nr):
            f = m[i][c]
            if f:
                f = f / prc
                row = m[i]
                for j in range(c, nc):
                    row[j] -= f * pr[j]
        r += 1
    return r


def winding_image_rank_linalg(bond_ids: list[int]) -> int:
    """Winding-image rank without any spanning tree.

    For edge subset S the winding map is w : Q^S -> Q^2, w(c) = sum of the
    lift displacements of the 1-chain c; the winding image of the subgraph is
    w(ker B_S) where B_S is the vertex-by-edge incidence matrix (any
    orientation).  Its rank is

        dim ker B_S - dim ker [B_S; W]  =  rank [B_S; W] - rank B_S,

    and both ranks are computed by exact rational elimination.
    """
    if not bond_ids:
        return 0
    cols = len(bond_ids)
    incidence: list[list[int]] = [[0] * cols for _ in range(NV)]
    for j, i in enumerate(bond_ids):
        u, v = BOND_UV[i]
        incidence[u][j] += 1
        incidence[v][j] -= 1
    winding = [[BOND_STEP[i][0] for i in bond_ids],
               [BOND_STEP[i][1] for i in bond_ids]]
    rank_b = _gauss_rank([[Fraction(x) for x in row] for row in incidence])
    rank_bw = _gauss_rank([[Fraction(x) for x in row]
                           for row in incidence + winding])
    return rank_bw - rank_b


# ------------------------------------------------------------------ census

def _edges_of(occ: list[int], order) -> list[tuple[int, int, int, int]]:
    return [(BOND_UV[i][0], BOND_UV[i][1], BOND_STEP[i][0], BOND_STEP[i][1])
            for i in order if occ[i]]


def _dual_edges_of(occ: list[int], geometric: bool = False
                   ) -> list[tuple[int, int, int, int]]:
    """Occupied dual bonds transported to primal coordinates.

    Two conventions, both audited:

    - ``geometric=False`` (PR #628's code path): iterate i ascending, map the
      dual of pair i to its primal index pi = DUAL_TO_PRIMAL[i], keep pi when
      *primal pi is vacant*.  Because DUAL_TO_PRIMAL is a fixed-point-free
      map on bond indices, this enumerates ``{pi : occ[pi] == 0}`` -- the
      plain bit-complement of the occupied set, NOT the transported geometric
      dual.  This is the convention error #42's note already flagged; it is
      reproduced here verbatim so the causal chain to the published numbers
      is closed.
    - ``geometric=True`` (the correct crossing-dual transport): dual bond i
      is occupied iff primal bond i is vacant, and the occupied dual bond is
      transported to the primal edge ``DUAL_TO_PRIMAL[i]``.  This is
      ``{DUAL_TO_PRIMAL[i] : occ[i] == 0}``, which differs from the
      complement on every configuration (DUAL_TO_PRIMAL is a derangement but
      NOT an involution on L=3).
    """
    out = []
    if geometric:
        for i in range(NB):
            if not occ[i]:
                pi = DUAL_TO_PRIMAL[i]
                out.append((BOND_UV[pi][0], BOND_UV[pi][1],
                            BOND_STEP[pi][0], BOND_STEP[pi][1]))
    else:
        for i in range(NB):
            pi = DUAL_TO_PRIMAL[i]
            if not occ[pi]:
                out.append((BOND_UV[pi][0], BOND_UV[pi][1],
                            BOND_STEP[pi][0], BOND_STEP[pi][1]))
    return out


_DESC = list(range(NB - 1, -1, -1))


def census_chunk(args: tuple[str, str, list[int]]) -> dict:
    which, dual_mode, masks = args
    rank_fn = {"e628": bond_ambient_rank_e628,
               "fixed": bond_ambient_rank_fixed}[which]
    geometric = (dual_mode == "geometric")
    comp: dict[tuple[int, int], list[int]] = {}
    pair_counts: dict[tuple[int, int], int] = {}
    dual_fail = 0
    order_unstable = 0
    for mask in masks:
        occ = [(mask >> i) & 1 for i in range(NB)]
        edges = _edges_of(occ, range(NB))
        dedges = _dual_edges_of(occ, geometric=geometric)
        rb = rank_fn(edges)
        rw = rank_fn(dedges)
        if which == "e628" and not geometric:
            # order-instability probe on the same function: descending order
            if rank_fn(_edges_of(occ, _DESC)) != rb:
                order_unstable += 1
        if rb + rw != 2:
            dual_fail += 1
        pair_counts[(rb, rw)] = pair_counts.get((rb, rw), 0) + 1
        k = sum(occ)
        comp.setdefault((rb, rw), [0] * (NB + 1))[k] += 1
    return {"comp": comp, "pair_counts": pair_counts,
            "dual_fail": dual_fail, "order_unstable": order_unstable,
            "n": len(masks)}


def _merge(target: dict, chunk: dict) -> None:
    target["dual_fail"] += chunk["dual_fail"]
    target["order_unstable"] += chunk["order_unstable"]
    for key, val in chunk["pair_counts"].items():
        target["pair_counts"][key] = target["pair_counts"].get(key, 0) + val
    for key, poly in chunk["comp"].items():
        tp = target["comp"].setdefault(key, [0] * (NB + 1))
        for k in range(NB + 1):
            tp[k] += poly[k]


def run_census(which: str, workers: int, dual_mode: str = "e628-complement"
               ) -> dict:
    out = {"comp": {}, "pair_counts": {}, "dual_fail": 0,
           "order_unstable": 0}
    chunk_size = 2048
    chunks = [(which, dual_mode,
               list(range(a, min(a + chunk_size, NCONF))))
              for a in range(0, NCONF, chunk_size)]
    with mp.get_context("fork").Pool(workers) as pool:
        for chunk in pool.imap_unordered(census_chunk, chunks):
            _merge(out, chunk)
    total = sum(out["pair_counts"].values())
    assert total == NCONF, (which, total)
    return out


# ------------------------------------------------------- independent sampling

def audit_sample(rng: random.Random, n_sample: int) -> dict:
    masks = rng.sample(range(NCONF), n_sample)
    extra = list(range(1024))  # all small-|S| structured configs
    bad = []
    checked = 0
    for mask in masks + extra:
        occ = [(mask >> i) & 1 for i in range(NB)]
        primal_ids = [i for i in range(NB) if occ[i]]
        dual_ids = [DUAL_TO_PRIMAL[i] for i in range(NB) if not occ[DUAL_TO_PRIMAL[i]]]
        # dual_ids as ascending primal indices, unsorted per e628 scheme:
        dual_ids_iter = [DUAL_TO_PRIMAL[i] for i in range(NB) if not occ[i]]
        for ids, expect in ((primal_ids, None), (dual_ids_iter, None)):
            edges = [(BOND_UV[i][0], BOND_UV[i][1],
                      BOND_STEP[i][0], BOND_STEP[i][1]) for i in ids]
            got = bond_ambient_rank_fixed(edges)
            want = winding_image_rank_linalg(list(dict.fromkeys(ids)))
            checked += 1
            if got != want:
                bad.append({"mask": mask, "side": "primal" if expect is None
                            and ids is primal_ids else "dual",
                            "fixed": got, "linalg": want})
    return {"masks_sampled": len(masks), "extra_structured": len(extra),
            "rank_comparisons": checked, "disagreements": bad[:20],
            "n_disagreements": len(bad)}


def main() -> None:
    workers = max(2, (mp.cpu_count() or 4))
    rng = random.Random(635635)

    # Counterexample of record: L=3 mask 63 dual graph.
    occ63 = [(63 >> i) & 1 for i in range(NB)]
    d63 = _dual_edges_of(occ63)
    counterexample = {
        "mask": 63, "side": "dual (transported)",
        "e628_ascending": bond_ambient_rank_e628(d63),
        "e628_descending": bond_ambient_rank_e628(d63[::-1]),
        "fixed": bond_ambient_rank_fixed(d63),
        "linalg": winding_image_rank_linalg(
            [DUAL_TO_PRIMAL[i] for i in range(NB) if not occ63[i]]),
    }
    assert counterexample["e628_ascending"] != counterexample["e628_descending"]
    assert counterexample["e628_ascending"] == 2  # published-function artifact
    assert counterexample["fixed"] == counterexample["linalg"] == 1

    audit = audit_sample(rng, 20000)
    assert audit["n_disagreements"] == 0, audit

    census_e628 = run_census("e628", workers)
    got_counts = {f"{a},{b}": n for (a, b), n in census_e628["pair_counts"].items()}
    assert census_e628["dual_fail"] == PUBLISHED_E628["dual_fail"], \
        (census_e628["dual_fail"],)
    assert got_counts == PUBLISHED_E628["pair_counts"], got_counts

    census_fixed = run_census("fixed", workers, dual_mode="e628-complement")
    census_geo_fixed = run_census("fixed", workers, dual_mode="geometric")
    census_geo_e628 = run_census("e628", workers, dual_mode="geometric")

    # 2x2 attribution: the published dual_fail decomposes over the two
    # independent defects (dual-occupation convention, rank winding sign).
    fail_complement_e628 = census_e628["dual_fail"]     # published: 118133
    fail_complement_fixed = census_fixed["dual_fail"]
    fail_geometric_e628 = census_geo_e628["dual_fail"]
    fail_geometric_fixed = census_geo_fixed["dual_fail"]
    # Exact duality law on the corrected path: r_b + r_w = 2 must hold on
    # every configuration, and per-k the three rank pairs must partition
    # C(18, k); (0,2)<->(2,0) mirror under k -> 18-k, (1,1) symmetric.
    assert fail_geometric_fixed == 0, fail_geometric_fixed
    from math import comb
    kmasses = {pair: dict(enumerate(poly))
               for pair, poly in census_geo_fixed["comp"].items()}
    for k in range(NB + 1):
        s = sum(kmasses.get(p, {}).get(k, 0) for p in kmasses)
        assert s == comb(NB, k), (k, s, comb(NB, k))
    m02 = kmasses.get((0, 2), {})
    m20 = kmasses.get((2, 0), {})
    m11 = kmasses.get((1, 1), {})
    assert all(m02[k] == m20.get(NB - k, 0) for k in m02)
    assert all(m11[k] == m11[NB - k] for k in m11)

    def polys(cens: dict) -> dict:
        comp = cens["comp"]
        get = lambda name: comp.get(name, [0] * (NB + 1))
        P11, P20, P02 = get((1, 1)), get((2, 0)), get((0, 2))
        M = [P20[k] - P02[k] for k in range(NB + 1)]
        F = [P20[k] + P11[k] / 2 for k in range(NB + 1)]
        return {"P11": P11, "P20": P20, "P02": P02, "M": M, "F": F}

    P_fixed = polys(census_fixed)
    M = P_fixed["M"]
    # Exact checks on the corrected census.
    assert M[0] == -1 and M[NB] == 1, (M[0], M[NB])
    assert all(M[k] + M[NB - k] == 0 for k in range(NB + 1)), M
    assert M[NB // 2] == 0
    sym_ok = all(census_fixed["pair_counts"].get((a, b), 0)
                 == census_fixed["pair_counts"].get((b, a), 0)
                 for a in range(3) for b in range(3))
    assert sym_ok

    counts_fixed = {f"{a},{b}": n for (a, b), n
                    in census_fixed["pair_counts"].items()}
    counts_geo_fixed = {f"{a},{b}": n for (a, b), n
                        in census_geo_fixed["pair_counts"].items()}
    counts_geo_e628 = {f"{a},{b}": n for (a, b), n
                       in census_geo_e628["pair_counts"].items()}
    diff_pairs = {k: (PUBLISHED_E628["pair_counts"].get(k, 0),
                      counts_fixed.get(k, 0))
                  for k in sorted(set(PUBLISHED_E628["pair_counts"])
                                  | set(counts_fixed))
                  if PUBLISHED_E628["pair_counts"].get(k, 0)
                  != counts_fixed.get(k, 0)}

    # M coefficients of the corrected geometric-dual census (the one that
    # satisfies exact duality): P20 - P02, exact, antisymmetric.
    comp_geo = census_geo_fixed["comp"]
    getp = lambda name: comp_geo.get(name, [0] * (NB + 1))
    M_geo = [getp((2, 0))[k] - getp((0, 2))[k] for k in range(NB + 1)]
    assert all(M_geo[k] + M_geo[NB - k] == 0 for k in range(NB + 1))
    assert M_geo[0] == -1 and M_geo[NB] == 1

    out = {
        "schema": "matching-one.probe635.bond-rank-audit.v2",
        "scope": "L=3 square bond torus, 2^18 configs, exact integer ranks",
        "suspicion": ("exact_rank_census.bond_ambient_rank adds the edge "
                      "displacement to the tree-lift displacement of the "
                      "fundamental cycle (ax+dx instead of ax-dx); the "
                      "result then depends on edge order, so it is not a "
                      "topological invariant of the subgraph"),
        "counterexample_mask63": counterexample,
        "independent_check": audit,
        "census_e628_verbatim": {
            "dual_fail": census_e628["dual_fail"],
            "pair_counts": got_counts,
            "order_unstable_configs": census_e628["order_unstable"],
            "matches_published": True,
        },
        "census_fixed": {
            "dual_fail": census_fixed["dual_fail"],
            "pair_counts": counts_fixed,
            "components": {f"{a},{b}": [str(x) for x in P_fixed[name]]
                           for (a, b), name in
                           (((1, 1), "P11"), ((2, 0), "P20"), ((0, 2), "P02"))},
            "M_coeffs": [str(x) for x in M],
            "M_quarter": str(sum(Fraction(c) * Fraction(1, 4) ** k
                                 * Fraction(3, 4) ** (NB - k)
                                 for k, c in enumerate(M))),
            "M_half": str(sum(Fraction(c) * Fraction(1, 2) ** NB
                              for c in M)),
            "caveat": ("same complement dual-occupation convention as #628; "
                       "isolates the rank-sign bug only"),
        },
        "census_geometric_dual": {
            "with_fixed_rank": {
                "dual_fail": fail_geometric_fixed,
                "pair_counts": counts_geo_fixed,
                "k_masses": {f"{a},{b}": [str(x) for x in poly]
                             for (a, b), poly in census_geo_fixed["comp"].items()},
                "M_coeffs": [str(x) for x in M_geo],
                "M_half": str(sum(Fraction(c) * Fraction(1, 2) ** NB
                                  for c in M_geo)),
            },
            "with_e628_rank": {
                "dual_fail": fail_geometric_e628,
                "pair_counts": counts_geo_e628,
            },
        },
        "attribution_2x2": {
            "complement + e628 rank": fail_complement_e628,
            "complement + fixed rank": fail_complement_fixed,
            "geometric + e628 rank": fail_geometric_e628,
            "geometric + fixed rank": fail_geometric_fixed,
            "reading": ("the published 118133 duality failures are a "
                        "compound artifact of two independent defects: the "
                        "complement dual-occupation convention (#42's "
                        "known convention error) and the rank winding-sign "
                        "bug.  Under the geometric crossing-dual transport "
                        "with the corrected rank, r_b + r_w = 2 holds on "
                        "all 262144 configurations; exact bond duality "
                        "does NOT break on the L=3 square bond torus"),
        },
        "pair_count_differences_e628_vs_fixed": {k: {"e628": a, "fixed": b}
                                                 for k, (a, b) in diff_pairs.items()},
    }
    dest = ROOT / "results" / "probe635-bond-rank-audit" / "latest.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out, indent=1))
    print(json.dumps({k: out[k] for k in
                      ("counterexample_mask63", "census_e628_verbatim",
                       "census_fixed", "pair_count_differences_e628_vs_fixed")},
                     indent=1)[:6000])


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Exact rank census for issue #622: site L=3,4 and square-bond L=3.

Imports the committed L=3,4 site enumerator from PR #606
(``scripts/homological_balance/exact_torus_enum.py``) as a library, verbatim,
re-exported here as ``_exact_torus_enum_imported``.  Nothing about its
algorithm or its published numbers is re-derived as a deliverable; what is new
here is only what #622 asks for:

* per-rank-pair decomposition ``(r_b, r_w)`` of ``M_L`` at the exact level
  (the #606 ledger aggregates over rank pairs; the Z question needs them
  separately);
* the *bond* analogue on the L=3 square torus: 2 L^2 = 18 primal bonds, the
  crossing-dual pairing convention of ``scripts/square_bond_kappa3.py``
  (``square_bond_pairs``), rank of the occupied-bond graph's H_1 image in
  H_1(T^2; Q) for both the occupied primal graph and the occupied dual graph
  (dual bonds expressed in primal coordinates via the crossing-dual
  identification), with ``r_b + r_w = 2`` checked configuration-wise
  (2^18 = 262144 configs, in budget per the probe);
* from each census: F_L = (1 + M_L)/2, its exact quantiles at a fixed level
  grid, and the anchored shape Z on that grid.

Conventions follow #606: black = 4-connected NN (occupied), white = 8-connected
NN+NNN (vacant), rank = rank of the winding image in H_1(T^2; Q), M_L(p) =
E[r_b] - 1 so that M_L(0) = -1 and M_L(1) = +1.  For bonds, black = occupied
primal bonds, white = vacant *dual* bonds transported to primal coordinates by
the crossing-dual identification (this is the geometric dual transport T of
``notes/square-bond-duality-tiny-torus.md``, not the bit complement — that note
proves the complement is a convention error).

Every L=3/L=4 site number here is a Fraction.  The bond L=3 census is likewise
exact (counts are integers); only the *quantile bisection* on the bond side
runs in floats to 1e-14, which the probe allows.

Rank of a bond subgraph on the torus: by the same winding argument as #606
(a 2-cell contributes only a contractible face boundary, winding 0), the
ambient H_1 rank of a bond subgraph G' is the rank over Q of the winding
vectors of its cycles.  I compute it as rank of the image of

    H_1(G'; Q) -> H_1(T^2; Q) ~= Q^2,

which for a spanning forest + non-tree edges is the rank of the integer
windings of the fundamental cycles.  Equivalently (and this is what the code
does, via the lift): build the graph on the lift Z^2, take a spanning forest
of the *lift* connectivity within one period window plus wrap edges, and
record lift displacement of each fundamental cycle.  The implementation below
uses the same fundamental-cycle/winding scheme as #606, adapted to edges.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "scripts"))

import _exact_torus_enum_imported as enum  # noqa: E402  (PR #606, verbatim)

from square_bond_kappa3 import square_bond_pairs  # noqa: E402

#: Anchors of the anchored shape Z, fixed before any of these numbers were
#: computed (probe text: "fixed anchors 0 < a < b < 1").  1/4, 3/4 is chosen
#: because it is symmetric under p -> 1-p; with a symmetric anchor pair the
#: odd part of Z is read directly off Z(u) + Z(1-u) - 1.
ANCHOR_A = Fraction(1, 4)
ANCHOR_B = Fraction(3, 4)

#: Level grid for Z, fixed before the runs (symmetric under u -> 1-u).
Z_LEVELS = (Fraction(1, 10), Fraction(1, 5), Fraction(3, 10), Fraction(2, 5),
            Fraction(1, 2), Fraction(3, 5), Fraction(7, 10), Fraction(4, 5),
            Fraction(9, 10))

TOL = Fraction(1, 10 ** 14)


# ------------------------------------------------------------- rank machinery

def bond_ambient_rank(edges: list[tuple[int, int, int, int]]) -> int:
    """Ambient H_1 rank of a bond subgraph on the LxL square torus.

    ``edges`` are (u, v, dx, dy) with (u,v) vertex ids in [0, L^2) and
    (dx, dy) the lift displacement of the directed edge v - u in Z^2 (in
    {-1, 0, 1}^2 for nearest-neighbour bonds, modulo the wrap convention).
    Rank of the winding image, same scheme as #606 but edges carry the lift
    step instead of vertices.
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
    from collections import defaultdict, deque
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


def site_rank_pair_components(L: int) -> dict[str, list[Fraction]]:
    """Coefficients of P11, P20, P02 (see module docstring) for site census.

    M_L = P20 - P02;  F_L = (1 + M_L)/2 = P20 + P11/2.
    """
    N = L * L
    P11 = [Fraction(0)] * (N + 1)
    P20 = [Fraction(0)] * (N + 1)
    P02 = [Fraction(0)] * (N + 1)
    for mask in range(1 << N):
        black = [(k // L, k % L) for k in range(N) if (mask >> k) & 1]
        rb = enum.ambient_rank(black, L, diagonal=False)
        rw = enum.ambient_rank([(k // L, k % L) for k in range(N)
                                if not ((mask >> k) & 1)], L, diagonal=True)
        k = len(black)
        if (rb, rw) == (1, 1):
            P11[k] += 1
        elif (rb, rw) == (2, 0):
            P20[k] += 1
        elif (rb, rw) == (0, 2):
            P02[k] += 1
        else:
            raise AssertionError((rb, rw))
    return {"P11": P11, "P20": P20, "P02": P02}


# ---------------------------------------------------------------- bond census

def bond_census(L: int) -> dict[str, object]:
    """Exact (r_b, r_w) decomposition on the LxL square bond torus.

    r_b: rank of the occupied *primal* bond graph.
    r_w: rank of the occupied *dual* bond graph, with every dual bond
         transported to primal coordinates by the crossing-dual pairing
         (geometric dual transport T; each dual edge equals a unique primal
         edge under the identification of ``square_bond_pairs``).
    Returns per-rank-pair count polynomials in p over |B| = number of
    occupied primal bonds, plus duality failures (must be 0).
    """
    pairs = square_bond_pairs(L)
    # Primal bond index and dual->primal map (crossing-dual identification).
    primal_index: dict[tuple[int, int], int] = {}
    dual_to_primal: dict[tuple[int, int], int] = {}
    for idx, bp in enumerate(pairs):
        (u, v, dx, dy) = bp.primal
        key = (min(u, v), max(u, v))
        assert key not in primal_index
        primal_index[key] = idx
        (du, dv) = bp.dual[0], bp.dual[1]
        dkey = (min(du, dv), max(du, dv))
        dual_to_primal[dkey] = -1  # placeholder, resolved after all primals registered
    assert len(primal_index) == 2 * L * L
    for dkey in list(dual_to_primal):
        assert dkey in primal_index, dkey
        dual_to_primal[dkey] = primal_index[dkey]

    # Lift displacement of each primal bond (as an undirected edge with a
    # canonical orientation: the one from square_bond_pairs, dx,dy in {-1,0,1}
    # with wrap already resolved by the vertex-id arithmetic of that function).
    lift_step: dict[int, tuple[int, int]] = {}
    for idx, bp in enumerate(pairs):
        (u, v, dx, dy) = bp.primal
        lift_step[idx] = (dx, dy)

    NB = 2 * L * L
    # Decompose by (r_b, r_w): count of configurations with |B| = k bonds.
    comp: dict[tuple[int, int], list[Fraction]] = {}
    dual_fail = 0
    pair_counts: dict[tuple[int, int], int] = {}
    total_masks = 1 << NB
    for mask in range(total_masks):
        occ = [(mask >> i) & 1 for i in range(NB)]
        # r_b: occupied primal graph.
        edges = [(pairs[i].primal[0], pairs[i].primal[1], *lift_step[i])
                 for i in range(NB) if occ[i]]
        rb = bond_ambient_rank(edges)
        # r_w: occupied dual graph, transported to primal indices.
        dedges = []
        for i in range(NB):
            (du, dv) = pairs[i].dual[0], pairs[i].dual[1]
            dkey = (min(du, dv), max(du, dv))
            pi = dual_to_primal[dkey]
            if not occ[pi]:
                dedges.append((pairs[pi].primal[0], pairs[pi].primal[1],
                               *lift_step[pi]))
        rw = bond_ambient_rank(dedges)
        if rb + rw != 2:
            dual_fail += 1
        pair_counts[(rb, rw)] = pair_counts.get((rb, rw), 0) + 1
        k = sum(occ)
        comp.setdefault((rb, rw), [Fraction(0)] * (NB + 1))[k] += 1

    return {
        "L": L,
        "bonds": NB,
        "configs": total_masks,
        "dual_fail": dual_fail,
        "pair_counts": pair_counts,
        "components": comp,
    }


# ------------------------------------------------------------ F, M, Q, Z math

def eval_poly_at(coeffs: list[Fraction], p: Fraction) -> Fraction:
    """sum_k C_k p^k (1-p)^{N-k}, exact."""
    N = len(coeffs) - 1
    total = Fraction(0)
    for k, c in enumerate(coeffs):
        total += c * p ** k * (1 - p) ** (N - k)
    return total


def M_and_F_from_components(components: dict[str, list[Fraction]]
                            ) -> tuple[list[Fraction], list[Fraction]]:
    """M = P20 - P02, F = P20 + P11/2 (both degree-N coefficient lists)."""
    P11, P20, P02 = components["P11"], components["P20"], components["P02"]
    N = len(P11) - 1
    M = [P20[k] - P02[k] for k in range(N + 1)]
    F = [P20[k] + P11[k] / 2 for k in range(N + 1)]
    return M, F


def eval_F(F: list[Fraction], p: Fraction) -> Fraction:
    return eval_poly_at(F, p)


def eval_F_float(F: list[Fraction], p: float) -> float:
    """F(p) in floats, log-space anchored, for the bisection quantiles."""
    N = len(F) - 1
    if p <= 0.0:
        return 0.0
    if p >= 1.0:
        return 1.0
    lp, lq = math.log(p), math.log1p(-p)
    lf = math.lgamma(N + 1)
    mode = min(N, max(0, int((N + 1) * p)))
    peak = lf - math.lgamma(mode + 1) - math.lgamma(N - mode + 1) \
        + mode * lp + (N - mode) * lq
    total = 0.0
    for k in range(N + 1):
        logw = lf - math.lgamma(k + 1) - math.lgamma(N - k + 1) \
            + k * lp + (N - k) * lq
        if peak - logw > 700.0:
            continue
        total += float(F[k]) * math.exp(logw - peak)
    return total


def quantile_exact(F: list[Fraction], u: Fraction,
                   lo: Fraction = Fraction(0), hi: Fraction = Fraction(1)
                   ) -> Fraction:
    """Exact bisection to 1e-14 on the strictly increasing polynomial F."""
    assert eval_F(F, lo) <= u <= eval_F(F, hi)
    while hi - lo > TOL:
        mid = (lo + hi) / 2
        if eval_F(F, mid) < u:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def z_from_quantiles(qa: Fraction, qb: Fraction, qs: list[Fraction]
                     ) -> list[Fraction]:
    return [(q - qa) / (qb - qa) for q in qs]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check-606", action="store_true",
                    help="assert margins against the published #606 ledger")
    ap.add_argument("--bond", action="store_true", help="also run bond L=3")
    args = ap.parse_args()

    out: dict[str, object] = {
        "schema": "matching-one.probe-invariant-shape.census.v1",
        "anchors": [str(ANCHOR_A), str(ANCHOR_B)],
        "z_levels": [str(u) for u in Z_LEVELS],
        "site": {},
        "bond": {},
        "checks": {},
    }

    published = {"3": {"0,2": 259, "1,1": 162, "2,0": 91},
                 "4": {"0,2": 36559, "1,1": 19932, "2,0": 9045}}

    comp_names = {(0, 2): "P02", (1, 1): "P11", (2, 0): "P20"}
    for L in (3, 4):
        comps = site_rank_pair_components(L)
        M, F = M_and_F_from_components(comps)
        N = L * L
        # Sanity against #606 cited facts (cited, not re-proven: M(0)=-1,
        # M(1)=+1, and the exact M(1/2) values from the ledger).
        m0 = eval_poly_at(M, Fraction(0))
        m1 = eval_poly_at(M, Fraction(1))
        m_half = eval_poly_at(M, Fraction(1, 2))
        assert m0 == -1 and m1 == 1
        if args.check_606:
            counts = {(rb, rw): sum(comps[name][k] for k in range(N + 1))
                      for (rb, rw), name in comp_names.items()}
            for key, val in published[str(L)].items():
                rb, rw = map(int, key.split(","))
                assert counts[(rb, rw)] == val, (L, key, counts)
        qa = quantile_exact(F, ANCHOR_A)
        qb = quantile_exact(F, ANCHOR_B)
        qs = [quantile_exact(F, u) for u in Z_LEVELS]
        zs = z_from_quantiles(qa, qb, qs)
        half = [float(z) for z in zs]
        odd_sum = [float(zs[i] + zs[-1 - i]) - 1.0 for i in range(4)] \
            + [0.0] + [float(zs[i] + zs[-1 - i]) - 1.0 for i in range(5, 9)]
        out["site"][str(L)] = {  # type: ignore[index]
            "N": N,
            "M_half": str(m_half),
            "M_half_float": float(m_half),
            "p_L_H": float(quantile_exact(
                [Fraction(0)] + [(M[k] + M[0]) * 0 for k in range(N)],  # dummy
                Fraction(0))) if False else None,
            "Q_quarter": str(qa),
            "Q_threequarters": str(qb),
            "Q_quarter_float": float(qa),
            "Q_threequarters_float": float(qb),
            "Z_levels_float": half,
            "Z_odd_part_sum_float": odd_sum,
            "rank_pair_counts": {f"{rb},{rw}": int(sum(comps[name][k]
                                  for k in range(N + 1)))
                                 for (rb, rw), name in comp_names.items()},
        }
        # p_L^H is #606's cited number; do not recompute it here.
        out["site"][str(L)]["p_L_H_cited"] = (  # type: ignore[index]
            0.5865114551126757 if L == 3 else 0.5906721123310283)

    if args.bond:
        cens = bond_census(3)
        comp = cens["components"]
        NB = cens["bonds"]
        P11 = comp.get((1, 1), [Fraction(0)] * (NB + 1))
        P20 = comp.get((2, 0), [Fraction(0)] * (NB + 1))
        P02 = comp.get((0, 2), [Fraction(0)] * (NB + 1))
        M = [P20[k] - P02[k] for k in range(NB + 1)]
        F = [P20[k] + P11[k] / 2 for k in range(NB + 1)]
        # float quantile bisection to 1e-14 (allowed for bond L=3)
        def qfloat(u: float) -> float:
            lo, hi = 0.0, 1.0
            for _ in range(200):
                mid = 0.5 * (lo + hi)
                if eval_F_float(F, mid) < u:
                    lo = mid
                else:
                    hi = mid
            return 0.5 * (lo + hi)
        qa, qb = qfloat(0.25), qfloat(0.75)
        zs = [(qfloat(float(u)) - qa) / (qb - qa) for u in
              (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)]
        m_half = eval_poly_at(M, Fraction(1, 2))
        out["bond"]["3"] = {  # type: ignore[index]
            "bonds": NB, "configs": cens["configs"],
            "dual_fail": cens["dual_fail"],
            "rank_pair_counts": {f"{rb},{rw}": int(n) for (rb, rw), n
                                 in cens["pair_counts"].items()},
            "M_half": str(m_half), "M_half_float": float(m_half),
            "Q_quarter": qa, "Q_threequarters": qb,
            "Z_levels_float": zs,
            "self_dual_p_half": 0.5,
        }

    dest = ROOT / "results" / "probe-invariant-shape" / "census-exact.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1)[:4000])


if __name__ == "__main__":
    main()

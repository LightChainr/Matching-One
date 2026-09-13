#!/usr/bin/env python3
"""Exact bounded search for a summary-compression witness.

Frozen before search (do not adapt):
  Summary = (S(z) coefficients, n, H2, b2, terminal-local r-neighborhood).
  Successor-hazard moments are NOT in the summary.
  Experiments = ordinary continuation plus shared prefix p<=2 with two
  independent continuations of length c<=2; observation is L-R disconnection.
  Delayed-fork (p=1,c=1) is computed but is not used as the first reported
  distinguisher when a deeper frozen experiment also splits.

Enumerated class G (explicit, bounded, planar-by-construction except the
n<=5 exhaustive plane filter):
  G0 n=0 unit edge used only as a composer;
  exhaustive connected two-terminal simple graphs, n<=5, no L-R edge,
    L-R path exists, G+{LR} planar;
  two-terminal series-parallel graphs n<=12;
  Wheatstone seed and SP compositions with it, n<=12;
  multi-path (theta) graphs n<=12;
  2xk and 3xk ladders/grids n<=12;
  path-hidden copies of n<=5 cores (left/right corridors of length 2 or 3).

Arithmetic is exact (fractions.Fraction / integer counts). No Monte Carlo.
"""

from __future__ import annotations

import json
import sys
from collections import defaultdict, deque
from fractions import Fraction
from itertools import combinations, permutations
from math import comb
from pathlib import Path
from typing import Iterable, Optional

SCHEMA = "matching-one/bounded-summary-search/v1"
_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parents[1] if _HERE.name == "summary_search" else _HERE
OUT_JSON = _ROOT / "artifacts" / "bounded_summary_search.json"
OUT_MD = _ROOT / "artifacts" / "bounded_summary_search.md"


class Net:
    __slots__ = ("n", "adj_L", "adj_R", "adj", "lr_edge", "name")

    def __init__(
        self,
        n: int,
        adj_L: int,
        adj_R: int,
        adj: tuple[int, ...],
        lr_edge: bool = False,
        name: str = "",
    ):
        self.n = n
        self.adj_L = adj_L
        self.adj_R = adj_R
        self.adj = adj
        self.lr_edge = lr_edge
        self.name = name

    def copy_named(self, name: str) -> "Net":
        return Net(self.n, self.adj_L, self.adj_R, self.adj, self.lr_edge, name)


def _validate(net: Net) -> None:
    n = net.n
    full = (1 << n) - 1
    assert 0 <= net.adj_L <= full and 0 <= net.adj_R <= full
    assert len(net.adj) == n
    for i, a in enumerate(net.adj):
        assert 0 <= a <= full
        assert ((a >> i) & 1) == 0
        for j in range(n):
            if (a >> j) & 1:
                assert (net.adj[j] >> i) & 1


def edge_count(net: Net) -> int:
    internal = sum(a.bit_count() for a in net.adj) // 2
    return internal + net.adj_L.bit_count() + net.adj_R.bit_count() + int(net.lr_edge)


def connected_carrier(net: Net) -> bool:
    """L, R and every switchable vertex lie in one connected component of the carrier."""
    n = net.n
    if n == 0:
        return net.lr_edge
    start = net.adj_L
    if not start and not net.lr_edge:
        return False
    q = start
    seen = start
    while q:
        v = (q & -q).bit_length() - 1
        q ^= 1 << v
        extra = net.adj[v] & ~seen
        seen |= extra
        q |= extra
    reached_R = bool(seen & net.adj_R) or net.lr_edge
    if not reached_R:
        return False
    return seen == (1 << n) - 1


def lr_path_exists(net: Net) -> bool:
    if net.lr_edge:
        return True
    start = net.adj_L
    seen = start
    q = start
    while q:
        v = (q & -q).bit_length() - 1
        q ^= 1 << v
        if (net.adj_R >> v) & 1:
            return True
        extra = net.adj[v] & ~seen
        seen |= extra
        q |= extra
    return False


def _full_adj_matrix(net: Net, plus_lr: bool) -> list[int]:
    """Vertices: 0..n-1 switchable, n = L, n+1 = R. Return neighbour bitmasks."""
    n = net.n
    L, R = n, n + 1
    vtot = n + 2
    adj = [0] * vtot
    for i, a in enumerate(net.adj):
        adj[i] |= a
    for i in range(n):
        if (net.adj_L >> i) & 1:
            adj[i] |= 1 << L
            adj[L] |= 1 << i
        if (net.adj_R >> i) & 1:
            adj[i] |= 1 << R
            adj[R] |= 1 << i
    if plus_lr or net.lr_edge:
        adj[L] |= 1 << R
        adj[R] |= 1 << L
    return adj


def _has_k5_subgraph(adj: list[int]) -> bool:
    v = len(adj)
    for comb5 in combinations(range(v), 5):
        ok = True
        for i, a in enumerate(comb5):
            for b in comb5[i + 1 :]:
                if ((adj[a] >> b) & 1) == 0:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            return True
    return False


def _has_k33_subgraph(adj: list[int]) -> bool:
    v = len(adj)
    if v < 6:
        return False
    verts = range(v)
    for A in combinations(verts, 3):
        rest = [x for x in verts if x not in A]
        for B in combinations(rest, 3):
            ok = True
            for a in A:
                for b in B:
                    if ((adj[a] >> b) & 1) == 0:
                        ok = False
                        break
                if not ok:
                    break
            if ok:
                return True
    return False


def _contract_edge(adj: list[int], u: int, w: int) -> list[int]:
    """Contract w into u, delete w, relabel >w down by 1."""
    v = len(adj)
    merged = (adj[u] | adj[w]) & ~((1 << u) | (1 << w))
    new = []
    for i in range(v):
        if i == w:
            continue
        row = adj[i]
        if i == u:
            row = merged
        else:
            if (row >> w) & 1:
                row |= 1 << u
            row &= ~(1 << w)
        low = row & ((1 << w) - 1)
        high = row >> (w + 1)
        row = low | (high << w)
        new.append(row)
    return new


def _nonplanar_forbidden(adj: list[int]) -> bool:
    v = len(adj)
    e = sum(a.bit_count() for a in adj) // 2
    if v <= 4:
        return False
    if v >= 3 and e > 3 * v - 6:
        return True
    if v == 5:
        return e >= 10
    if _has_k5_subgraph(adj) or _has_k33_subgraph(adj):
        return True
    if v <= 6:
        return False
    if v > 8:
        return False
    for u in range(v):
        nb = adj[u]
        while nb:
            w = (nb & -nb).bit_length() - 1
            nb ^= 1 << w
            if w < u:
                continue
            H = _contract_edge(adj, u, w)
            if _nonplanar_forbidden(H):
                return True
    return False


def is_plane_two_terminal(net: Net) -> bool:
    """True iff G ∪ {L-R} is planar (L,R can sit on a common face)."""
    adj = _full_adj_matrix(net, plus_lr=True)
    return not _nonplanar_forbidden(adj)


def _encoding(n: int, adj_L: int, adj_R: int, adj: tuple[int, ...]) -> tuple:
    upper = []
    for i in range(n):
        for j in range(i + 1, n):
            upper.append((adj[i] >> j) & 1)
    return (n, adj_L, adj_R, tuple(upper))


def _apply_perm(net: Net, perm: tuple[int, ...]) -> tuple[int, int, tuple[int, ...]]:
    """perm[old] = new index."""
    n = net.n
    adj_L = 0
    adj_R = 0
    for old in range(n):
        new = perm[old]
        if (net.adj_L >> old) & 1:
            adj_L |= 1 << new
        if (net.adj_R >> old) & 1:
            adj_R |= 1 << new
    adj = [0] * n
    for old in range(n):
        a = net.adj[old]
        no = perm[old]
        b = 0
        m = a
        while m:
            j = (m & -m).bit_length() - 1
            m ^= 1 << j
            b |= 1 << perm[j]
        adj[no] = b
    return adj_L, adj_R, tuple(adj)


def _wl_colors(net: Net) -> tuple[int, ...]:
    n = net.n
    color = []
    for i in range(n):
        tL = (net.adj_L >> i) & 1
        tR = (net.adj_R >> i) & 1
        deg = net.adj[i].bit_count() + tL + tR
        color.append((tL, tR, deg))
    for _ in range(max(n, 1)):
        nxt = []
        for i in range(n):
            nb = []
            a = net.adj[i]
            m = a
            while m:
                j = (m & -m).bit_length() - 1
                m ^= 1 << j
                nb.append(color[j])
            nb.sort()
            nxt.append((color[i], tuple(nb)))
        mapping = {}
        for k in sorted(set(nxt)):
            mapping[k] = len(mapping)
        new_color = tuple(mapping[k] for k in nxt)
        if isinstance(color[0], int) and new_color == tuple(color):
            return new_color
        color = list(new_color)
        if len(set(color)) == n:
            return tuple(color)
    return tuple(color)


def canonical_key(net: Net) -> tuple:
    n = net.n
    if n == 0:
        return (0, int(net.lr_edge))
    colors = _wl_colors(net)
    cells: dict[int, list[int]] = defaultdict(list)
    for i, c in enumerate(colors):
        cells[c].append(i)
    cell_list = [cells[c] for c in sorted(cells)]
    best: Optional[tuple] = None
    sizes = [len(c) for c in cell_list]
    if max(sizes, default=1) <= 8 and _prod(sizes) <= 40320:
        for assignment in _product_perms(cell_list):
            perm = tuple(assignment)
            adj_L, adj_R, adj = _apply_perm(net, perm)
            enc = _encoding(n, adj_L, adj_R, adj)
            if best is None or enc < best:
                best = enc
        assert best is not None
        return best
    return _canon_backtrack(net, colors)


def _prod(xs: list[int]) -> int:
    p = 1
    for x in xs:
        p *= x
        if p > 10**9:
            return p
    return p


def _product_perms(cell_list: list[list[int]]):
    n = sum(len(c) for c in cell_list)
    blocks = []
    start = 0
    for cell in cell_list:
        blocks.append((cell, list(range(start, start + len(cell)))))
        start += len(cell)

    def rec(i, perm):
        if i == len(blocks):
            yield tuple(perm)
            return
        olds, news = blocks[i]
        for seq in permutations(news):
            for o, nv in zip(olds, seq):
                perm[o] = nv
            yield from rec(i + 1, perm)

    yield from rec(0, [0] * n)


def _canon_backtrack(net: Net, colors: tuple[int, ...]) -> tuple:
    n = net.n
    best: Optional[tuple] = None

    def refine(col: list[int]) -> list[int]:
        nxt = []
        for i in range(n):
            nb = []
            a = net.adj[i]
            m = a
            while m:
                j = (m & -m).bit_length() - 1
                m ^= 1 << j
                nb.append(col[j])
            nb.sort()
            nxt.append((col[i], tuple(nb), (net.adj_L >> i) & 1, (net.adj_R >> i) & 1))
        mapping = {k: i for i, k in enumerate(sorted(set(nxt)))}
        return [mapping[k] for k in nxt]

    def search(col: list[int]) -> None:
        nonlocal best
        cells: dict[int, list[int]] = defaultdict(list)
        for i, c in enumerate(col):
            cells[c].append(i)
        clist = [cells[c] for c in sorted(cells)]
        if all(len(c) == 1 for c in clist):
            perm = [0] * n
            order = sorted(range(n), key=lambda i: col[i])
            for new, old in enumerate(order):
                perm[old] = new
            adj_L, adj_R, adj = _apply_perm(net, tuple(perm))
            enc = _encoding(n, adj_L, adj_R, adj)
            if best is None or enc < best:
                best = enc
            return
        cell = max(clist, key=len)
        if len(cell) > 9:
            perm = [0] * n
            order = sorted(range(n), key=lambda i: (col[i], i))
            for new, old in enumerate(order):
                perm[old] = new
            adj_L, adj_R, adj = _apply_perm(net, tuple(perm))
            enc = _encoding(n, adj_L, adj_R, adj)
            if best is None or enc < best:
                best = enc
            return
        for v in cell:
            ncol = col[:]
            ncol[v] = -1
            ncol = refine(ncol)
            search(ncol)

    search(list(colors))
    assert best is not None
    return best


def apply_canon(net: Net) -> Net:
    n = net.n
    if n == 0:
        return net
    target = canonical_key(net)
    colors = _wl_colors(net)
    cells: dict[int, list[int]] = defaultdict(list)
    for i, c in enumerate(colors):
        cells[c].append(i)
    cell_list = [cells[c] for c in sorted(cells)]
    sizes = [len(c) for c in cell_list]
    if max(sizes, default=1) <= 8 and _prod(sizes) <= 40320:
        for assignment in _product_perms(cell_list):
            adj_L, adj_R, adj = _apply_perm(net, tuple(assignment))
            if _encoding(n, adj_L, adj_R, adj) == target:
                return Net(n, adj_L, adj_R, adj, net.lr_edge, net.name)
    if canonical_key(net) == _encoding(n, net.adj_L, net.adj_R, net.adj):
        return net
    if n <= 8:
        for perm in permutations(range(n)):
            adj_L, adj_R, adj = _apply_perm(net, perm)
            if _encoding(n, adj_L, adj_R, adj) == target:
                return Net(n, adj_L, adj_R, adj, net.lr_edge, net.name)
    return net


def safe_table(net: Net) -> list[bool]:
    """safe[mask] True iff occupied mask does NOT connect L-R."""
    n = net.n
    N = 1 << n
    table = [False] * N
    if net.lr_edge:
        return table
    adj = net.adj
    adj_L = net.adj_L
    adj_R = net.adj_R
    for mask in range(N):
        start = mask & adj_L
        if not start:
            table[mask] = True
            continue
        seen = start
        q = start
        absorbed = False
        while q:
            v = (q & -q).bit_length() - 1
            q ^= 1 << v
            if (adj_R >> v) & 1:
                absorbed = True
                break
            extra = adj[v] & mask & ~seen
            seen |= extra
            q |= extra
        table[mask] = not absorbed
    return table


def S_coeffs(table: list[bool], n: int) -> tuple[int, ...]:
    s = [0] * (n + 1)
    for mask, safe in enumerate(table):
        if safe:
            s[mask.bit_count()] += 1
    return tuple(s)


def frac(num: int, den: int) -> Fraction:
    if den == 0:
        return Fraction(0, 1)
    return Fraction(num, den)


def experiments(table: list[bool], n: int) -> dict[str, Fraction]:
    """Exact probabilities in the frozen depth-2 language."""
    full = (1 << n) - 1
    s = S_coeffs(table, n)
    out: dict[str, Fraction] = {}
    for k in range(n + 1):
        out[f"ord_{k}"] = frac(s[k], comb(n, k) if k <= n else 1)

    def p_surv(occupied: int, c: int) -> Fraction:
        alive_bits = []
        m = (~occupied) & full
        tmp = m
        while tmp:
            w = (tmp & -tmp).bit_length() - 1
            tmp ^= 1 << w
            alive_bits.append(w)
        mcount = len(alive_bits)
        if c == 0:
            return Fraction(1, 1) if table[occupied] else Fraction(0, 1)
        if c > mcount:
            return Fraction(1, 1) if table[occupied | m] else Fraction(0, 1)
        if c == 1:
            good = 0
            for w in alive_bits:
                if table[occupied | (1 << w)]:
                    good += 1
            return frac(good, mcount)
        if c == 2:
            tot = comb(mcount, 2)
            good = 0
            for i in range(mcount):
                a = 1 << alive_bits[i]
                for j in range(i + 1, mcount):
                    if table[occupied | a | (1 << alive_bits[j])]:
                        good += 1
            return frac(good, tot)
        raise ValueError("c<=2")

    def prefix_both(p: int, c1: int, c2: int) -> Fraction:
        if p == 0:
            return p_surv(0, c1) * p_surv(0, c2)
        verts = list(range(n))
        num = Fraction(0, 1)
        den = 0
        if p == 1:
            for v in verts:
                den += 1
                occ = 1 << v
                if not table[occ]:
                    continue
                num += p_surv(occ, c1) * p_surv(occ, c2)
            return num / den if den else Fraction(0, 1)
        if p == 2:
            for i in verts:
                for j in verts:
                    if i == j:
                        continue
                    den += 1
                    occ = (1 << i) | (1 << j)
                    if not table[occ]:
                        continue
                    num += p_surv(occ, c1) * p_surv(occ, c2)
            return num / den if den else Fraction(0, 1)
        raise ValueError("p<=2")

    out["E0_c1"] = prefix_both(0, 1, 1)
    out["E0_c2"] = prefix_both(0, 2, 2)
    out["E0_mix"] = prefix_both(0, 1, 2)
    out["E1_c1"] = prefix_both(1, 1, 1)
    out["E1_c2"] = prefix_both(1, 2, 2)
    out["E1_mix"] = prefix_both(1, 1, 2)
    out["E2_c1"] = prefix_both(2, 1, 1)
    out["E2_c2"] = prefix_both(2, 2, 2)
    out["E2_mix"] = prefix_both(2, 1, 2)
    return out


def experiments_fast(table: list[bool], n: int) -> dict[str, Fraction]:
    """Integer-arithmetic equivalent of experiments(); same values, less Fraction churn."""
    full = (1 << n) - 1
    s = S_coeffs(table, n)
    out: dict[str, Fraction] = {}
    for k in range(n + 1):
        out[f"ord_{k}"] = frac(s[k], comb(n, k) if k <= n else 1)

    def good_c(occupied: int, c: int) -> tuple[int, int]:
        alive = []
        m = (~occupied) & full
        tmp = m
        while tmp:
            w = (tmp & -tmp).bit_length() - 1
            tmp ^= 1 << w
            alive.append(w)
        mc = len(alive)
        if c == 1:
            good = 0
            for w in alive:
                if table[occupied | (1 << w)]:
                    good += 1
            return good, (mc if mc else 1)
        if c == 2:
            if mc < 2:
                return (1 if table[occupied | m] else 0), 1
            good = 0
            for i in range(mc):
                a = 1 << alive[i]
                for j in range(i + 1, mc):
                    if table[occupied | a | (1 << alive[j])]:
                        good += 1
            return good, comb(mc, 2)
        raise ValueError("c<=2")

    def prefix_both_int(p: int, c1: int, c2: int) -> Fraction:
        if p == 0:
            g1, d1 = good_c(0, c1)
            g2, d2 = good_c(0, c2)
            return Fraction(g1 * g2, d1 * d2)
        if p == 1:
            if n - 1 >= max(c1, c2) and n - 1 >= 2:
                d1 = (n - 1) if c1 == 1 else comb(n - 1, 2)
                d2 = (n - 1) if c2 == 1 else comb(n - 1, 2)
                acc = 0
                for v in range(n):
                    occ = 1 << v
                    if not table[occ]:
                        continue
                    g1, _ = good_c(occ, c1)
                    g2, _ = good_c(occ, c2)
                    acc += g1 * g2
                den = n * d1 * d2
                return Fraction(acc, den) if den else Fraction(0, 1)
            num = Fraction(0, 1)
            for v in range(n):
                occ = 1 << v
                if not table[occ]:
                    continue
                g1, d1 = good_c(occ, c1)
                g2, d2 = good_c(occ, c2)
                num += Fraction(g1 * g2, d1 * d2)
            return num / n if n else Fraction(0, 1)
        if p == 2:
            nden = n * (n - 1)
            if n - 2 >= max(c1, c2) and n - 2 >= 2:
                d1 = (n - 2) if c1 == 1 else comb(n - 2, 2)
                d2 = (n - 2) if c2 == 1 else comb(n - 2, 2)
                acc = 0
                for i in range(n):
                    for j in range(i + 1, n):
                        occ = (1 << i) | (1 << j)
                        if not table[occ]:
                            continue
                        g1, _ = good_c(occ, c1)
                        g2, _ = good_c(occ, c2)
                        acc += 2 * g1 * g2
                den = nden * d1 * d2
                return Fraction(acc, den) if den else Fraction(0, 1)
            num = Fraction(0, 1)
            den = 0
            for i in range(n):
                for j in range(n):
                    if i == j:
                        continue
                    den += 1
                    occ = (1 << i) | (1 << j)
                    if not table[occ]:
                        continue
                    g1, d1 = good_c(occ, c1)
                    g2, d2 = good_c(occ, c2)
                    num += Fraction(g1 * g2, d1 * d2)
            return num / den if den else Fraction(0, 1)
        raise ValueError("p<=2")

    out["E0_c1"] = prefix_both_int(0, 1, 1)
    out["E0_c2"] = prefix_both_int(0, 2, 2)
    out["E0_mix"] = prefix_both_int(0, 1, 2)
    out["E1_c1"] = prefix_both_int(1, 1, 1)
    out["E1_c2"] = prefix_both_int(1, 2, 2)
    out["E1_mix"] = prefix_both_int(1, 1, 2)
    out["E2_c1"] = prefix_both_int(2, 1, 1)
    out["E2_c2"] = prefix_both_int(2, 2, 2)
    out["E2_mix"] = prefix_both_int(2, 1, 2)
    return out


FROZEN_EXPERIMENT_ORDER = [
    "E1_c2",
    "E1_mix",
    "E2_c1",
    "E2_c2",
    "E2_mix",
    "E0_c2",
    "E0_mix",
    "E0_c1",
    "E1_c1",
]


def successor_h2_moments(table: list[bool], n: int) -> tuple[Fraction, Fraction]:
    """Diagnostic only; NOT in the frozen summary."""
    if n == 0:
        return Fraction(0, 1), Fraction(0, 1)
    xs = []
    for v in range(n):
        occ = 1 << v
        if not table[occ]:
            continue
        h2 = 0
        for w in range(n):
            if w == v:
                continue
            if not table[occ | (1 << w)]:
                h2 += 1
        xs.append(h2)
    if not xs:
        return Fraction(0, 1), Fraction(0, 1)
    m1 = Fraction(sum(xs), len(xs))
    m2 = Fraction(sum(x * x for x in xs), len(xs))
    return m1, m2


def distances(net: Net) -> tuple[list[int], list[int]]:
    n = net.n
    INF = 99

    def bfs(start_mask: int) -> list[int]:
        dist = [INF] * n
        q = deque()
        m = start_mask
        while m:
            v = (m & -m).bit_length() - 1
            m ^= 1 << v
            dist[v] = 1
            q.append(v)
        while q:
            v = q.popleft()
            a = net.adj[v]
            while a:
                w = (a & -a).bit_length() - 1
                a ^= 1 << w
                if dist[w] > dist[v] + 1:
                    dist[w] = dist[v] + 1
                    q.append(w)
        return dist

    return bfs(net.adj_L), bfs(net.adj_R)


def neighborhood_key(net: Net, radius: int) -> tuple:
    """Rooted typed induced ball around {L,R} of radius r, L and R labeled."""
    n = net.n
    dL, dR = distances(net)
    ball = [i for i in range(n) if min(dL[i], dR[i]) <= radius]
    k = len(ball)
    if k == 0:
        return (0, radius, int(net.lr_edge))
    idx = {old: new for new, old in enumerate(ball)}
    types = []
    adj_local = [0] * k
    tL = 0
    tR = 0
    for new, old in enumerate(ball):
        types.append((dL[old], dR[old], (net.adj_L >> old) & 1, (net.adj_R >> old) & 1))
        if (net.adj_L >> old) & 1:
            tL |= 1 << new
        if (net.adj_R >> old) & 1:
            tR |= 1 << new
        a = net.adj[old]
        loc = 0
        mm = a
        while mm:
            w = (mm & -mm).bit_length() - 1
            mm ^= 1 << w
            if w in idx:
                loc |= 1 << idx[w]
        adj_local[new] = loc
    cells: dict[tuple, list[int]] = defaultdict(list)
    for i, t in enumerate(types):
        cells[t].append(i)
    cell_list = [cells[t] for t in sorted(cells)]
    sizes = [len(c) for c in cell_list]
    best = None
    if k <= 9 and _prod(sizes) <= 40320:

        def rec(bi, perm):
            nonlocal best
            if bi == len(cell_list):
                adj_L_n = 0
                adj_R_n = 0
                adj_n = [0] * k
                for old in range(k):
                    nv = perm[old]
                    if (tL >> old) & 1:
                        adj_L_n |= 1 << nv
                    if (tR >> old) & 1:
                        adj_R_n |= 1 << nv
                    b = 0
                    a = adj_local[old]
                    while a:
                        j = (a & -a).bit_length() - 1
                        a ^= 1 << j
                        b |= 1 << perm[j]
                    adj_n[nv] = b
                types_n = [None] * k
                for old in range(k):
                    types_n[perm[old]] = types[old]
                upper = []
                for i in range(k):
                    for j in range(i + 1, k):
                        upper.append((adj_n[i] >> j) & 1)
                enc = (k, tuple(types_n), adj_L_n, adj_R_n, tuple(upper))
                if best is None or enc < best:
                    best = enc
                return
            olds = cell_list[bi]
            for seq in permutations(olds):
                for a, b in zip(olds, seq):
                    perm[a] = b
                rec(bi + 1, perm)

        rec(0, list(range(k)))
        assert best is not None
        return (radius, best)
    order = sorted(range(k), key=lambda i: (types[i], i))
    perm = [0] * k
    for new, old in enumerate(order):
        perm[old] = new
    adj_L_n = 0
    adj_R_n = 0
    adj_n = [0] * k
    types_n = [None] * k
    for old in range(k):
        nv = perm[old]
        types_n[nv] = types[old]
        if (tL >> old) & 1:
            adj_L_n |= 1 << nv
        if (tR >> old) & 1:
            adj_R_n |= 1 << nv
        b = 0
        a = adj_local[old]
        while a:
            j = (a & -a).bit_length() - 1
            a ^= 1 << j
            b |= 1 << perm[j]
        adj_n[nv] = b
    upper = tuple((adj_n[i] >> j) & 1 for i in range(k) for j in range(i + 1, k))
    return (radius, (k, tuple(types_n), adj_L_n, adj_R_n, upper))


def path_net(k: int) -> Net:
    """L-0-1-...-(k-1)-R, k>=1. k==0 is the unit edge."""
    if k == 0:
        return Net(0, 0, 0, (), True, "edge")
    adj = [0] * k
    for i in range(k - 1):
        adj[i] |= 1 << (i + 1)
        adj[i + 1] |= 1 << i
    adj_L = 1 << 0
    adj_R = 1 << (k - 1)
    return Net(k, adj_L, adj_R, tuple(adj), False, f"path{k}")


def wheatstone() -> Net:
    """n=2: L-a, L-b, a-R, b-R, a-b. Planar, not series-parallel."""
    adj_L = 0b11
    adj_R = 0b11
    adj = (0b10, 0b01)
    return Net(2, adj_L, adj_R, adj, False, "wheatstone")


def series(A: Net, B: Net, name: str = "") -> Net:
    """Identify A.R with B.L as a new switchable vertex z."""
    n1, n2 = A.n, B.n
    n = n1 + n2 + 1
    z = n1 + n2
    adj = [0] * n
    for i in range(n1):
        adj[i] = A.adj[i]
    for i in range(n2):
        a = B.adj[i]
        shifted = 0
        tmp = a
        while tmp:
            j = (tmp & -tmp).bit_length() - 1
            tmp ^= 1 << j
            shifted |= 1 << (n1 + j)
        adj[n1 + i] = shifted
    zmask = 0
    tmp = A.adj_R
    while tmp:
        j = (tmp & -tmp).bit_length() - 1
        tmp ^= 1 << j
        zmask |= 1 << j
        adj[j] |= 1 << z
    tmp = B.adj_L
    while tmp:
        j = (tmp & -tmp).bit_length() - 1
        tmp ^= 1 << j
        idx = n1 + j
        zmask |= 1 << idx
        adj[idx] |= 1 << z
    adj[z] = zmask
    adj_L = A.adj_L
    adj_R = 0
    tmp = B.adj_R
    while tmp:
        j = (tmp & -tmp).bit_length() - 1
        tmp ^= 1 << j
        adj_R |= 1 << (n1 + j)
    if A.lr_edge:
        adj_L |= 1 << z
    if B.lr_edge:
        adj_R |= 1 << z
    return Net(n, adj_L, adj_R, tuple(adj), False, name or f"series({A.name}|{B.name})")


def parallel(A: Net, B: Net, name: str = "") -> Net:
    n1, n2 = A.n, B.n
    n = n1 + n2
    adj = [0] * n
    for i in range(n1):
        adj[i] = A.adj[i]
    for i in range(n2):
        a = B.adj[i]
        shifted = 0
        tmp = a
        while tmp:
            j = (tmp & -tmp).bit_length() - 1
            tmp ^= 1 << j
            shifted |= 1 << (n1 + j)
        adj[n1 + i] = shifted
    adj_L = A.adj_L
    adj_R = A.adj_R
    tmp = B.adj_L
    while tmp:
        j = (tmp & -tmp).bit_length() - 1
        tmp ^= 1 << j
        adj_L |= 1 << (n1 + j)
    tmp = B.adj_R
    while tmp:
        j = (tmp & -tmp).bit_length() - 1
        tmp ^= 1 << j
        adj_R |= 1 << (n1 + j)
    lr = A.lr_edge or B.lr_edge
    return Net(n, adj_L, adj_R, tuple(adj), lr, name or f"par({A.name}|{B.name})")


def hide_corridors(core: Net, left: int, right: int) -> Net:
    g = core
    if left:
        g = series(path_net(left), g, "")
    if right:
        g = series(g, path_net(right), "")
    g.name = f"hide(L{left},{core.name},R{right})"
    return g


def ladder(rows: int, cols: int) -> Net:
    n = rows * cols
    adj = [0] * n

    def vid(r, c):
        return r * cols + c

    for r in range(rows):
        for c in range(cols):
            u = vid(r, c)
            if c + 1 < cols:
                v = vid(r, c + 1)
                adj[u] |= 1 << v
                adj[v] |= 1 << u
            if r + 1 < rows:
                v = vid(r + 1, c)
                adj[u] |= 1 << v
                adj[v] |= 1 << u
    adj_L = 0
    adj_R = 0
    for r in range(rows):
        adj_L |= 1 << vid(r, 0)
        adj_R |= 1 << vid(r, cols - 1)
    return Net(n, adj_L, adj_R, tuple(adj), False, f"ladder{rows}x{cols}")


def multipath(lengths: tuple[int, ...]) -> Net:
    g = path_net(lengths[0])
    g.name = f"mp{lengths}"
    for L in lengths[1:]:
        g = parallel(g, path_net(L), f"mp{lengths}")
    return g


def _edge_slots(n: int) -> list[tuple[str, int, int]]:
    slots = []
    for i in range(n):
        slots.append(("L", i, -1))
    for i in range(n):
        slots.append(("R", i, -1))
    for i in range(n):
        for j in range(i + 1, n):
            slots.append(("S", i, j))
    return slots


def decode_mask(n: int, mask: int, slots: list[tuple[str, int, int]]) -> Net:
    adj_L = 0
    adj_R = 0
    adj = [0] * n
    for b, (kind, i, j) in enumerate(slots):
        if (mask >> b) & 1:
            if kind == "L":
                adj_L |= 1 << i
            elif kind == "R":
                adj_R |= 1 << i
            else:
                adj[i] |= 1 << j
                adj[j] |= 1 << i
    return Net(n, adj_L, adj_R, tuple(adj), False, f"exh{n}_{mask}")


def exhaustive_n(n: int) -> list[Net]:
    slots = _edge_slots(n)
    m = len(slots)
    seen = set()
    out = []
    vtot = n + 2
    emax = 3 * vtot - 6
    for mask in range(1 << m):
        adj_L = 0
        adj_R = 0
        adj = [0] * n
        e = 0
        bit = 1
        for kind, i, j in slots:
            if mask & bit:
                e += 1
                if kind == "L":
                    adj_L |= 1 << i
                elif kind == "R":
                    adj_R |= 1 << i
                else:
                    adj[i] |= 1 << j
                    adj[j] |= 1 << i
            bit <<= 1
        if e + 1 > emax:
            continue
        if adj_L == 0 or adj_R == 0:
            continue
        net = Net(n, adj_L, adj_R, tuple(adj), False, "")
        if not connected_carrier(net):
            continue
        key = canonical_key(net)
        if key in seen:
            continue
        seen.add(key)
        net.name = f"exh_n{n}_{len(out)}"
        out.append(net)
    planar = []
    for net in out:
        if is_plane_two_terminal(net):
            net.name = f"exh_n{net.n}_{len(planar)}"
            planar.append(net)
    return planar


def generate_spsp(max_n: int) -> list[Net]:
    by_n: dict[int, list[Net]] = {0: [path_net(0)]}
    seen: dict[int, set] = {0: {canonical_key(path_net(0))}}
    p1 = path_net(1)
    by_n[1] = [p1]
    seen[1] = {canonical_key(p1)}

    def add(net: Net) -> None:
        n = net.n
        if n > max_n or n < 0:
            return
        if net.lr_edge and n > 0:
            return
        if n > 0 and not connected_carrier(net):
            return
        key = canonical_key(net)
        bucket = seen.setdefault(n, set())
        if key in bucket:
            return
        bucket.add(key)
        by_n.setdefault(n, []).append(net)

    for n in range(2, max_n + 1):
        for n1 in range(0, n):
            n2 = n - 1 - n1
            if n2 < 0:
                continue
            for A in by_n.get(n1, []):
                for B in by_n.get(n2, []):
                    add(series(A, B, f"ser_n{n}"))
        for n1 in range(1, n):
            n2 = n - n1
            if n1 > n2:
                break
            for A in by_n.get(n1, []):
                for B in by_n.get(n2, []):
                    add(parallel(A, B, f"par_n{n}"))
    out = []
    for n, lst in sorted(by_n.items()):
        if n == 0:
            continue
        out.extend(lst)
    return out


def generate_wheatstone_family(max_n: int, sp_graphs: list[Net]) -> list[Net]:
    out = []
    seen = set()
    W = wheatstone()
    key = canonical_key(W)
    seen.add(key)
    out.append(W)
    sp_by_n: dict[int, list[Net]] = defaultdict(list)
    for g in sp_graphs:
        sp_by_n[g.n].append(g)
    for g in list(out):
        for n2 in range(0, max_n - g.n + 1):
            partners = [path_net(0)] if n2 == 0 else sp_by_n.get(n2, [])
            for H in partners:
                if n2 == 0:
                    continue
                for composer, tag in ((series, "s"), (parallel, "p")):
                    cands = []
                    if composer is series:
                        if g.n + H.n + 1 <= max_n:
                            cands.append(series(g, H, f"W{tag}{H.name}"))
                            cands.append(series(H, g, f"{H.name}{tag}W"))
                    else:
                        if g.n + H.n <= max_n:
                            cands.append(parallel(g, H, f"Wpar{H.name}"))
                    for c in cands:
                        k = canonical_key(c)
                        if k not in seen and connected_carrier(c) and not c.lr_edge:
                            seen.add(k)
                            out.append(c)
    return out


def generate_multipaths(max_n: int) -> list[Net]:
    out = []
    seen = set()

    def rec(remaining, acc):
        if remaining == 0 and acc:
            g = multipath(tuple(acc))
            k = canonical_key(g)
            if k not in seen:
                seen.add(k)
                out.append(g)
            return
        start = acc[-1] if acc else 1
        for L in range(start, remaining + 1):
            rec(remaining - L, acc + [L])

    for n in range(1, max_n + 1):
        rec(n, [])
    return out


def generate_grids(max_n: int) -> list[Net]:
    out = []
    seen = set()
    for rows in (2, 3, 4):
        for cols in range(2, max_n + 1):
            if rows * cols > max_n:
                break
            g = ladder(rows, cols)
            k = canonical_key(g)
            if k not in seen and connected_carrier(g):
                seen.add(k)
                out.append(g)
    return out


def generate_hidden(cores: list[Net], max_n: int) -> list[Net]:
    out = []
    seen = set()
    for core in cores:
        for left in (2, 3):
            for right in (2, 3):
                if left + right + core.n + 2 > max_n:
                    continue
                g = core
                if left:
                    g = series(path_net(max(left - 1, 0) if left else 0), g, "")
                if right:
                    g = series(g, path_net(max(right - 1, 0)), "")
                if g.n > max_n:
                    continue
                g.name = f"hide(L{left},{core.name},R{right})"
                if not connected_carrier(g) or g.lr_edge:
                    continue
                k = canonical_key(g)
                if k not in seen:
                    seen.add(k)
                    out.append(g)
    for core in cores:
        for L in range(1, 4):
            for R in range(1, 4):
                ntot = L + core.n + R + 2
                if ntot > max_n:
                    continue
                g = series(path_net(L), series(core, path_net(R), ""), f"P{L}-{core.name}-P{R}")
                if g.n > max_n or not connected_carrier(g):
                    continue
                k = canonical_key(g)
                if k not in seen:
                    seen.add(k)
                    out.append(g)
    return out


def analyze(net: Net) -> dict:
    table = safe_table(net)
    n = net.n
    s = S_coeffs(table, n)
    h2 = n - s[1] if n >= 1 else 0
    b2 = (comb(n, 2) - s[2]) if n >= 2 else 0
    exps = experiments(table, n)
    n1 = neighborhood_key(net, 1)
    n2 = neighborhood_key(net, 2)
    m1, m2 = successor_h2_moments(table, n)
    return {
        "n": n,
        "S": s,
        "H2": h2,
        "b2": b2,
        "neigh1": n1,
        "neigh2": n2,
        "exps": exps,
        "succ_m1": m1,
        "succ_m2": m2,
        "edges": edge_count(net),
        "name": net.name,
        "adj_L": net.adj_L,
        "adj_R": net.adj_R,
        "adj": net.adj,
        "lr_edge": net.lr_edge,
        "canon": canonical_key(net),
    }


def summary_tuple(rec: dict, radius: int) -> tuple:
    neigh = rec["neigh1"] if radius == 1 else rec["neigh2"]
    return (rec["n"], rec["S"], rec["H2"], rec["b2"], neigh)


def behavior_tuple(rec: dict) -> tuple:
    e = rec["exps"]
    return tuple(e[k] for k in FROZEN_EXPERIMENT_ORDER)


def first_split(a: dict, b: dict) -> Optional[str]:
    ea, eb = a["exps"], b["exps"]
    for k in FROZEN_EXPERIMENT_ORDER:
        if ea[k] != eb[k]:
            return k
    for k in ea:
        if k.startswith("ord_") and ea[k] != eb[k]:
            return k
    return None


def incidence_repr(net: Net) -> dict:
    n = net.n
    edges = []
    for i in range(n):
        if (net.adj_L >> i) & 1:
            edges.append(["L", i])
        if (net.adj_R >> i) & 1:
            edges.append(["R", i])
        for j in range(i + 1, n):
            if (net.adj[i] >> j) & 1:
                edges.append([i, j])
    if net.lr_edge:
        edges.append(["L", "R"])
    return {
        "switchable": n,
        "terminals": ["L", "R"],
        "edges": edges,
        "name": net.name,
        "adj_L": net.adj_L,
        "adj_R": net.adj_R,
        "adj": list(net.adj),
    }


def frac_text(x: Fraction) -> str:
    return f"{x.numerator}/{x.denominator}"


def self_checks() -> None:
    p1 = path_net(1)
    t = safe_table(p1)
    assert S_coeffs(t, 1) == (1, 0)
    assert connected_carrier(p1) and is_plane_two_terminal(p1)
    p2 = path_net(2)
    t = safe_table(p2)
    assert S_coeffs(t, 2) == (1, 2, 0), S_coeffs(t, 2)
    two = parallel(path_net(1), path_net(1))
    t = safe_table(two)
    assert S_coeffs(t, 2) == (1, 0, 0), S_coeffs(t, 2)
    W = wheatstone()
    t = safe_table(W)
    assert S_coeffs(t, 2) == (1, 0, 0)
    assert neighborhood_key(W, 1) != neighborhood_key(two, 1)
    s2 = series(path_net(1), path_net(1))
    assert s2.n == 3
    t = safe_table(s2)
    assert S_coeffs(t, 3) == (1, 3, 3, 0), S_coeffs(t, 3)
    ex = experiments(safe_table(p2), 2)
    assert ex["E1_c1"] == 0
    assert ex["E0_c1"] == 1
    # fast == slow
    for g in (p1, p2, two, W, s2):
        tbl = safe_table(g)
        a, b = experiments(tbl, g.n), experiments_fast(tbl, g.n)
        for k in FROZEN_EXPERIMENT_ORDER:
            assert a[k] == b[k], (g.name, k, a[k], b[k])
    print("self-checks passed", flush=True)


def main() -> int:
    self_checks()
    print("library restored; use run_full_search.py for the hunt", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())



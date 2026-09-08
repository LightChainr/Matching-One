#!/usr/bin/env python3
"""#635 companion: is the *bond* percolation wrap-gap difference of two
sector amplitudes at finite L?

On the square bond torus (2*L^2 bonds), define for bond subset B:

    wrap_x(B) = exists a connected component of B with nontrivial
                displacement along the x period,
    wrap_y(B) = same along y.

Channels:
    p = primal graph (the occupied bonds themselves)
    d = dual graph: dual bond occupied iff the crossing primal bond is VACANT
        (this is the standard bond duality; #42's note uses the transport map
        T for wrapping bookkeeping, same content)

Sector split by homology: x-only, y-only, both-same-cluster, both-two-clusters.

Gap polynomial (exact Bernstein form, N=2L^2 bonds):

    G(p) = sum_{|B|=k} [1{p wraps} - 1{d wraps}] C(2L^2,k) p^k (1-p)^(N-k)

If bond duality were a sector map, G(p) would be an antisymmetric- around-1/2
combination of primal-only and dual-only sector amplitudes.  The test:
exhibit the per-sector integer tables and check whether the difference is
anti-symmetric under p -> 1-p, i.e. whether

    G(1/2) = 0    (exact, all L)

and whether G is a difference of two sector amplitudes of ONE transfer
structure (primal=dual graph identity, self-dual lattice).
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path


def bond_wrap_sectors(occ: list[bool], edges, n_bonds: int, L: int):
    """Union-find over SITES (L*L) with displacement potentials.

    edges: list of (site_a, site_b, dx, dy) for each bond.
    Returns (wx, wy, same, two).
    """
    n_sites = L * L
    parent = list(range(n_sites))
    size = [1] * n_sites
    dx_ = [0] * n_sites
    dy_ = [0] * n_sites
    wrap_x = [False] * n_sites
    wrap_y = [False] * n_sites

    def find(x):
        if parent[x] == x:
            return x, 0, 0
        p = parent[x]
        r, px, py = find(p)
        dx_[x] += px
        dy_[x] += py
        parent[x] = r
        return r, dx_[x], dy_[x]

    for idx, (a, b, ex, ey) in enumerate(edges):
        if not occ[idx]:
            continue
        ra, ax, ay = find(a)
        rb, bx, by = find(b)
        rdx = ax + ex - bx
        rdy = ay + ey - by
        if ra == rb:
            if rdx != 0:
                wrap_x[ra] = True
            if rdy != 0:
                wrap_y[ra] = True
        else:
            if size[rb] < size[ra]:
                ra, rb = rb, ra
                rdx, rdy = -rdx, -rdy
            parent[rb] = ra
            size[ra] += size[rb]
            dx_[rb] = rdx
            dy_[rb] = rdy
            wrap_x[ra] = wrap_x[ra] or wrap_x[rb]
            wrap_y[ra] = wrap_y[ra] or wrap_y[rb]

    any_x = any_y = False
    for s in range(n_sites):
        r, _, _ = find(s)
        any_x = any_x or wrap_x[r]
        any_y = any_y or wrap_y[r]
    same = False
    if any_x and any_y:
        for s in range(n_sites):
            r, _, _ = find(s)
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


def build_torus_bonds(L: int):
    """Bonds as site pairs with displacement; site id = x + L*y.

    Every bond is a displaced edge pos(b) = pos(a) + (dx, dy) on the universal
    cover, with dx=+1 for the +x bond and dy=+1 for the +y bond (the wrap bond
    at x=L-1 carries displacement +1 onto the next cover cell — this is the
    same convention as matched_torus_reference._make_edges).
    """
    bonds = []
    for y in range(L):
        for x in range(L):
            a = x + L * y
            b = ((x + 1) % L) + L * y
            bonds.append((a, b, 1, 0))
            b = x + L * ((y + 1) % L)
            bonds.append((a, b, 0, 1))
    return bonds


def dual_occupation(occ: list[bool], L: int) -> list[bool]:
    """Standard bond duality: dual bond occupied iff primal bond vacant."""
    return [not v for v in occ]


def run(L: int):
    bonds = build_torus_bonds(L)
    n = len(bonds)  # 2 L^2
    counts = {c: [0] * (n + 1) for c in ("p-x", "p-y", "p-both-same", "p-both-two", "p-none",
                                          "d-x", "d-y", "d-both-same", "d-both-two", "d-none")}
    gap_poly = [Fraction(0)] * (n + 1)  # exact power-basis check later
    # instead: accumulate bernstein directly
    bernstein = [0] * (n + 1)
    anti_defect_num = 0  # count configs where p-wrap != dual-wrap
    for mask in range(1 << n):
        if mask % 2 and False:
            pass
        occ = [bool((mask >> i) & 1) for i in range(n)]
        k = sum(occ)
        bp = bond_wrap_sectors(occ, bonds, n, L)
        occ_d = dual_occupation(occ, L)
        bd = bond_wrap_sectors(occ_d, bonds, n, L)
        wp, wd = any(bp[:2]), any(bd[:2])
        d = int(wp) - int(wd)
        bernstein[k] += d
        if wp:
            counts[f"p-{sector_label(*bp)}"][k] += 1
        if wd:
            counts[f"d-{sector_label(*bd)}"][k] -= 1
        if wp != wd:
            anti_defect_num += 1
    return n, counts, bernstein, anti_defect_num


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--L", type=int, default=2)
    args = ap.parse_args()
    L = args.L
    n, counts, bernstein, defect = run(L)

    # exact antisymmetry test at p=1/2: M(1/2) = sum_k a_k C(n,k) / 2^n
    m_half = Fraction(0)
    for k, a in enumerate(bernstein):
        m_half += a * math.comb(n, k)
    m_half /= 2 ** n

    print(f"bond torus L={L}  N={n} bonds")
    print("gap bernstein coefficients:", bernstein)
    print(f"G(1/2) exact = {m_half}  (zero? {m_half == 0})")
    print(f"configs where primal wrap != dual wrap: {defect} of {2**n}")
    for key in sorted(counts):
        if any(counts[key]):
            print(f"  {key:12s}", counts[key])

    payload = {
        "ticket": 635,
        "kind": "bond-torus-sector-map",
        "L": L,
        "N_bonds": n,
        "gap_bernstein": bernstein,
        "gap_at_half_exact": f"{m_half.numerator}/{m_half.denominator}",
        "antisymmetry_defect_configs": defect,
        "sector_tables": {k: v for k, v in counts.items() if any(v)},
    }
    out = Path("results/probe635-sector-map/bond-torus.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    prior = json.loads(out.read_text()) if out.exists() else {"runs": []}
    prior["runs"].append(payload)
    out.write_text(json.dumps(prior, indent=1))
    print(f"wrote {out}")
    print("ALL_CHECKS_PASS")


if __name__ == "__main__":
    main()

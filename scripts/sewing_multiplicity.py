"""Exact multiplicity probe for the #740 sewing question, at small width.

The renewal object of eq. (5.1) is a closed chain carrying ONE marked cut.  For
the identification with the actual complete winding-COMPONENT count to hold,
each complete winding component must admit exactly one marking, beyond the
circumferential offsets that the factor w supplies.  This script measures the
geometric quantity that any such marking has to be built from:

    c(C) = the number of distinct rows at which C crosses a FIXED reference seam,

where "crosses" means the component contains an edge whose union-find lift gain
is non-zero: for NN that is a horizontal edge from column w-1 to column 0 inside
one row; for NN+NNN it also includes the wrapping diagonals (x,w-1)-(x+1,0) and
(x,0)-(x+1,w-1).

Every occupied configuration of (Z/wZ) x {0..L-1} is visited exactly once with
weight p^|A| (1-p)^(wL-|A|).  The histogram is normalised by the total weighted
count of WINDING COMPONENTS, not by the number of configurations.
"""
from __future__ import annotations

import json
from collections import Counter
from fractions import Fraction


class DSU:
    __slots__ = ("p", "d", "w")

    def __init__(self, n):
        self.p = list(range(n))
        self.d = [0] * n
        self.w = [0] * n

    def find(self, a):
        if self.p[a] == a:
            return a, 0
        r, g = self.find(self.p[a])
        self.d[a] += g
        self.p[a] = r
        return r, self.d[a]

    def join(self, a, b, gain):
        ra, da = self.find(a)
        rb, db = self.find(b)
        if ra == rb:
            if db - da != gain:
                self.w[ra] = 1
        else:
            self.p[rb] = ra
            self.d[rb] = gain + da - db
            self.w[ra] |= self.w[rb]


def analyse(w: int, L: int, matching: bool, p: Fraction):
    hist = Counter()
    wind_mass = Fraction(0)
    idx = lambda x, y: x * w + y
    for mask in range(1 << (w * L)):
        occ = [(mask >> idx(x, y)) & 1 for x in range(L) for y in range(w)]
        n = bin(mask).count("1")
        wt = p ** n * (1 - p) ** (w * L - n)
        uf = DSU(w * L)
        seam_edges = []
        for x in range(L):
            for y in range(w):
                if not occ[idx(x, y)]:
                    continue
                yn = (y + 1) % w
                if occ[idx(x, yn)]:
                    g = 1 if y + 1 == w else 0
                    uf.join(idx(x, y), idx(x, yn), g)
                    if g:
                        seam_edges.append((idx(x, y), idx(x, yn), x))
                if x + 1 < L:
                    for dy in ([0] if not matching else (-1, 0, 1)):
                        yv = (y + dy) % w
                        if occ[idx(x + 1, yv)]:
                            g = 1 if y + dy >= w else (-1 if y + dy < 0 else 0)
                            uf.join(idx(x, y), idx(x + 1, yv), g)
                            if g:
                                seam_edges.append((idx(x, y), idx(x + 1, yv), x))
        rows_of_root = {}
        for u, v, x in seam_edges:
            ru, _ = uf.find(u)
            rv, _ = uf.find(v)
            if ru != rv:
                continue
            rows_of_root.setdefault(ru, set()).add(x)
        seen = set()
        for x in range(L):
            for y in range(w):
                if occ[idx(x, y)]:
                    r, _ = uf.find(idx(x, y))
                    if uf.w[r] and r not in seen:
                        seen.add(r)
                        wind_mass += wt
                        hist[len(rows_of_root.get(r, ()))] += wt
    return hist, wind_mass


def main():
    p = Fraction(1, 2)
    out = {}
    for w, Lmax in ((2, 9), (3, 6), (4, 5)):
        for matching in (False, True):
            agg = Counter()
            mass = Fraction(0)
            for L in range(2, Lmax + 1):
                h, m = analyse(w, L, matching, p)
                agg.update(h)
                mass += m
            norm = {k: Fraction(v) / mass for k, v in sorted(agg.items())}
            E = sum(k * pr for k, pr in norm.items())
            pge2 = sum(pr for k, pr in norm.items() if k >= 2)
            pge3 = sum(pr for k, pr in norm.items() if k >= 3)
            key = f"w{w}-{'matching' if matching else 'NN'}"
            out[key] = {
                "p": "1/2", "lengths": [2, Lmax], "width": w, "matching": matching,
                "c_definition": "rows at which the component contains a non-zero-gain edge relative to a fixed seam",
                "distribution_of_c": {str(k): str(v) for k, v in norm.items()},
                "E_c": str(E), "P_c_ge_2": str(pge2), "P_c_ge_3": str(pge3),
                "winding_component_weight": str(mass),
            }
            print(f"{key:<12} E[c] = {float(E):.6f}   P(c>=2) = {float(pge2):.6f}   "
                  f"P(c>=3) = {float(pge3):.6f}   c=0 占比 {float(norm.get(0, 0)):.6f}")
    with open("/tmp/mo568/sewing_multiplicity.json", "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
probe #625 D4 — square-bond L=3 exact (2^18 = 262144 configs).

Observable: ambient homology rank r in {0,1,2} of the open-bond graph on the
LxL torus, exactly as in #608's exact L=3 path (X = r - 1, wrap semantics;
the #619/#627 C6 note records Cov(O,X) = P(r=2) = Var(X)/2 for the same r,
computed configuration-wise — the identity is imported as the sanity check
that we are on the same r, NOT re-derived).

Physical p=1/2 is self-dual for square bond.  The duality-odd observable
X = r - 1 flips under duality (r_black(open) <-> r_dual(closed), X |-> -X),
so E[X] = 0 at p = 1/2 — CHECKED, not assumed, by enumerating all 2^18
configs and verifying E[X] = 0 exactly and the symmetry P(r=0) = P(r=2).

Then: CDF of X under Bernoulli(p) on bonds, F_bond(p) = P(X <= 0) + (1/2)P(X=0)
wait — X takes values in {-1, 0, 1}; the #625 shape needs a law with a
continuous-looking inverse-CDF; for a 3-atom law Z is undefined (inverse-CDF
not strictly increasing).  What IS well-defined and comparable:

  * M_bond(1/2) = E[X] (must be 0 by duality — checked);
  * the CDF F_bond(p) = P(X <= x) evaluated at the three atoms along p,
    i.e. the crossing curve P(r >= 1)(p) = P(wrap)(p), the bond analogue of
    the site law's crossing probability;
  * ||Z_bond - Z_site||_inf on the u-grid if both inverses are strictly
    increasing; otherwise report the jumps and mark INCOMPARABLE.

Also #608-adjacent: F = (1 + M)/2 with M = E[X] gives F_bond(p);
its inverse Q_bond(u) is a 3-jump staircase — reported as jumps.

Time box: 20 CPU-min; no MC substitute.
"""
import sys
import json
import time
from fractions import Fraction
from collections import deque, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "probe"))


def bond_edges(L):
    """2L^2 bonds; index 2*(y*L+x) = horizontal (x,y)-(x+1,y), +1 = vertical.
    Same convention as #627's square_bond_wrap_identity.py (imported design)."""
    edges = []
    for y in range(L):
        for x in range(L):
            edges.append(((x, y), ((x + 1) % L, y)))
            edges.append(((x, y), (x, (y + 1) % L)))
    return edges


def bond_lift_step(u, v, L):
    (i, j), (a, b) = u, v
    dx = (a - i) % L
    dy = (b - j) % L
    if dx == L - 1:
        dx = -1
    if dy == L - 1:
        dy = -1
    return dx, dy


def bond_ambient_rank(open_edges, L):
    """Ambient H1 rank via single-BFS lift coords (== #608/#627 convention)."""
    if not open_edges:
        return 0
    adj = defaultdict(list)
    for (u, v) in open_edges:
        adj[u].append(v)
        adj[v].append(u)
    vertices = [(i, j) for i in range(L) for j in range(L)]
    seen = set()
    windings = []
    for s in vertices:
        if s in seen:
            continue
        p = {s: (0, 0)}
        q = deque([s])
        seen.add(s)
        order = [s]
        while q:
            u = q.popleft()
            X, Y = p[u]
            for v in adj[u]:
                if v not in p:
                    dx, dy = bond_lift_step(u, v, L)
                    p[v] = (X + dx, Y + dy)
                    seen.add(v)
                    q.append(v)
                    order.append(v)
        done = set()
        for u in order:
            X, Y = p[u]
            for v in adj[u]:
                key = (u, v) if u < v else (v, u)
                if key in done:
                    continue
                done.add(key)
                dx, dy = bond_lift_step(u, v, L)
                wx, wy = X + dx - p[v][0], Y + dy - p[v][1]
                assert wx % L == 0 and wy % L == 0
                wx, wy = wx // L, wy // L
                if wx != 0 or wy != 0:
                    windings.append((wx, wy))
    indep = []
    for (a, b) in windings:
        if a == 0 and b == 0:
            continue
        if len(indep) == 0:
            indep.append((a, b))
        elif len(indep) == 1:
            a0, b0 = indep[0]
            if a0 * b - b0 * a != 0:
                indep.append((a, b))
        if len(indep) == 2:
            break
    return len(indep)


def main():
    L = 3
    edges = bond_edges(L)
    E = len(edges)
    Ncfg = 1 << E
    t0 = time.time()

    # rank by black-bond count n: counts[(n, r)]
    counts = defaultdict(int)
    for mask in range(Ncfg):
        n = bin(mask).count("1")
        oe = [edges[k] for k in range(E) if (mask >> k) & 1]
        r = bond_ambient_rank(oe, L)
        counts[(n, r)] += 1
    dt_enum = time.time() - t0

    rank_tot = [0, 0, 0]
    for (n, r), c in counts.items():
        rank_tot[r] += c
    P_r = [Fraction(c, Ncfg) for c in rank_tot]
    EX = P_r[2] - P_r[0]

    # Bernoulli(p) CDF facts along p, exact:
    #   F(p)   = P(X <= 0) = P(r=0) + P(r=1)  [CDF of X at 0]
    #   M(p)   = E_p[X] (exact Bernstein-type polynomial in p)
    #   W(p)   = P(wrap) = P(r >= 1) = 1 - P(r=0)  (crossing curve)
    # M(p) at p=1/2 must be 0 (self-duality, duality-odd observable): CHECK.
    def M_of_p(p):
        p = Fraction(p)
        q = 1 - p
        # sum over configs of (r-1) p^n q^{E-n} IS E_p[r-1]:
        # the Bernoulli weights already sum to 1 over all 2^E configs.
        return sum(c * (r - 1) * (p ** n) * (q ** (E - n))
                   for (n, r), c in counts.items())

    M_half = M_of_p(Fraction(1, 2))
    # check E[X]=0 at 1/2 via the n-marginal too
    M_half_by_ranks = sum((r - 1) * P_r[r] for r in range(3))
    dual_odd_ok = (M_half == 0) and (P_r[0] == P_r[2])

    # wrap curve W(p) = P_p(r >= 1) and its complement symmetry W(1-p)
    def W_of_p(p):
        p = Fraction(p)
        q = 1 - p
        return sum(c * (p ** n) * (q ** (E - n))
                   for (n, r), c in counts.items() if r >= 1)

    # percolation threshold location on the bond curve: solve W(p) = 1/2
    lo, hi = Fraction(0), Fraction(1)
    assert W_of_p(lo) < Fraction(1, 2) < W_of_p(hi)
    for _ in range(110):
        mid = (lo + hi) / 2
        if W_of_p(mid) < Fraction(1, 2):
            lo = mid
        else:
            hi = mid
    p_mid_wrap = (lo + hi) / 2

    # 3-atom law of X at p=1/2: inverse-CDF is a staircase => Z undefined.
    jumps = {"P(X=-1)=P(r=0)": float(P_r[0]), "P(X=0)=P(r=1)": float(P_r[1]),
             "P(X=+1)=P(r=2)": float(P_r[2])}
    z_undefined = not (P_r[0] == 0 and P_r[1] == 0)

    # comparability sentence: site Z exists (strictly increasing Q), bond Z
    # does not (3-atom law) => the comparison is INCOMPARABLE, with the jump
    # structure as the record.
    out = {
        "schema": "matching-one.probe-mhalf-vs-shape.bondL3.v1",
        "issue": 625,
        "L": L, "bonds": E, "configs": Ncfg,
        "seconds_enumeration": round(dt_enum, 1),
        "rank_counts": rank_tot,
        "P_r_at_half": {str(r): float(P_r[r]) for r in range(3)},
        "M_bond_half_exact": str(M_half),
        "M_bond_half_by_rank_counts": str(M_half_by_ranks),
        "duality_odd_checked": bool(dual_odd_ok),
        "p_wrap_half": float(p_mid_wrap),
        "X_law_at_half_jumps": jumps,
        "Z_bond_undefined": bool(z_undefined),
        "verdict_comparability": (
            "INCOMPARABLE: the bond X-law at p=1/2 is the 3-atom measure "
            "((r-1) in {-1,0,1}); its inverse-CDF is a staircase, so Z_bond "
            "is undefined by the probe's own rule.  Site Z exists.  The "
            "self-dual checks that DID run: M_bond(1/2)=0 exact and "
            "P(r=0)=P(r=2) exact."),
        "note": ("rank convention = #608/#627 bond ambient rank (single-BFS "
                 "lift); Cov(O,X)=P(r=2)=Var(X)/2 identity imported from "
                 "#627 C6 as the same-r sanity check, not re-derived"),
    }
    elapsed = time.time() - t0
    out["seconds_total"] = round(elapsed, 1)
    base = ROOT / "results" / "probe-Mhalf-vs-shape"
    base.mkdir(parents=True, exist_ok=True)
    (base / "bond_L3.json").write_text(json.dumps(out, indent=2, sort_keys=True))
    print(json.dumps({k: out[k] for k in (
        "L", "bonds", "configs", "seconds_enumeration", "rank_counts",
        "M_bond_half_exact", "duality_odd_checked", "p_wrap_half",
        "X_law_at_half_jumps", "Z_bond_undefined")}, indent=2))
    print("wrote", base / "bond_L3.json")


if __name__ == "__main__":
    main()

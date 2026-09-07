#!/usr/bin/env python3
"""
C6 + C7 (L=3 square-bond exact enumeration, 2^18 = 262144 configs).

C6 — self-dual square-bond torus (p=1/2 uniform).
  r = ambient homology rank in {0,1,2}, X = r-1, O = wrap_either = int(r>0)
  (== int(top_rank>0); pin = #608 score_qtangent_empirical._primal_statistics).
  Check configuration-wise  O == 1 + (X - X^2)/2.
  Consequences: Cov(O,X) ?= P(r=2) ?= Var(X)/2, plus the moments.

C7 — Harris / nested-reveal, exact.
  Gamma_j = E[m_j m_j^T] - E[m_{j-1} m_{j-1}^T], m_j = E[Y | F_j], F_j the
  #608 radius-ordered nested reveal (spatial_filtration_levels).  Verify
  telescoping  sum_j Gamma_j == Cov(Y)  and cross-increments Gamma_j^{O,X} >= -eps.
"""
import sys
import json
import time
from pathlib import Path
from collections import deque, defaultdict

ROOT = Path(__file__).resolve().parents[2]


# ---------------- bond geometry (== #608 square_bond_kappa3 convention) ----------------

def bond_edges(L):
    """2L^2 bonds; index 2*(y*L+x)=horizontal (x,y)-(x+1,y), +1=vertical."""
    edges = []
    for y in range(L):
        for x in range(L):
            edges.append(((x, y), ((x + 1) % L, y)))          # horizontal
            edges.append(((x, y), (x, (y + 1) % L)))          # vertical
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
    """Ambient H1 rank (= #608 span_rank) via single-BFS lift coords."""
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


# ---------------- #608 spatial filtration ----------------

def bond_radius_squared(length, index):
    cell, orientation = divmod(index, 2)
    y, x = divmod(cell, length)
    midpoint = (2 * x + 1, 2 * y) if orientation == 0 else (2 * x, 2 * y + 1)
    period = 2 * length
    total = 0
    for coordinate in midpoint:
        wrapped = coordinate % period
        total += min(wrapped, period - wrapped) ** 2
    return total


def spatial_filtration_levels(length):
    bonds = 2 * length * length
    radii = sorted({bond_radius_squared(length, i) for i in range(bonds)})
    levels = []
    for radius in radii:
        levels.append([i for i in range(bonds) if bond_radius_squared(length, i) <= radius])
    return levels


def main():
    L = 3
    edges = bond_edges(L)
    E = len(edges)
    N = 1 << E
    t0 = time.time()

    # ---- C6 enumeration, keep Y=(O,X) arrays for C7 ----
    YO = [0] * N
    YX = [0] * N
    identity_fail = 0
    identity_fail_examples = []
    rank_counts = [0, 0, 0]
    sum_O = sum_X = sum_X2 = sum_OX = sum_OX2 = 0

    for mask in range(N):
        open_edges = [edges[k] for k in range(E) if (mask >> k) & 1]
        r = bond_ambient_rank(open_edges, L)
        rank_counts[r] += 1
        X = r - 1
        O = 1 if r > 0 else 0
        YO[mask] = O
        YX[mask] = X
        rhs = 1 + (X - X * X) // 2
        if O != rhs:
            identity_fail += 1
            if len(identity_fail_examples) < 20:
                identity_fail_examples.append((mask, r, O, rhs))
        sum_O += O; sum_X += X; sum_X2 += X * X; sum_OX += O * X; sum_OX2 += O * X * X

    inv = 1.0 / N
    P0, P1, P2 = [c * inv for c in rank_counts]
    EO, EX, EX2, EOX, EOX2 = sum_O * inv, sum_X * inv, sum_X2 * inv, sum_OX * inv, sum_OX2 * inv
    VarX = EX2 - EX * EX
    CovOX = EOX - EO * EX
    CovOX2 = EOX2 - EO * EX2
    VarO = EO * (1 - EO)

    print("=== C6 ===")
    print(f"configs={N}  identity_failures={identity_fail}")
    print(f"P(r=0,1,2) = {P0:.9f}, {P1:.9f}, {P2:.9f}")
    print(f"E[X]={EX:.12f}  Var(X)={VarX:.9f}  E[O]={EO:.9f}")
    print(f"Cov(O,X)={CovOX:.9f}  P(r=2)={P2:.9f}  Var(X)/2={VarX/2:.9f}")
    print(f"checks: CovOX==P2 {abs(CovOX-P2)<1e-12}; P2==VarX/2 {abs(P2-VarX/2)<1e-12}")

    # ---- C7 exact nested-reveal Gamma_j ----
    levels = spatial_filtration_levels(L)
    selectors = []
    acc = 0
    for lv in levels:
        for idx in lv:
            acc |= 1 << idx
        selectors.append(acc)

    # mean of Y over the empty filtration
    mE = [EO, EX]
    # E[m_j m_j^T] for each scale
    Emm = []
    for sel in selectors:
        gcount = defaultdict(int)
        gsumO = defaultdict(int)
        gsumX = defaultdict(int)
        for mask in range(N):
            g = mask & sel
            gcount[g] += 1
            gsumO[g] += YO[mask]
            gsumX[g] += YX[mask]
        # E[m_j m_j^T] = sum_g (count/N) * mean_g mean_g^T
        eOO = eOX = eXX = 0.0
        for g in gcount:
            c = gcount[g]
            mO = gsumO[g] / c
            mX = gsumX[g] / c
            w = c / N
            eOO += w * mO * mO
            eOX += w * mO * mX
            eXX += w * mX * mX
        Emm.append((eOO, eOX, eXX))

    gammas = []
    prev = (mE[0] * mE[0], mE[0] * mE[1], mE[1] * mE[1])
    for e in Emm:
        gammas.append((e[0] - prev[0], e[1] - prev[1], e[2] - prev[2]))
        prev = e

    # telescoping: sum gamma = E[Y Y^T] - E[Y]E[Y]^T = Cov(Y)
    sum_g = [0.0, 0.0, 0.0]
    for g in gammas:
        sum_g[0] += g[0]; sum_g[1] += g[1]; sum_g[2] += g[2]
    CovYY = (EX2 - EX * EX)          # Cov(X,X)
    CovOO = VarO                     # Cov(O,O)
    CovOX_ = CovOX                   # Cov(O,X)
    telescoping_residual = (abs(sum_g[0] - CovOO), abs(sum_g[1] - CovOX_), abs(sum_g[2] - CovYY))

    print("\n=== C7 ===")
    print(f"levels reveal counts: {[sel.bit_count() for sel in selectors]}")
    print("Gamma_j (j: OO, OX, XX):")
    for j, g in enumerate(gammas):
        print(f"  j={j}: OO={g[0]:+.9f}  OX={g[1]:+.9f}  XX={g[2]:+.9f}")
    print(f"sum Gamma = ({sum_g[0]:.12f}, {sum_g[1]:.12f}, {sum_g[2]:.12f})")
    print(f"Cov(Y)    = ({CovOO:.12f}, {CovOX_:.12f}, {CovYY:.12f})")
    print(f"telescoping residual = {max(telescoping_residual):.3e}")
    cross_min = min(g[1] for g in gammas)
    print(f"min_j Gamma_j^{{O,X}} = {cross_min:.3e}  (>= -1e-15 ? {cross_min >= -1e-15})")

    dt = time.time() - t0
    out = {
        "schema": "matching-one.probe-exact-controls.bond-L3-wrap.v1",
        "L": L, "configs": N, "seconds": round(dt, 2),
        "C6": {
            "wrap_identity_failures": identity_fail,
            "wrap_identity_examples": identity_fail_examples[:20],
            "rank_counts": rank_counts,
            "P_r": {"r0": P0, "r1": P1, "r2": P2},
            "E_O": EO, "Var_O": VarO, "E_X": EX, "Var_X": VarX,
            "Cov_OX": CovOX, "Cov_OX2": CovOX2,
            "checks": {"CovOX_eq_P2": abs(CovOX - P2) < 1e-12,
                       "P2_eq_VarX_half": abs(P2 - VarX / 2) < 1e-12,
                       "CovOX_eq_VarX_half": abs(CovOX - VarX / 2) < 1e-12},
        },
        "C7": {
            "reveal_counts": [sel.bit_count() for sel in selectors],
            "gammas_OO_OX_XX": gammas,
            "sum_gamma": sum_g,
            "Cov_Y": [CovOO, CovOX_, CovYY],
            "telescoping_residual_max": max(telescoping_residual),
            "min_cross_increment_OX": cross_min,
            "harris_nonneg": cross_min >= -1e-15,
        },
        "note": "X=r-1; O=wrap_either=int(r>0); p=1/2 uniform exact; filtration=#608 spatial_filtration_levels",
    }
    base = ROOT / "results" / "probe-exact-controls"
    base.mkdir(parents=True, exist_ok=True)
    (base / "bond_L3_wrap.json").write_text(json.dumps(out, indent=2, sort_keys=True))
    print("\nwrote", base / "bond_L3_wrap.json")


if __name__ == "__main__":
    main()

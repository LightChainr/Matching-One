#!/usr/bin/env python3
"""verify802src -- C4: why does the staggered row field have a first-order
response only on ODD torus length?

Everything is exact full enumeration on small tori.  For each configuration we
compute
  * K_y = #black sites in row y,  H_s = sum_y s_y K_y  with s_y = (+1)^y
    (the team's / note's staggering; sign +1 on even rows, -1 on odd rows),
  * the ambient homology rank r in {0,1,2} of the black 4-neighbour subgraph,
  * the Bernoulli weight w0 = p^K (1-p)^(N-K),
and then report
  * <H_s>_{w0}                       (the global mean)
  * E[H_s | r] for r=0,1,2           (the rank-conditional means)
  * b_g = (E[H_s|0]-E[H_s|2])/2      (the conditional source scoring used)
  * the "seam defect"  H_s(T cfg) + H_s(cfg)   under the y-shift T,
    which measures how far H_s is from being shift-ODD.

Outputs /workspace/v802src/out/v2_c4.json
"""
from __future__ import annotations
import json
import math
import os
from collections import defaultdict

OUT = "/workspace/v802src/out"


def rank_of(omega, w, m):
    N = w * m
    occupied = [i for i in range(N) if (omega >> i) & 1]
    occ = set(occupied)
    pot = {}
    cyc = []
    for start in occupied:
        if start in pot:
            continue
        pot[start] = (0, 0)
        stack = [start]
        while stack:
            u = stack.pop()
            i, y = u % w, u // w
            pu = pot[u]
            for (j, yy, di, dy) in (((i + 1) % w, y, 1, 0),
                                    ((i - 1) % w, y, -1, 0),
                                    (i, (y + 1) % m, 0, 1),
                                    (i, (y - 1) % m, 0, -1)):
                v = yy * w + j
                if v not in occ:
                    continue
                cand = (pu[0] + di, pu[1] + dy)
                if v not in pot:
                    pot[v] = cand
                    stack.append(v)
                else:
                    dv = pot[v]
                    cyc.append((cand[0] - dv[0], cand[1] - dv[1]))
    nz = [c for c in cyc if c != (0, 0)]
    if not nz:
        return 0
    a, b = nz[0]
    for c, d in nz[1:]:
        if a * d - b * c != 0:
            return 2
    return 1


def shift_y(omega, w, m):
    """Cyclic shift of the configuration by one row (y -> y+1)."""
    out = 0
    for y in range(m):
        for i in range(w):
            if (omega >> (y * w + i)) & 1:
                out |= 1 << (((y + 1) % m) * w + i)
    return out


def Hk(omega, w, m):
    """K_y per row (list) and H_s = sum_y (-1)^y K_y."""
    Ky = []
    for y in range(m):
        Ky.append(sum((omega >> (y * w + i)) & 1 for i in range(w)))
    Hs = sum((1 if y % 2 == 0 else -1) * Ky[y] for y in range(m))
    return Ky, Hs


def run_case(w, m, p):
    N = w * m
    n = 1 << N
    Z = defaultdict(float)          # rank -> sum w0
    SH = defaultdict(float)         # rank -> sum w0 * H_s
    SH2 = defaultdict(float)        # rank -> sum w0 * H_s^2
    tot = 0.0
    seam_nonzero = 0                # configs with H_s(T cfg) != -H_s(cfg)
    seam_absmax = 0
    shift_invariance_viol = 0       # w0(T cfg) != w0(cfg)?
    for omega in range(n):
        K = bin(omega).count("1")
        w0 = p ** K * (1 - p) ** (N - K)
        r = rank_of(omega, w, m)
        _, Hs = Hk(omega, w, m)
        Z[r] += w0
        SH[r] += w0 * Hs
        SH2[r] += w0 * Hs * Hs
        tot += w0
        o2 = shift_y(omega, w, m)
        K2 = bin(o2).count("1")
        _, Hs2 = Hk(o2, w, m)
        if K2 != K:
            shift_invariance_viol += 1
        d = Hs2 + Hs
        if d != 0:
            seam_nonzero += 1
        seam_absmax = max(seam_absmax, abs(d))
    P = {r: Z[r] / tot for r in Z}
    EH = {r: (SH[r] / Z[r]) if Z[r] > 0 else float("nan") for r in Z}
    VarH = {r: (SH2[r] / Z[r] - EH[r] ** 2) if Z[r] > 0 else float("nan") for r in Z}
    Hbar = sum(SH[r] for r in SH) / tot
    # theory: <H_s> = w p * sum_y s_y  = w p * (1 if m odd else 0)
    summ = sum(1 if y % 2 == 0 else -1 for y in range(m))
    rec = {
        "w": w, "m": m, "p": p, "N": N, "n_configs": n,
        "period_parity": "odd" if m % 2 else "even",
        "sum_y_pm1": summ,
        "theory_Hbar": w * p * summ,
        "Hbar": Hbar,
        "Hbar_minus_theory": Hbar - w * p * summ,
        "P_r": {str(r): P[r] for r in sorted(P)},
        "EH_by_rank": {str(r): EH[r] for r in sorted(EH)},
        "VarH_by_rank": {str(r): VarH[r] for r in sorted(VarH)},
        "b_g_cond_score": 0.5 * (EH.get(0, float("nan")) - EH.get(2, float("nan"))),
        "base_translation_invariance_violations": shift_invariance_viol,
        "n_configs_where_Hs_T_is_not_minus_Hs": seam_nonzero,
        "max_abs_seam_defect_Hs_T_plus_Hs": seam_absmax,
    }
    print("C4 w=%d m=%d (%s) p=%.4f: <H_s>=%.10f (theory %.10f)  "
          "E[H|0]=%.6f E[H|2]=%.6f  b_g=%+.6e  seam_nonzero=%d/%d  "
          "w0_T_inv_viol=%d"
          % (w, m, rec["period_parity"], p, Hbar, w * p * summ,
             EH.get(0, float("nan")), EH.get(2, float("nan")),
             rec["b_g_cond_score"], seam_nonzero, n, shift_invariance_viol),
          flush=True)
    return rec


def main():
    os.makedirs(OUT, exist_ok=True)
    cases = []
    for (w, m) in [(3, 3), (4, 4), (3, 4), (4, 3), (3, 5)]:
        for p in [0.4, 0.59274605079210]:
            cases.append(run_case(w, m, p))
    with open(OUT + "/v2_c4.json", "w") as f:
        json.dump({"cases": cases}, f, indent=1)
    print("-> %s/v2_c4.json" % OUT)


if __name__ == "__main__":
    main()

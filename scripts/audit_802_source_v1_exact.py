#!/usr/bin/env python3
"""verify802src -- C1 / C2 / C3 independent verification (exact + full precision).

Nothing here re-runs the w=4..8 transfer engines; C2 only post-processes the
team's own printed JSON values (its stated error class).  All exact claims use
Fraction arithmetic and full enumeration of small tori.

Outputs: /workspace/v802src/out/v1_exact.json
"""
from __future__ import annotations
import json
import math
import itertools
from fractions import Fraction as F

IN = "/workspace/v802src/in"
OUT = "/workspace/v802src/out"
PC_S = "0.5927460507921"


# ------------------------------------------------------------------ helpers
def torus_configs(w, m):
    """Yield (K, H_B, H_W) for every configuration of a w x m torus.

    H_B / H_W = total number of horizontally (cyclic, within-row) adjacent
    BLACK / WHITE occupied pairs, summed over rows.  Same convention as the
    team's s4_torus_bruteforce.py.
    """
    N = w * m
    for omega in range(1 << N):
        K = 0
        HB = 0
        HW = 0
        for y in range(m):
            base = y * w
            row = [(omega >> (base + i)) & 1 for i in range(w)]
            K += sum(row)
            for i in range(w):
                a = row[i]
                b = row[(i + 1) % w]
                if a and b:
                    HB += 1
                elif (not a) and (not b):
                    HW += 1
        yield K, HB, HW


def row_marks(w):
    """All single-row masks with (K, H_B, H_W) for the row of length w."""
    for mask in range(1 << w):
        K = bin(mask).count("1")
        HB = 0
        HW = 0
        for j in range(w):
            a = (mask >> j) & 1
            b = (mask >> ((j + 1) % w)) & 1
            if a and b:
                HB += 1
            elif (not a) and (not b):
                HW += 1
        yield K, HB, HW


def logit_shift_peff(p, Y):
    """Exact p_eff with logit(p_eff) = logit(p) - 2*log(Y)  (Y = exp(g)).

    p_eff = p / (p + (1-p) Y^2).
    """
    return p / (p + (1 - p) * Y * Y)


# --------------------------------------------------------------- part C1a
def c1_identity():
    out = {}
    for (w, m) in [(3, 3), (4, 4), (3, 4)]:
        N = w * m
        viol = 0
        worst = 0
        n = 0
        for K, HB, HW in torus_configs(w, m):
            n += 1
            rhs = N - 2 * K + HB
            d = abs(HW - rhs)
            if d != 0:
                viol += 1
                worst = max(worst, d)
        out["%dx%d" % (w, m)] = {"n_configs": n, "N": N,
                                 "identity_HW_eq_N_minus_2K_plus_HB":
                                 "violations=%d max_abs_dev=%d" % (viol, worst),
                                 "violations": viol}
        print("C1a %dx%d: %d configs, violations=%d" % (w, m, n, viol), flush=True)
    # per-row proof-level table (all row masks, w=3,4,5).  Note cyclic pairs.
    rowtbl = {}
    for w in (3, 4, 5, 6):
        bad = 0
        for K, HB, HW in row_marks(w):
            if HW != w - 2 * K + HB:
                bad += 1
        rowtbl[str(w)] = {"n_row_masks": 1 << w, "violations": bad}
        print("C1a row w=%d: %d masks, violations=%d" % (w, 1 << w, bad), flush=True)
    out["per_row_all_masks"] = rowtbl
    return out


# --------------------------------------------------------------- part C1b
def c1_measure_equivalence():
    """Exact check of (11): mu_(p,g H_W) == mu_(p_eff,g H_B), logit shift -2g.

    Take Y = exp(g) rational, p rational.  Claim: for every configuration the
    ratio

        W(cfg)/B(cfg) = [p^K(1-p)^(N-K) Y^H_W] / [p_eff^K(1-p_eff)^(N-K) Y^H_B]

    is independent of cfg (so the two normalised measures coincide), with
    p_eff/(1-p_eff) = (p/(1-p)) Y^{-2}.
    """
    cases = []
    for (w, m) in [(3, 3), (4, 3)]:
        N = w * m
        for p in [F(1, 2), F(3, 5), F(2, 3)]:
            for Y in [F(2), F(3, 2), F(5)]:
                pe = logit_shift_peff(p, Y)
                # exact logit identity check
                lhs = pe / (1 - pe)
                rhs = (p / (1 - p)) / (Y * Y)
                logit_ok = (lhs == rhs)
                ratios = set()
                allpos = True
                for K, HB, HW in torus_configs(w, m):
                    Wn = p ** K * (1 - p) ** (N - K) * Y ** HW
                    Bn = pe ** K * (1 - pe) ** (N - K) * Y ** HB
                    if Bn == 0:
                        allpos = False
                        continue
                    ratios.add(Wn / Bn)
                cases.append({
                    "w": w, "m": m, "p": str(p), "Y_exp_g": str(Y),
                    "p_eff": str(pe),
                    "logit_identity_exact": logit_ok,
                    "distinct_ratios": len(ratios),
                    "common_ratio": str(sorted(ratios)[0]) if ratios else None,
                    "measure_equiv_exact": (len(ratios) == 1 and allpos and logit_ok),
                })
                print("C1b w=%d m=%d p=%s Y=%s: ratios=%d logit_ok=%s"
                      % (w, m, p, Y, len(ratios), logit_ok), flush=True)
    return {"cases": cases}


# --------------------------------------------------------------- part C1c
def c1_T_identity():
    with open(IN + "/root-response.json") as f:
        d = json.load(f)
    pc = float(PC_S)
    pred = 2 * pc * (1 - pc)
    rows = []
    for w in ["4", "5", "6", "7", "8"]:
        tb = d["sources"]["black"]["at_pc"][w]["T"]
        tw = d["sources"]["white"]["at_pc"][w]["T"]
        diff = tw - tb
        rows.append({"w": int(w), "T_black": tb, "T_white": tw,
                     "T_W_minus_T_B": diff, "two_p_one_minus_p": pred,
                     "abs_dev": abs(diff - pred),
                     "rel_dev": abs(diff - pred) / abs(pred)})
        print("C1c w=%s T_W-T_B=%.17g  2p(1-p)=%.17g  dev=%.3e"
              % (w, diff, pred, abs(diff - pred)), flush=True)
    # the exact decimal subtraction quoted in the task spec
    spec_meas = 0.17325528 - (-0.30954106)
    return {"rows": rows, "two_p_one_minus_p_full": pred,
            "spec_8dec_measurement": spec_meas}


# --------------------------------------------------------------- part C1d
def c1d_finite_torus_TN(w=3, m=3, p0=0.4, h=1e-5):
    """Independent finite-torus check of (12): N_W == N_B and T_W - T_B = 2p(1-p).

    Uses ONLY full enumeration + the rank function: builds the rank law (P0,P1,P2)
    as a function of (p, gamma) for both sources (black pairs / white pairs) and
    differentiates numerically.  No transfer engine, no asymptotic identification.
    """
    N = w * m
    rankcache = {}
    HB = {}
    HW = {}
    Kc = {}
    for omega in range(1 << N):
        Kc[omega] = bin(omega).count("1")
        rankcache[omega] = _rank(None, omega, w, m)
        hb = hw = 0
        for y in range(m):
            for i in range(w):
                a = (omega >> (y * w + i)) & 1
                b = (omega >> (y * w + (i + 1) % w)) & 1
                if a and b:
                    hb += 1
                elif (not a) and (not b):
                    hw += 1
        HB[omega] = hb
        HW[omega] = hw
    # index configs by rank once
    byrank = {0: [], 1: [], 2: []}
    for omega in range(1 << N):
        byrank[rankcache[omega]].append(omega)

    def law(p, g, src):
        Z = {0: 0.0, 1: 0.0, 2: 0.0}
        for r in (0, 1, 2):
            s = 0.0
            for omega in byrank[r]:
                K = Kc[omega]
                Hv = HB[omega] if src == "B" else HW[omega]
                s += math.exp(K * math.log(p) + (N - K) * math.log(1 - p) + g * Hv)
            Z[r] = s
        return Z

    def bd(p, g, src):
        Z = law(p, g, src)
        bb = 0.5 * math.log(Z[0] / Z[2])
        dd = math.log(Z[1] / math.sqrt(Z[0] * Z[2]))
        return bb, dd

    out = {"w": w, "m": m, "p0": p0, "h": h, "sources": {}}
    for src in ("B", "W"):
        bgp = (bd(p0 + h, 0.0, src)[0] - bd(p0 - h, 0.0, src)[0]) / (2 * h)
        dgp = (bd(p0 + h, 0.0, src)[1] - bd(p0 - h, 0.0, src)[1]) / (2 * h)
        bgg = (bd(p0, h, src)[0] - bd(p0, -h, src)[0]) / (2 * h)
        dgg = (bd(p0, h, src)[1] - bd(p0, -h, src)[1]) / (2 * h)
        T = -bgg / bgp
        Nn = dgg - (dgp / bgp) * bgg
        out["sources"][src] = {"b_p": bgp, "d_p": dgp, "b_g": bgg, "d_g": dgg,
                               "T": T, "N": Nn}
        print("C1d src=%s T=%.12f N=%.6e" % (src, T, Nn), flush=True)
    TW = out["sources"]["W"]["T"]
    TB = out["sources"]["B"]["T"]
    out["T_W_minus_T_B"] = TW - TB
    out["two_p_one_minus_p"] = 2 * p0 * (1 - p0)
    out["T_identity_residual"] = (TW - TB) - 2 * p0 * (1 - p0)
    out["N_W_minus_N_B"] = out["sources"]["W"]["N"] - out["sources"]["B"]["N"]
    print("C1d: T_W-T_B=%.12f (2p(1-p)=%.12f, resid=%.2e)  N_W-N_B=%.2e"
          % (out["T_W_minus_T_B"], out["two_p_one_minus_p"],
             out["T_identity_residual"], out["N_W_minus_N_B"]), flush=True)
    return out


# --------------------------------------------------------------- part C2f
def c2_f_identities():
    """Exact (Fraction) single-row partition function facts, w = 4..8.

    Z_row(p,g) = sum_m p^K (1-p)^(w-K) exp(g H(m)),  f = log Z_row.
    Verify at g=0:  Z_row = 1,  f_p = 0,
                    d_g f = E[H] = w p^2 (black pairs) / w (1-p)^2 (white pairs).
    """
    res = {}
    for w in (4, 5, 6, 7, 8):
        marks = list(row_marks(w))
        per = {}
        for p in [F(1, 2), F(3, 5), F(7, 13)]:
            Z = sum((p ** K) * ((1 - p) ** (w - K)) for (K, HB, HW) in marks)
            # f_p at g=0 : sum (K/p - (w-K)/(1-p)) * weight  == 0
            fp = sum((F(K)/p - F(w - K) / (1 - p)) * (p ** K) * ((1 - p) ** (w - K))
                     for (K, HB, HW) in marks)
            # d_g f at g=0 = E[H] = sum H * weight   (f_g, not raw normalised)
            fgB = sum(HB * (p ** K) * ((1 - p) ** (w - K)) for (K, HB, HW) in marks)
            fgW = sum(HW * (p ** K) * ((1 - p) ** (w - K)) for (K, HB, HW) in marks)
            per[str(p)] = {
                "Z_row_at_g0": str(Z),
                "Z_row_is_one": (Z == 1),
                "f_p_at_g0": str(fp),
                "f_p_is_zero": (fp == 0),
                "f_g_black": str(fgB),
                "f_g_black_eq_w_p2": (fgB == w * p * p),
                "f_g_white": str(fgW),
                "f_g_white_eq_w_1mp2": (fgW == w * (1 - p) ** 2),
                "w_p2": str(w * p * p),
                "w_1mp2": str(w * (1 - p) ** 2),
            }
            print("C2f w=%d p=%s: Z=1? %s  f_p=0? %s  f_gB=wp^2? %s  f_gW=w(1-p)^2? %s"
                  % (w, p, Z == 1, fp == 0, fgB == w * p * p,
                     fgW == w * (1 - p) ** 2), flush=True)
        res[str(w)] = per
    return res


# --------------------------------------------------------------- part C2
def c2_recompute():
    with open(IN + "/root-response.json") as f:
        d = json.load(f)
    pc = float(PC_S)
    fb = lambda w: w * pc * pc
    fw = lambda w: w * (1 - pc) ** 2
    rows = []
    diffs = []
    wNhat = []
    for w in [4, 5, 6, 7, 8]:
        s = d["sources"]
        nb = s["black"]["at_pc"][str(w)]["Nhat_per_row"]
        nw = s["white"]["at_pc"][str(w)]["Nhat_per_row"]
        cb = nb + fb(w)
        cw = nw + fw(w)
        rows.append({"w": w, "raw_black": nb, "raw_white": nw,
                     "f_g_black": fb(w), "f_g_white": fw(w),
                     "corr_black": cb, "corr_white": cw,
                     "diff_corr": cb - cw,
                     "raw_white_minus_black": nw - nb,
                     "w_times_2p_minus_1": w * (2 * pc - 1),
                     "w_times_corr_black": w * cb,
                     "w_times_corr_white": w * cw})
        diffs.append(abs(cb - cw))
        wNhat.append(w * cb)
        print("C2 w=%d rawB=%.11f rawW=%.11f corrB=%.11f corrW=%.11f "
              "diff=%.3e w*corr=%.8f"
              % (w, nb, nw, cb, cw, cb - cw, w * cb), flush=True)
    ref = {4: 0.05534176915, 5: 0.04284515731, 6: 0.03483750063,
           7: 0.02943453316, 8: 0.02550667528}
    refw = [0.22137, 0.21423, 0.20903, 0.20604, 0.20405]
    cmp_rows = []
    for i, w in enumerate([4, 5, 6, 7, 8]):
        cmp_rows.append({"w": w, "my_corr_black": rows[i]["corr_black"],
                         "analysis_table": ref[w],
                         "abs_dev_vs_analysis": abs(rows[i]["corr_black"] - ref[w]),
                         "my_w_times_corr": wNhat[i],
                         "analysis_w_times": refw[i],
                         "abs_dev_w_times": abs(wNhat[i] - refw[i])})
        print("   vs analysis: dev=%.3e (w*Nhat dev=%.3e)"
              % (abs(rows[i]["corr_black"] - ref[w]), abs(wNhat[i] - refw[i])), flush=True)
    return {"rows": rows, "max_abs_diff_corrected": max(diffs),
            "max_rel_diff_corrected": max(diffs) / max(abs(r["corr_black"]) for r in rows),
            "w_times_Nhat": wNhat, "vs_analysis": cmp_rows,
            "analysis_claim_1p5e-15": 1.5e-15}


# --------------------------------------------------------------- part C3
def c3_counterexample(w=3, m=3):
    """Exact: popcount source exp(g K) is a pure logit shift.

    It is a source whose rank-law tangent is parallel to d/dp (so N == 0) yet
    whose T equals -p(1-p), not -1.  Verified by exact evaluation of the
    finite-torus b(p) curve together with the exact identification
    b(p, g) = b_0(p_eff(p,g)).
    """
    N = w * m
    # rank of a config on the torus (ambient homology rank of the black set)
    def rank_of(omega):
        return _rank(None, omega, w, m)

    # collect P0,P2 as exact Fraction polynomials in p for every rank class
    # (weights are p^K (1-p)^(N-K); accumulate exact integer coefficient lists)
    # coefficient[K] = #configs of rank r with K black sites
    from collections import defaultdict
    coeff = {0: defaultdict(int), 1: defaultdict(int), 2: defaultdict(int)}
    for omega in range(1 << N):
        K = bin(omega).count("1")
        coeff[rank_of(omega)][K] += 1

    def P(r, p):
        return sum(c * (p ** K) * ((1 - p) ** (N - K)) for K, c in coeff[r].items())

    def coords(p):
        P0, P1, P2 = P(0, p), P(1, p), P(2, p)
        b = F(0) if (P0 == 0 or P2 == 0) else (F(1) / 2) * _log_ratio(P0, P2)
        d = None
        return b, P0, P1, P2

    # b is an irrational (log) of a rational -> evaluate numerically in high
    # precision using Fraction -> float is enough for a counterexample of O(1);
    # we also give the exact p-derivative identity at the level of the measure.
    def bn(pf):
        p = float(pf)
        P0 = sum(c * p ** K * (1 - p) ** (N - K) for K, c in coeff[0].items())
        P2 = sum(c * p ** K * (1 - p) ** (N - K) for K, c in coeff[2].items())
        return 0.5 * math.log(P0 / P2)

    def dn(pf):
        p = float(pf)
        P0 = sum(c * p ** K * (1 - p) ** (N - K) for K, c in coeff[0].items())
        P1 = sum(c * p ** K * (1 - p) ** (N - K) for K, c in coeff[1].items())
        P2 = sum(c * p ** K * (1 - p) ** (N - K) for K, c in coeff[2].items())
        return math.log(P1 / math.sqrt(P0 * P2))

    p0 = 0.4
    out = {"w": w, "m": m, "p0": p0, "rank_class_counts_by_K":
           {str(r): dict(sorted(coeff[r].items())) for r in (0, 1, 2)}}

    # (i) exact measure identification: weight of config under popcount source
    #     at coupling g equals Bernoulli weight at p_eff with
    #     logit(p_eff) = logit(p) + g  (exact for rational Y=exp(g))
    tests = []
    for pr in [F(1, 3), F(2, 5)]:
        for Y in [F(2), F(3, 2)]:
            r0 = pr / (1 - pr) * Y          # odds ratio of p_eff
            pe = r0 / (1 + r0)
            ok = True
            eK = defaultdict(int)
            for omega in range(1 << N):
                K = bin(omega).count("1")
                eK[K] += 1
            # ratio of full weights per config with same K must be constant
            ratios = set()
            for K in range(N + 1):
                Wn = pr ** K * (1 - pr) ** (N - K) * Y ** K
                Bn = pe ** K * (1 - pe) ** (N - K)
                if Bn == 0:
                    if Wn != 0:
                        ok = False
                    continue
                ratios.add(Wn / Bn)
            tests.append({"p": str(pr), "Y_exp_g": str(Y), "p_eff": str(pe),
                          "distinct_ratios": len(ratios), "exact": ok and len(ratios) == 1})
    out["popcount_is_exact_logit_shift"] = tests

    # (ii) numeric: b(g) from the source reweighting  == b_0(p_eff(g))
    def b_of_g(g, p):
        # explicit sum over configs (small torus)
        Z0 = Z1 = Z2 = 0.0
        for omega in range(1 << N):
            K = bin(omega).count("1")
            w0 = math.exp(K * math.log(p) + (N - K) * math.log(1 - p))
            r = rank_of(omega)
            t = w0 * math.exp(g * K)
            if r == 0:
                Z0 += t
            elif r == 1:
                Z1 += t
            else:
                Z2 += t
        bb = 0.5 * math.log(Z0 / Z2)
        dd = math.log(Z1 / math.sqrt(Z0 * Z2))
        return bb, dd

    h = 1e-5
    gs = [-2 * h, -h, h, 2 * h]
    bg = [b_of_g(g, p0) for g in gs]
    db_dg = (-bg[3][0] + 8 * bg[2][0] - 8 * bg[1][0] + bg[0][0]) / (12 * h)
    dd_dg = (-bg[3][1] + 8 * bg[2][1] - 8 * bg[1][1] + bg[0][1]) / (12 * h)
    # p-derivatives of b,d at p0
    ph = 1e-5
    def bd(p):
        return bn(p), dn(p)
    bp = (bd(p0 + ph)[0] - bd(p0 - ph)[0]) / (2 * ph)
    dp = (bd(p0 + ph)[1] - bd(p0 - ph)[1]) / (2 * ph)
    T = -db_dg / bp
    N = dd_dg - (dp / bp) * db_dg
    out["numeric"] = {"p0": p0, "b_g": db_dg, "d_g": dd_dg, "b_p": bp, "d_p": dp,
                      "T": T, "N": N,
                      "p(1-p)": p0 * (1 - p0),
                      "T_equals_minus_p_1mp": abs(T + p0 * (1 - p0)),
                      "N_is_zero_abs": abs(N),
                      "minus_one_is_the_truth_abs": abs(T + 1.0)}
    print("C3 counterexample (popcount source, %dx%d, p0=%.3f): T=%.10f N=%.3e "
          "(p(1-p)=%.6f)" % (w, m, p0, T, N, p0 * (1 - p0)), flush=True)
    return out


def _log_ratio(a, b):
    # exact log not needed; only used for display of exact Fractions
    return math.log(float(a) / float(b))


def _rank(_, omega, w, m):
    """Ambient homology rank in {0,1,2} of the black set (4-neighbour torus)."""
    N = w * m
    occupied = set(i for i in range(N) if (omega >> i) & 1)
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
                if v not in occupied:
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


def main():
    import os
    os.makedirs(OUT, exist_ok=True)
    res = {}
    res["C1a_identity"] = c1_identity()
    res["C1b_measure_equivalence"] = c1_measure_equivalence()
    res["C1c_T_identity"] = c1_T_identity()
    res["C1d_finite_torus"] = c1d_finite_torus_TN()
    res["C2f_row_partition_exact"] = c2_f_identities()
    res["C2_recompute"] = c2_recompute()
    res["C3_counterexample"] = c3_counterexample()
    with open(OUT + "/v1_exact.json", "w") as f:
        json.dump(res, f, indent=1, default=str)
    print("-> %s/v1_exact.json" % OUT)


if __name__ == "__main__":
    main()

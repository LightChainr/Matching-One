#!/usr/bin/env python3
"""
probe #625 D2 — physical Z_L(u; p) tables, exact.

F_L(p) = (1 + M_L(p))/2 (exact Bernstein polynomial from the joint);
Q_L(u) = F^-1(u) by exact Fraction bisection (<=1e-14 after rounding);
Z_L(u; a, b) = (Q(u) - Q(a)) / (Q(b) - Q(a)), anchors (0.2, 0.8) and (0.1, 0.9).

New tables reported (not a recap of p_L^H):
  * Z_L(u) at p = p_L^H and at p = 1/2 — read directly off the same Q_L,
    since Q is parameter-free once F is fixed; "at p" selects which rows of
    the (p, u) sheet are compared: the physical sheet's Z at the two
    p-choices 1/2 and p_L^H is compared as the Z of the two p-slices of the
    joint (n, r) histogram reweighted by Bernoulli(p).  Concretely we emit
    Z^{(p)} computed from F_p = (1+M_p)/2 with M_p the SAME polynomial —
    the Z rows below are the quantile-shape rows at the fixed thresholds
    p=1/2 and p=p_L^H evaluated as quantile levels u_p with F_L(p)=... —
    the operative comparison is: Z_L(u) from F_L, and the same from the
    p-restricted sheets, both reported.

Also: Q(u) + Q(1-u) - 1 at p=1/2 and p=p_L^H (rows added to #619's table).
Also: ||Z_3 - Z_4||_inf on {0.1..0.9} at each p-choice (first W1 datum).
"""
import sys
import json
import time
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "probe"))
sys.path.insert(0, str(ROOT / "scripts" / "homological_balance"))

from mhalf_common import (  # noqa: E402
    stream_joint, joint_to_bernstein, eval_bernstein, invert_F_bisect,
)

UGRID = [Fraction(k, 10) for k in range(1, 10)]
ANCHORS = (Fraction(2, 10), Fraction(8, 10))


def physical_quantiles(counts, N):
    coeff = joint_to_bernstein(counts, N)
    F = lambda p: (1 + eval_bernstein(coeff, Fraction(p))) / 2
    Q = {u: invert_F_bisect(F, u) for u in UGRID}
    pLH = Q[Fraction(1, 2)]
    return Q, pLH, coeff


def main():
    t0 = time.time()
    out = {"schema": "matching-one.probe-mhalf-vs-shape.physicalZ.v1",
           "issue": 625, "bisection": "exact Fractions, 110 iterations (~2^-110)",
           "anchors": [0.2, 0.8], "u_grid": [0.1, 0.2, 0.3, 0.4, 0.5,
                                             0.6, 0.7, 0.8, 0.9]}

    info = {}
    for L in (3, 4):
        counts, N = stream_joint(L)
        Q, pLH, coeff = physical_quantiles(counts, N)
        Mhalf = 2 * (1 + Fraction(0)) / 2  # placeholder; real value below
        # F(1/2) and M(1/2) exact
        Fh = (1 + eval_bernstein(coeff, Fraction(1, 2))) / 2
        Mh = 2 * Fh - 1
        Z = {u: (Q[u] - Q[ANCHORS[0]]) / (Q[ANCHORS[1]] - Q[ANCHORS[0]]) for u in UGRID}
        sym = {u: Q[u] + Q[1 - u] - 1 for u in UGRID if u <= Fraction(4, 10)}
        info[L] = {
            "N": N,
            "pLH": pLH,
            "F_half": Fh,
            "M_half": Mh,
            "Q": Q,
            "Z": Z,
            "sym": sym,
            "counts": counts,
            "coeff": coeff,
        }
        print(f"L={L}: p_L^H={float(pLH):.12f}  M(1/2)={Mh}  "
              f"F(1/2)={Fh} ({float(Fh):.9f})")
        print(f"  Z(u; 0.2/0.8): " +
              " ".join(f"{float(u)}:{float(Z[u]):.9f}" for u in UGRID))
        print(f"  Q(u)+Q(1-u)-1 (u<=0.4): " +
              " ".join(f"{float(u)}:{float(sym[u]):.9f}" for u in sym))

    # p-choice reading: "Z at p" = Z at the u with Q(u) = p.
    #   p = p_L^H  -> u = 1/2 (definition of p_L^H)
    #   p = 1/2    -> u = F_L(1/2) = (1+M(1/2))/2  (differs from 1/2 since M(1/2)!=0)
    pchoice = {}
    for L in (3, 4):
        F = lambda p, L=L: (1 + eval_bernstein(info[L]["coeff"], Fraction(p))) / 2
        u_p12 = info[L]["F_half"]
        Q_ph = info[L]["Q"][Fraction(1, 2)]                # Q(1/2) = p_L^H
        Q_p12 = invert_F_bisect(F, u_p12)                  # Q(F(1/2))
        a, b = ANCHORS
        Qa, Qb = info[L]["Q"][a], info[L]["Q"][b]
        Z_ph = (Q_ph - Qa) / (Qb - Qa)
        Z_p12 = (Q_p12 - Qa) / (Qb - Qa)
        # Q(u)+Q(1-u)-1 at the two p-choices:
        #   at p=p_L^H: u=1/2 -> 2*p_L^H - 1
        #   at p=1/2:   u=F(1/2), partner 1-F(1/2) -> Q(1-F(1/2)) via bisection
        Q_partner = invert_F_bisect(F, 1 - u_p12)
        sym_p12 = Q_p12 + Q_partner - 1
        sym_ph = 2 * Q_ph - 1
        pchoice[L] = {
            "Z_at_u_half (=Z at p_L^H)": Z_ph,
            "u_at_p_1_2": u_p12,
            "Z_at_u_Fhalf (=Z at p=1/2)": Z_p12,
            "sym_at_pLH": sym_ph,
            "sym_at_p12": sym_p12,
        }
        print(f"  p-choice L={L}: Z(u=1/2)={float(Z_ph):.9f}  "
              f"u(p=1/2)=F(1/2)={float(u_p12):.9f}  Z(that u)={float(Z_p12):.9f}")
        print(f"  selfsym: 2p_L^H-1={float(sym_ph):.9f}  "
              f"Q(F(1/2))+Q(1-F(1/2))-1={float(sym_p12):.9f}")

    dz_pLH = abs(pchoice[3]["Z_at_u_half (=Z at p_L^H)"]
                 - pchoice[4]["Z_at_u_half (=Z at p_L^H)"])
    dz_p12 = abs(pchoice[3]["Z_at_u_Fhalf (=Z at p=1/2)"]
                 - pchoice[4]["Z_at_u_Fhalf (=Z at p=1/2)"])
    print(f"  |Z3-Z4| at p=p_L^H (u=1/2): {float(dz_pLH):.9f}")
    print(f"  |Z3-Z4| at p=1/2 (u=F(1/2)): {float(dz_p12):.9f}")

    # ||Z_3 - Z_4||_inf on the u-grid (same anchors) — the W1 datum
    dz = max(abs(info[3]["Z"][u] - info[4]["Z"][u]) for u in UGRID)
    # self-symmetry row deltas (max over grid, both L)
    sym_max = {L: max(abs(v) for v in info[L]["sym"].values()) for L in (3, 4)}
    print(f"||Z_3 - Z_4||_inf = {float(dz):.9f}  ({dz})")
    out["L3"] = {
        "pLH": float(info[3]["pLH"]), "pLH_frac": str(info[3]["pLH"]),
        "M_half": str(info[3]["M_half"]), "F_half": str(info[3]["F_half"]),
        "Q": {str(u): float(info[3]["Q"][u]) for u in UGRID},
        "Z": {str(u): float(info[3]["Z"][u]) for u in UGRID},
        "Q_plus_Q1minus_1": {str(u): float(info[3]["sym"][u]) for u in info[3]["sym"]},
    }
    out["L4"] = {
        "pLH": float(info[4]["pLH"]), "pLH_frac": str(info[4]["pLH"]),
        "M_half": str(info[4]["M_half"]), "F_half": str(info[4]["F_half"]),
        "Q": {str(u): float(info[4]["Q"][u]) for u in UGRID},
        "Z": {str(u): float(info[4]["Z"][u]) for u in UGRID},
        "Q_plus_Q1minus_1": {str(u): float(info[4]["sym"][u]) for u in info[4]["sym"]},
    }
    out["Z3_minus_Z4_sup"] = float(dz)
    out["Z3_minus_Z4_sup_frac"] = str(dz)
    out["p_choice_table"] = {
        str(L): {
            "Z_at_u_half": float(pchoice[L]["Z_at_u_half (=Z at p_L^H)"]),
            "u_at_p_1_2": float(pchoice[L]["u_at_p_1_2"]),
            "Z_at_u_Fhalf": float(pchoice[L]["Z_at_u_Fhalf (=Z at p=1/2)"]),
            "sym_at_pLH": float(pchoice[L]["sym_at_pLH"]),
            "sym_at_p12": float(pchoice[L]["sym_at_p12"]),
        } for L in (3, 4)
    }
    out["Z3_minus_Z4_at_pchoices"] = {
        "at_pLH_u_half": float(dz_pLH), "at_p12_u_Fhalf": float(dz_p12)}
    out["selfsym_max"] = {str(L): float(sym_max[L]) for L in (3, 4)}
    out["p_choice_note"] = ("Z is computed from the physical F_L alone; the "
                            "p-choices enter as F(1/2) vs F(pLH) row levels "
                            "(M(1/2) != 0 means the two rows differ).")

    elapsed = time.time() - t0
    out["seconds"] = round(elapsed, 1)
    base = ROOT / "results" / "probe-Mhalf-vs-shape"
    base.mkdir(parents=True, exist_ok=True)
    (base / "physicalZ.json").write_text(json.dumps(out, indent=2, sort_keys=True))
    print("wrote", base / "physicalZ.json")


if __name__ == "__main__":
    main()

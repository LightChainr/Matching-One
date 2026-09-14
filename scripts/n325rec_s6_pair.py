#!/usr/bin/env python3
"""n325rec step 6 -- the small-scale rehearsal that IS affordable.

N=325's three-orientation rehearsal is blocked (step 4/5: every orientation
exceeds the 200k-state gate, and ell=sqrt(325) is already the SMALLEST
circumference with >=3 same-ell D4 classes).  The smallest same-ell
multi-orientation set is

    M = ell^2 = 25  (ell = 5):   u=(0,1) with n=5   and   u=(3,4) with n=1

-- two D4-distinct orientations at ONE fixed ell, both tiny automata
(G4/G8 = 45/147 states each).

This script runs BOTH implementations on both orientations, forms the exact
design matrix A=[1,cos4], fits the same-ell harmonic coefficients
Q = A^{-1} Omega for each implementation, propagates the measured oblique
dispersion into error bars, and quantifies the fixed-ell confound.
Self-contained (does not depend on the certificate file).
"""
from __future__ import annotations

import json
import math
import sys
import time
from fractions import Fraction

import numpy as np

sys.path.insert(0, "/workspace/dpfloor/scripts")
sys.path.insert(0, "/workspace/n325rec/scripts")

import fl_path_oblique as B            # noqa: E402
import fl_common as C                  # noqa: E402
import oblique_indep as N              # noqa: E402
import s1_orient as S                  # noqa: E402

PC = 0.5927460507921
M3 = {"P0": +0.00279, "P4": +0.29043, "B0": -0.18453, "B4": +0.45446}
PAIR = [("n25_0_1_n5", (0, 1), 5), ("n25_3_4", (3, 4), 1)]


def main():
    outp = sys.argv[1] if len(sys.argv) > 1 else \
        "/workspace/n325rec/out/s6_pair.json"
    res = {"schema": "n325rec.pair-rehearsal.v1", "p_c": PC,
           "pair_M": 25, "pair_ell": 5.0, "orientations": {}}
    for tag, u, n in PAIR:
        t0 = time.time()
        b4 = B.ObliqueSafeTransfer(n, u, False)
        b8 = B.ObliqueSafeTransfer(n, u, True)
        mode = "dense" if (b4.n <= 2500 and b8.n <= 2500) else "arpack"

        def eq(p):
            l4, _ = b4.lambda0(p, mode=mode, arpack_tol=1e-15)
            l8, _ = b8.lambda0(1.0 - p, mode=mode, arpack_tol=1e-15)
            return math.log(l4) - math.log(l8)

        rb, _w, _nf, fb = C.solve_root(eq, 0.5, 0.8, xtol=1e-16)
        n4 = N.IndependentOblique(n, u, False)
        n8 = N.IndependentOblique(n, u, True)
        rn, fn, bracket = N.solve_p_root(n4, n8, PC, dtype=np.longdouble)
        ell = n * math.hypot(*u)
        res["orientations"][tag] = {
            "u": list(u), "n": n, "ell": ell,
            "n_safe_G4_B": b4.n, "n_safe_G4_N": n4.n,
            "n_safe_G8_B": b8.n, "n_safe_G8_N": n8.n,
            "counts_match": (b4.n == n4.n and b8.n == n8.n),
            "p_root_B": rb, "p_root_N": repr(rn),
            "p_root_diff": float(rn) - rb,
            "Delta_at_root_B": float(fb), "Delta_at_root_N": float(fn),
            "Omega_B": float(-(rb - PC) * ell ** 4),
            "Omega_N": float(-(float(rn) - PC) * ell ** 4),
            "seconds": round(time.time() - t0, 2),
        }
        print("%-12s counts %s (%d/%d)  p_root_diff=%+.3e"
              % (tag, "OK" if res["orientations"][tag]["counts_match"]
                 else "MISMATCH", b4.n, b8.n,
                 res["orientations"][tag]["p_root_diff"]), flush=True)

    ell4 = 25 ** 2
    rows = [u for _t, u, _n in PAIR]
    A = [[Fraction(1), S.cos4m(a, b, 25, 1)] for a, b in rows]
    n_ = 2
    Abar = [[A[j][i] for j in range(n_)] for i in range(n_)]
    Ainv = []
    for t in range(n_):
        e = [Fraction(1 if k == t else 0) for k in range(n_)]
        w, _r, ok = S.solve_exact([r[:] for r in Abar], e)
        if not ok:
            raise SystemExit("singular design matrix")
        Ainv.append(w)
    Ainv = [[Ainv[j][i] for j in range(n_)] for i in range(n_)]
    amps = [math.sqrt(sum(float(Ainv[m][i]) ** 2 for i in range(n_)))
            for m in range(n_)]
    disp = max(abs(res["orientations"][t]["p_root_diff"]) for t, _u, _n in PAIR)
    F_oblique = disp
    fits = {}
    for path in ("B", "N"):
        o = [res["orientations"][t]["Omega_" + path] for t, _u, _n in PAIR]
        Q = [sum(float(Ainv[m][i]) * o[i] for i in range(n_)) for m in range(n_)]
        fits[path] = {"Q": Q, "sigma_from_measured_dispersion":
                      [F_oblique * ell4 * amps[m] for m in range(n_)]}
    e2 = 1.0 / 25
    res["design"] = {
        "A_exact": [[str(x) for x in r] for r in A],
        "rank_exact": S.rank_exact(A),
        "A_inverse": [[float(x) for x in r] for r in Ainv],
        "row_l2_norm": amps,
        "ell4_exact": ell4,
    }
    res["fit"] = fits
    res["measured_oblique_dispersion"] = F_oblique
    res["confound"] = {
        "B0_over_ell2": M3["B0"] * e2,
        "B4_over_ell2": M3["B4"] * e2,
        "B0_minus_B4_over_ell2": (M3["B0"] - M3["B4"]) * e2,
        "P0_minus_P4_imported": M3["P0"] - M3["P4"],
        "ratio_confound_over_signal": abs((M3["B0"] - M3["B4"]) * e2)
        / abs(M3["P0"] - M3["P4"]),
    }
    print("A rank=%d  Omega_B=%s" % (res["design"]["rank_exact"],
          ["%.12f" % res["orientations"][t]["Omega_B"] for t, _u, _n in PAIR]))
    print("Q(B)=%s  Q(N)=%s  dQ=%s"
          % (["%.9f" % x for x in fits["B"]["Q"]],
             ["%.9f" % x for x in fits["N"]["Q"]],
             ["%.2e" % (fits["N"]["Q"][m] - fits["B"]["Q"][m])
              for m in range(n_)]))
    print("sigma from measured dispersion %.2e: %s"
          % (F_oblique, ["%.2e" % x
                         for x in fits["B"]["sigma_from_measured_dispersion"]]))
    with open(outp, "w") as fh:
        json.dump(res, fh, indent=1)
    print("->", outp)


if __name__ == "__main__":
    main()

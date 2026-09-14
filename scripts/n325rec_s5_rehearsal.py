#!/usr/bin/env python3
"""n325rec step 5 -- the rehearsal algebra that does NOT need the automaton.

Three things:
 (A) SCAN: for which circumferences ell=sqrt(M) do >=3 D4-distinct orientations
     exist at the SAME ell?  (same-ell means n^2*(a^2+b^2) equal for all points.)
     For each candidate set the frame memory / frontier size is printed, so
     "change the modulus" can be answered with a cost estimate.
 (B) N=325 design: A = [[1,cos4,cos8]] exact, its inverse, the noise
     propagation from the dpfloor floor F=4.68e-16 into the fitted harmonic
     coefficients, and N=1105 for comparison.
 (C) The confound: at a FIXED ell the 3-orientation fit determines
     Q_m = P_m + B_m/ell^2 + ..., never P_m alone.  Using n1105mix's own M3
     global fit (imported, clearly marked) we quantify how large the ell^-2
     H0/H8 split is compared with the intrinsic H0/H8 split, i.e. whether a
     fixed-ell rehearsal can test P0=P8 (or B0=B8) even at zero noise.
"""
from __future__ import annotations

import json
import math
import sys
from fractions import Fraction

sys.path.insert(0, "/workspace/n325rec/scripts")
import oblique_indep as N            # noqa: E402
import s1_orient as S                # noqa: E402  (exact helpers)

FLOOR_F = 4.68e-16                   # dpfloor's adopted implementation floor
PC = 0.5927460507921

# n1105mix's own 10-point M3 global fit (imported for the confound estimate
# ONLY -- clearly marked as an external model, not measured here):
M3 = {"P0": +0.00279, "P4": +0.29043, "P8": +0.00257,
      "dP0": 0.00109, "dP4": 0.00231, "dP8": 0.00129,
      "B0": -0.18453, "B4": +0.45446, "B8": -0.04185,
      "dB0": 0.06773, "dB4": 0.14384, "dB8": 0.08362}
C_AXIAL = 3.54                       # the ell^-4 coefficient n1105mix fitted


def scan(Nmax=900, ellmax=40.0):
    out = []
    for M in range(1, Nmax + 1):
        ell = math.sqrt(M)
        if ell > ellmax:
            break
        found = {}
        for a in range(0, int(ell) + 1):
            for b in range(0, a + 1):
                if a == 0 and b == 0:
                    continue
                if math.gcd(a, b) != 1:
                    continue
                L2 = a * a + b * b
                for n in range(1, int(ellmax / math.sqrt(L2)) + 1):
                    if n * n * L2 == M:
                        key = (a, b) if a >= b else (b, a)
                        found.setdefault(key, []).append(n)
        if len(found) >= 3:
            ents = []
            for (a, b), ns in sorted(found.items()):
                u = (a, b)
                _, e4, m4 = N.frame(u, False)
                _, e8, m8 = N.frame(u, True)
                ents.append({"u": [a, b], "n": ns, "cos4": float(
                    S.cos4m(a, b, M, 1) if a * b else Fraction(
                        (a ** 4 - 6 * a * a * b * b + b ** 4),
                        (a * a + b * b) ** 2)),
                    "memory_G4": m4, "memory_G8": m8,
                    "frontier_G4": m4 * ns[-1], "frontier_G8": m8 * ns[-1]})
            out.append({"M": M, "ell": ell, "n_classes": len(found),
                        "orientations": ents})
    return out


def design_and_noise(M, rows):
    """rows: list of (a,b).  Exact A, exact/inverse, noise propagation."""
    A = [[Fraction(1)] + [S.cos4m(a, b, M, m) for m in (1, 2)] for a, b in rows]
    n = len(rows)
    if n != 3:
        return None
    Mmat = [[A[j][i] for j in range(n)] for i in range(n)]   # A^T
    Ainv = []
    for t in range(n):
        e = [Fraction(1 if k == t else 0) for k in range(n)]
        w, _r, ok = S.solve_exact([r[:] for r in Mmat], e)
        if not ok:
            return None
        Ainv.append(w)                                       # = A^{-T}
    Ainv = [[Ainv[j][i] for j in range(n)] for i in range(n)]  # A^{-1}
    amps = [math.sqrt(sum(float(Ainv[m][i]) ** 2 for i in range(n)))
            for m in range(n)]
    ell4 = Fraction(M) ** 2
    sigma_omega = FLOOR_F * float(ell4)
    return {"M": M, "ell": math.sqrt(M), "ell4_exact": str(ell4),
            "A": [[str(x) for x in r] for r in A],
            "A_inverse": [[float(x) for x in r] for r in Ainv],
            "row_l2_norm_of_Ainv": amps,
            "sigma_Omega_from_F": sigma_omega,
            "sigma_coefficient": [sigma_omega * a for a in amps],
            "sigma_coefficient_exact_ell4": [FLOOR_F * a for a in amps]}


def main():
    outp = sys.argv[1] if len(sys.argv) > 1 else \
        "/workspace/n325rec/out/s5_rehearsal.json"
    res = {"schema": "n325rec.rehearsal-algebra.v1", "p_c": PC,
           "floor_F": FLOOR_F, "imported_model_M3": M3,
           "imported_C_axial": C_AXIAL}
    # ---- (A) which ell admit >=3 same-ell D4 classes
    sc = scan()
    res["same_ell_scan"] = sc
    print("(A) circumferences with >=3 D4-distinct orientations (ell<=40):")
    for e in sc[:8]:
        print("    M=%-5d ell=%6.3f  classes=%d  %s"
              % (e["M"], e["ell"], e["n_classes"],
                 ", ".join("%s n=%s mem=%d/%d"
                           % (o["u"], o["n"], o["memory_G4"], o["memory_G8"])
                           for o in e["orientations"])))
    print("    total candidates up to ell=40: %d" % len(sc))
    # ---- (B) N=325 and N=1105 designs
    res["design_325"] = design_and_noise(325, [(1, 18), (6, 17), (10, 15)])
    res["design_1105"] = design_and_noise(1105, [(4, 33), (9, 32), (12, 31),
                                                 (23, 24)])
    print("(B) N=325 sigma(coeff) from F=4.68e-16:",
          ["%.3e" % x for x in res["design_325"]["sigma_coefficient"]])
    # ---- (C) the fixed-ell confound
    conf = {}
    for tag, M in (("N=325", 325), ("N=1105", 1105)):
        e2 = 1.0 / M
        conf[tag] = {
            "ell": math.sqrt(M),
            "B0_over_ell2": M3["B0"] * e2,
            "B4_over_ell2": M3["B4"] * e2,
            "B8_over_ell2": M3["B8"] * e2,
            "B0_minus_B8_over_ell2": (M3["B0"] - M3["B8"]) * e2,
            "P0_minus_P8_imported": M3["P0"] - M3["P8"],
            "P0_minus_P8_imported_sigma":
                math.hypot(M3["dP0"], M3["dP8"]),
            "C_over_ell4": C_AXIAL / M ** 2,
            "confound_over_intrinsic_signal":
                abs((M3["B0"] - M3["B8"]) * e2)
                / max(abs(M3["P0"] - M3["P8"]),
                      math.hypot(M3["dP0"], M3["dP8"])),
            "sigma_Omega_from_F": FLOOR_F * M ** 2,
            "noise_vs_confound": (FLOOR_F * M ** 2)
                / abs((M3["B0"] - M3["B8"]) * e2),
        }
    res["fixed_ell_confound"] = conf
    print("(C) fixed-ell confound (imported M3):")
    for k, v in conf.items():
        print("    %-7s ell=%6.3f  (B0-B8)/ell^2=%+.3e  P0-P8=%+.3e+-%.1e  "
              "confound/signal=%.2f  noise/confound=%.2e"
              % (k, v["ell"], v["B0_minus_B8_over_ell2"],
                 v["P0_minus_P8_imported"], v["P0_minus_P8_imported_sigma"],
                 v["confound_over_intrinsic_signal"], v["noise_vs_confound"]))
    with open(outp, "w") as fh:
        json.dump(res, fh, indent=1)
    print("->", outp)


if __name__ == "__main__":
    main()

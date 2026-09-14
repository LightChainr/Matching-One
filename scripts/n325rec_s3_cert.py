#!/usr/bin/env python3
"""n325rec step 3 -- CERTIFICATE: does the independent oblique automaton build
the same object as PATH B?

Three levels, weakest to strongest:
  (i)   frame identity: same Bezout complement and same transformed edge set;
  (ii)  automaton identity: same safe-state count, same row memory, and (for
        matrices small enough) the same sorted real spectrum of the safe
        transfer R(p) at a probe density -- a permutation-invariant fingerprint;
  (iii) p_root identity: Delta(p)=0 solved by each path's own machinery in the
        same p_c convention (0.5927460507921).
"""
from __future__ import annotations

import json
import math
import sys
import time

import numpy as np

sys.path.insert(0, "/workspace/dpfloor/scripts")
sys.path.insert(0, "/workspace/n325rec/scripts")

import fl_path_oblique as B            # PATH B (reference)
import fl_common as C                  # shared p_c / root helper
import oblique_indep as N              # my independent implementation

PC = 0.5927460507921
PROBE_P = 0.6
SPECTRUM_CAP = 2600

GEOMS = [
    ("diag_n4",     (1, 1), 4),
    ("diag_n5",     (1, 1), 5),
    ("slope21_n3",  (2, 1), 3),
    ("slope21_n4",  (2, 1), 4),
    ("slope31_n3",  (3, 1), 3),
    ("slope32_n2",  (3, 2), 2),
    ("slope52_n2",  (5, 2), 2),
    ("axis_n4",     (1, 0), 4),
    ("axis_n6",     (1, 0), 6),
    ("axis_n8",     (1, 0), 8),
    ("n25_3_4",     (3, 4), 1),
    ("n25_0_1_n5",  (0, 1), 5),
]


def b_root(u, n, pc, xtol=1e-16):
    """PATH B, tight configuration, exactly as fl_path_oblique.run does it."""
    g4 = B.ObliqueSafeTransfer(n, u, False)
    g8 = B.ObliqueSafeTransfer(n, u, True)
    mode = "dense" if (g4.n <= 2500 and g8.n <= 2500) else "arpack"

    def eq(p):
        l4, _ = g4.lambda0(p, mode=mode, arpack_tol=1e-15)
        l8, _ = g8.lambda0(1.0 - p, mode=mode, arpack_tol=1e-15)
        return math.log(l4) - math.log(l8)

    root, width, _nf, fres = C.solve_root(eq, 0.5, 0.8, xtol=xtol)
    return root, fres, g4, g8


def main():
    outp = sys.argv[1] if len(sys.argv) > 1 else \
        "/workspace/n325rec/out/s3_cert.json"
    only = sys.argv[2].split(",") if len(sys.argv) > 2 else None
    res = {"schema": "n325rec.certificate.v1", "p_c": PC,
           "probe_p": PROBE_P, "records": {}}
    for tag, u, n in GEOMS:
        if only and tag not in only:
            continue
        t0 = time.time()
        rec = {"tag": tag, "direction": list(u), "n": n,
               "ell": n * math.hypot(*u)}
        # ---- (i) frame identity
        frame_ok = {}
        for sect, matching in (("G4", False), ("G8", True)):
            cb, eb, mb = B.transformed_edges(u, matching)
            cn, en, mn = N.frame(u, matching)
            frame_ok[sect] = {
                "comp_B": list(cb), "comp_N": list(cn),
                "edges_B": [list(e) for e in eb],
                "edges_N": [list(e) for e in en],
                "memory_B": mb, "memory_N": mn,
                "identical": (list(cb) == list(cn)
                              and [list(e) for e in eb] == [list(e) for e in en]
                              and mb == mn),
            }
        rec["frame"] = frame_ok
        # ---- (ii) automaton identity
        aut = {}
        for sect, matching in (("G4", False), ("G8", True)):
            tb = B.ObliqueSafeTransfer(n, u, matching)
            tn = N.IndependentOblique(n, u, matching)
            entry = {"n_safe_B": tb.n, "n_safe_N": tn.n,
                     "memory_B": tb.memory, "memory_N": tn.memory,
                     "edges_B": [list(e) for e in tb.edges],
                     "edges_N": [list(e) for e in tn.edges],
                     "nnz_B": int(len(tb.counts)), "nnz_N": int(tn.nnz)}
            if max(tb.n, tn.n) <= SPECTRUM_CAP and tb.n == tn.n:
                Rb = tb.dense(PROBE_P)
                Rn = tn.matrix_f64(PROBE_P)
                sb = np.sort(np.linalg.eigvals(Rb).real)
                sn = np.sort(np.linalg.eigvals(Rn).real)
                entry["spectrum_max_abs_dev"] = float(np.max(np.abs(sb - sn)))
                entry["spectrum_scale"] = float(np.max(np.abs(sb)))
                entry["row_sum_max_abs_dev"] = float(
                    np.max(np.abs(np.sort(Rb.sum(axis=1))
                                  - np.sort(Rn.sum(axis=1)))))
            aut[sect] = entry
        rec["automaton"] = aut
        # ---- (iii) p_root identity
        try:
            rb, fresb, b4, b8 = b_root(u, n, PC)
        except Exception as exc:                            # noqa: BLE001
            rec["b_root_error"] = "%s: %s" % (type(exc).__name__, exc)
            rb = None
        try:
            t4 = N.IndependentOblique(n, u, False)
            t8 = N.IndependentOblique(n, u, True)
            rn_ld, fresn, bracket = N.solve_p_root(t4, t8, PC,
                                                   dtype=np.longdouble)
            rn = float(rn_ld)
        except Exception as exc:                            # noqa: BLE001
            rec["n_root_error"] = "%s: %s" % (type(exc).__name__, exc)
            rn = None
        if rb is not None and rn is not None:
            rec["p_root_B_tight"] = rb
            rec["p_root_N_longdouble"] = repr(rn_ld)
            rec["p_root_N_float"] = rn
            rec["p_root_diff"] = float(rn - rb)
            rec["Delta_at_root_B"] = float(fresb)
            rec["Delta_at_root_N"] = float(fresn)
            rec["bracket_width_N"] = float(bracket)
            rec["Omega_diff"] = -(rn - rb) * (n * math.hypot(*u)) ** 4
        rec["seconds"] = round(time.time() - t0, 2)
        res["records"][tag] = rec
        print("CERT %-12s n_safe=%d/%d  p_root_diff=%.3e  spec_dev=%s  %.1fs"
              % (tag, aut["G4"]["n_safe_B"], aut["G8"]["n_safe_B"],
                 (rec.get("p_root_diff") or float("nan")),
                 aut["G4"].get("spectrum_max_abs_dev"), rec["seconds"]),
              flush=True)
        with open(outp, "w") as fh:
            json.dump(res, fh, indent=1)
    print("->", outp)


if __name__ == "__main__":
    main()

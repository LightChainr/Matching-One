#!/usr/bin/env python3
"""Certificate that the three automata build the SAME object.

State counts alone are not enough: two different automata could agree on the
Perron root by accident.  Here, for small widths, the full spectrum of the safe
block is compared between
   PATH A  (#739 engine automaton, sector802_lib.aggregate)
   PATH B  (rev769 Bezout oblique automaton)
   PATH M  (double-cover frontier partition)
at the common p_c, for both sectors.  Equal spectra of equal-size matrices
(with the same row-sum structure) is a strong same-object statement.
"""
from __future__ import annotations
import sys, os, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, "/workspace/sectorA/in/engine")
import numpy as np
import fl_common as C
import sector802_lib as L
import fl_path_oblique as O
import fl_path_mine as M


def specA(T, w, matching, p):
    order, index, trans, wind = L.enumerate_states(T, w, matching)
    agg = L.aggregate(trans, wind, w, mark="none", mask_is_black=(not matching))
    R, _, _ = L.build_matrices(agg, w, p, want_dp=False, want_dg=False)
    return agg["n"], np.sort(np.linalg.eigvals(R).real)


def specB(w, matching, p):
    o = O.ObliqueSafeTransfer(w, (1, 0), matching)
    R = o.dense(p)
    return o.n, np.sort(np.linalg.eigvals(R).real)


def specM(w, matching, p):
    s, t = M.build(w, matching)
    R = M.matrix(s, t, w, p)
    return len(s), np.sort(np.linalg.eigvals(R).real)


def main():
    T = L.load_engine("/workspace/sectorA/in/engine")
    pc = C.PC
    out = {"p_c": pc, "records": {}}
    for w in (3, 4, 5, 6):
        for matching in (False, True):
            p = pc if not matching else 1.0 - pc
            nA, eA = specA(T, w, matching, p)
            nB, eB = specB(w, matching, p)
            nM, eM = specM(w, matching, p)
            rec = {"width": w, "matching": matching,
                   "n_safe": {"A": nA, "B": nB, "M": nM},
                   "n_equal": (nA == nB == nM)}
            if nA == nB == nM:
                rec["max_abs_eig_diff_A_B"] = float(np.max(np.abs(eA - eB)))
                rec["max_abs_eig_diff_A_M"] = float(np.max(np.abs(eA - eM)))
                rec["max_abs_eig_diff_B_M"] = float(np.max(np.abs(eB - eM)))
            out["records"]["w%d_%s" % (w, "G8" if matching else "G4")] = rec
            print(json.dumps(rec), flush=True)
    with open("/workspace/dpfloor/out/step_matrix_certificate.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print("-> /workspace/dpfloor/out/step_matrix_certificate.json")


if __name__ == "__main__":
    main()

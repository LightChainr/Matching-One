#!/usr/bin/env python3
"""PATH A ("lib") -- the pipeline that produced BOTH shipped files.

sector802_lib.enumerate_states (BFS over the pinned #739 engine's `step`,
tagged=False) -> aggregate(safe only) -> build_matrices (dense R) ->
perron_fh (scipy.linalg.eig) -> charge-coexistence root by bisection.

This is the *same* library that produced
  sector802 : out/root-response-raw.json          (file A)
  sector802b: out/closure-amplitude-raw.json      (file B)
so A and B are NOT independent of each other -- that is one of the findings.

Only axis directions are available here (a plain cylinder), so this path is
run on the axis geometries and, for the record, on the diag_n4 geometry with
width = n (which is the *oblique* (1,1) torus and is NOT what this cylinder
code implements -- included only to show the geometry differs).
"""
from __future__ import annotations
import sys, os, json, math, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, "/workspace/sectorA/in/engine")

import numpy as np
import fl_common as C
import sector802_lib as L


def build(T, w, matching):
    order, index, trans, is_wind = L.enumerate_states(T, w, matching)
    agg = L.aggregate(trans, is_wind, w, mark="none", mask_is_black=(not matching))
    return agg, len(order), agg["n_safe"]


def run(engine_dir, widths_matching, pc=C.PC, sink=None):
    T = L.load_engine(engine_dir)
    out = {"path": "A_lib", "engine": engine_dir, "p_c": pc, "records": {}}
    for w, matching, tag in widths_matching:
        t0 = time.time()
        aggN, ntot4, nsafe4 = build(T, w, False)
        aggQ, ntot8, nsafe8 = build(T, w, True)

        # I0 and diagnostics at a given black density p (G8 evaluated at 1-p)
        def I0_G4(p):
            return L.I0_and_derivs(aggN, w, p, want_dg=False)

        def I0_G8(p):
            return L.I0_and_derivs(aggQ, w, 1.0 - p, want_dg=False)

        def f(p):
            return I0_G4(p)["I0"] - I0_G8(p)["I0"]

        root, width, nfev, fres = C.solve_root(f, 0.5, 0.8, xtol=1e-16)
        a = I0_G4(root)
        b = I0_G8(root)
        rec = {
            "tag": tag, "width": w, "matching_pair": bool(matching),
            "n_states_total_G4": ntot4, "n_safe_G4": nsafe4,
            "n_states_total_G8": ntot8, "n_safe_G8": nsafe8,
            "p_root": root, "bracket_width": width, "nfev": nfev,
            "Delta_at_root": fres,
            "I0_G4": a["I0"], "I0_G8": b["I0"],
            "rho_G4": a["rho"], "rho_G8": b["rho"],
            "residual_G4": a["residual"], "residual_G8": b["residual"],
            "cw_width_G4": a["cw_width"], "cw_width_G8": b["cw_width"],
            "gap_ratio_G4": a["gap_ratio"], "gap_ratio_G8": b["gap_ratio"],
            "Omega": C.omega(root, tag, pc),
            "root_minus_pc": root - pc,
            "seconds": round(time.time() - t0, 2),
        }
        out["records"][tag] = rec
        print("A %-11s w=%d  p_root=%.17g  (bracket %.1e)  Omega=%.12e  %.1fs"
              % (tag, w, root, width, rec["Omega"], rec["seconds"]), flush=True)
        if sink:
            with open(sink, "w") as fh:
                json.dump(out, fh, indent=2)
    return out


if __name__ == "__main__":
    eng = sys.argv[1] if len(sys.argv) > 1 else "/workspace/sectorA/in/engine"
    pc = float(sys.argv[2]) if len(sys.argv) > 2 else C.PC
    ws = [int(x) for x in (sys.argv[3].split(",") if len(sys.argv) > 3
                           else ["4", "5", "6", "7", "8"])]
    out_path = sys.argv[4] if len(sys.argv) > 4 else "/workspace/dpfloor/out/pathA_lib.json"
    axis = [(w, False, "axis_n%d" % w) for w in ws]
    res = run(eng, axis, pc, sink=out_path)
    res["note"] = ("sector802_lib is byte-identical to the one shipped by both "
                   "sector802 and sector802b (md5 4646de953442c8cc71ae255f6cc3eec2); "
                   "this path therefore reproduces file A and file B exactly and "
                   "is NOT independent of them.")
    with open(out_path, "w") as fh:
        json.dump(res, fh, indent=2)
    print("->", out_path)

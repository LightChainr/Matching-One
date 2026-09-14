#!/usr/bin/env python3
"""n325rec step 9b -- the largest oblique geometry (slope52_n2, G8=131677) with
the float64 branch of the independent solver, compared against PATH B's
recorded tight value (that geometry is n>2500 on BOTH sides, so B is on ARPACK).
Cost gate: this is the only geometry whose longdouble root was too slow."""
from __future__ import annotations
import json, math, sys, time
import numpy as np
sys.path.insert(0, "/workspace/n325rec/scripts")
import oblique_indep as N

PC = 0.5927460507921
REF = json.load(open("/workspace/dpfloor/out/pathB_oblique_tight2.json"))["records"]["slope52_n2"]
u, n = (5, 2), 2

def main():
    outp = sys.argv[1] if len(sys.argv) > 1 else "/workspace/n325rec/out/s9b_big.json"
    t0 = time.time()
    t4 = N.IndependentOblique(n, u, False)
    t8 = N.IndependentOblique(n, u, True)
    build = time.time() - t0
    out = {"tag": "slope52_n2", "direction": list(u), "n": n,
           "n_safe_G4_indep": t4.n, "n_safe_G8_indep": t8.n,
           "n_safe_G4_ref": REF["n_safe_G4"], "n_safe_G8_ref": REF["n_safe_G8"],
           "counts_match": (t4.n == REF["n_safe_G4"] and t8.n == REF["n_safe_G8"]),
           "memory_indep": [t4.memory, t8.memory],
           "memory_ref": [REF["row_memory_G4"], REF["row_memory_G8"]],
           "nnz": [t4.nnz, t8.nnz], "build_seconds": round(build, 1),
           "p_root_ref": REF["p_root"]}
    t1 = time.time()
    rf, resid, bracket = None, None, None
    for tol, xtol in ((1e-14, 1e-16), (1e-12, 3e-16)):
        try:
            rf, resid, bracket = N.solve_p_root(t4, t8, PC, dtype=np.float64,
                                                xtol=xtol, tol=tol)
            out["solver_tol_used"] = tol
            break
        except Exception as exc:                                # noqa: BLE001
            out["error_%g" % tol] = "%s: %s" % (type(exc).__name__, exc)
    if rf is not None:
        ell = n * math.hypot(*u)
        out.update({"p_root_indep_float64": rf, "root_minus_pc": rf - PC,
                    "Delta_at_root": float(resid),
                    "bracket_width": float(bracket),
                    "p_root_diff": float(rf) - REF["p_root"],
                    "Omega": -(rf - PC) * ell ** 4,
                    "solve_seconds": round(time.time() - t1, 1)})
    out["seconds"] = round(time.time() - t0, 1)
    with open(outp, "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps({k: v for k, v in out.items()}, indent=1), flush=True)
    print("->", outp)

main()

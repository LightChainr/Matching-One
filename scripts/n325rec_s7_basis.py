#!/usr/bin/env python3
"""n325rec step 7 -- basis-invariance of the frame (the SHARED layer S1).

The Bezout complement is NOT unique.  Sliding it by the period,
    v  ->  v + n*u ,   i.e. (c,d) -> (c + n*a, d + n*b),
leaves det(u,v)=1 and leaves the quotient lattice <n*u, v> IDENTICAL (because
n*u is already in the lattice), but it changes the edge set (ds,dt), the row
memory and therefore the whole automaton.  So this is a genuinely different
transfer matrix for the SAME physical cylinder: the two p_root values must
agree.  This tests the frame layer that two same-convention paths cannot test.

Run with the independent implementation only (B's class hard-codes the
extended-gcd complement).
"""
from __future__ import annotations

import json
import math
import sys
import time

import numpy as np

sys.path.insert(0, "/workspace/n325rec/scripts")
import oblique_indep as N            # noqa: E402

PC = 0.5927460507921
GEOMS = [
    ("diag_n4",    (1, 1), 4),
    ("diag_n5",    (1, 1), 5),
    ("slope21_n3", (2, 1), 3),
    ("slope32_n2", (3, 2), 2),
    ("axis_n8",    (1, 0), 8),
    ("n25_3_4",    (3, 4), 1),
]


def main():
    outp = sys.argv[1] if len(sys.argv) > 1 else \
        "/workspace/n325rec/out/s7_basis.json"
    res = {"schema": "n325rec.basis-invariance.v1", "p_c": PC, "records": {}}
    for tag, u, n in GEOMS:
        a, b = u
        c0, d0 = N.bezout_complement_iter(a, b)
        alt = (c0 + n * a, d0 + n * b)
        rec = {"tag": tag, "direction": list(u), "n": n,
               "comp_base": [c0, d0], "comp_shifted": list(alt),
               "det_check_base": a * d0 - b * c0,
               "det_check_shifted": a * alt[1] - b * alt[0],
               "same_sublattice": bool(n * a * alt[1] - n * b * alt[0]
                                       == n * a * d0 - n * b * c0)}
        try:
            out = []
            for name, comp in (("base", (c0, d0)), ("shifted", alt)):
                t0 = time.time()
                t4 = N.IndependentOblique(n, u, False, complement=comp)
                t8 = N.IndependentOblique(n, u, True, complement=comp)
                root, resid, width = N.solve_p_root(t4, t8, PC,
                                                    dtype=np.longdouble)
                out.append({"name": name, "comp": list(comp),
                            "edges_G4": [list(e) for e in t4.edges],
                            "edges_G8": [list(e) for e in t8.edges],
                            "memory_G4": t4.memory, "memory_G8": t8.memory,
                            "n_safe_G4": t4.n, "n_safe_G8": t8.n,
                            "p_root": repr(root), "Delta": float(resid),
                            "seconds": round(time.time() - t0, 2)})
            rec["runs"] = out
            rec["p_root_diff_shifted_minus_base"] = \
                float(np.longdouble(out[1]["p_root"])
                      - np.longdouble(out[0]["p_root"]))
            rec["ok"] = True
        except Exception as exc:                            # noqa: BLE001
            rec["error"] = "%s: %s" % (type(exc).__name__, exc)
            rec["ok"] = False
        res["records"][tag] = rec
        if rec.get("ok"):
            print("BASIS %-12s base %s n_safe=%d/%d -> shifted n_safe=%d/%d "
                  "diff=%.3e"
                  % (tag, rec["comp_base"], rec["runs"][0]["n_safe_G4"],
                     rec["runs"][0]["n_safe_G8"], rec["runs"][1]["n_safe_G4"],
                     rec["runs"][1]["n_safe_G8"],
                     rec["p_root_diff_shifted_minus_base"]), flush=True)
        else:
            print("BASIS %-12s ERROR %s" % (tag, rec["error"]), flush=True)
        with open(outp, "w") as fh:
            json.dump(res, fh, indent=1)
    print("->", outp)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""n325rec step 8 -- does the "same ell" quantity depend on the frame gauge?

The Bezout complement v is defined only up to v -> v + k*u.  For n | k the
QUOTIENT LATTICE <n*u, v> is unchanged (step 7: p_root identical to the last
bit).  For n NOT dividing k the sublattice -- hence the physical n-site torus
on the square lattice -- is DIFFERENT while det(n*u,v)=n is the same, so the
circumference ell = n*|u| is the same as well.

If those give different p_root, then "Omega(theta, ell)" is not a function of
(theta, ell) alone: it carries an extra discrete argument (which n-site torus /
which transversal), and no same-ell multi-orientation fit can be interpreted
without fixing it.  This is a systematic effect that no amount of path
duplication can detect, so it must be measured separately.

Only the two smallest oblique geometries are probed (cost).
"""
from __future__ import annotations

import json
import math
import sys

import numpy as np

sys.path.insert(0, "/workspace/n325rec/scripts")
import oblique_indep as N            # noqa: E402

PC = 0.5927460507921

CASES = [
    ("diag_n4",   (1, 1), 4, [(0, 1), (1, 2), (2, 3), (3, 4)]),
    ("slope21_n3", (2, 1), 3, [(-1, 0), (1, 1), (3, 2)]),
    ("diag_n5",   (1, 1), 5, [(0, 1), (1, 2), (2, 3)]),
]


def main():
    outp = sys.argv[1] if len(sys.argv) > 1 else \
        "/workspace/n325rec/out/s8_lattice.json"
    res = {"schema": "n325rec.frame-gauge.v1", "p_c": PC,
           "note": ("for each lattice basis the automaton is rebuilt from "
                    "scratch; det(n*u,v)=n and |n*u| are identical by "
                    "construction"), "records": {}}
    for tag, u, n, comps in CASES:
        a, b = u
        rec = {"tag": tag, "direction": list(u), "n": n,
               "ell": n * math.hypot(a, b), "runs": []}
        for comp in comps:
            c, d = comp
            if a * d - b * c != 1:
                rec["runs"].append({"comp": list(comp), "error": "det != 1"})
                continue
            # gauge class: v -> v + k*u, k = ((c-c0)*... ) -- just report k
            try:
                t4 = N.IndependentOblique(n, u, False, complement=comp)
                t8 = N.IndependentOblique(n, u, True, complement=comp)
                root, resid, width = N.solve_p_root(t4, t8, PC,
                                                    dtype=np.longdouble)
                rec["runs"].append({
                    "comp": list(comp),
                    "torus_sites": abs(n * a * d - n * b * c),
                    "n_safe_G4": t4.n, "n_safe_G8": t8.n,
                    "memory_G4": t4.memory, "memory_G8": t8.memory,
                    "p_root": repr(root),
                    "p_root_float": float(root),
                    "Omega": float(-(root - PC) * (n * math.hypot(a, b)) ** 4),
                    "Delta": float(resid)})
            except Exception as exc:                        # noqa: BLE001
                rec["runs"].append({"comp": list(comp),
                                    "error": "%s: %s" % (type(exc).__name__,
                                                         exc)})
        ok = [r for r in rec["runs"] if "p_root_float" in r]
        if ok:
            base = ok[0]["p_root_float"]
            rec["spread_p_root"] = max(abs(r["p_root_float"] - base)
                                       for r in ok)
            rec["spread_Omega"] = max(abs(r["Omega"] - ok[0]["Omega"])
                                      for r in ok)
        res["records"][tag] = rec
        print("GAUGE %-12s ell=%.4f" % (tag, rec["ell"]), flush=True)
        for r in rec["runs"]:
            if "p_root_float" in r:
                print("      v=%-8s states=%d/%d  p_root=%.17f  Omega=%+.12e"
                      % (r["comp"], r["n_safe_G4"], r["n_safe_G8"],
                         r["p_root_float"], r["Omega"]), flush=True)
            else:
                print("      v=%-8s %s" % (r["comp"], r.get("error")),
                      flush=True)
        if "spread_p_root" in rec:
            print("      -> spread(p_root)=%.3e   spread(Omega)=%.3e"
                  % (rec["spread_p_root"], rec["spread_Omega"]), flush=True)
        with open(outp, "w") as fh:
            json.dump(res, fh, indent=1)
    print("->", outp)


if __name__ == "__main__":
    main()

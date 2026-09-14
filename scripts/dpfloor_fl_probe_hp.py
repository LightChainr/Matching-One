#!/usr/bin/env python3
"""Certificate for PATH H: is the longdouble power iteration actually converged?

At the PATH-H root of axis_n8 (and axis_n4) compare, at the SAME p:
  (i)  float64 LAPACK dense Perron root of R,
  (ii) longdouble power iteration, cold start, forced 400 iterations (no early stop),
  (iii) longdouble power iteration with the production stopping rule.
Also report Delta(root) reached by PATH H and the longdouble residual.
"""
from __future__ import annotations
import sys, os, json, math
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, "/workspace/sectorA/in/engine")
import numpy as np
import fl_common as C
import sector802_lib as L
import fl_hp as H

LD = np.longdouble


def main():
    T = L.load_engine("/workspace/sectorA/in/engine")
    hp = json.load(open("/workspace/dpfloor/out/pathH_longdouble.json"))
    out = {"records": {}}
    for tag in ("axis_n4", "axis_n6", "axis_n8"):
        if tag not in hp["records"]:
            continue
        w = hp["records"][tag]["width"]
        root = LD(hp["records"][tag]["p_root"])
        rec = {"width": w, "p_root_H": str(root)}
        for name, matching in (("G4", False), ("G8", True)):
            x = root if name == "G4" else (LD(1) - root)
            agg = H.safe_counts(T, w, matching)
            R, _ = H.assemble(agg, w, x)
            rho_f = float(np.max(np.linalg.eigvals(np.asarray(R, dtype=float)).real))
            lam_c, _, it_c = H.perron_power(R, None, tol=LD(0), maxiter=400)
            lam_p, _, it_p = H.perron_power(R, None, tol=LD(1e-30), maxiter=20000)
            res_c = float(np.linalg.norm(R @ (R @ np.ones(R.shape[0], dtype=LD))
                                         / LD(R.shape[0]) - lam_c
                                         * (R @ np.ones(R.shape[0], dtype=LD)) / LD(R.shape[0])))
            rec[name] = {"rho_lapack_float64": rho_f,
                         "rho_power_ld_forced400": str(lam_c), "iters_forced": it_c,
                         "rho_power_ld_production": str(lam_p), "iters_prod": it_p,
                         "rho_ld_minus_float64": str(LD(lam_c) - LD(rho_f))}
        out["records"][tag] = rec
        print(tag, json.dumps(rec, indent=1), flush=True)
    # PATH H's own Delta at its root
    print("PATH H Delta_at_root:", {t: hp["records"][t]["Delta_at_root"]
                                    for t in out["records"]})
    with open("/workspace/dpfloor/out/step_hp_certificate.json", "w") as fh:
        json.dump(out, fh, indent=2)


if __name__ == "__main__":
    main()

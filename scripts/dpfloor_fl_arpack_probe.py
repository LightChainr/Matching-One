#!/usr/bin/env python3
"""Configuration probe: how much of p_root is set by the eigensolver tolerance?

On ONE fixed geometry (axis_n8, i.e. the l=8 case of the shipped comparison)
compute the safe Perron root at the common p_c by
  * dense LAPACK (float64),
  * ARPACK eigs at tol = 1e-6, 1e-8, 1e-10, 1e-12, 1e-14, twice each,
and translate every difference into an equivalent error in p_root and in Omega.

The shipped oblique file was produced with ARPACK tol=1e-10 and
scipy brentq xtol=3e-11, so this is exactly the configuration cost of the file
whose 2.1e-12 difference the task is about.
"""
from __future__ import annotations
import sys, os, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import fl_common as C
import fl_path_oblique as O

from scipy.sparse.linalg import eigs as sp_eigs


def main():
    w = 8
    pc = C.PC
    out = {"geometry": "axis_n8", "width": w, "p_c": pc, "sectors": {}}
    # Delta'_w at the root (from closure-amplitude-raw.json) for error conversion
    dp = 2.023597286707526
    g4 = O.ObliqueSafeTransfer(w, (1, 0), False)
    g8 = O.ObliqueSafeTransfer(w, (1, 0), True)
    print("n_safe_G4=%d  n_safe_G8=%d  (dense LAPACK used: %s)"
          % (g4.n, g8.n, g4.n < 200), flush=True)
    for nm, obj, x in (("G4", g4, pc), ("G8", g8, 1.0 - pc)):
        R = obj.sparse(x)
        dense = float(np.max(np.linalg.eigvals(obj.dense(x)).real))
        rec = {"dense_float64": dense, "arpack": {}}
        for tol in (1e-6, 1e-8, 1e-10, 1e-12, 1e-14):
            vals = []
            for _ in range(3):
                v = float(sp_eigs(R, k=1, which="LM", tol=tol,
                                  maxiter=500000, return_eigenvectors=False)[0].real)
                vals.append(v)
            rec["arpack"]["tol_%g" % tol] = {
                "values": vals,
                "repeat_spread": max(vals) - min(vals),
                "minus_dense": vals[0] - dense,
                "abs_err_rel": abs(vals[0] - dense) / dense,
            }
        out["sectors"][nm] = rec
        print(nm, json.dumps(rec, indent=1), flush=True)
    # translate the worst |lam0 - dense|/lam0 into a p_root error
    worst = 0.0
    for nm, rec in out["sectors"].items():
        for tol, r in rec["arpack"].items():
            worst = max(worst, abs(r["minus_dense"]) / rec["dense_float64"])
    out["worst_relative_lambda_error"] = worst
    out["implied_p_root_error"] = worst / dp
    out["implied_Omega_error_at_ell8"] = worst / dp * 8.0 ** 4
    out["conversion"] = {"Delta_prime_w8": dp,
                         "rule": "dp_root = dlambda/lambda / Delta' ; dOmega = dp_root * ell^4"}
    print(json.dumps({k: out[k] for k in ("worst_relative_lambda_error",
                                          "implied_p_root_error",
                                          "implied_Omega_error_at_ell8")}, indent=2))
    with open("/workspace/dpfloor/out/step_arpack_probe.json", "w") as fh:
        json.dump(out, fh, indent=2)
    print("-> /workspace/dpfloor/out/step_arpack_probe.json")


if __name__ == "__main__":
    main()

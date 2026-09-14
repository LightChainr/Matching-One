#!/usr/bin/env python3
"""PATH S ("solvers") -- SAME matrix (from sector802_lib / PATH A), different
linear algebra AND different root-finder.

This path is explicitly NOT independent of PATH A for the automaton: it asks a
different question -- how much of the p_root value is set by the eigensolver and
the root finder rather than by the matrix.  Solver variants:

  scipy_eig      scipy.linalg.eig on dense R, argmax Re
  numpy_eigvals  numpy.linalg.eigvals on dense R, max Re
  scipy_eigvals  scipy.linalg.eigvals on dense R, max Re
  arpack_1e10    scipy.sparse.linalg.eigs (ARPACK), tol=1e-10
  arpack_1e14    ARPACK, tol=1e-14
  power          power iteration on dense R to a 1e-16 eigen-residual
  lib2_sectorgap lib_sectorgap.perron (numpy.eig right + numpy.eig(R.T) left)

Root finders: plain bisection to 1e-16, and scipy brentq at its default
xtol=2e-12 / rtol=8.9e-16.
"""
from __future__ import annotations
import sys, os, json, math, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, "/workspace/sectorA/in/engine")

import numpy as np
import fl_common as C
import sector802_lib as L

try:
    import scipy.linalg as sla
    from scipy.optimize import brentq
    from scipy.sparse import csr_matrix
    from scipy.sparse.linalg import eigs as sp_eigs
    HAVE_SCIPY = True
except Exception:
    HAVE_SCIPY = False


def power_perron(R, tol=1e-16, maxiter=200000):
    n = R.shape[0]
    v = np.ones(n)
    v /= v.sum()
    rho = 0.0
    lam_prev = 0.0
    for it in range(maxiter):
        w = R @ v
        s = w.sum()
        if s <= 0:
            break
        lam = s / v.sum()
        v = w / s
        if it % 5 == 4:
            if abs(lam - lam_prev) <= tol * abs(lam):
                rho = lam
                break
        lam_prev = lam
    else:
        rho = lam_prev
    # Rayleigh quotient with the converged (right) vector
    rho = float(v @ (R @ v) / (v @ v))
    return rho, it + 1


def build_R(agg, w, p):
    R, _, _ = L.build_matrices(agg, w, p, want_dp=False, want_dg=False)
    return R


SOLVERS = ["scipy_eig", "numpy_eigvals", "scipy_eigvals", "arpack_1e10",
           "arpack_1e14", "power"]


def rho_of(R, name, sp=None):
    if name == "scipy_eig":
        ev = sla.eig(R, left=False, right=False)
        return float(np.max(ev.real))
    if name == "numpy_eigvals":
        return float(np.max(np.linalg.eigvals(R).real))
    if name == "scipy_eigvals":
        return float(np.max(sla.eigvals(R).real))
    if name == "arpack_1e10":
        val = sp_eigs(sp, k=1, which="LM", tol=1e-10, maxiter=500000,
                      return_eigenvectors=False)[0]
        return float(val.real)
    if name == "arpack_1e14":
        val = sp_eigs(sp, k=1, which="LM", tol=1e-14, maxiter=500000,
                      return_eigenvectors=False)[0]
        return float(val.real)
    if name == "power":
        return power_perron(R)[0]
    raise ValueError(name)


def run(widths, pc, solvers=SOLVERS, sink=None):
    T = L.load_engine("/workspace/sectorA/in/engine")
    out = {"path": "S_solvers", "p_c": pc, "records": {}}
    for w in widths:
        t0 = time.time()
        order4, index4, trans4, wind4 = L.enumerate_states(T, w, False)
        aggN = L.aggregate(trans4, wind4, w, mark="none", mask_is_black=True)
        order8, index8, trans8, wind8 = L.enumerate_states(T, w, True)
        aggQ = L.aggregate(trans8, wind8, w, mark="none", mask_is_black=False)
        tag = "axis_n%d" % w
        rec = {"tag": tag, "width": w, "n_safe_G4": aggN["n"], "n_safe_G8": aggQ["n"],
               "variants": {}}
        for name in solvers:
            def f(p, name=name):
                R4 = build_R(aggN, w, p)
                sp4 = csr_matrix(R4) if name.startswith("arpack") else None
                l4 = rho_of(R4, name, sp4)
                R8 = build_R(aggQ, w, 1.0 - p)
                sp8 = csr_matrix(R8) if name.startswith("arpack") else None
                l8 = rho_of(R8, name, sp8)
                return -math.log(l4) + math.log(l8)

            root, width, nfev, fres = C.solve_root(f, 0.5, 0.8, xtol=1e-16)
            vrec = {"p_root_bisect": root, "bracket_width": width,
                    "Delta_at_root": fres,
                    "Omega_bisect": C.omega(root, tag, pc)}
            if HAVE_SCIPY:
                rb = brentq(f, 0.5, 0.8)
                vrec["p_root_brentq_default"] = rb
                vrec["brentq_minus_bisect"] = rb - root
            rec["variants"][name] = vrec
            print("S %-9s %-13s p_root(bisect)=%.17g  brentq-bisect=%+.2e"
                  % (tag, name, root, vrec.get("brentq_minus_bisect", float("nan"))),
                  flush=True)
        roots = [v["p_root_bisect"] for v in rec["variants"].values()]
        rec["spread_bisect"] = max(roots) - min(roots)
        rec["seconds"] = round(time.time() - t0, 2)
        out["records"][tag] = rec
        print("  -> %s solver spread (p) = %.3e   %.1fs" % (tag, rec["spread_bisect"], rec["seconds"]),
              flush=True)
        if sink:
            with open(sink, "w") as fh:
                json.dump(out, fh, indent=2)
    out["note"] = ("same automaton and same matrix as PATH A; isolates the "
                   "eigensolver + root-finder contribution.  NOT an independent "
                   "automaton path.")
    return out


if __name__ == "__main__":
    pc = float(sys.argv[1]) if len(sys.argv) > 1 else C.PC
    ws = [int(x) for x in (sys.argv[2].split(",") if len(sys.argv) > 2 else ["4", "6", "8"])]
    outp = sys.argv[3] if len(sys.argv) > 3 else "/workspace/dpfloor/out/pathS_solvers.json"
    solvers = sys.argv[4].split(",") if len(sys.argv) > 4 else SOLVERS
    res = run(ws, pc, solvers, sink=outp)
    with open(outp, "w") as fh:
        json.dump(res, fh, indent=2)
    print("->", outp)

#!/usr/bin/env python3
"""PATH H ("high precision") -- same geometry and same definition as PATH A,
but the Perron root and the root solve are done in extended precision.

This is an ARITHMETIC-level independent path: the automaton is integer-exact
(it only produces (i, j, k) multiplicities), so the only float step in PATH A
is (i) the assembly of R and (ii) the eigensolver.  Here R is assembled in
numpy.longdouble directly from the integer multiplicities and the Perron root
is obtained by longdouble power iteration (no LAPACK), with a warm start
carried along the root search; the charge-coexistence root is then found by
bisection in longdouble.

If PATH H reproduces the float64 paths to ~1e-25 in p, the float64 pipeline is
NOT the limiter.  If it does not, the floor is arithmetic, not algorithmic.
"""
from __future__ import annotations
import sys, os, json, math, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, "/workspace/sectorA/in/engine")

import numpy as np
import fl_common as C
import sector802_lib as L

LD = np.longdouble


def safe_counts(T, w, matching):
    order, index, trans, is_wind = L.enumerate_states(T, w, matching)
    agg = L.aggregate(trans, is_wind, w, mark="none", mask_is_black=(not matching))
    return agg


def assemble(agg, w, p):
    """R in longdouble from the integer multiplicities; p is longdouble."""
    n = agg["n"]
    R = np.zeros((n, n), dtype=LD)
    dR = np.zeros((n, n), dtype=LD)
    p = LD(p)
    one = LD(1)
    wt = [p ** k * (one - p) ** (w - k) for k in range(w + 1)]
    dwt = []
    for k in range(w + 1):
        if k == 0:
            dwt.append(-w * (one - p) ** (w - 1))
        elif k == w:
            dwt.append(w * p ** (w - 1))
        else:
            dwt.append((k - w * p) * p ** (k - 1) * (one - p) ** (w - k - 1))
    for i in range(n):
        for (j, k, h), c in agg["agg"][i].items():
            R[i, j] += LD(c) * wt[k]
            dR[i, j] += LD(c) * dwt[k]
    return R, dR


def perron_power(R, v0=None, tol=LD(1e-30), maxiter=20000):
    """Perron root by power iteration in longdouble, warm-started."""
    n = R.shape[0]
    if v0 is None:
        v = np.ones(n, dtype=LD)
    else:
        v = np.asarray(v0, dtype=LD).copy()
    s = v.sum()
    if s == 0:
        v = np.ones(n, dtype=LD)
        s = v.sum()
    v /= s
    lam = LD(0)
    for it in range(maxiter):
        wv = R @ v
        s = wv.sum()
        if s <= 0:
            break
        newlam = s / v.sum()
        v = wv / s
        if it and abs(newlam - lam) <= tol * abs(newlam):
            lam = newlam
            break
        lam = newlam
    lam = LD(v @ (R @ v)) / LD(v @ v)
    return lam, v, it + 1


def left_vector(R, v):
    """left Perron vector via float64 eig of R.T (only used for the derivative)."""
    ev, VL = np.linalg.eig(np.asarray(R, dtype=float).T)
    k = int(np.argmax(np.abs(ev)))
    l = np.real(VL[:, k])
    if l.sum() < 0:
        l = -l
    return l


def run(widths, pc, seeds=None):
    print("longdouble:", np.finfo(LD).dtype, "eps=%.3e" % np.finfo(LD).eps, flush=True)
    T = L.load_engine("/workspace/sectorA/in/engine")
    seeds = seeds or {}
    out = {"path": "H_longdouble", "p_c": pc,
           "longdouble_dtype": str(np.finfo(LD).dtype),
           "longdouble_eps": float(np.finfo(LD).eps), "records": {}}
    for w in widths:
        t0 = time.time()
        aggN = safe_counts(T, w, False)
        aggQ = safe_counts(T, w, True)
        tag = "axis_n%d" % w
        cache = {}

        def rho(mat_agg, x, key):
            R, dR = assemble(mat_agg, w, x)
            key0 = (key, "v")
            lam, v, it = perron_power(R, cache.get(key0))
            cache[key0] = v
            return lam, R, dR, it

        def Dlt(p):
            l4, R4, dR4, _ = rho(aggN, p, "g4")
            l8, R8, dR8, _ = rho(aggQ, 1.0 - LD(p), "g8")
            # np.log keeps longdouble (math.log would silently go back to float64)
            return -np.log(l4) + np.log(l8)

        # secant in longdouble, started from the float64 root; ~6 evaluations.
        p0 = LD(seeds.get(tag, 0.5927))
        p1 = p0 + LD("1e-9")
        f0, f1 = Dlt(p0), Dlt(p1)
        it = 0
        for it in range(40):
            if f1 == f0:
                break
            p2 = p1 - f1 * (p1 - p0) / (f1 - f0)
            if p2 <= 0.4 or p2 >= 0.9:
                p2 = (p0 + p1) / 2
            f2 = Dlt(p2)
            if abs(p2 - p1) <= LD("1e-18") or abs(f2) <= LD("1e-26"):
                p0, p1, f0, f1 = p1, p2, f1, f2
                break
            p0, p1, f0, f1 = p1, p2, f1, f2
        root = p1
        l4, R4, dR4, it4 = rho(aggN, root, "g4")
        l8, R8, dR8, it8 = rho(aggQ, 1.0 - LD(root), "g8")
        rec = {"tag": tag, "width": w, "n_safe_G4": aggN["n"], "n_safe_G8": aggQ["n"],
               "p_root": root, "secant_last_step": float(abs(p1 - p0)),
               "secant_iters": it, "seed": float(seeds.get(tag, 0.5927)),
               "Delta_at_root": Dlt(root),
               "rho_G4": l4, "rho_G8": l8,
               "power_iters": [it4, it8],
               "Omega": C.omega(root, tag, pc),
               "seconds": round(time.time() - t0, 2)}
        out["records"][tag] = rec
        print("H %-9s w=%d  p_root=%.25g  Omega=%.18e  iters=%s  %.1fs"
              % (tag, w, root, rec["Omega"], rec["power_iters"], rec["seconds"]),
              flush=True)
    out["note"] = ("arithmetic-independent path: integer multiplicities -> "
                   "longdouble assembly -> power iteration -> longdouble root "
                   "search.  No LAPACK eigenvalue solver at all.")
    return out


if __name__ == "__main__":
    pc = float(sys.argv[1]) if len(sys.argv) > 1 else C.PC
    ws = [int(x) for x in (sys.argv[2].split(",") if len(sys.argv) > 2 else ["4", "6"])]
    outp = sys.argv[3] if len(sys.argv) > 3 else "/workspace/dpfloor/out/pathH_longdouble.json"
    seed_path = sys.argv[4] if len(sys.argv) > 4 else "/workspace/dpfloor/out/pathA_lib.json"
    seeds = {}
    if os.path.exists(seed_path):
        for tag, r in json.load(open(seed_path)).get("records", {}).items():
            seeds[tag] = r["p_root"]
        print("seeds from", seed_path, flush=True)
    res = run(ws, pc, seeds)
    with open(outp, "w") as fh:
        json.dump(res, fh, indent=2, default=str)
    print("->", outp)

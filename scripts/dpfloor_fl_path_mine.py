#!/usr/bin/env python3
"""PATH M ("mine") -- third automaton, written from scratch, different state
representation and a *different winding test*.

State: the frontier row of the cylinder represented in its **double cover**:
positions X = 0..2w-1, where X and X+w are the two lifts of physical column
X mod w.  The state is the canonical partition of the occupied lifted
positions by cluster connectivity.

Winding test: a cluster wraps the cylinder **iff its lift identifies X with
X+w** for some occupied lifted X.  A transition is rejected when the newly
formed frontier partition identifies any such pair.  (A wrapping cluster must
occupy the just-added row, so testing inside old-row + new-row is exact; and
the BFS is restricted to non-winding states, so the winding identification
never has to be carried forward.)

This is the same *definition* as PATH A but a different implementation AND a
different formulation of "winding" (partition/double cover vs. a per-component
deck-gain flag).  Axis cylinder only.
"""
from __future__ import annotations
import sys, os, json, math, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import fl_common as C


class DSU:
    __slots__ = ("p",)

    def __init__(self, n):
        self.p = list(range(n))

    def find(self, a):
        p = self.p
        while p[a] != a:
            p[a] = p[p[a]]
            a = p[a]
        return a

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[rb] = ra


def canonical(labels):
    """relabel to first-occurrence order (stable canonical form)."""
    remap = {}
    out = []
    for x in labels:
        if x < 0:
            out.append(-1)
            continue
        if x not in remap:
            remap[x] = len(remap)
        out.append(remap[x])
    return tuple(out)


def frontier_after(state2w, mask, w, dxs):
    """Return (new_state, winding_flag)."""
    W2 = 2 * w
    nn = 2 * W2                              # old nodes 0..W2-1, new W2..2W2-1
    d = DSU(nn)
    # union old nodes sharing a label
    first = {}
    for X in range(W2):
        lab = state2w[X]
        if lab < 0:
            continue
        if lab in first:
            d.union(first[lab], X)
        else:
            first[lab] = X
    newocc = [bool((mask >> (X % w)) & 1) for X in range(W2)]
    # horizontal (within new row) edges, lifted
    for X in range(W2):
        if newocc[X] and newocc[(X + 1) % W2]:
            d.union(W2 + X, W2 + ((X + 1) % W2))
    # vertical edges new -> old
    for X in range(W2):
        if not newocc[X]:
            continue
        for dx in dxs:
            Y = (X + dx) % W2
            if state2w[Y] >= 0:
                d.union(W2 + X, Y)
    # winding test on the new row
    for X in range(w):
        if newocc[X] and newocc[X + w]:
            if d.find(W2 + X) == d.find(W2 + X + w):
                return None, True
    labels = [-1] * W2
    for X in range(W2):
        if newocc[X]:
            labels[X] = d.find(W2 + X)
    return canonical(labels), False


def build(w, matching, state_cap=200000):
    dxs = (-1, 0, 1) if matching else (0,)
    start = canonical([-1] * (2 * w))
    index = {start: 0}
    states = [start]
    trans = []
    i = 0
    while i < len(states):
        s = states[i]
        row = []
        for mask in range(1 << w):
            tgt, wind = frontier_after(s, mask, w, dxs)
            if wind or tgt is None:
                row.append(-1)
                continue
            j = index.get(tgt)
            if j is None:
                if len(states) >= state_cap:
                    raise RuntimeError("state cap")
                j = len(states)
                index[tgt] = j
                states.append(tgt)
            row.append(j)
        trans.append(row)
        i += 1
    return states, trans


def matrix(states, trans, w, p):
    n = len(states)
    R = np.zeros((n, n))
    for i, row in enumerate(trans):
        for mask, j in enumerate(row):
            if j >= 0:
                k = bin(mask).count("1")
                R[i, j] += p ** k * (1.0 - p) ** (w - k)
    return R


def lambda0(states, trans, w, p):
    R = matrix(states, trans, w, p)
    return float(np.max(np.linalg.eigvals(R).real))


def run(widths, pc):
    out = {"path": "M_mine_doublecover", "p_c": pc, "records": {}}
    for w in widths:
        t0 = time.time()
        sN, tN = build(w, False)
        sQ, tQ = build(w, True)

        def f(p):
            return -math.log(lambda0(sN, tN, w, p)) + math.log(lambda0(sQ, tQ, w, 1.0 - p))

        root, width, nfev, fres = C.solve_root(f, 0.5, 0.8, xtol=1e-16)
        tag = "axis_n%d" % w
        rec = {"tag": tag, "width": w,
               "n_safe_G4": len(sN), "n_safe_G8": len(sQ),
               "p_root": root, "bracket_width": width, "nfev": nfev,
               "Delta_at_root": fres,
               "rho_G4": lambda0(sN, tN, w, root),
               "rho_G8": lambda0(sQ, tQ, w, 1.0 - root),
               "Omega": C.omega(root, tag, pc), "root_minus_pc": root - pc,
               "seconds": round(time.time() - t0, 2)}
        out["records"][tag] = rec
        print("M %-9s w=%d  n_safe=%d/%d  p_root=%.17g  Omega=%.12e  %.1fs"
              % (tag, w, len(sN), len(sQ), root, rec["Omega"], rec["seconds"]),
              flush=True)
    out["note"] = ("third automaton, double-cover frontier partition; winding = "
                   "lift identifies X and X+w.  Independent of both the #739 "
                   "engine (PATH A) and the Bezout/deck-gain code (PATH B).")
    return out


if __name__ == "__main__":
    pc = float(sys.argv[1]) if len(sys.argv) > 1 else C.PC
    ws = [int(x) for x in (sys.argv[2].split(",") if len(sys.argv) > 2 else ["3", "4", "5", "6"])]
    outp = sys.argv[3] if len(sys.argv) > 3 else "/workspace/dpfloor/out/pathM_mine.json"
    res = run(ws, pc)
    with open(outp, "w") as fh:
        json.dump(res, fh, indent=2)
    print("->", outp)

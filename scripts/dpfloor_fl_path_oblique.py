#!/usr/bin/env python3
"""PATH B ("oblique") -- independent code base, SL(2,Z) oblique basis.

Algorithm (different from PATH A / sector802_lib):
  * complete the primitive direction u=(a,b) to an SL(2,Z) basis (u,v) with
    Bezout;  physical NN/matching edges become a finite edge set (ds,dt) in the
    (s,t) coordinates;  rows advance in t with a fixed frontier memory
    memory = max dt.
  * frontier state = (labels, gains) ONLY (no explicit winding tuple): a cycle
    with non-zero deck gain is detected as a contradictory union and the whole
    row is discarded (`dsu.bad`) -- an "eager rejection" formulation, whereas
    PATH A records a per-component winding flag and rejects later.
  * Perron root: dense numpy eigvals for small n, ARPACK eigs for n >= 200.
  * root: brentq (default xtol=3e-11) or bisection (tight).

This is the algorithm of rev769 scripts/oblique_charge_transfer.py, which
produced the shipped oblique-spin4-controls.json (file B').

Two configurations are run per geometry:
  "loose"  : ARPACK tol=1e-10  + brentq xtol=3e-11   (the shipped settings)
  "tight"  : ARPACK tol=1e-15  + bisection 1e-16     (matched to PATH A)
so that configuration noise and implementation noise can be separated.
"""
from __future__ import annotations
import sys, os, json, math, time
from collections import Counter, deque
from dataclasses import dataclass

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import numpy as np
import fl_common as C

try:
    from scipy.optimize import brentq
    from scipy.sparse import coo_matrix
    from scipy.sparse.linalg import eigs as sp_eigs
    HAVE_SCIPY = True
except Exception:                                    # pragma: no cover
    HAVE_SCIPY = False


@dataclass(frozen=True)
class State:
    labels: tuple
    gains: tuple


class DSU:
    def __init__(self, size):
        self.parent = list(range(size))
        self.delta = [0] * size
        self.bad = False

    def find(self, a):
        if self.parent[a] != a:
            root, gain = self.find(self.parent[a])
            self.delta[a] += gain
            self.parent[a] = root
        return self.parent[a], self.delta[a]

    def join(self, a, b, g):
        ra, da = self.find(a)
        rb, db = self.find(b)
        if ra == rb:
            if db - da != g:
                self.bad = True
            return
        self.parent[rb] = ra
        self.delta[rb] = g + da - db


def extended_gcd(a, b):
    if b == 0:
        return abs(a), 1 if a >= 0 else -1, 0
    g, x1, y1 = extended_gcd(b, a % b)
    return g, y1, x1 - (a // b) * y1


def bezout_complement(a, b):
    g, x, y = extended_gcd(a, b)
    if g != 1:
        raise ValueError("direction must be primitive")
    d = x
    c = -y
    assert a * d - b * c == 1
    return c, d


def transformed_edges(direction, matching):
    a, b = direction
    c, d = bezout_complement(a, b)
    gens = [(1, 0), (0, 1)]
    if matching:
        gens += [(1, 1), (1, -1)]
    edges = set()
    for dx, dy in gens:
        ds = d * dx - c * dy
        dt = -b * dx + a * dy
        if dt < 0 or (dt == 0 and ds < 0):
            ds, dt = -ds, -dt
        edges.add((ds, dt))
    ordered = sorted(edges, key=lambda e: (e[1], e[0]))
    return (c, d), ordered, max(dt for _, dt in ordered)


def empty_state(width, memory):
    return State((-1,) * (width * memory), (0,) * (width * memory))


def step(state, mask, width, edges, memory):
    old_size = memory * width
    new_base = old_size
    dsu = DSU((memory + 1) * width)
    reps = {}
    for idx, lab in enumerate(state.labels):
        if lab < 0:
            continue
        if lab in reps:
            dsu.join(reps[lab], idx, state.gains[idx])
        else:
            reps[lab] = idx
    for ds, dt in edges:
        if dt == 0:
            for i in range(width):
                if not (mask >> i) & 1:
                    continue
                j = (i + ds) % width
                if (mask >> j) & 1:
                    dsu.join(new_base + i, new_base + j, (i + ds) // width)
            continue
        row_position = memory - dt
        for i in range(width):
            oi = row_position * width + i
            if state.labels[oi] < 0:
                continue
            j = (i + ds) % width
            if (mask >> j) & 1:
                dsu.join(oi, new_base + j, (i + ds) // width)
    if dsu.bad:
        return None
    retained = list(range(width, old_size)) + [new_base + i for i in range(width)]
    labels = [-1] * old_size
    gains = [0] * old_size
    canon = {}
    nxt = 0
    for out_i, idx in enumerate(retained):
        occ = state.labels[idx] >= 0 if idx < old_size else bool(mask >> (idx - new_base) & 1)
        if not occ:
            continue
        root, gain = dsu.find(idx)
        if root not in canon:
            canon[root] = (nxt, gain)
            nxt += 1
        lab, base = canon[root]
        labels[out_i] = lab
        gains[out_i] = gain - base
    return State(tuple(labels), tuple(gains))


def build(width, direction, matching, state_cap=400000):
    comp, edges, memory = transformed_edges(direction, matching)
    start = empty_state(width, memory)
    states = [start]
    index = {start: 0}
    queue = deque([start])
    transitions = []
    while queue:
        state = queue.popleft()
        row = []
        for mask in range(1 << width):
            nxt = step(state, mask, width, edges, memory)
            if nxt is None:
                row.append(-1)
                continue
            if nxt not in index:
                if len(states) >= state_cap:
                    raise RuntimeError("state cap reached")
                index[nxt] = len(states)
                states.append(nxt)
                queue.append(nxt)
            row.append(index[nxt])
        transitions.append(row)
    return states, transitions, comp, edges, memory


class ObliqueSafeTransfer:
    def __init__(self, width, direction, matching, state_cap=400000):
        self.width = width
        self.direction = direction
        self.states, transitions, self.complement, self.edges, self.memory = \
            build(width, direction, matching, state_cap)
        counter = Counter()
        for src, row in enumerate(transitions):
            for mask, dst in enumerate(row):
                if dst >= 0:
                    counter[(src, dst, mask.bit_count())] += 1
        keys = list(counter)
        self.src = np.array([k[0] for k in keys], dtype=np.int32)
        self.dst = np.array([k[1] for k in keys], dtype=np.int32)
        self.occ = np.array([k[2] for k in keys], dtype=np.int16)
        self.counts = np.array([counter[k] for k in keys], dtype=float)
        self.n = len(self.states)

    def dense(self, p):
        q = 1.0 - p
        R = np.zeros((self.n, self.n))
        data = self.counts * p ** self.occ * q ** (self.width - self.occ)
        np.add.at(R, (self.src, self.dst), data)
        return R

    def sparse(self, p):
        q = 1.0 - p
        data = self.counts * p ** self.occ * q ** (self.width - self.occ)
        return coo_matrix((data, (self.src, self.dst)),
                          shape=(self.n, self.n)).tocsr()

    def lambda0(self, p, mode="auto", arpack_tol=1e-10):
        if mode == "auto":
            mode = "dense" if self.n < 200 else "arpack"
        if mode == "dense":
            return float(np.max(np.linalg.eigvals(self.dense(p)).real)), "dense"
        sp = self.sparse(p)
        val = sp_eigs(sp, k=1, which="LM", tol=arpack_tol, maxiter=500000,
                      return_eigenvectors=False)[0]
        return float(val.real), "arpack(%g)" % arpack_tol


def run(engine_note, spec, pc, config, sink=None):
    """spec: list of (tag, direction, n).  config: 'loose' | 'tight'."""
    out = {"path": "B_oblique", "config": config, "p_c": pc, "records": {}}
    for tag, direction, n in spec:
        t0 = time.time()
        g4 = ObliqueSafeTransfer(n, direction, False)
        g8 = ObliqueSafeTransfer(n, direction, True)
        tol = 1e-10 if config == "loose" else 1e-15
        # tight config: dense (full-accuracy LAPACK) whenever the matrix is
        # small enough for a dense eig; ARPACK with tol=1e-14 above that
        # (empirically ~1e-16 in p, same as dense -- see the ARPACK probe).
        MODE_CAP = 2500
        if config == "loose":
            mode = "auto"
        else:
            mode = "dense" if g4.n <= MODE_CAP and g8.n <= MODE_CAP else "arpack"

        def equation(p):
            l4, _ = g4.lambda0(p, mode=mode, arpack_tol=tol)
            l8, _ = g8.lambda0(1.0 - p, mode=mode, arpack_tol=tol)
            return math.log(l4) - math.log(l8)

        if config == "loose":
            if HAVE_SCIPY:
                root = brentq(equation, 0.5, 0.8, xtol=3e-11)
                width = float("nan")
            else:
                root, width, _, _ = C.solve_root(equation, 0.5, 0.8, xtol=1e-16)
        else:
            root, width, _, _ = C.solve_root(equation, 0.5, 0.8, xtol=1e-16)
        l4, m4 = g4.lambda0(root, mode=mode, arpack_tol=tol)
        l8, m8 = g8.lambda0(1.0 - root, mode=mode, arpack_tol=tol)
        ell = n * math.hypot(*direction)
        rec = {"tag": tag, "direction": list(direction), "n": n,
               "ell": ell, "cos4": C.cos4(direction),
               "n_safe_G4": g4.n, "n_safe_G8": g8.n,
               "row_memory_G4": g4.memory, "row_memory_G8": g8.memory,
               "edges_G4": [list(e) for e in g4.edges],
               "edges_G8": [list(e) for e in g8.edges],
               "lambda0_mode_G4": m4, "lambda0_mode_G8": m8,
               "p_root": float(root), "bracket_width": width,
               "Delta_at_root": math.log(l4) - math.log(l8),
               "rho_G4": l4, "rho_G8": l8,
               "Omega_ell_common": -(root - pc) * ell ** 4,
               "root_minus_pc": root - pc,
               "seconds": round(time.time() - t0, 2)}
        out["records"][tag] = rec
        print("B/%-5s %-11s n=%d p_root=%.17g  Omega=%.12e  modes=%s/%s  %.1fs"
              % (config, tag, n, root, rec["Omega_ell_common"], m4, m8,
                 rec["seconds"]), flush=True)
        if sink:
            with open(sink, "w") as fh:
                json.dump(out, fh, indent=2)
    out["note"] = engine_note
    return out


if __name__ == "__main__":
    pc = float(sys.argv[1]) if len(sys.argv) > 1 else C.PC
    cfg = sys.argv[2] if len(sys.argv) > 2 else "tight"
    out_path = sys.argv[3] if len(sys.argv) > 3 else \
        "/workspace/dpfloor/out/pathB_oblique_%s.json" % cfg
    only = sys.argv[4].split(",") if len(sys.argv) > 4 else None
    spec = [(t, u, n) for (t, u, n, _e) in C.GEOMS]
    if only:
        spec = [s for s in spec if s[0] in only]
        print("restricted to", [s[0] for s in spec], flush=True)
    res = run("implementation of rev769 scripts/oblique_charge_transfer.py",
              spec, pc, cfg, sink=out_path)
    with open(out_path, "w") as fh:
        json.dump(res, fh, indent=2)
    print("->", out_path)

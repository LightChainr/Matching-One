#!/usr/bin/env python3
"""n325rec -- SECOND, independent oblique safe-charge-transfer implementation.

Object (identical to PATH B by construction; see the SHARED layer below):
  primitive direction u=(a,b), repeat count n (circumference ell = n*|u|),
  frontier = `memory` rows of `width`=n cells in the SL(2,Z) frame
  x = s*u + t*v.  G4 = black NN graph at density p, G8 = white NN+matching graph
  at density 1-p.  The SAFE transfer keeps only frontier states with no
  cylinder-wrapping cluster;  Delta(p) = log lambda0_G4(p) - log lambda0_G8(1-p)
  and p_root is its zero.

INDEPENDENCE -- where this file differs from B_oblique
  (B = /workspace/dpfloor/scripts/fl_path_oblique.py, md5 c402cdf6a0d92d317864d793b6eece44,
   a verbatim copy of rev769 scripts/oblique_charge_transfer.py,
   md5 d5c0e5a4896bc6a244be8460e6929dd3):

  LAYER 1 -- state encoding.  B stores (labels, gains) per frontier cell with
  labels=-1 meaning "empty": the occupancy is only implicit and the gains live
  inside a DSU's parent/delta arrays.  Here a state is an explicit 3-tuple
  (occupancy bitmask, label tuple, lift tuple): a real bitmask, cells ordered
  NEWEST ROW FIRST (B orders oldest row first), and an explicit integer lift per
  cell.  Canonical form: labels by first occurrence, lifts relative to each
  component's first cell.

  LAYER 2 -- winding test.  B rejects a transition the instant a DSU union
  closes a cycle with a non-zero deck gain (eager, incremental, inside the DSU).
  Here the constraint graph of the candidate transition (new-row cells + old
  components as super-nodes, weighted edges carrying the deck displacement) is
  materialised FIRST and every edge is checked while the spanning forest is
  grown; reject iff any edge closes a non-zero fundamental cycle.
  Same MATHEMATICAL predicate, different algorithm and a different failure
  moment (see SHARED layer S2).

  LAYER 3 -- eigensolver / root finder / arithmetic.  B uses dense LAPACK
  eigvals (n<2500) or ARPACK tol<=1e-15 (n>=2500) and float64 brentq/bisection.
  Here lambda0 comes from an explicit power iteration on the sparse count-list
  matrix, warm-started between successive p, evaluated in np.longdouble (IEEE
  binary128 on this aarch64 host, ~34 decimal digits), and the root is found by
  longdouble bisection with xtol=1e-25.  No LAPACK, no ARPACK, no scipy.

  SHARED LAYER -- NOT covered by this independence (stated explicitly because a
  second path cannot certify it):
    (S1) the choice of Bezout complement v=(c,d), a*d-b*c=1.  This is a
         CONVENTION, not a theorem: v -> v+k*u gives the same cylinder only when
         n|k.  Both paths take the extended-gcd complement, so a wrong
         convention would be wrong in both.
    (S2) the predicate "wrapping <=> the frontier constraint graph admits no
         consistent integer potential".
    (S3) row memory = max dt and the frontier bookkeeping itself.
"""
from __future__ import annotations

import json
import math
import sys
import time
from collections import deque

import numpy as np

# ------------------------------------------------------------------ the frame
def egcd_iter(a: int, b: int):
    """Iterative extended gcd (no recursion).  a*x + b*y = g."""
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r != 0:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    if old_r < 0:
        old_r, old_s, old_t = -old_r, -old_s, -old_t
    return old_r, old_s, old_t


def bezout_complement_iter(a: int, b: int):
    """(c,d) with a*d - b*c = 1, same normalisation as the reference code."""
    g, x, y = egcd_iter(a, b)
    if g != 1:
        raise ValueError("direction (%d,%d) is not primitive" % (a, b))
    d, c = x, -y
    if a * d - b * c != 1:
        raise AssertionError("bezout failure for (%d,%d)" % (a, b))
    return c, d


def frame(direction, matching, complement=None):
    """Return (comp, edges, memory) for the SL(2,Z) frame of `direction`."""
    a, b = direction
    c, d = complement if complement is not None else bezout_complement_iter(a, b)
    if a * d - b * c != 1:
        raise ValueError("complement must satisfy det(u,v)=1")
    gens = [(1, 0), (0, 1)]
    if matching:
        gens = gens + [(1, 1), (1, -1)]
    edges = set()
    for dx, dy in gens:
        ds = d * dx - c * dy
        dt = -b * dx + a * dy
        if dt < 0 or (dt == 0 and ds < 0):
            ds, dt = -ds, -dt
        edges.add((ds, dt))
    ordered = sorted(edges, key=lambda e: (e[1], e[0]))
    return (c, d), ordered, max(dt for _, dt in ordered)


# ------------------------------------------------------------- the automaton
class IndependentOblique:
    """Safe transfer automaton for (width=n, direction=u, matching)."""

    def __init__(self, width, direction, matching, state_cap=300000,
                 complement=None):
        self.width = width
        self.direction = tuple(direction)
        self.matching = bool(matching)
        self.comp, self.edges, self.memory = frame(direction, matching,
                                                  complement=complement)
        t0 = time.time()
        self._build(state_cap)
        self.build_seconds = time.time() - t0
        self._warm = None

    # ---- one transition -----------------------------------------------------
    def _step(self, state, mask):
        W = self.width
        M = self.memory
        occ, lab, lif = state
        parent = list(range(W))
        delta = [0] * W
        comp_of_node = {}

        def find(x):
            root, acc = x, 0
            while parent[root] != root:
                acc += delta[root]
                root = parent[root]
            cur, g = x, 0
            while parent[cur] != cur:
                nxt = parent[cur]
                dcur = delta[cur]
                parent[cur] = root
                delta[cur] = acc - g
                g += dcur
                cur = nxt
            return root, acc

        def join(u, v, g):
            """enforce lift[v] - lift[u] = g; False on contradiction."""
            ru, gu = find(u)
            rv, gv = find(v)
            if ru == rv:
                return (gv - gu) == g
            parent[rv] = ru
            delta[rv] = g + gu - gv
            return True

        def oldnode(lbl):
            node = comp_of_node.get(lbl)
            if node is None:
                node = len(parent)
                comp_of_node[lbl] = node
                parent.append(node)
                delta.append(0)
            return node

        # --- constraints from the frame edges
        for ds, dt in self.edges:
            if dt == 0:
                for i in range(W):
                    if not (mask >> i) & 1:
                        continue
                    j = (i + ds) % W
                    if not (mask >> j) & 1:
                        continue
                    if not join(i, j, ds):
                        return None
            else:
                row = dt - 1                    # 0 == newest retained row
                for i in range(W):
                    idx = row * W + i
                    if not (occ >> idx) & 1:
                        continue
                    j = (i + ds) % W
                    if not (mask >> j) & 1:
                        continue
                    node = oldnode(lab[idx])
                    if not join(node, j, lif[idx] + ds):
                        return None

        # --- new frontier occupancy: new row on top, old rows shifted down
        new_occ = mask & ((1 << W) - 1)
        for r in range(M - 1):
            base = r * W
            if (occ >> base) & ((1 << W) - 1):
                new_occ |= (((occ >> base) & ((1 << W) - 1)) << ((r + 1) * W))

        newlab = [-1] * (M * W)
        newlif = [0] * (M * W)
        root2lab = {}
        for r in range(M):
            for c in range(W):
                nidx = r * W + c
                if not (new_occ >> nidx) & 1:
                    continue
                if r == 0:
                    node, base = c, 0
                else:
                    oidx = (r - 1) * W + c
                    node, base = oldnode(lab[oidx]), lif[oidx]
                root, g = find(node)
                if root not in root2lab:
                    root2lab[root] = (len(root2lab), g + base)
                nl, b0 = root2lab[root]
                newlab[nidx] = nl
                newlif[nidx] = (g + base) - b0
        return (new_occ, tuple(newlab), tuple(newlif))

    def _build(self, state_cap):
        W = self.width
        M = self.memory
        start = (0, tuple([-1] * (M * W)), tuple([0] * (M * W)))
        states = [start]
        index = {start: 0}
        queue = deque([0])
        trans = []
        while queue:
            si = queue.popleft()
            st = states[si]
            row = []
            for mask in range(1 << W):
                nxt = self._step(st, mask)
                if nxt is None:
                    row.append(-1)
                    continue
                j = index.get(nxt)
                if j is None:
                    if len(states) >= state_cap:
                        raise RuntimeError(
                            "state cap %d reached (%s n=%d matching=%s)"
                            % (state_cap, self.direction, W, self.matching))
                    j = len(states)
                    index[nxt] = j
                    states.append(nxt)
                    queue.append(j)
                row.append(j)
            trans.append(row)
        if len(trans) != len(states):
            raise AssertionError("BFS bookkeeping")
        self.states = states
        self.trans = trans
        self.n = len(states)
        agg = {}
        for src, row in enumerate(trans):
            for mask, dst in enumerate(row):
                if dst < 0:
                    continue
                k = (src, dst, bin(mask).count("1"))
                agg[k] = agg.get(k, 0) + 1
        keys = sorted(agg)
        self.src = np.array([k[0] for k in keys], dtype=np.int64)
        self.dst = np.array([k[1] for k in keys], dtype=np.int64)
        self.occ = np.array([k[2] for k in keys], dtype=np.int64)
        self.cnt = np.array([agg[k] for k in keys], dtype=float)
        self.nnz = len(keys)

    # ---- lambda0 ------------------------------------------------------------
    def matrix_f64(self, p):
        w = self.width
        pw = p ** self.occ
        qw = (1.0 - p) ** (w - self.occ)
        R = np.zeros((self.n, self.n))
        np.add.at(R, (self.src, self.dst), self.cnt * pw * qw)
        return R

    def matvec(self, p, v, dtype):
        d = dtype(p)
        q = dtype(1) - d
        w = self.width
        pw = np.empty(w + 1, dtype=dtype)
        qw = np.empty(w + 1, dtype=dtype)
        pw[0] = dtype(1)
        qw[0] = dtype(1)
        for i in range(1, w + 1):
            pw[i] = pw[i - 1] * d
            qw[i] = qw[i - 1] * q
        coeff = self.cnt.astype(dtype) * pw[self.occ] * qw[w - self.occ]
        out = np.zeros(self.n, dtype=dtype)
        np.add.at(out, self.dst, coeff * v[self.src])
        return out

    def lambda0(self, p, dtype=np.longdouble, tol=1e-25, maxit=1000,
                warm=True):
        """Perron root by power iteration (1-norm growth ratio).

        Returns (lambda, iters, rel_residual).  The ratio ||R v||_1/||v||_1
        converges to rho(R) as v -> the Perron vector.
        """
        dty = dtype
        if warm and self._warm is not None:
            v = self._warm
        else:
            v = np.ones(self.n, dtype=dty)
        v = v / v.max()
        prev = None
        resid = float("nan")
        it = 0
        for it in range(1, maxit + 1):
            w = self.matvec(p, v, dty)
            s = w.sum()
            vs = v.sum()
            if not (s > 0):
                raise RuntimeError("zero transfer matrix")
            lam = s / vs
            v = w / w.max()
            if prev is not None:
                resid = float(abs(lam - prev) / abs(lam))
                if resid <= tol:
                    self._warm = v
                    return float(lam), it, resid
            prev = lam
        self._warm = v
        return float(prev), it, resid


# ------------------------------------------------------------------ geometry
def solve_p_root(t4, t8, pc, dtype=np.longdouble, xtol=None, lo=0.50, hi=0.80,
                 maxiter=400, tol=1e-25):
    """Delta(p) = log l4(p) - log l8(1-p) ; zero by longdouble bisection."""
    if xtol is None:
        xtol = 1e-25 if dtype is np.longdouble else 1e-16
    dt = dtype
    l4 = lambda p: t4.lambda0(p, dtype=dt, tol=tol)[0]
    l8 = lambda p: t8.lambda0(p, dtype=dt, tol=tol)[0]

    def f(p):
        return math.log(l4(p)) - math.log(l8(dt(1) - dt(p)))

    a, b = dt(lo), dt(hi)
    fa, fb = f(a), f(b)
    if (fa > 0) == (fb > 0):
        raise ValueError("no sign change on [%s,%s]: f=%r,%r" % (lo, hi, fa, fb))
    for _ in range(maxiter):
        if b - a < dt(xtol):
            break
        m = (a + b) / 2
        fm = f(m)
        if fm == 0:
            a = b = m
            break
        if (fm > 0) == (fa > 0):
            a, fa = m, fm
        else:
            b, fb = m, fm
    root = (a + b) / 2
    return root, f(root), abs(b - a)


GEOMS_REF = {
    "diag_n4": ((1, 1), 4),
    "diag_n5": ((1, 1), 5),
    "slope21_n3": ((2, 1), 3),
    "slope21_n4": ((2, 1), 4),
    "slope31_n3": ((3, 1), 3),
    "slope32_n2": ((3, 2), 2),
    "slope52_n2": ((5, 2), 2),
    "axis_n4": ((1, 0), 4),
    "axis_n6": ((1, 0), 6),
    "axis_n8": ((1, 0), 8),
    "n325_1_18": ((1, 18), 1),
    "n325_6_17": ((6, 17), 1),
    "n325_2_3_n5": ((2, 3), 5),
    "n25_3_4": ((3, 4), 1),
    "n25_0_1_n5": ((0, 1), 5),
}


def run(spec, pc, state_cap=300000, dtype=np.longdouble, xtol=None,
        complement=None, sink=None, tol=1e-25):
    out = {"path": "n325rec_independent_oblique", "p_c": pc,
           "dtype": str(np.dtype(dtype)), "tol": tol, "records": {}}
    for tag, u, n in spec:
        t0 = time.time()
        try:
            t4 = IndependentOblique(n, u, False, state_cap=state_cap,
                                    complement=complement)
            t8 = IndependentOblique(n, u, True, state_cap=state_cap,
                                    complement=complement)
            build_t = time.time() - t0
            root, resid, width = solve_p_root(t4, t8, pc, dtype=dtype,
                                              xtol=xtol, tol=tol)
            ell = n * math.hypot(u[0], u[1])
            a, b = u
            cos4 = (a ** 4 - 6 * a * a * b * b + b ** 4) / (a * a + b * b) ** 2
            rec = {"tag": tag, "direction": list(u), "n": n, "ell": ell,
                   "cos4": cos4,
                   "comp": list(t4.comp),
                   "edges_G4": [list(e) for e in t4.edges],
                   "edges_G8": [list(e) for e in t8.edges],
                   "memory_G4": t4.memory, "memory_G8": t8.memory,
                   "nnz_G4": t4.nnz, "nnz_G8": t8.nnz,
                   "n_safe_G4": t4.n, "n_safe_G8": t8.n,
                   "build_seconds": round(build_t, 2),
                   "p_root_ld": repr(root),
                   "p_root_float": float(root),
                   "root_minus_pc": float(root - pc),
                   "Delta_at_root": float(resid),
                   "bracket_width": float(width),
                   "Omega": float(-(root - pc) * ell ** 4),
                   "seconds": round(time.time() - t0, 2)}
        except RuntimeError as exc:
            rec = {"tag": tag, "direction": list(u), "n": n,
                   "error": str(exc), "seconds": round(time.time() - t0, 2)}
        out["records"][tag] = rec
        if "error" in rec:
            print("N %-13s ABORT %s" % (tag, rec["error"]), flush=True)
        else:
            print("N %-13s n_safe=%d/%d  p_root=%s  Delta=%.3e  %.1fs"
                  % (tag, rec["n_safe_G4"], rec["n_safe_G8"],
                     rec["p_root_ld"], rec["Delta_at_root"], rec["seconds"]),
                  flush=True)
        if sink:
            with open(sink, "w") as fh:
                json.dump(out, fh, indent=1)
    return out


if __name__ == "__main__":
    pc = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5927460507921
    outp = sys.argv[2] if len(sys.argv) > 2 else \
        "/workspace/n325rec/out/s3_indep.json"
    only = sys.argv[3].split(",") if len(sys.argv) > 3 else None
    spec = [(t, u, n) for t, (u, n) in GEOMS_REF.items()]
    if only:
        spec = [s for s in spec if s[0] in only]
        print("restricted to", [s[0] for s in spec], flush=True)
    res = run(spec, pc, sink=outp)
    with open(outp, "w") as fh:
        json.dump(res, fh, indent=1)
    print("->", outp)

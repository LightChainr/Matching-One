#!/usr/bin/env python3
"""C10 cheap tests T1 and T2, packaged as one reproducible script.

T1 (P398 even-language dynamic resolution): fingerprints of the even
readout dictionaries D0/D2 over a time grid; two states in one C2 orbit can
never be separated by an even language, so the resolution is bounded by the
orbit count; measure how close it gets (saturates at n_orbit on w4-8).
T2 (C7 hidden coupling through a marked channel): on the two-copy C2 family
G(c) of c7_no_go, even tasks are c-independent while an odd-source/odd-readout
task sees c.

Deterministic, no sampling.  Writes two result JSONs.
"""

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from p398core import (Generator, orbit_labels, readout_functions,  # noqa: E402
                      reflection_perm, state_permutation)
from phase_c_quotient_factor import SparseOp, expv_cache  # noqa: E402
from c7_no_go import build  # noqa: E402

WS = Path(__file__).resolve().parents[2]


# ----------------------------------------------------------------- T1
def dynamic_resolution(w: int, readout_names, times) -> int:
    g = Generator(w)
    fns = readout_functions(w)
    pi = state_permutation(g, reflection_perm(w, w - 1))
    labels = orbit_labels(g, pi)
    blocks = {}
    for i, b in enumerate(labels):
        blocks.setdefault(b, []).append(i)
    rows = g.rows(g.rate_baseline())
    op = SparseOp(rows, g.size)
    vals = []
    for t in times:
        for name in readout_names:
            f = np.array([fns[name](s) for s in g.states], dtype=float)
            f = f - f.mean()
            vals.append(expv_cache(op, f, (float(t),))[float(t)])
    F = np.column_stack(vals)
    reps = np.array([F[blocks[b][0]] for b in sorted(blocks)])
    tol = 1e-6
    classes = []
    for r in reps:
        if not any(np.max(np.abs(r - c)) < tol for c in classes):
            classes.append(r.copy())
    return len(classes)


def run_t1() -> dict:
    D0 = ("blocks", "singletons", "wrap")
    D2 = ("blocks", "singletons", "wrap", "max_block", "linked_pairs",
          "boundary_span", "halves_linked", "covering_depth")
    out = {}
    for w in range(4, 9):
        g = Generator(w)
        pi = state_permutation(g, reflection_perm(w, w - 1))
        row = {"n": g.size, "orbits": len(set(orbit_labels(g, pi)))}
        for name, D in (("D0", D0), ("D2", D2)):
            for T, tag in ((0.5, "short"), (2.0, "mid"), (6.0, "long")):
                row[f"{name}_{tag}"] = dynamic_resolution(w, D, np.linspace(0.01, T, 8))
        out[str(w)] = row
    dest = WS / "results" / "probe-c10-t1-dynamic-signatures" / "latest.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out, indent=1))
    return out


# ----------------------------------------------------------------- T2
def run_t2() -> dict:
    times = [0.25, 0.5, 1.0, 2.0, 4.0]
    out = {}
    for m in (2, 3):
        src_e = np.zeros(2 * m)
        src_e[0] = src_e[1] = 0.5
        rd_e = np.array([1.0, 1.0, -1.0, -1.0] + [0.0, 0.0] * (m - 2))
        src_o = np.zeros(2 * m)
        rd_o = np.zeros(2 * m)
        for i in range(m):
            src_o[2 * i] = 0.5
            src_o[2 * i + 1] = -0.5
            rd_o[2 * i] = 1.0
            rd_o[2 * i + 1] = -1.0

        def responses(M, s, r):
            w_, V = np.linalg.eig(M)
            Vi = np.linalg.inv(V)
            return [float(r @ ((V * np.exp(w_ * t)) @ Vi @ s)) for t in times]

        rows = []
        for c in (0.0, 0.25, 0.5, 1.0):
            M = build(m, c)
            rows.append({"c": c, "even_response": responses(M, src_e, rd_e),
                         "odd_response": responses(M, src_o, rd_o)})
        out[str(m)] = rows
    dest = WS / "results" / "probe-c10-t2-marked-channel" / "latest.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    t1 = run_t1()
    print("T1 dynamic resolutions (w8 D0_mid, D2_short):",
          t1["8"]["D0_mid"], t1["8"]["D2_short"])
    t2 = run_t2()
    eq = all(np.allclose(t2["2"][0]["even_response"], r["even_response"], atol=1e-12)
             for r in t2["2"])
    print("T2 even-invariance check m=2 (allclose 1e-12):", eq)
    print("wrote T1 and T2 JSONs")

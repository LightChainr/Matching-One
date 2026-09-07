#!/usr/bin/env python3
"""C5 numerics: inner-product invariance and the isotypic-fibre realisation.

Two things #603 did not do:

Part A (P398, widths 4..8).  The round-1 Phase C verification used the
counting inner product.  #598/#600 left open whether the factorization is
specific to that convention.  Here the whole full-space vs orbit-quotient
task-spectrum comparison is repeated in the stationary L2(pi) inner product
(pi = invariant measure of G0, reflection-invariant, hence orbit-constant).
Claim verified: full-space and quotient spectra still agree (machine-ish
precision) in the pi product, while the pi-spectrum itself differs from the
counting spectrum - i.e. the factorization is convention-robust and the
numbers are convention-dependent.

Part B (C3 synthetic object, Theorem 2(c) of the round-2 note).  For the
six-state C3 example with V = 2 triv + 2 W, the task data of a W-sector task
computed on the full space equals the task data of the same task computed on
the *isotypic fibre system* (G restricted to the 4-dimensional W-isotypic
subspace), and the joint spectrum of a two-sector task is the union of the
per-sector spectra (block-diagonality in the isotypic decomposition).

Deterministic, no sampling.
"""

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from p398core import (Generator, orbit_labels, readout_functions,  # noqa: E402
                      reflection_perm, source_vectors, state_permutation)
from phase_c_quotient_factor import SparseOp, spectrum, transpose_op  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
DEST = ROOT / "results" / "probe-c5-l2pi-fibre" / "latest.json"
T = 4.0
GL = 24


# ---------------------------------------------------------------------------
# part A: L2(pi) factorization on P398
# ---------------------------------------------------------------------------


def stationary_measure(g: Generator, iters: int = 30000, tol: float = 1e-13):
    """pi by power iteration on (I + G^T / lam), normalized."""
    rows = g.rows(g.rate_baseline())
    src, dst, val = [], [], []
    for s, row in enumerate(rows):
        for c, v in row:
            src.append(s)
            dst.append(c)
            val.append(v)
    src = np.array(src, dtype=np.int64)
    dst = np.array(dst, dtype=np.int64)
    val = np.array(val, dtype=np.float64)
    abs_sums = np.zeros(g.size)
    np.add.at(abs_sums, src, np.abs(val))
    rate = max(2.0 * float(abs_sums.max()), 1e-9)

    def mvT(x):
        y = np.zeros(g.size)
        np.add.at(y, dst, val * x[src])  # G^T x
        return y

    v = np.ones(g.size) / g.size
    for _ in range(iters):
        v = v + mvT(v) / rate
        v /= v.sum()
    resid = float(np.max(np.abs(mvT(v))))
    if resid > 1e-6:
        raise RuntimeError(f"stationary solve did not converge: resid {resid:.2e}")
    return v, rate


def transform_rows(rows, weights: np.ndarray) -> list:
    """Rows in sqrt-weight coordinates: A'_{ij} = A_ij * w_i / w_j."""
    out = []
    for s, row in enumerate(rows):
        out.append([(c, val * weights[s] / weights[c]) for c, val in row])
    return out


def l2pi_compare(w: int) -> dict:
    g = Generator(w)
    n = g.size
    pi, _ = stationary_measure(g)
    sq = np.sqrt(pi)

    # readouts: pi-centred then sqrt(pi)-scaled; restrict to the R-even
    # readouts (the protected dictionary).  At odd widths halves_linked is
    # R-odd and must be excluded exactly as in #603, otherwise the quotient
    # comparison is a bug report, not a bug.
    fns = readout_functions(w)
    readout_vecs = {name: np.array([fns[name](s) for s in g.states]) for name in fns}
    centre = {name: float(np.dot(pi, v)) for name, v in readout_vecs.items()}
    pi_state = state_permutation(g, reflection_perm(w, w - 1))
    even_names = []
    for name, v in readout_vecs.items():
        u = v - centre[name]
        if np.linalg.norm(u - u[pi_state]) < 1e-9:
            even_names.append(name)
    reads = [sq * (readout_vecs[name] - centre[name]) for name in even_names]
    # sources
    src = source_vectors(g)
    names = ["delta_all_singletons", "delta_single_block", "delta_wrapped_pair", "uniform"]
    srcs = [sq * np.array(src[name]) for name in names]

    rows = g.rows(g.rate_baseline())
    op = SparseOp(transform_rows(rows, sq), n)
    spec_full = spectrum(op, srcs, reads)

    # quotient in pi-weighted orbit coordinates
    pi_state = state_permutation(g, reflection_perm(w, w - 1))
    labels = orbit_labels(g, pi_state)
    blocks = {}
    for i, b in enumerate(labels):
        blocks.setdefault(b, []).append(i)
    order = sorted(blocks)
    nq = len(order)
    orb_pi = np.array([pi[blocks[b][0]] for b in order])   # orbit-constant pi
    orb_size = np.array([len(blocks[b]) for b in order])
    W = np.sqrt(orb_size * orb_pi)

    # lumped generator rows
    base = g.rows(g.rate_baseline())
    qrows = []
    for b in order:
        rep = blocks[b][0]
        totals = {}
        for c, v in base[rep]:
            totals[labels[c]] = totals.get(labels[c], 0.0) + v
        qrows.append(sorted(totals.items()))
    raw_bound = max(sum(abs(val) for _c, val in row) for row in qrows)
    qop = SparseOp(transform_rows(qrows, W), nq, rate_override=max(2.0 * raw_bound, 1e-9))

    def to_orbit(vec):
        return np.array([vec[blocks[b][0]] for b in order])

    # quotient vectors: v_o = sqrt(|o| pi_o) * f_o  where f_o is the raw
    # (orbit-constant) value of the object; the readout is centred in the raw
    # space first, exactly as on the full side.
    q_srcs = [W * to_orbit(np.array(src[name])) for name in names]
    raw_centred = {name: readout_vecs[name] - centre[name] for name in readout_vecs}
    q_reads = [W * to_orbit(raw_centred[name]) for name in even_names]
    spec_quot = spectrum(qop, q_srcs, q_reads)

    k = min(len(spec_full["singular_values"]), len(spec_quot["singular_values"]))
    a = np.array(spec_full["singular_values"][:k])
    b = np.array(spec_quot["singular_values"][:k])
    rel = float(np.linalg.norm(a - b) / np.linalg.norm(a)) if k else 0.0
    return {
        "width": w, "n": n, "n_orbit": nq,
        "full_sv_head": [round(float(x), 6) for x in spec_full["singular_values"][:4]],
        "quot_sv_head": [round(float(x), 6) for x in spec_quot["singular_values"][:4]],
        "rel_diff_full_vs_quot_L2pi": rel,
    }


# ---------------------------------------------------------------------------
# part B: isotypic fibre realisation on the C3 object
# ---------------------------------------------------------------------------


def build_c3():
    n = 6
    idx = lambda o, i: o * 3 + i  # noqa: E731
    G = np.zeros((n, n))
    for i in range(3):
        G[idx(0, i), idx(0, (i + 1) % 3)] = 1.0
        G[idx(0, i), idx(1, i)] = 2.0
        G[idx(1, i), idx(1, (i - 1) % 3)] = 1.0
        G[idx(1, i), idx(0, i)] = 1.0
    for i in range(n):
        G[i, i] = -G[i].sum()
    import math
    ang = 2 * math.pi / 3
    W = np.zeros((n, 4))
    for i in range(3):
        W[idx(0, i), 0] = math.sqrt(2 / 3) * math.cos(i * ang)
        W[idx(0, i), 1] = math.sqrt(2 / 3) * math.sin(i * ang)
        W[idx(1, i), 2] = math.sqrt(2 / 3) * math.cos(i * ang)
        W[idx(1, i), 3] = math.sqrt(2 / 3) * math.sin(i * ang)
    return G, W, n


def c3_fibre_compare() -> dict:
    G, Wbasis, n = build_c3()
    # task vectors
    mu_triv = np.ones(n) / n
    f_triv = np.array([1.0, 1, 1, -1, -1, -1])
    f_triv = f_triv / np.linalg.norm(f_triv)
    mu_W = Wbasis[:, 0] + Wbasis[:, 2]
    mu_W = mu_W / np.linalg.norm(mu_W)
    f_W = Wbasis[:, 1] + Wbasis[:, 3]
    f_W = f_W / np.linalg.norm(f_W)

    rows = [[(c, G[s, c]) for c in range(n) if G[s, c] != 0] for s in range(n)]
    op_full = SparseOp(rows, n)

    def dense_rows(A):
        m = A.shape[0]
        return [[(c, A[s, c]) for c in range(m) if A[s, c] != 0] for s in range(m)]

    out = {}
    # (i) single-sector W task, full space vs fibre system
    spec_full_w = spectrum(op_full, [mu_W], [f_W])
    A_W = Wbasis.T @ G @ Wbasis          # 4x4 fibre operator (two copies of W)
    op_fibre = SparseOp(dense_rows(A_W), 4)
    mu_w = Wbasis.T @ mu_W
    f_w = Wbasis.T @ f_W
    spec_fibre = spectrum(op_fibre, [mu_w], [f_w])
    k = min(len(spec_full_w["singular_values"]), len(spec_fibre["singular_values"]))
    a = np.array(spec_full_w["singular_values"][:k])
    b = np.array(spec_fibre["singular_values"][:k])
    out["W_task_full_vs_fibre_rel"] = float(np.linalg.norm(a - b) / np.linalg.norm(a)) if k else 0.0
    out["W_task_full_sv"] = [round(float(x), 6) for x in spec_full_w["singular_values"][:4]]
    out["W_task_fibre_sv"] = [round(float(x), 6) for x in spec_fibre["singular_values"][:4]]
    # (ii) triv-task (full = triv sector, which is itself the orbit quotient)
    spec_full_t = spectrum(op_full, [mu_triv], [f_triv])
    out["triv_task_full_sv"] = [round(float(x), 6) for x in spec_full_t["singular_values"][:4]]
    # (iii) joint two-sector task spectrum vs union of per-sector spectra
    spec_joint = spectrum(op_full, [mu_triv, mu_W], [f_triv, f_W])
    out["joint_sv"] = [round(float(x), 6) for x in spec_joint["singular_values"]]
    union = sorted(
        list(spec_full_t["singular_values"]) + list(spec_full_w["singular_values"]),
        reverse=True)
    out["union_of_sector_sv"] = [round(float(x), 6) for x in union]
    return out


def main() -> None:
    part_a = {}
    for w in range(4, 9):
        t0 = time.time()
        part_a[str(w)] = l2pi_compare(w)
        part_a[str(w)]["seconds"] = round(time.time() - t0, 1)
        print(f"w={w}: rel diff (full vs quotient, L2(pi)) = "
              f"{part_a[str(w)]['rel_diff_full_vs_quot_L2pi']:.2e}", flush=True)
    part_b = c3_fibre_compare()
    print("C3 fibre:")
    print(" ", part_b)
    payload = {
        "schema": "matching-one.probe-c5-l2pi-fibre.v1",
        "probe_issue": 599,
        "part_a": {"inner_product": "L2(pi)", "widths": part_a},
        "part_b": part_b,
    }
    DEST.parent.mkdir(parents=True, exist_ok=True)
    DEST.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print("wrote", DEST)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""#598 Phase C: does the task balanced realization factor through the quotient?

Question (open half of #598; also #599 program 3 and high-value outcome 4):

    factor  microscopic dynamics
             -> C2 reflection isotypic split (even sector V^+ = orbit-quotient
                function space, odd sector V^-)
             -> finite-horizon task balanced reduction.

Expected theorem-level statement (tested numerically here):

    * with an R-equivariant generator, R-even sources B and R-even readouts C,
      and any inner product that makes the reflection an isometry, the finite-
      horizon task Hankel singular data computed on the full state space equals
      the one computed on the even sector / orbit quotient, and the anti-
      invariant (odd) sector contributes exactly zero;
    * the identification of ``even sector`` with ``orbit-quotient functions``
      uses the *orbit-weighted* inner product; using the plain (unweighted)
      Euclidean product on the quotient silently changes the singular values,
      which is precisely the inner-product convention that #598 flags as open.

Two tasks are compared at every width 4..8:

    EVEN task   : 4 sources (all R-even) x even readouts.
                  (halves_linked excluded at odd widths; all readouts are even
                  at even widths.)
    MARKED task : the same plus halves_linked at odd widths, where it carries a
                  nonzero odd component -> the anti-invariant sector becomes
                  visible.  A quotient realization built from the even part of
                  that readout *cannot* reproduce the extra singular values.

Conventions (declared, all deterministic):

    * Euclidean (counting-measure) inner product on the full space;
    * readouts centred by subtracting their counting mean;
    * the constant mode needs no explicit removal: centred readouts make the
      observability directions orthogonal to the constant function, and the
      generator's row sums kill the constant direction under e^{tG}1 = 1;
    * finite horizon T = 4 with Gauss-Legendre quadrature.

The response kernel is K_{a,j}(t) = v_a^T e^{tG} u_j  (sources v_a, centred
readouts u_j).  Task Hankel singular values are the singular values of Z^T C
where C holds columns sqrt(w_l) e^{tau_l G} v_a  (controllability factors) and
Z holds columns sqrt(w_l) e^{tau_l G^T} u_j (observability factors).

For the quotient, all objects are pulled back to the orbit-function space and
expressed in the sqrt(orbit-weight) coordinates so the inner product agrees
with the counting measure on the full space.
"""

from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from p398core import (  # noqa: E402
    Generator,
    orbit_labels,
    readout_functions,
    reflection_perm,
    source_vectors,
    state_permutation,
)

HORIZON = 4.0
GL_NODES = 32
TOL = 1e-11


# --------------------------------------------------------------------------
# sparse-matrix exponentials by uniformization


class SparseOp:
    """Row-oriented operator as numpy arrays.

    Rows are stored per-state as (col, value) so that ``rows[i]`` means
    ``(G f)_i = sum_j G_ij f_j``.  ``matvec`` therefore indexes the output by
    the row index ``i`` and pulls ``f`` at the column index ``j``.
    """

    def __init__(self, rows, size, rate_override: float = None) -> None:
        src, dst, val = [], [], []
        for s, row in enumerate(rows):
            for c, v in row:
                src.append(s)
                dst.append(c)
                val.append(v)
        self.src = np.array(src, dtype=np.int64)
        self.dst = np.array(dst, dtype=np.int64)
        self.val = np.array(val, dtype=np.float64)
        self.size = size
        # safe uniformization rate: lam >= 2 * spectral radius.  A sufficient
        # bound is 2 * max_i sum_j |A_ij| (Gershgorin row bound).
        abs_sums = np.zeros(size)
        np.add.at(abs_sums, self.src, np.abs(self.val))
        bound = float(abs_sums.max()) if size else 0.0
        if rate_override is not None:
            self.rate = float(rate_override)
        else:
            self.rate = max(2.0 * bound, 1e-9)

    def matvec(self, v: np.ndarray) -> np.ndarray:
        """A v: output indexed by the row index."""
        out = np.zeros(self.size)
        np.add.at(out, self.src, self.val * v[self.dst])
        return out

    def matvec_T(self, v: np.ndarray) -> np.ndarray:
        """A^T v: output indexed by the column index."""
        out = np.zeros(self.size)
        np.add.at(out, self.dst, self.val * v[self.src])
        return out


def poisson_terms(rate: float, t: float, tail: float = 1e-14, cap: int = 3000) -> int:
    mean = rate * t
    terms, mass, weight = 1, math.exp(-mean), math.exp(-mean)
    while mass < 1.0 - tail and terms < cap:
        weight *= mean / terms
        mass += weight
        terms += 1
    return terms + 1


def expv_cache(op: SparseOp, v: np.ndarray, times: Sequence[float]):
    """Return dict t -> e^{tA} v with one cached (I + A/lam)^k sequence."""
    rate = op.rate
    max_terms = max(poisson_terms(rate, t) for t in times) if times else 1
    seq = [np.array(v, dtype=np.float64)]
    cur = np.array(v, dtype=np.float64)
    for _ in range(max_terms - 1):
        cur = cur + op.matvec(cur) / rate
        seq.append(cur)
    out = {}
    for t in times:
        if t == 0.0:
            out[t] = seq[0].copy()
            continue
        mean = rate * t
        acc = np.zeros(op.size)
        weight = math.exp(-mean)
        k = 0
        while k < len(seq) and weight > 0.0:
            acc += weight * seq[k]
            k += 1
            weight *= mean / k
        out[t] = acc
    return out


# --------------------------------------------------------------------------
# quotient construction


def build_quotient(g: Generator, labels: Sequence[int]) -> Dict:
    """Orbit-quotient generator, orbit weights, and the embedding map.

    ``labels[i]`` is the orbit (block) of state i.  The lumped generator acts on
    orbit-constant functions; rows are aggregated exactly because the orbit
    partition is a strong lumping (verified in verify_gate1.py).
    """
    n = g.size
    blocks: Dict[int, List[int]] = {}
    for i, b in enumerate(labels):
        blocks.setdefault(b, []).append(i)
    order = sorted(blocks)
    nq = len(order)
    member: Dict[int, List[int]] = {b: blocks[b] for b in order}
    weights = np.array([len(member[b]) for b in order], dtype=np.float64)

    # aggregate the baseline generator G0 rows: rate from orbit o to orbit o'
    base = g.rows(g.rate_baseline())
    qrows: List[List[Tuple[int, float]]] = []
    for b in order:
        rep = member[b][0]
        totals: Dict[int, float] = {}
        for c, v in base[rep]:
            totals[labels[c]] = totals.get(labels[c], 0.0) + v
        qrows.append(sorted(totals.items()))
    # embed: map an orbit-constant vector (as function on orbits) to full space
    def lift(vo: np.ndarray) -> np.ndarray:
        vf = np.zeros(n)
        for idx, b in enumerate(order):
            vf[member[b]] = vo[idx]
        return vf

    return {"nq": nq, "order": order, "member": member, "weights": weights,
            "qrows": qrows, "lift": lift}


# --------------------------------------------------------------------------
# task singular spectra


def centred(readout_vecs: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
    out = {}
    for name, v in readout_vecs.items():
        u = v - float(v.mean())
        if np.linalg.norm(u) > 1e-14:
            u = u / np.linalg.norm(u)
        out[name] = u
    return out


def gl_quad(horizon: float, n: int) -> Tuple[np.ndarray, np.ndarray]:
    x, w = np.polynomial.legendre.leggauss(n)
    return horizon / 2 * (x + 1), horizon / 2 * w


def factor_columns(op: SparseOp, dirs: List[np.ndarray], times: np.ndarray,
                   weights: np.ndarray) -> np.ndarray:
    """Columns sqrt(w_l) e^{t_l A} d for every direction d (A = op as a matrix)."""
    cols = []
    for d in dirs:
        ev = expv_cache(op, d, tuple(float(t) for t in times))
        for t, w in zip(times, weights):
            cols.append(math.sqrt(float(w)) * ev[float(t)])
    return np.column_stack(cols) if cols else np.zeros((op.size, 0))


def transpose_op(op: SparseOp) -> SparseOp:
    tr = [[] for _ in range(op.size)]
    for s, c, v in zip(op.src, op.dst, op.val):
        tr[c].append((s, float(v)))
    return SparseOp([sorted(row) for row in tr], op.size)


def spectrum(op: SparseOp, source_vecs: List[np.ndarray],
             readout_vecs: List[np.ndarray]) -> Dict:
    """Full-space finite-horizon task singular values (descending)."""
    times, wts = gl_quad(HORIZON, GL_NODES)
    C = factor_columns(op, source_vecs, times, wts)
    Z = factor_columns(transpose_op(op), readout_vecs, times, wts)
    G = Z.T @ C
    sv = np.linalg.svd(G, compute_uv=False)
    return {"singular_values": [float(x) for x in sv if x > 1e-13],
            "matrix_shape": [int(x) for x in G.shape]}


def quotient_spectrum(g: Generator, labels: Sequence[int], source_vecs_full: List[np.ndarray],
                      readout_vecs_full: List[np.ndarray],
                      weighted: bool = True) -> Dict:
    """Same task computed on the orbit quotient.

    ``weighted=True`` uses orbit-weight coordinates (== counting measure on the
    full space, so the spectra must equal the full-space ones).  ``weighted=
    False`` uses the plain Euclidean product on the quotient -- the *wrong*
    convention, kept as a diagnostic of #598's open inner-product question.
    """
    q = build_quotient(g, labels)
    nq = q["nq"]
    wts_orbit = q["weights"]
    W = np.sqrt(wts_orbit) if weighted else np.ones(nq)
    lift = q["lift"]

    # represent a full-space vector (function) on orbits: it must be orbit-constant
    def to_orbits(vf: np.ndarray) -> np.ndarray:
        vo = np.zeros(nq)
        for idx, b in enumerate(q["order"]):
            val = vf[q["member"][b][0]]
            vo[idx] = val
        return vo

    # similarity-transform the lumped generator so that the W-inner product is
    # the plain Euclidean product in transformed coordinates.  The spectral
    # radius is unchanged by the similarity, so a uniformization rate derived
    # from the *untransformed* lumped rows is safe.
    raw_bound = max(
        (sum(abs(val) for _c, val in row) for row in q["qrows"]), default=0.0
    )
    qop_rows = []
    for s, row in enumerate(q["qrows"]):
        nrow = []
        for c, val in row:
            nrow.append((c, val * W[s] / W[c]))
        qop_rows.append(nrow)
    qop = SparseOp(qop_rows, nq, rate_override=max(2.0 * raw_bound, 1e-9))

    sv_sources = [to_orbits(v) * W for v in source_vecs_full]
    sv_reads = [to_orbits(v) * W for v in readout_vecs_full]

    times, wts = gl_quad(HORIZON, GL_NODES)
    C = factor_columns(qop, sv_sources, times, wts)
    Z = factor_columns(transpose_op(qop), sv_reads, times, wts)
    G = Z.T @ C
    sv = np.linalg.svd(G, compute_uv=False)
    return {"singular_values": [float(x) for x in sv if x > 1e-13],
            "matrix_shape": [int(x) for x in G.shape]}


# --------------------------------------------------------------------------


def run_width(w: int) -> Dict:
    g = Generator(w)
    n = g.size
    pi = state_permutation(g, reflection_perm(w, w - 1))
    labels = orbit_labels(g, pi)
    n_orb = len(set(labels))

    fns = readout_functions(w)
    srcs = source_vectors(g)
    source_vecs = [np.array(srcs[name], dtype=np.float64) for name in
                   ("delta_all_singletons", "delta_single_block", "delta_wrapped_pair", "uniform")]
    readout_vecs = {name: np.array([fns[name](s) for s in g.states], dtype=np.float64)
                    for name in fns}
    centred_all = centred(readout_vecs)

    # EVEN dictionary: all readouts that are even under R
    even_names = [name for name, v in centred_all.items()
                  if np.linalg.norm(v - v[pi]) < 1e-9]
    odd_names = [name for name, v in centred_all.items()
                 if np.linalg.norm(v - v[pi]) >= 1e-9]
    even_vecs = [centred_all[name] for name in even_names]

    rows = g.rows(g.rate_baseline())
    op = SparseOp(rows, n)

    spec_full_even = spectrum(op, source_vecs, even_vecs)

    # quotient: use even readouts only (odd readouts have no orbit-constant lift)
    q_weighted = quotient_spectrum(g, labels, source_vecs, even_vecs, weighted=True)
    q_plain = quotient_spectrum(g, labels, source_vecs, even_vecs, weighted=False)

    # a marked odd source channel: signed delta on one 2-cycle orbit {i, pi(i)}
    odd_orbit = next((i, pi[i]) for i in range(n) if pi[i] != i)
    marked_src = np.zeros(n)
    marked_src[odd_orbit[0]] = 1.0
    marked_src[odd_orbit[1]] = -1.0
    marked_src /= np.linalg.norm(marked_src)

    # spectrum with one odd input channel and the odd readout(s) included
    spec_marked = spectrum(op, source_vecs + [marked_src], list(centred_all.values()))

    # direct response checks on the odd readout(s)
    times, wts = gl_quad(HORIZON, 24)
    even_src_mat = np.column_stack(source_vecs)
    odd_response_peak_even_src = {}
    odd_response_peak_odd_src = {}
    for name in odd_names:
        uo = centred_all[name] - centred_all[name][pi]  # odd part only
        if np.linalg.norm(uo) < 1e-12:
            continue
        peak_even, peak_odd = 0.0, 0.0
        ev = expv_cache(transpose_op(op), uo, tuple(float(t) for t in times))
        for t in times:
            zt = ev[float(t)]
            peak_even = max(peak_even, float(np.max(np.abs(even_src_mat.T @ zt))))
            peak_odd = max(peak_odd, float(abs(marked_src @ zt)))
        odd_response_peak_even_src[name] = peak_even
        odd_response_peak_odd_src[name] = peak_odd

    return {
        "width": w,
        "n": n,
        "n_orbit": n_orb,
        "even_readouts": even_names,
        "odd_readouts": odd_names,
        "even_task": {"full_sv": spec_full_even["singular_values"],
                      "quotient_weighted_sv": q_weighted["singular_values"],
                      "quotient_plain_sv": q_plain["singular_values"]},
        "marked_odd_input_sv": spec_marked["singular_values"],
        "odd_readout_peak_response_from_even_sources": odd_response_peak_even_src,
        "odd_readout_peak_response_from_odd_source": odd_response_peak_odd_src,
        "odd_orbit_used": [int(odd_orbit[0]), int(odd_orbit[1])],
    }


def main() -> None:
    out = {}
    for w in range(4, 9):
        t0 = time.time()
        out[str(w)] = run_width(w)
        out[str(w)]["seconds"] = round(time.time() - t0, 1)
        print(f"w={w} done in {out[str(w)]['seconds']}s", flush=True)
    payload = {
        "schema": "matching-one.probe-p398-quotient-balanced.v1",
        "issue": 598,
        "probe_issue": 599,
        "status": "exact finite computation; no sampling, no new width, no fit",
        "conventions": {
            "inner_product": "counting measure on full space; orbit-weighted on quotient",
            "horizon": 4.0,
            "gl_nodes": 32,
            "readout_centring": "counting mean subtracted",
            "constant_mode": "removed exactly on output side via centring",
        },
        "by_width": out,
    }
    dest = Path(__file__).resolve().parents[2] / "results" / "probe-p398-quotient-balanced" / "latest.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open(dest, "w") as f:
        json.dump(payload, f, indent=2, sort_keys=True)
    # summary
    print()
    for w in range(4, 9):
        r = out[str(w)]
        ev = r["even_task"]
        a = np.array(ev["full_sv"])
        b = np.array(ev["quotient_weighted_sv"])
        k = min(len(a), len(b))
        rel = np.linalg.norm(a[:k] - b[:k]) / np.linalg.norm(a[:k]) if k else 0.0
        tail_diff = abs(float(np.sum(a[k:])) - float(np.sum(b[k:])))
        cp = np.array(ev["quotient_plain_sv"])
        cp = np.pad(cp, (0, max(0, k - len(cp))))[:k]
        relp = np.linalg.norm(cp - a[:k]) / np.linalg.norm(a[:k]) if k else 0.0
        m = np.array(r["marked_odd_input_sv"])
        extra = len(m) - len(a)
        print(f"w={w} n={r['n']} n_orbit={r['n_orbit']} odd_readouts={r['odd_readouts']}")
        print(f"   full sv[{len(a)}]  head = {np.round(a[:6],5)}")
        print(f"   quot-weighted[{len(b)}] rel diff(head {k}) = {rel:.2e}  tail-sum-diff = {tail_diff:.1e}")
        print(f"   quot-PLAIN rel diff = {relp:.2e}   (inner-product diagnosis)")
        if extra > 0 and len(a) < len(m):
            print(f"   marked odd input: +{extra} extra nonzero sv; first extra = {m[len(a)]:.3e}")
        print(f"   odd readout response: even-src peak {r['odd_readout_peak_response_from_even_sources']}, "
              f"odd-src peak {r['odd_readout_peak_response_from_odd_source']}")


if __name__ == "__main__":
    main()

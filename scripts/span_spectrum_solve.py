#!/usr/bin/env python3
"""span_solve.py — stationary analysis of the span-binned retirement chain (v2).

v2: vectorised sparse construction (chunked COO aggregation), power-iteration
path for large chains, adaptive integer scale for the exact int64 certificate.

Precision policy
  n <= 300            : exact Fraction stationary solve, exact everything.
  300 < n             : float64 stationary solve (splu for moderate n, power
                        iteration for large), pi_hat = round(pi * 2^k) with
                        k chosen so the exact int64 residual fits, exact
                        residual ||pi_hat K - pi_hat||_1 in int64, certificate
                        |pi_hat.g - nu| <= ||g||_inf * ||pi_hat K - pi_hat||_1/delta.
                        If the rounding mass is too large at this k the
                        certificate is reported as None (float residual only).
Closure check (strongest validation): sum(d_h) + tail_bin == nu_w EXACTLY,
because the depth-clamped chain is an exact lumping: bins 1..D_MAX are the
exact span spectrum and the tail bin is the exact mass of span >= D_MAX+1.
"""
import json
import struct
import sys
from fractions import Fraction

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu


def load_header(path):
    dt = np.dtype([("next", "<u4"), ("nret", "u1"), ("bins", "u1", (8,)),
                   ("cnts", "u1", (8,)), ("pad", "u1", (3,))])
    f = open(path, "rb")
    W, matching, d_max, n = struct.unpack("<4i", f.read(16))
    return W, bool(matching), d_max, n, f, dt


def _popcounts(masks):
    m = np.arange(masks, dtype=np.int64)
    pc = np.zeros(masks, dtype=np.int64)
    for _ in range(64):
        if m.max() == 0:
            break
        pc += m & 1
        m >>= 1
    return pc


def build_sparse(W, d_max, n, f, dt, p_num, p_den, chunk_states=50000, out_dtype=np.int64):
    """K (n x n row-stochastic) and G (n x d_max+2).

    Streams the binary table with sequential read() in state-chunks: peak RSS is
    one chunk (~100 MB), never the 15 GB table, and no page cache accumulates.
    """
    q_num = p_den - p_num
    P = p_den
    masks = 1 << W
    pc = _popcounts(masks)
    mw = (p_num ** pc.astype(object) * q_num ** (W - pc).astype(object)).astype(np.int64)
    assert mw.max() * masks < (1 << 62)
    K_blocks, G_blocks = [], []
    for s0 in range(0, n, chunk_states):
        s1 = min(s0 + chunk_states, n)
        m = s1 - s0
        buf = f.read(m * masks * 24)
        rec = np.frombuffer(buf, dtype=dt)
        nx_ = rec["next"].astype(np.int64)
        nr_ = rec["nret"]
        bins_ch = rec["bins"]
        cnts_ch = rec["cnts"]
        rows_all = np.repeat(np.arange(m, dtype=np.int64), masks)
        mw_all = np.tile(mw, m)
        Kb = sparse.coo_matrix((mw_all, (rows_all, nx_)), shape=(m, n), dtype=np.int64).tocsr().astype(out_dtype)
        K_blocks.append(Kb)
        Gd = np.zeros((m, d_max + 2), dtype=np.int64)
        for q in range(8):
            sel = nr_ > q
            if not sel.any():
                continue
            rows_local = np.nonzero(sel)[0] // masks
            bcol = bins_ch[:, q][sel]
            cval = mw_all[sel] * cnts_ch[:, q][sel].astype(np.int64)
            np.add.at(Gd, (rows_local, bcol), cval)
        G_blocks.append(Gd)
    K = K_blocks[0] if len(K_blocks) == 1 else sparse.vstack(K_blocks, format="csr", dtype=out_dtype)
    G = np.vstack(G_blocks)
    assert (np.asarray(K.sum(axis=1)).ravel() == P ** W).all(), "rows not stochastic"
    return K.tocsr(), G


def stationary(K_float, mode="auto"):
    n = K_float.shape[0]
    if mode == "splu" or (mode == "auto" and n <= 50000):
        A = K_float.T.tocsr() - sparse.eye(n, dtype=np.float64)
        A = A.tolil()
        A[n - 1, :] = 1.0
        A = A.tocsr()
        rhs = np.zeros(n); rhs[n - 1] = 1.0
        lu = splu(A.tocsc())
        x = lu.solve(rhs)
        x = x - lu.solve(A @ x - rhs)          # one refinement
        return x
    # power iteration: pi_{t+1} = pi_t K   (row vector)
    pi = np.full(n, 1.0 / n)
    KT = K_float.T.tocsr()
    tol = 1e-13 if K_float.dtype == np.float64 else 1e-7
    for it in range(20000):
        nxt_pi = KT @ pi
        d = np.abs(nxt_pi - pi).sum()
        pi = nxt_pi
        if d < tol:
            break
    return pi


def analyse(path, p_num, p_den, nu_ref=None, label="", light=False, d_max_override=None):
    W, matching, d_max, n, f, dt = load_header(path)
    if d_max_override is not None:
        d_max = d_max_override
    K_int, Gd = build_sparse(W, d_max, n, f, dt, p_num, p_den,
                             out_dtype=(np.float32 if light else np.int64))
    f.close()
    P = p_den
    P_W = P ** W
    delta = Fraction(P - p_num, P) ** W
    eq9 = Fraction(W) * (1 - (1 - Fraction(p_num, p_den)) ** W) ** d_max
    ginf_num = int(np.abs(Gd).max()) if Gd.size else 0

    exact = (n <= 300) and not light
    if exact:
        Kf = K_int.toarray().astype(object)
        A = [[Fraction(int(Kf[j, i]), P_W) - (1 if i == j else 0) for j in range(n)] for i in range(n)]
        A[-1] = [Fraction(1)] * n
        M = [A[i][:] + [Fraction(0)] for i in range(n)]
        M[n - 1][n] = Fraction(1)
        for c in range(n):
            piv = next(r for r in range(c, n) if M[r][c] != 0)
            M[c], M[piv] = M[piv], M[c]
            pv = M[c][c]
            M[c] = [x / pv for x in M[c]]
            for r_ in range(n):
                if r_ != c and M[r_][c] != 0:
                    f = M[r_][c]
                    M[r_] = [x - f * y for x, y in zip(M[r_], M[c])]
        pi_f = [M[i][n] for i in range(n)]
        nu = [sum(pi_f[s] * int(Gd[s, h]) for s in range(n)) / P_W for h in range(d_max + 2)]
        resid = sum(abs(sum(pi_f[s] * Fraction(int(Kf[s, j]), P_W) for s in range(n)) - pi_f[j])
                    for j in range(n))
        ginf = max(Fraction(int(Gd[s, h]), P_W) for s in range(n) for h in range(d_max + 2)) if Gd.size else Fraction(0)
        cert = str(ginf * resid / delta)
        mode = "exact-rational"
        dh = nu[1:d_max + 1]
        tail = nu[d_max + 1]
        dh_f = [float(x) for x in dh]
        exact_sum = sum(dh)
        exact_tail = tail
        tail_f = float(nu[d_max + 1])
        total_f = float(sum(nu[1:d_max + 1]))
        nu_total_f = float(sum(nu))
        fres = None
    else:
        if light:
            # memory-light: float64 chain only, no int64 certificate.
            # Validation rests on the exact closure against the independently
            # certified nu_w of the winding_build engine (#741).
            Kf = K_int
            pi = stationary(Kf)
            pi = np.maximum(pi, 0.0)
            pi /= pi.sum()
            nu_f = (pi @ Gd) / P_W
            dh_f = [float(x) for x in nu_f[1:d_max + 1]]
            tail_f = float(nu_f[d_max + 1])
            total_f = float(sum(dh_f))
            nu_total_f = total_f + tail_f
            fres = float(np.abs(Kf.T @ pi - pi).sum())
            cert = "skipped(light; closure vs certified nu_w)"
            mode = "float64(power, light)"
            out = {
                "file": path, "label": label, "width": W, "matching": matching,
                "p": f"{p_num}/{p_den}", "states": n, "d_max": d_max, "mode": mode,
                "d_h_float": dh_f, "tail_bin_float": tail_f,
                "sum_dh_float": total_f, "nu_total_float": nu_total_f,
                "certificate_bound": cert, "float_residual_l1": fres,
                "tail_bound_eq9": str(eq9), "tail_bound_eq9_float": float(eq9),
                "delta": str(delta),
            }
            if nu_ref is not None:
                ref = Fraction(nu_ref)
                out["nu_ref"] = nu_ref
                out["closure_exact"] = False
                out["gap_float"] = float(ref) - nu_total_f
            return out
        Kf = K_int.astype(np.float64) / P_W
        pi = stationary(Kf)
        pi = np.maximum(pi, 0.0)
        pi /= pi.sum()
        nu_f = (pi @ Gd) / P_W
        dh_f = [float(x) for x in nu_f[1:d_max + 1]]
        tail_f = float(nu_f[d_max + 1])
        total_f = float(sum(dh_f))
        nu_total_f = total_f + tail_f
        fres = float(np.abs(Kf.T @ pi - pi).sum())
        cert = None
        # exact int64 certificate where the scale fits
        try:
            k = min(40, max(20, 62 - int(P_W).bit_length() - 4))
            a = np.round(pi * (1 << k)).astype(np.int64)
            a = np.maximum(a, 0)
            a[int(np.argmax(a))] -= int(a.sum()) - (1 << k)
            assert (a >= 0).all() and int(a.sum()) == (1 << k)
            r_num = np.zeros(n, dtype=np.int64)
            cbits = max(4, 58 - int(n).bit_length() - int(P_W).bit_length())
            for t in range(0, k, cbits):
                ch = ((a >> t) & ((1 << cbits) - 1)).astype(np.int64)
                if not ch.any():
                    continue
                c = K_int.T @ ch
                if np.abs(c).max() >= (1 << (60 - t)):
                    raise OverflowError
                r_num += c << t
            r_num -= a * P_W
            assert (np.abs(r_num) < (1 << 62)).all()
            l1 = int(np.abs(r_num).sum())
            bound = Fraction(ginf_num, P_W) * Fraction(l1, (1 << k) * P_W) / delta
            cert = str(bound) if bound < Fraction(1, 10 ** 6) else f"WEAK({float(bound):.2e})"
        except (OverflowError, AssertionError):
            cert = None
        mode = f"float64({'splu' if n <= 50000 else 'power'})"

    out = {
        "file": path, "label": label, "width": W, "matching": matching,
        "p": f"{p_num}/{p_den}", "states": n, "d_max": d_max, "mode": mode,
        "d_h": ([str(x) for x in dh] if exact else None),
        "tail_bin": (str(tail) if exact else None),
        "sum_dh": (str(exact_sum) if exact else None),
        "d_h_float": dh_f,
        "tail_bin_float": tail_f,
        "sum_dh_float": total_f,
        "nu_total_float": nu_total_f,
        "certificate_bound": cert,
        "tail_bound_eq9": str(eq9), "tail_bound_eq9_float": float(eq9),
        "delta": str(delta),
    }
    if fres is not None:
        out["float_residual_l1"] = fres
    if nu_ref is not None:
        ref = Fraction(nu_ref)
        out["nu_ref"] = nu_ref
        if exact:
            got = exact_sum + exact_tail
            out["closure_exact"] = bool(ref == got)
            out["gap_to_ref"] = str(ref - got)
            out["gap_float"] = float(ref) - (total_f + tail_f)
        else:
            out["closure_exact"] = False
            out["gap_float"] = float(ref) - (total_f + tail_f)
    return out



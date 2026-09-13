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


CENSOR_TOL = 1e-6


def spectrum_moments(dh, tail_bin, d_max):
    """Binned moments of the span spectrum d_1..d_{d_max} plus the tail bin.

    The tail bin collects span >= d_max+1 and is EXCLUDED from the moments, so
    E[L] and Var/E^2 are censored readouts. `censored` flags any configuration
    whose tail fraction is not negligible; at such a width the moments are a
    lower bound and must be quoted as censored (see the 2026-09-13 erratum).
    """
    n = len(dh)
    tot = float(sum(dh))
    m1 = float(sum((i + 1) * dh[i] for i in range(n)))
    m2 = float(sum(float(i + 1) ** 2 * dh[i] for i in range(n)))
    E = m1 / tot if tot else float("nan")
    E2 = m2 / tot if tot else float("nan")
    V = E2 - E * E
    tail = float(tail_bin)
    return {
        "d_max": d_max,
        "sum_dh": tot,
        "tail_bin": tail,
        "nu_total": tot + tail,
        "tail_fraction": tail / (tot + tail) if (tot + tail) else 0.0,
        "E_L": E,
        "E_L2": E2,
        "var_over_E2": V / (E * E) if E else float("nan"),
        "censored": bool(tail / (tot + tail) > CENSOR_TOL if (tot + tail) else False),
    }


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
            #
            # FIX (erratum 2026-09-13, PR #739 issuecomment-5653178247): this
            # branch used to hand the RAW integer-weight float32 chain to
            # stationary(). The stationary distribution is scale invariant, but
            # neither solver is: with K unscaled, (K^T - I) is nonsingular so the
            # splu path returns a non-stationary vector whose residual is O(p_den**W)
            # (measured: residual 255.0, nu total x84 wrong at w=4), and the power
            # iteration grows by p_den**W per sweep and overflows to inf/nan
            # (measured: nan at w=5 NN p=1/4). Divide exactly as the normal branch.
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
            cert = "skipped(light; closure vs certified nu_w)"
            mode = f"float64(light/{'splu' if n <= 50000 else 'power'})"
            out = {
                "file": path, "label": label, "width": W, "matching": matching,
                "p": f"{p_num}/{p_den}", "states": n, "d_max": d_max, "mode": mode,
                "d_h_float": dh_f, "tail_bin_float": tail_f,
                "sum_dh_float": total_f, "nu_total_float": nu_total_f,
                "certificate_bound": cert, "float_residual_l1": fres,
                "tail_bound_eq9": str(eq9), "tail_bound_eq9_float": float(eq9),
                "delta": str(delta),
                "ginf_num": ginf_num,
                "moments": spectrum_moments(dh_f, tail_f, d_max),
                "cert_candidate": None,
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
        pi_cert = None
        k = None
        # The precision claim has two separate parts (erratum 2026-09-13):
        #   bound_cert  bounds |observable(pi_cert) - observable(pi_exact)|, where
        #               pi_cert = a/2^k is the RATIONAL candidate, not the float pi;
        #   bound_float bounds |observable(pi_float) - observable(pi_exact)| by the
        #               same lemma, using the reported float residual.
        # Reporting bound_cert alone certifies the candidate, not the printed value;
        # we therefore emit the certified candidate's observables and the float-side
        # bound as well, so no gap is left implicit.
        bound_float = float(Fraction(ginf_num, P_W) * Fraction(fres) / delta)
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
            pi_cert = a.astype(np.float64) / float(1 << k)
        except (OverflowError, AssertionError):
            cert = None
            pi_cert = None
            k = None
        run_bound_cert = None
        if pi_cert is not None:
            nu_c = (pi_cert @ Gd) / P_W
            dh_c = [float(x) for x in nu_c[1:d_max + 1]]
            tail_c = float(nu_c[d_max + 1])
            mc = spectrum_moments(dh_c, tail_c, d_max)
            mf = spectrum_moments(dh_f, tail_f, d_max)
            run_bound_cert = float(cert.split("(")[1].rstrip(")")) if cert.startswith("WEAK") else float(
                Fraction(cert) if cert else 0.0)
            _cert_block = {
                "k": k,
                "d_h_float": dh_c,
                "sum_dh_float": mc["sum_dh"],
                "tail_bin_float": tail_c,
                "moments": mc,
                "max_abs_diff_dh_vs_float": max(abs(dh_c[i] - dh_f[i]) for i in range(len(dh_f))),
                "rel_diff_sum_dh_vs_float": abs(mc["sum_dh"] - mf["sum_dh"]) / mf["sum_dh"],
                "observable_bound_cert": run_bound_cert,
                "observable_bound_float": bound_float,
                "observable_bound_sum": (
                    run_bound_cert + bound_float if run_bound_cert is not None else None),
            }
        else:
            _cert_block = None
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
        "ginf_num": ginf_num,
        "moments": spectrum_moments(dh_f, tail_f, d_max),
        "cert_candidate": (_cert_block if not exact else None),
    }
    if exact:
        out["moments_exact_denominator"] = str(Fraction(P_W) ** 2)
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


def run_jobs(spec_path, out_path):
    """Rerun a list of configurations from a small JSON job spec.

    spec: [{"table":..., "p":"1/4", "d_max":48, "light":false,
            "nu_ref":"...", "label":"NN w4 p=1/4"}, ...]
    Written so that any delivered spectrum can be regenerated from the committed
    scripts without the ad-hoc driver used on 2026-09-13 (erratum note, §5).
    """
    import json
    with open(spec_path) as fh:
        jobs = json.load(fh)
    rows = []
    for j in jobs:
        pn, pd = j["p"].split("/")
        r = analyse(j["table"], int(pn), int(pd), nu_ref=j.get("nu_ref"),
                    label=j.get("label", ""), light=bool(j.get("light", False)),
                    d_max_override=j.get("d_max"))
        rows.append(r)
        m = r["moments"]
        print("%-24s W=%d %-6s n=%-8d mode=%-28s E[L]=%.6f var/E^2=%.6f tail=%.2e%s"
              % (r["label"], r["width"], r["p"], r["states"], r["mode"], m["E_L"],
                 m["var_over_E2"], m["tail_fraction"], "  CENSORED" if m["censored"] else ""))
    with open(out_path, "w") as fh:
        json.dump(rows, fh, indent=1)
    return rows


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="span-spectrum solver / moment reporter")
    ap.add_argument("spec", help="job-spec JSON (list of configurations)")
    ap.add_argument("out", help="where to write the result JSON")
    args = ap.parse_args()
    run_jobs(args.spec, args.out)



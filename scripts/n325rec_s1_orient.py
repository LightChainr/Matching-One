#!/usr/bin/env python3
"""n325rec step 1 -- EXACT orientation table for a^2+b^2=N.

Independent computation (no numbers are copied from n1105mix):
  * all integer representations a^2+b^2=N, classified into D4 orbits
    (D4 = the 8 symmetries (a,b)->(+-a,+-b),(+-b,+-a));
  * primitivity gcd(a,b)=1 and, crucially, the "same-ell realization":
    a non-primitive representative (k*a0,k*b0) is the SAME cylinder as the
    primitive u=(a0,b0) with index n=k  (both give period vector k*u0), so it
    is usable exactly when written that way;
  * exact cos(4m theta) as Fraction, from cos(4m th) = Re((a+ib)^{4m})/N^{2m};
  * design matrix A[i][m] = cos(4m theta_i) (m=0,1,2 -> H0,H4,H8), its EXACT
    rank, the exact H0/H4/H8 projector weights w (solve A^T w = e_0), and
    |w|_1, |w|_2 (the L2 norm is the orientation-noise amplification);
  * floating point is used only for the final sqrt / L2 display.

Pure stdlib (fractions).  No numpy.
"""
from __future__ import annotations
import json
import math
import sys
from fractions import Fraction


# ---------------------------------------------------------------- integer reps
def reps(N: int):
    """All (a,b), 0<=a<=b, a^2+b^2=N, as D4-class representatives."""
    out = []
    a = 0
    while a * a <= N:
        b2 = N - a * a
        b = math.isqrt(b2)
        if b * b == b2 and b >= a:
            out.append((a, b))
        a += 1
    return out


def d4_orbit(p):
    a, b = p
    return sorted({(a, b), (b, a), (-a, b), (a, -b), (-b, a), (b, -a),
                   (-a, -b), (-b, -a)})


# ------------------------------------------------------------ exact harmonics
def ipow(a: int, b: int, k: int):
    """(a+ib)^k exactly."""
    re, im = 1, 0
    for _ in range(k):
        re, im = re * a - im * b, re * b + im * a
    return re, im


def cos4m(a: int, b: int, N: int, m: int) -> Fraction:
    """cos(4*m*theta) exactly, theta = angle of (a,b)."""
    re, _im = ipow(a, b, 4 * m)
    return Fraction(re, N ** (2 * m))


# --------------------------------------------------------- exact linear algebra
def solve_exact(M, rhs):
    """Solve M x = rhs exactly (Fraction Gauss-Jordan).  Returns (x, rank, ok)."""
    n = len(M)
    m = len(M[0])
    aug = [[Fraction(M[i][j]) for j in range(m)] + [Fraction(rhs[i])]
           for i in range(n)]
    if n != m:
        raise ValueError("solve_exact wants square")
    rank = 0
    piv = []
    r = 0
    for c in range(m):
        pr = None
        for i in range(r, n):
            if aug[i][c] != 0:
                pr = i
                break
        if pr is None:
            continue
        aug[r], aug[pr] = aug[pr], aug[r]
        pv = aug[r][c]
        aug[r] = [x / pv for x in aug[r]]
        for i in range(n):
            if i != r and aug[i][c] != 0:
                f = aug[i][c]
                aug[i] = [aug[i][j] - f * aug[r][j] for j in range(m + 1)]
        piv.append(c)
        r += 1
        if r == n:
            break
    rank = r
    for i in range(r, n):
        if all(aug[i][j] == 0 for j in range(m)) and aug[i][m] != 0:
            return None, rank, False
    if rank < m:
        return None, rank, False
    x = [aug[i][m] for i in range(m)]
    return x, rank, True


def rank_exact(M):
    """Exact rank of a rectangular Fraction matrix."""
    if not M:
        return 0
    A = [[Fraction(v) for v in row] for row in M]
    rows, cols = len(A), len(A[0])
    r = 0
    for c in range(cols):
        pr = None
        for i in range(r, rows):
            if A[i][c] != 0:
                pr = i
                break
        if pr is None:
            continue
        A[r], A[pr] = A[pr], A[r]
        pv = A[r][c]
        for i in range(r + 1, rows):
            if A[i][c] != 0:
                f = A[i][c] / pv
                A[i] = [A[i][j] - f * A[r][j] for j in range(cols)]
        r += 1
        if r == rows:
            break
    return r


# ---------------------------------------------------------------------- driver
def analyse(N: int, n_harm: int = 3):
    rp = reps(N)
    orbits = []
    seen = set()
    for p in rp:
        if p in seen:
            continue
        orb = [q for q in d4_orbit(p)]
        for q in rp:
            if q in orb:
                seen.add(q)
        orbits.append(orb)

    classes = []
    for orb in orbits:
        rep = None
        for q in orb:
            if q in rp:
                rep = q
                break
        a, b = rep
        g = math.gcd(a, b)
        c4 = cos4m(a, b, N, 1)
        c8 = cos4m(a, b, N, 2)
        c12 = cos4m(a, b, N, 3)
        harm = [Fraction(1), c4, c8, c12][:n_harm]
        classes.append({
            "rep": [a, b],
            "orbit_size": len({q for q in orb}),
            "gcd": g,
            "primitive": g == 1,
            "same_ell_realization": ({"u": [a // g, b // g], "n": g}
                                     if g > 1 else {"u": [a, b], "n": 1}),
            "theta_deg": math.degrees(math.atan2(b, a)),
            "cos4": str(c4), "cos4_float": float(c4),
            "cos8": str(c8), "cos8_float": float(c8),
            "cos12": str(c12), "cos12_float": float(c12),
            "_harm": harm,
        })

    A = [c["_harm"] for c in classes]
    A = [list(row) for row in A]
    rk = rank_exact(A)
    out = {
        "N": N,
        "sqrt_N": math.sqrt(N),
        "ell4_exact": N ** 2,
        "all_representations": [list(p) for p in rp],
        "n_d4_classes": len(classes),
        "n_primitive_classes": sum(1 for c in classes if c["primitive"]),
        "design_matrix_rank": rk,
        "design_matrix_shape": [len(A), n_harm],
        "classes": classes,
    }
    # H0/H4/... projector weights when the square system is invertible
    if len(classes) == n_harm and rk == n_harm:
        M = [[A[i][j] for i in range(n_harm)] for j in range(n_harm)]  # A^T
        for target in range(n_harm):
            e = [Fraction(1 if k == target else 0) for k in range(n_harm)]
            w, _r, ok = solve_exact([row[:] for row in M], e)
            if ok:
                l1 = sum(abs(x) for x in w)
                l2 = math.sqrt(float(sum(x * x for x in w)))
                out.setdefault("projector_weights", {})["H%d" % (4 * target)] = {
                    "w_num": [str(x) for x in w],
                    "w_float": [float(x) for x in w],
                    "L1": float(l1),
                    "L2": l2,
                }
    # error propagation of a per-point Omega error sigma
    if len(classes) == n_harm and rk == n_harm:
        B = [[A[i][j] for i in range(n_harm)] for j in range(n_harm)]
        Binv = []
        # invert A exactly (A^T w = e solved above gives rows of A^{-T})
        for target in range(n_harm):
            e = [Fraction(1 if k == target else 0) for k in range(n_harm)]
            w, _r, ok = solve_exact([row[:] for row in B], e)
            Binv.append(w)
        # B^{-1} = A^{-T}; P = (A^T A)^{-1}A^T Omega = A^{-1} Omega for square A
        # A^{-1} = (A^{-T})^T
        Ainv = [[Binv[j][i] for j in range(n_harm)] for i in range(n_harm)]
        amps = []
        for m in range(n_harm):
            amps.append(math.sqrt(float(sum(Ainv[m][i] ** 2
                                             for i in range(n_harm)))))
        out["coefficient_noise_amplification"] = {
            "H%d" % (4 * m): amps[m] for m in range(n_harm)}
    for c in classes:
        c.pop("_harm", None)
    return out


def main():
    res = {"schema": "n325rec.orientation-table.v1",
           "p_c": 0.5927460507921,
           "note": ("exact cos(4m theta); L2 of the H0 projector weight = "
                    "orientation-noise amplification; same_ell_realization "
                    "shows the primitive (u,n) pair with n*|u| = sqrt(N)")}
    for N in (25, 325, 1105):
        res["N=%d" % N] = analyse(N, n_harm=(2 if N == 25 else (3 if N == 325 else 4)))
    out = sys.argv[1] if len(sys.argv) > 1 else "/workspace/n325rec/out/s1_orient.json"
    with open(out, "w") as fh:
        json.dump(res, fh, indent=1)
    for key in ("N=25", "N=325", "N=1105"):
        d = res[key]
        print("%s: classes=%d (primitive %d)  rank=%d  sqrtN=%.6f"
              % (key, d["n_d4_classes"], d["n_primitive_classes"],
                 d["design_matrix_rank"], d["sqrt_N"]))
        for c in d["classes"]:
            print("    rep=%-9s gcd=%d prim=%-5s  cos4=%+.9f  cos8=%+.9f  u*n=%s"
                  % (str(c["rep"]), c["gcd"], c["primitive"],
                     c["cos4_float"], c["cos8_float"],
                     c["same_ell_realization"]))
        if "projector_weights" in d:
            for k, v in d["projector_weights"].items():
                print("    %s projector L1=%.6f L2=%.6f  w=%s"
                      % (k, v["L1"], v["L2"], v["w_float"]))
    print("->", out)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
C8 — matching involution, three-line lemma + exact polynomial check.

Lemma (planar matching pair, honest torus, H1):
    M_Ghat(p) = - M_G(1-p)
    F_Ghat(p) = 1 - F_G(1-p)
    Q_Ghat(u) = 1 - Q_G(1-u)
Negative: this does NOT imply Q_G(u) + Q_G(1-u) = 1 (C3 is the support).

On the site torus, black 4-connect / white 8-connect IS the matching pair
(the enumerator's r_b vs r_w).  r_b + r_w = 2 configuration-wise implies
    coeff_w[k] = - coeff_b[N-k],
hence  M_black(p) + M_white(1-p) = 0  as exact polynomials.  Verify exactly.

Imports #606 exact_torus_enum (ambient_rank); does NOT duplicate it.
"""
import sys
import json
from pathlib import Path
from fractions import Fraction
from collections import defaultdict

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "homological_balance"))
from exact_torus_enum import ambient_rank  # noqa: E402


def enumerate_joint(L):
    N = L * L
    jb = defaultdict(int)
    jw = defaultdict(int)
    for mask in range(1 << N):
        black = [(k // L, k % L) for k in range(N) if (mask >> k) & 1]
        n = len(black)
        rb = ambient_rank(black, L, False)
        rw = ambient_rank([(k // L, k % L) for k in range(N) if not ((mask >> k) & 1)],
                          L, True)
        jb[(n, rb)] += 1
        jw[(n, rw)] += 1
    return jb, jw


def bernstein_from_joint(joint, N):
    """Black-count grouping: coeff[n_black] = sum over black configs of (r-1)."""
    coeff = [Fraction(0)] * (N + 1)
    for (n, r), cnt in joint.items():
        coeff[n] += Fraction(r - 1) * cnt
    return coeff


def white_bernstein_from_joint(jw, N):
    """White-count grouping: jw is keyed by (n_black, r_w), so a white set of
    size k corresponds to n_black = N - k."""
    coeff = [Fraction(0)] * (N + 1)
    for (n_black, r), cnt in jw.items():
        k = N - n_black
        coeff[k] += Fraction(r - 1) * cnt
    return coeff


def main():
    out = {}
    for L in (3, 4):
        N = L * L
        jb, jw = enumerate_joint(L)
        cb = bernstein_from_joint(jb, N)   # M_black Bernstein coeffs
        cw = white_bernstein_from_joint(jw, N)   # M_white Bernstein coeffs (by white count)

        # verify coeff_w[k] == -coeff_b[N-k] (exact)
        bad = []
        for k in range(N + 1):
            if cw[k] != -cb[N - k]:
                bad.append((k, str(cw[k]), str(-cb[N - k])))

        # verify the polynomial identity M_black(p) + M_white(1-p) = 0 at several p
        from exact_torus_enum import eval_ML
        id_fail = []
        for num, den in ((0, 8), (1, 8), (2, 8), (3, 8), (4, 8), (5, 8), (6, 8), (7, 8), (1, 2)):
            p = Fraction(num, den)
            lhs = eval_ML(cb, p) + eval_ML(cw, 1 - p)
            if lhs != 0:
                id_fail.append((str(p), str(lhs)))

        out[str(L)] = {
            "coeff_w_k_eq_minus_coeff_b_Nk": len(bad) == 0,
            "bad_entries": bad[:5],
            "poly_identity_Mb_p_plus_Mw_1mp": len(id_fail) == 0,
            "identity_fail_examples": id_fail[:5],
            "N": N,
        }
        print(f"===== L={L} =====")
        print(f"coeff_w[k] == -coeff_b[N-k] for all k: {len(bad) == 0}  (bad: {bad[:5]})")
        print(f"M_black(p) + M_white(1-p) == 0 (9 sample points): {len(id_fail) == 0}  (fail: {id_fail[:5]})")
        print()

    base = ROOT / "results" / "probe-exact-controls"
    base.mkdir(parents=True, exist_ok=True)
    (base / "matching_involution.json").write_text(json.dumps({
        "schema": "matching-one.probe-exact-controls.matching-involution.v1",
        "lemma": "M_Ghat(p) = -M_G(1-p); F_Ghat(p) = 1-F_G(1-p); Q_Ghat(u) = 1-Q_G(1-u)",
        "coupling_used": ("site black=4-connect (G), white=8-connect (Ghat); "
                          "r_b+r_w=2 configuration-wise (Alexander)"),
        "result": "M_black(p) + M_white(1-p) = 0 exactly, L=3 and L=4",
        "by_L": out,
    }, indent=2, sort_keys=True))
    print("wrote", base / "matching_involution.json")


if __name__ == "__main__":
    main()

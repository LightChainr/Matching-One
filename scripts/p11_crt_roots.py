#!/usr/bin/env python3
"""Combine CRT passes from p11_exact_dp.cpp into exact A(k), run sanity
checks, compute p_med / p_cell roots with mpmath, and compare digit-by-digit
with the published Mertens 2022 table (data/mertens_2022_square_site_estimators.csv).

Exactness contract: A(k) are reconstructed as exact integers (Chinese
remainder over 30-bit primes); no float enters the integer path.  Roots are
computed at 60 significant digits; the reported digit match counts the common
prefix with the published string (30 decimals "exact as printed" per the
source).

Usage:
  p11_crt_roots.py --crt file1.json [file2.json ...] --n N [--cell-crt f2 ...]
      --out out.json
  p11_crt_roots.py --f128-audit --n N --exact file.json --digits-file root.json
"""

from __future__ import annotations

import argparse
import json
import sys
from fractions import Fraction
from pathlib import Path

from mpmath import mp, mpf

ROOT = Path(__file__).resolve().parents[1]


def crt_pair(r1, m1, r2, m2):
    g = pow(m1, -1, m2)
    k = ((r2 - r1) * g) % m2
    return (r1 + m1 * k) % (m1 * m2), m1 * m2


def combine(passes):
    rs, ms = [], []
    for p in passes:
        rs.append(int(p["prime"]))
        ms.append(p["A"])
    nm = len(ms[0]) - 1
    A = [0] * (nm + 1)
    for k in range(nm + 1):
        r, m = ms[0][k], rs[0]
        for j in range(1, len(rs)):
            r, m = crt_pair(r, m, ms[j][k], rs[j])
        A[k] = r
    return A


def sanity(n, m, A):
    assert A[m] == n, f"A[{m}]={A[m]} != {n}"
    from math import comb
    for k, a in enumerate(A):
        assert 0 <= a <= comb(n * m, k), f"A[{k}] out of range"
    parity = sum(c * (-1) ** k for k, c in enumerate(A) if c)
    return parity


def R_of_poly(A, p, nm):
    total = mpf(0)
    q = 1 - p
    pk = mpf(1)
    for k, c in enumerate(A):
        if c:
            total += c * pk * q ** (nm - k)
            pk *= p
        else:
            pk *= p
    return total


def bisect_root(f, lo, hi, iters=170):
    flo = f(lo)
    assert flo < 0, "f(lo) must be negative"
    for _ in range(iters):
        mid = (lo + hi) / 2
        fm = f(mid)
        if (fm < 0) == (flo < 0):
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def p_med_root(A):
    mp.dps = 60
    nm = len(A) - 1
    f = lambda p: R_of_poly(A, p, nm) - mpf(1) / 2
    return bisect_root(f, mpf("0.4"), mpf("0.75"))


def p_cell_root(A_n, A_nm1):
    mp.dps = 60
    f = lambda p: (R_of_poly(A_n, p, len(A_n) - 1)
                   - R_of_poly(A_nm1, p, len(A_nm1) - 1))
    return bisect_root(f, mpf("0.4"), mpf("0.75"))


def published():
    path = ROOT / "data" / "mertens_2022_square_site_estimators.csv"
    rows = {}
    lines = path.read_text().strip().splitlines()
    hdr = lines[0].split(",")
    for ln in lines[1:]:
        parts = ln.split(",")
        d = dict(zip(hdr, parts))
        n = int(d["n"])
        rows[n] = (d["p_med"], d["p_cell"])
    return rows


def match_digits(mine_str, pub_str):
    if not pub_str:
        return None
    a, b = mine_str, pub_str
    common = 0
    for x, y in zip(a, b):
        if x == y:
            common += 1
        else:
            break
    return common


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--crt", nargs="+", required=True,
                    help="CRT json file(s) for width n (n-1 files may be "
                         "given via --cell-crt or inferred)")
    ap.add_argument("--cell-crt", nargs="*", default=None,
                    help="CRT json file(s) for width n-1 (for p_cell)")
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--m", type=int, default=None)
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()
    m = args.m or args.n

    def load_and_combine(files):
        passes = []
        for fp in files:
            d = json.loads(Path(fp).read_text())
            passes.extend(d["passes"])
        return combine(passes)

    A = load_and_combine(args.crt)
    parity = sanity(args.n, m, A)
    mp.dps = 60
    med = p_med_root(A)

    A_prev = None
    if args.cell_crt:
        A_prev = load_and_combine(args.cell_crt)
    elif args.n > 1:
        # try sibling file convention: same dir, name with n-1
        cands = []
        for fp in args.crt:
            p = Path(fp)
            s = str(p)
            for tag in (f"n{args.n}", f"_{args.n}.", f"-{args.n}."):
                if tag in s:
                    cands.append(Path(s.replace(tag, tag.replace(
                        str(args.n), str(args.n - 1)))))
        for c in cands:
            if c.exists():
                A_prev = load_and_combine([str(c)])
                break

    out = {
        "n": args.n, "m": m,
        "A_exact": [str(x) for x in A],
        "F_minus1": str(parity),
        "p_med": mp.nstr(med, 42),
    }
    pub = published()
    pub_med, pub_cell = pub[args.n]
    out["published_p_med"] = pub_med
    out["p_med_digit_match"] = match_digits(out["p_med"], pub_med)
    if A_prev is not None:
        sanity(args.n - 1, args.n - 1, A_prev)
        cell = p_cell_root(A, A_prev)
        out["p_cell"] = mp.nstr(cell, 42)
        if pub_cell:
            out["published_p_cell"] = pub_cell
            out["p_cell_digit_match"] = match_digits(out["p_cell"], pub_cell)
    print(json.dumps(out, indent=2))
    if args.out:
        Path(args.out).write_text(json.dumps(out))


if __name__ == "__main__":
    main()

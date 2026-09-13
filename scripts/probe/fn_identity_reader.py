#!/usr/bin/env python3
"""
Issue #678 tiny reader — re-verify, exactly, that F_L = (1 + M_L)/2 is a
polynomial identity of the rank triple on the committed axis L=3,4 tables.

Reads ONLY committed artifacts (no enumeration):
  results/probe-exact-controls/joints.json        (per-(k, r_b) counts)
  results/probe-exact-controls/M_poly_L{3,4}.json (committed M, power basis)
  results/issue665-joint-rank-wrap/axis-L{3,4}.json (PR #668 five-cell rows)

All arithmetic is fractions.Fraction. No floats. Exits nonzero on any failure.

Checks per L:
  A  P_0 + P_1 + P_2 = 1                         (layer partition)
  B  P_2 - P_0 == committed M_coeffs             (coefficientwise)
  C  (P_1 + 2 P_2)/2 == (1 + M)/2                (the identity, coefficientwise)
  D  F(0)=0, F(1)=1, F(1/2) == (1+M(1/2))/2      (Fraction evaluation)
  E  a_k = #cross_k - #none_k == PR #668 collapsed_D / committed_bernstein
  F  five-cell: P_1 == x + y + spiral mass; P_0 == none; P_2 == exclusive cross
"""
import json
from fractions import Fraction
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "results"
COMMITTED_F_HALF = {3: Fraction(43, 128), 4: Fraction(19011, 65536)}


def counts_from_joints(L):
    doc = json.loads((RES / "probe-exact-controls" / "joints.json").read_text())
    d = doc["by_L"][str(L)]["joint_black_count_rank"]
    per = {}
    for k, v in d.items():
        n, r = map(int, k.split(","))
        per[(n, r)] = per.get((n, r), 0) + v
    return per


def counts_from_pr668(L):
    p = RES / "issue665-joint-rank-wrap" / f"axis-L{L}.json"
    if not p.exists():
        return None, None
    doc = json.loads(p.read_text())
    per = {}
    for row in doc["joint"]:
        key = (row["k"], row["r_black"])
        per[key] = per.get(key, 0) + row["count"]
    return per, doc


def counts_to_power(per, N):
    """Power-basis coeffs of P_r(p) = sum_k c_{k,r} p^k (1-p)^{N-k}, per rank."""
    A = [[Fraction(0)] * (N + 1) for _ in range(3)]
    for (k, r), c in per.items():
        for j in range(N - k + 1):
            sign = -1 if (j & 1) else 1
            A[r][k + j] += c * sign * comb(N - k, j)
    return A


def ev(a, x):
    r = Fraction(0)
    for c in reversed(a):
        r = r * x + c
    return r


def check(L):
    N = L * L
    per_main = counts_from_joints(L)
    per_668, doc668 = counts_from_pr668(L)
    mdoc = json.loads((RES / "probe-exact-controls" / f"M_poly_L{L}.json").read_text())
    M = [Fraction(e["c"]) for e in mdoc["M_coeffs"]]
    M_half = Fraction(mdoc["M_half"])

    results = []
    for tag, per in (("main", per_main), ("pr668", per_668)):
        if per is None:
            continue
        # A: layer partition
        layers_ok = all(
            sum(per.get((k, r), 0) for r in range(3)) == comb(N, k) for k in range(N + 1)
        )
        A = counts_to_power(per, N)
        one = [Fraction(1) if m == 0 else Fraction(0) for m in range(N + 1)]
        B = [A[2][m] - A[0][m] for m in range(N + 1)] == M
        F1 = [(A[1][m] + 2 * A[2][m]) / 2 for m in range(N + 1)]
        F2 = [(Fraction(1 if m == 0 else 0) + M[m]) / 2 for m in range(N + 1)]
        C = F1 == F2
        D = (
            ev(F1, Fraction(0)) == 0
            and ev(F1, Fraction(1)) == 1
            and ev(F1, Fraction(1, 2)) == (1 + M_half) / 2
        )
        results.append((tag, layers_ok, B, C, D))

    E = None
    F5 = None
    if doc668 is not None:
        ak = [per_668.get((k, 2), 0) - per_668.get((k, 0), 0) for k in range(N + 1)]
        E = ak == doc668["collapsed_D"] == doc668["committed_bernstein"]
        fn = doc668["five_name_by_rank_black"]
        rp = doc668["rank_pair_totals"]
        F5 = (
            fn["none"][0] == rp["0,2"]
            and fn["both-same"][2] == rp["2,0"]
            and fn["x"][1] + fn["y"][1] + fn["both-same"][1] == rp["1,1"]
            and fn["both-two"] == [0, 0, 0]
        )

    ok = all((a and b and c and d) for _, a, b, c, d in results) and E is not False and F5 is not False
    return ok, results, E, F5, ev(
        [(Fraction(1 if m == 0 else 0) + M[m]) / 2 for m in range(N + 1)],
        Fraction(1, 2),
    )


def main():
    all_ok = True
    for L in (3, 4):
        ok, results, E, F5, fhalf = check(L)
        all_ok &= ok
        status = "OK" if ok else "FAIL"
        print(f"L={L}: {status}  F_L(1/2) = {fhalf}")
        for tag, a, b, c, d in results:
            print(f"  [{tag}] A(layers)={a}  B(M==P2-P0)={b}  C(F==(1+M)/2)={c}  D(boundaries)={d}")
        if E is not None:
            print(f"  [pr668] E(a_k collapse)={E}  F(five-cell P1=x+y+spiral)={F5}")
        assert fhalf == COMMITTED_F_HALF[L], f"F_L(1/2) mismatch at L={L}"
    print("ALL CHECKS EXACT-OK" if all_ok else "FAILURES PRESENT")
    raise SystemExit(0 if all_ok else 1)


if __name__ == "__main__":
    main()

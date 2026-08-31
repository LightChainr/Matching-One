#!/usr/bin/env python3
"""One fixed uniform interval certificate from the locked N25 histograms."""
from collections import defaultdict
from fractions import Fraction as F
from math import comb, prod
from pathlib import Path
import csv
import hashlib
import json
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
N, DELTA = 25, F(1152, 625)
HASHES = ["2d23fecc98d276d9ad15ad1867199cd308f0570cb5040ef94eb6b923b4c53458",
          "225031e612929ed922ba75c55e76703d59990f5283e7ac39b94f022841798da5"]


def add(*polys):
    out = defaultdict(int)
    for p in polys:
        for j, a in p.items():
            out[j] += a
    return {j: a for j, a in out.items() if a}


def neg(p):
    return {j: -a for j, a in p.items()}


def mul(a, b):
    out = defaultdict(int)
    for i, ai in a.items():
        for j, bj in b.items():
            out[i+j] += ai*bj
    return {j: a for j, a in out.items() if a}


def norm(p, x, offset=0):
    assert all(j >= offset for j in p)
    return sum((abs(a)*x**(j-offset) for j, a in p.items()), F(0))


def packet(rows, shift):
    # Exact substitution h=1+x^shift (or h=1 for shift=None), AFTER h differentiation.
    out = {}
    for name, derivative in [("z", 0), ("zh", 1), ("zhh", 2),
                             ("q", 0), ("qh", 1), ("r", 0),
                             ("rh", 1), ("rhh", 2)]:
        p = defaultdict(int)
        for row in rows:
            k, e, q, count = (row[t] for t in ("k", "e", "q", "count"))
            if k < derivative:
                continue
            weight = q if name.startswith("q") else int(q == 0) if name.startswith("r") else 1
            coef = count*weight*prod(range(k-derivative+1, k+1))
            if not coef:
                continue
            for j in range(k-derivative+1) if shift else [0]:
                p[e+(shift or 0)*j] += coef*comb(k-derivative, j)
        out[name] = dict(p)
    return out


def root_poly(a, b):
    return add(mul(a["q"], b["z"]), mul(b["q"], a["z"]))


def angular_poly(a, b):
    pa = add(mul(a["rh"], a["z"]), neg(mul(a["r"], a["zh"])))
    pb = add(mul(b["rh"], b["z"]), neg(mul(b["r"], b["zh"])))
    return add(mul(pb, mul(a["z"], a["z"])),
               neg(mul(pa, mul(b["z"], b["z"]))))


def positive_bound(rows, x, h, derivative=0, rankone=False, excited=False, offset=0):
    total = F(0)
    for row in rows:
        k, e, count = row["k"], row["e"], row["count"]
        if k < derivative or (rankone and row["q"] != 0) or (excited and e == 0):
            continue
        assert e >= offset
        total += count*prod(range(k-derivative+1, k+1))*h**(k-derivative)*x**(e-offset)
    return total


def serial(x):
    if isinstance(x, F):
        return {"exact": str(x), "decimal": float(x)}
    if isinstance(x, dict):
        return {k: serial(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [serial(v) for v in x]
    return x


def certify(pair, law):
    # Rational upper enclosure for lambda^(1/25); no finite-coupling samples.
    x, gap, shift, residual_order, rorder, lead = (
        (F(1, 64), 2, 2, 6, 9, 11) if law == "star"
        else (F(17, 20), 52, None, 52, 210, 210))
    if law == "drop":
        assert x**25 >= F(1, 64)
    h = 1+2*x**gap
    rows_pair = [[dict(row, e=row["g"] if law == "star"
                       else 25*(row["g"]-row["q"]-1)+2*row["k"])
                  for row in rows] for rows in pair]
    at_one = [packet(rows, None) for rows in rows_pair]
    f_one = root_poly(*at_one)
    root_lead = f_one.get(gap, 0)
    root_remainder = norm({j: c for j, c in f_one.items() if j != gap}, x, gap)
    # F_h ground =100 h^49. Bound all excited products in q_a Z_b+q_b Z_a.
    z0, zh0, q0 = 1+h**25, 25*h**24, h**25-1
    ze = [positive_bound(rows, x, h, excited=True) for rows in rows_pair]
    zeh = [positive_bound(rows, x, h, derivative=1, excited=True) for rows in rows_pair]
    fh_error = sum(zh0*ze[j]+zeh[i]*z0+zeh[i]*ze[j]
                   +q0*zeh[j]+ze[i]*zh0+ze[i]*zeh[j]
                   for i, j in [(0, 1), (1, 0)])
    fh_lower = 100-fh_error
    root_numer = abs(root_lead)+root_remainder
    root_gate = root_lead+root_remainder < 0 and fh_lower > 0 and root_numer < 2*fh_lower
    if not root_gate:
        return serial({"certified": False, "obstacle": "implicit-root interval",
                       "root_lead": root_lead, "root_remainder": root_remainder,
                       "Fh_lower": fh_lower, "root_numerator_bound": root_numer})
    at_bar = [packet(rows, shift) for rows in rows_pair]
    residual = norm(root_poly(*at_bar), x, residual_order)
    root_error_coefficient = residual/fh_lower
    angular = angular_poly(*at_bar)
    expected_lead = -200 if law == "star" else 600
    assert angular.get(lead) == expected_lead and min(angular) == lead
    angular_remainder = norm({j: c for j, c in angular.items() if j != lead}, x, lead)
    # |d_h angular_numerator| / x^rorder, throughout the root tube.
    zb = [z0+v for v in ze]
    zhb = [zh0+v for v in zeh]
    zhhb = [positive_bound(rows, x, h, derivative=2) for rows in rows_pair]
    rb = [positive_bound(rows, x, h, rankone=True, offset=rorder) for rows in rows_pair]
    rhb = [positive_bound(rows, x, h, derivative=1, rankone=True, offset=rorder) for rows in rows_pair]
    rhhb = [positive_bound(rows, x, h, derivative=2, rankone=True, offset=rorder) for rows in rows_pair]
    pb = [rhb[i]*zb[i]+rb[i]*zhb[i] for i in range(2)]
    phb = [rhhb[i]*zb[i]+rb[i]*zhhb[i] for i in range(2)]
    nh_bound = sum(phb[i]*zb[j]**2+2*pb[i]*zb[j]*zhb[j] for i, j in [(0, 1), (1, 0)])
    motion_error = nh_bound*root_error_coefficient*x**(rorder+residual_order-lead)
    error = angular_remainder+motion_error
    num_interval = (F(expected_lead)-error, F(expected_lead)+error)
    # M_i=q_ih Z_i-q_i Z_ih, with ground 50 h^24.
    m_error = [2*zh0*ze[i]+2*h**25*zeh[i]+2*ze[i]*zeh[i] for i in range(2)]
    m_lower = [50-v for v in m_error]
    g_lower = 4*sum(m_lower)
    g_upper = sum((50*h**24+m_error[i])*zb[j]**2 for i, j in [(0, 1), (1, 0)])
    denominator_gate = min(m_lower) > 0
    ratios = [2*a/(DELTA*b) for a in num_interval for b in (g_lower, g_upper)]
    u_interval = (min(ratios), max(ratios))
    leading = F(-1, 1)/DELTA if law == "star" else F(3, 1)/DELTA
    remainder_bound = max(abs(v-leading) for v in u_interval)
    sign_gate = u_interval[1] < 0 if law == "star" else u_interval[0] > 0
    return serial({"certified": denominator_gate and sign_gate,
        "coordinate": "lambda=1/m" if law == "star" else "x=m^(-1/25)",
        "coordinate_upper_bound": x, "root_tube": f"1<=h_or_d<=1+2*x^{gap}",
        "root_existence_uniqueness_gate": root_gate, "F_at_one_leading": root_lead,
        "F_at_one_remainder_over_x_gap": root_remainder, "Fh_lower": fh_lower,
        "root_approximation": "1+lambda^2" if shift else "1",
        "root_error_order": residual_order, "root_error_coefficient": root_error_coefficient,
        "angular_leading_coefficient": expected_lead, "angular_power": lead,
        "angular_polynomial_remainder_scaled": angular_remainder,
        "angular_root_motion_remainder_scaled": motion_error,
        "angular_numerator_over_x_power": num_interval,
        "positive_denominator_G": [g_lower, g_upper],
        "U_over_A_over_lambda_power": u_interval,
        "lambda_power": str(F(lead, 1 if law == "star" else 25)),
        "leading_coefficient": leading,
        "uniform_absolute_remainder_coefficient_at_same_power": remainder_bound,
        "remainder_interpretation": "|U/A-leading*lambda^power| <= bound*lambda^power on the entire m>=64 interval"})


def main():
    start = time.perf_counter()
    pair, sources = [], []
    for name, expected in zip(("axis", "tilted"), HASHES):
        path = ROOT/f"results/p337-closed-source-finite-coupling/{name}.csv"
        data = path.read_bytes()
        assert hashlib.sha256(data).hexdigest() == expected
        with path.open() as stream:
            rows = [{k: int(v) for k, v in row.items()} for row in csv.DictReader(stream)]
        pair.append(rows)
        sources.append({"path": str(path.relative_to(ROOT)), "sha256": expected})
    result = {"schema": "p337.uniform-projection-tail.v1", "domain": "all real m>=64 at fixed N25",
              "N": N, "delta": str(DELTA), "laws": {law: certify(pair, law) for law in ("star", "drop")},
              "sources": sources, "new_enumerations": 0, "new_samples": 0,
              "coupling_point_evaluations": 0, "minimum_m_search": False,
              "code_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
              "elapsed_seconds": time.perf_counter()-start}
    out = ROOT/"results/p337-uniform-projection-tail/score.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2)+"\n")
    for law, row in result["laws"].items():
        print(law, row["certified"], row.get("U_over_A_over_lambda_power", row))
    print("seconds", result["elapsed_seconds"])


if __name__ == "__main__":
    main()

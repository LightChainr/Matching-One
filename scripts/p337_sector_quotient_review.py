#!/usr/bin/env python3
"""Exact finite-sector diagnostics; no sampling, optimizer, or network access.

Consumes two immutable N25 histograms. Rational interval bounds certify the
new moment calculations at a sign-bracketed pooled root. This is NOT an
implementation of an importance sampler or a global complexity lower bound.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from dataclasses import dataclass
from decimal import Decimal, localcontext, ROUND_FLOOR, ROUND_CEILING
from fractions import Fraction as F
from pathlib import Path
from typing import Callable, Dict, List, Tuple

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / "experiments/p337-sector-quotient-review-20260831"
N = 25
DELTA = F(1152, 625)
SOURCE_COMMIT = "cae9c8997b5994c218bfe060f75656137f745755"
HASHES = {
    "axis": "078be6411ac481b3381cf0df56901300de88d3b3",
    "tilted": "fb53bbffe3f15f242bfa8d49ffa129a3661db054",
}
Row = Tuple[int, int, int, int]  # K, defect action g, q=r-1, multiplicity


@dataclass(frozen=True)
class Interval:
    lo: F
    hi: F

    def __post_init__(self):
        lo, hi = F(self.lo), F(self.hi)
        if lo > hi:
            raise ValueError("reversed interval")
        # Directed exact dyadic rounding prevents denominator explosion.
        scale = 1 << 512
        low = (lo.numerator * scale) // lo.denominator
        high = -((-hi.numerator * scale) // hi.denominator)
        object.__setattr__(self, "lo", F(low, scale))
        object.__setattr__(self, "hi", F(high, scale))

    @staticmethod
    def point(x):
        return x if isinstance(x, Interval) else Interval(F(x), F(x))

    def __add__(self, other):
        b = self.point(other)
        return Interval(self.lo + b.lo, self.hi + b.hi)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-self.point(other))

    def __rsub__(self, other):
        return self.point(other) - self

    def __mul__(self, other):
        b = self.point(other)
        products = (self.lo*b.lo, self.lo*b.hi, self.hi*b.lo, self.hi*b.hi)
        return Interval(min(products), max(products))

    __rmul__ = __mul__

    def __truediv__(self, other):
        b = self.point(other)
        if b.lo <= 0 <= b.hi:
            raise ZeroDivisionError("interval denominator contains zero")
        return self * Interval(1 / b.hi, 1 / b.lo)

    def __rtruediv__(self, other):
        return self.point(other) / self

    def square(self):
        lo = F(0) if self.lo <= 0 <= self.hi else min(self.lo**2, self.hi**2)
        return Interval(lo, max(self.lo**2, self.hi**2))

    def overlaps(self, other):
        b = self.point(other)
        return max(self.lo, b.lo) <= min(self.hi, b.hi)

    def contains(self, x):
        return self.lo <= x <= self.hi


ZERO = Interval.point(0)
ONE = Interval.point(1)


def load_rows(path: Path, expected_blob: str) -> List[Row]:
    data = path.read_bytes()
    digest = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
    if digest != expected_blob:
        raise ValueError("input Git blob mismatch: " + str(path))
    reader = csv.reader(data.decode("utf-8").splitlines())
    if next(reader) != ["k", "g", "q", "count"]:
        raise ValueError("unexpected CSV schema")
    rows = [tuple(map(int, row)) for row in reader]
    if any(len(r) != 4 or not 0 <= r[0] <= N or r[1] < 0 or
           r[2] not in (-1, 0, 1) or r[3] <= 0 for r in rows):
        raise ValueError("invalid histogram row")
    if len(set(r[:3] for r in rows)) != len(rows):
        raise ValueError("duplicate histogram key")
    for k in range(N + 1):
        if sum(r[3] for r in rows if r[0] == k) != math.comb(N, k):
            raise ValueError("configuration multiplicity mismatch")
    return rows


def homogeneous_eval(coefficients: List[int], h: F) -> int:
    """Return denominator(h)^N * polynomial(h), using integer Horner."""
    if len(coefficients) != N + 1 or h <= 0:
        raise ValueError("require degree <= N and positive h")
    a, b = h.numerator, h.denominator
    value, bpower = coefficients[N], b
    for k in range(N - 1, -1, -1):
        value = value * a + coefficients[k] * bpower
        bpower *= b
    return value


def polynomial_bound(coefficients: List[int], h: Interval) -> Interval:
    # Split positive/negative coefficients: endpoint monotonicity is then exact.
    pos = [max(0, c) for c in coefficients]
    neg = [max(0, -c) for c in coefficients]
    def value(c, x):
        return F(homogeneous_eval(c, x), x.denominator**N)
    return Interval(value(pos, h.lo) - value(neg, h.hi),
                    value(pos, h.hi) - value(neg, h.lo))


class FiniteLaw:
    def __init__(self, rows: List[Row], m: int, drop: bool = False):
        if not isinstance(m, int) or m < 1:
            raise ValueError("m must be a positive integer")
        self.rows, self.m, self.drop = rows, m, drop
        ceiling = max(row[1] for row in rows)
        self.scaled = [(k, g, q, count * m**(ceiling-g+(q+1 if drop else 0)))
                       for k, g, q, count in rows]
        self.total = self.coefficients(lambda k, g, q: 1)
        self.qnumer = self.coefficients(lambda k, g, q: q)

    def coefficients(self, feature: Callable[[int, int, int], int]) -> List[int]:
        result = [0] * (N + 1)
        for k, g, q, weight in self.scaled:
            result[k] += weight * feature(k, g, q)
        return result

    def expectation(self, feature, h: Interval) -> Interval:
        return polynomial_bound(self.coefficients(feature), h) / polynomial_bound(self.total, h)


def bracket_root(laws: List[FiniteLaw], steps: int = 160) -> Interval:
    if len(laws) != 2 or steps < 32:
        raise ValueError("two geometries and at least 32 bisections required")
    def sign(x):
        z = [homogeneous_eval(law.total, x) for law in laws]
        q = [homogeneous_eval(law.qnumer, x) for law in laws]
        return q[0]*z[1] + q[1]*z[0]
    lo, hi = F(1, 10), F(3)
    if not sign(lo) < 0 < sign(hi):
        raise ValueError("pooled root not bracketed")
    for _ in range(steps):
        mid = (lo + hi) / 2
        v = sign(mid)
        if v == 0:
            return Interval(mid, mid)
        if v < 0:
            lo = mid
        else:
            hi = mid
    # For the actual law, uniqueness is inherited from the positive-source
    # root theorem, not inferred from bisection alone. A local slope is checked.
    return Interval(lo, hi)


def moments(law: FiniteLaw, h: Interval, thermal_shift: int = 0) -> Dict:
    E = lambda f: law.expectation(f, h)
    K = E(lambda k,g,q:k)
    K2 = E(lambda k,g,q:k*k)
    S = E(lambda k,g,q:-g + thermal_shift*k)
    KS = E(lambda k,g,q:k*(-g+thermal_shift*k))
    out = {"K": K, "K2": K2, "S": S, "KS": KS}
    for name, obs in (("q", lambda q:q), ("E", lambda q:q*q)):
        o = E(lambda k,g,q:obs(q))
        ok = E(lambda k,g,q:obs(q)*k)
        ok2 = E(lambda k,g,q:obs(q)*k*k)
        os = E(lambda k,g,q:obs(q)*(-g+thermal_shift*k))
        oks = E(lambda k,g,q:obs(q)*k*(-g+thermal_shift*k))
        out[name] = o
        out[name+"x"] = ok-o*K
        out[name+"xx"] = ok2-o*K2-2*K*out[name+"x"]
        out["J"+name] = os-o*S
        out["J"+name+"x"] = oks-os*K-ok*S-o*KS+2*o*K*S
    sectors = []
    for r in range(3):
        def er(f):
            return E(lambda k,g,q: f(k,g,q) if q == r-1 else 0)
        pr = er(lambda k,g,q:1)
        kr, sr = er(lambda k,g,q:k)/pr, er(lambda k,g,q:-g+thermal_shift*k)/pr
        sectors.append({"p":pr, "k":kr, "s":sr,
                        "vk":er(lambda k,g,q:k*k)/pr-kr.square(),
                        "csk":er(lambda k,g,q:k*(-g+thermal_shift*k))/pr-kr*sr})
    out["sectors"] = sectors
    for suffix, key in (("x","k"),("xx","vk"),("t","s"),("tx","csk")):
        out["eta_"+suffix] = (sectors[2][key]-sectors[0][key])/2
        out["xi_"+suffix] = sectors[1][key]-(sectors[0][key]+sectors[2][key])/2
    return out


def comoving_split(laws: List[FiniteLaw], h: Interval, thermal_shift: int = 0) -> Dict:
    f, s = [moments(law, h, thermal_shift) for law in laws]
    D, Qxx = (f["qx"]+s["qx"])/2, (f["qxx"]+s["qxx"])/2
    if D.lo <= 0:
        raise ValueError("positive local pooled slope not certified")
    Yx, Yxx = (f["Ex"]-s["Ex"])/DELTA, (f["Exx"]-s["Exx"])/DELTA
    Jq, Jqx = (f["Jq"]+s["Jq"])/2, (f["Jqx"]+s["Jqx"])/2
    w = Jq/D
    wx = Jqx/D-w*Qxx/D
    for row in (f, s):
        for mode in ("eta", "xi"):
            row[mode+"_fixedQ"] = row[mode+"_t"]-w*row[mode+"_x"]
            row[mode+"_fixedQ_x"] = row[mode+"_tx"]-wx*row[mode+"_x"]-w*row[mode+"_xx"]
        P1 = 1-row["E"]
        row["bias_x"] = (-row["Ex"]*row["q"]+P1*row["qx"])*row["eta_fixedQ"] + P1*row["q"]*row["eta_fixedQ_x"]
        row["middle_x"] = row["Ex"]*(2*row["E"]-1)*row["xi_fixedQ"]-P1*row["E"]*row["xi_fixedQ_x"]
    bias = (f["bias_x"]-s["bias_x"])/(DELTA*D)
    middle = (f["middle_x"]-s["middle_x"])/(DELTA*D)
    direct = ((f["JEx"]-s["JEx"])/DELTA-Yxx*w)/D - Yx*(Jqx-Qxx*w)/D.square()
    if not (bias+middle).overlaps(direct):
        raise AssertionError("independent direct and split formulas disagree")
    return {"h":h,"U_over_A":Yx/D,"V_over_A":direct,"bias_over_A":bias,
            "middle_over_A":middle,"clock_w":w,"geometry":{"axis":f,"tilted":s}}


def variance_budget(laws: List[FiniteLaw], h: Interval) -> Dict:
    geometry = {}
    for name, law in zip(("axis","tilted"), laws):
        E = lambda fn: law.expectation(fn, h)
        mu = E(lambda k,g,q:k)
        qbar = E(lambda k,g,q:q)
        D = E(lambda k,g,q:k*q)-mu*qbar
        gmin = min(g for k,g,q,n in law.rows if q == 0)
        def block(remove_minimal):
            def is_in(g,q):
                return q == 0 and (not remove_minimal or g > gmin)
            P = E(lambda k,g,q:int(is_in(g,q)))
            K = E(lambda k,g,q:k if is_in(g,q) else 0)
            K2 = E(lambda k,g,q:k*k if is_in(g,q) else 0)
            C = K-mu*P
            second = K2-2*mu*K+mu.square()*P
            threshold = mu.lo.numerator // mu.lo.denominator
            if not threshold < mu.lo <= mu.hi < threshold+1:
                raise ValueError("absolute-value sign not certified at integer K")
            plusP = E(lambda k,g,q:int(is_in(g,q) and k > threshold))
            plusK = E(lambda k,g,q:k if is_in(g,q) and k > threshold else 0)
            absolute = 2*(plusK-mu*plusP)-C
            return {"mass":P,"mean":C,"absolute":absolute,
                    "iid_variance":second-C.square(),
                    "one_proposal_floor":absolute.square()-C.square(),
                    "conditional_oracle_variance":P*K2-K.square()}
        full, remainder = block(False), block(True)
        geometry[name] = {"mu":mu,"D":D,"minimal_g":gmin,"full":full,
                          "remainder":remainder,"exact_control_mean":full["mean"]-remainder["mean"]}
    theta = geometry["axis"]["full"]["mean"]-geometry["tilted"]["full"]["mean"]
    if theta.lo <= 0 <= theta.hi:
        raise ValueError("signal sign not certified")
    def ncost(block,key):
        return 9*sum((g[block][key] for g in geometry.values()), ZERO)/theta.square()
    return {"h":h,"theta":theta,"geometry":geometry,
            "iid_n_each":ncost("full","iid_variance"),
            "one_proposal_n_each_floor":ncost("full","one_proposal_floor"),
            "rank1_conditional_n_each_oracle":ncost("full","conditional_oracle_variance"),
            "after_shell_n_each_floor":ncost("remainder","one_proposal_floor"),
            "after_shell_conditional_n_each_oracle":ncost("remainder","conditional_oracle_variance")}


def decimal_string(x: F, precision: int = 35) -> str:
    with localcontext() as ctx:
        ctx.prec = precision
        return str(Decimal(x.numerator)/Decimal(x.denominator))


def encode(value):
    if isinstance(value, Interval):
        with localcontext() as ctx:
            ctx.prec = 40
            ctx.rounding = ROUND_FLOOR
            lo = str(Decimal(value.lo.numerator)/Decimal(value.lo.denominator))
            ctx.rounding = ROUND_CEILING
            hi = str(Decimal(value.hi.numerator)/Decimal(value.hi.denominator))
        if F(lo) > value.lo or F(hi) < value.hi:
            raise AssertionError("decimal endpoints are not outward")
        return {"midpoint":decimal_string((value.lo+value.hi)/2),
                "lower":lo,"upper":hi}
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {k:encode(v) for k,v in value.items()}
    if isinstance(value, (tuple,list)):
        return [encode(v) for v in value]
    return value


def analyze(input_dir: Path) -> Dict:
    rows = [load_rows(input_dir/(name+".csv"),HASHES[name]) for name in ("axis","tilted")]
    laws1 = [FiniteLaw(r,1) for r in rows]
    root1 = bracket_root(laws1)
    split = comoving_split(laws1,root1)
    shifted = comoving_split(laws1,root1,thermal_shift=-1)
    for key in ("V_over_A","bias_over_A","middle_over_A"):
        if not split[key].overlaps(shifted[key]):
            raise AssertionError("thermal gauge check failed")
    laws64 = [FiniteLaw(r,64) for r in rows]
    budget = variance_budget(laws64,bracket_root(laws64))
    # Area factor is irrational; all decisive signs and costs are certified
    # without it. This separate display uses high precision, not an interval.
    with localcontext() as ctx:
        ctx.prec=60
        area=(Decimal(25).ln()*Decimal(13)/Decimal(8)).exp()/2
        displays={key: str(area*Decimal(decimal_string((split[key].lo+split[key].hi)/2,60)))
                  for key in ("U_over_A","V_over_A","bias_over_A","middle_over_A")}
    split["geometry"] = {
        name: {k: row[k] for k in ("q", "E", "qx", "Ex", "eta_x", "xi_x",
                 "eta_fixedQ", "xi_fixedQ", "eta_fixedQ_x", "xi_fixedQ_x")}
        for name, row in split["geometry"].items()}
    return {"schema":"matching-one/p337-sector-quotient-review/v1",
            "source_commit":SOURCE_COMMIT,"input_blobs":HASHES,
            "new_random_samples":0,"official_scorers_run":0,
            "priority":"P0 recommendation; not a production authorization or claim upgrade",
            "claim_boundary":"Finite identities and oracle variance comparison; no practical sampler, continuum field, or independent evidence block.",
            "Q1_comoving_split":encode(split),"area_restored_display_only":displays,
            "m64_star_oracle_budget":encode(budget)}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir",type=Path,default=PACKAGE/"inputs")
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    if args.output.exists():
        parser.error("refusing to overwrite an existing result")
    result=analyze(args.input_dir)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print("Certified finite identities and oracle moment budgets written; no sampler ran.")


if __name__ == "__main__":
    main()

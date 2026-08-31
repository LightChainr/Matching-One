#!/usr/bin/env python3
"""Exact closed Sstar/Bvac -> original-U transmission on the fixed N25 pair.

Compile/enumerate once, then evaluate integer profiles with Fraction intervals.
No Monte Carlo, data-trained coefficients, continuum fit or test campaign.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
from dataclasses import dataclass
from datetime import datetime, timezone
from decimal import Decimal, localcontext
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path
import platform
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis/p337_closed_source_n25_contract.json"
CPP = ROOT / "scripts/p337_closed_source_exact.cpp"
N = 25
DELTA = F(1152, 625)
FIELDS = ("count", "q", "e", "sstar", "bv", "qsstar", "esstar", "qbv", "ebv")


@dataclass(frozen=True)
class Interval:
    lo: F
    hi: F

    @staticmethod
    def of(x):
        return x if isinstance(x, Interval) else Interval(F(x), F(x))

    def __add__(self, x):
        x = self.of(x)
        return Interval(self.lo+x.lo, self.hi+x.hi)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, x):
        return self + -self.of(x)

    def __rsub__(self, x):
        return self.of(x) + -self

    def __mul__(self, x):
        x = self.of(x)
        bounds = (self.lo*x.lo, self.lo*x.hi, self.hi*x.lo, self.hi*x.hi)
        return Interval(min(bounds), max(bounds))

    __rmul__ = __mul__

    def __truediv__(self, x):
        x = self.of(x)
        if x.lo <= 0 <= x.hi:
            raise ArithmeticError("interval denominator contains zero")
        return self * Interval(1/x.hi, 1/x.lo)

    def __pow__(self, exponent):
        if exponent < 0:
            raise ValueError("nonnegative powers only")
        result = self.of(1)
        for _ in range(exponent):
            result = result*self
        return result


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def decimal_fraction(x):
    return Decimal(x.numerator)/Decimal(x.denominator)


def middle(x):
    return (x.lo+x.hi)/2


def interval_json(x):
    # Avoid unwieldy multi-thousand-digit interval endpoints. This is exact
    # outward rational rounding, never a float truncation or narrowed bound.
    scale = 10**40
    lo = F(math.floor(x.lo*scale), scale)
    hi = F(math.ceil(x.hi*scale), scale)
    return {
        "lower_fraction": str(lo), "upper_fraction": str(hi),
        "midpoint_approx": float(middle(x)), "width_approx": float(hi-lo),
        "excludes_zero": lo > 0 or hi < 0,
        "serialization": "exact_outward_rounding_to_1e_minus40_rational_grid",
    }


def profiles(path):
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != N+1:
        raise ValueError("the exact profile must contain all K=0..25")
    data = {f: [] for f in FIELDS}
    for k, row in enumerate(rows):
        if int(row["k"]) != k or int(row["count"]) != math.comb(N, k):
            raise ValueError("incomplete exact configuration population")
        for f in FIELDS:
            data[f].append(int(row[f if f == "count" else "sum_"+f]))
        # A direct count identity, not a second enumeration or statistical test.
        vacant_edges = 2*N*math.comb(N-2, k) if k <= N-2 else 0
        if data["bv"][-1] != vacant_edges:
            raise ValueError("vacant-edge count is inconsistent with the fixed graph")
    return data


def polynomial(coefficients, p):
    return sum(c*p**k*(1-p)**(N-k) for k, c in enumerate(coefficients))


def root_interval(pair, steps):
    coeff = [pair[0]["q"][k]+pair[1]["q"][k] for k in range(N+1)]
    lo, hi = F(11, 20), F(13, 20)
    if not polynomial(coeff, lo) < 0 < polynomial(coeff, hi):
        raise ValueError("predeclared matching root bracket fails")
    for _ in range(steps):
        mid = (lo+hi)/2
        val = polynomial(coeff, mid)
        if val == 0:
            return Interval(mid, mid)
        if val < 0:
            lo = mid
        else:
            hi = mid
    return Interval(lo, hi)


def moments(data, p):
    # Raw enumerated sums already include multiplicity: no extra binomial factor.
    pp = [p**k for k in range(N+1)]
    zz = [(1-p)**k for k in range(N+1)]
    result = {f: [Interval.of(0) for _ in range(3)] for f in FIELDS}
    for k in range(N+1):
        m = N-k
        w = pp[k]*zz[m]
        wp = (k*pp[k-1]*zz[m] if k else 0) - (m*pp[k]*zz[m-1] if m else 0)
        wpp = ((k*(k-1)*pp[k-2]*zz[m] if k >= 2 else 0)
               - (2*k*m*pp[k-1]*zz[m-1] if k and m else 0)
               + (m*(m-1)*pp[k]*zz[m-2] if m >= 2 else 0))
        for f in FIELDS:
            for order, weight in enumerate((w, wp, wpp)):
                result[f][order] += data[f][k]*weight
    return result


def score(pair, root):
    packet = [moments(data, root) for data in pair]
    D = (packet[0]["q"][1]+packet[1]["q"][1])/2
    B = (packet[0]["e"][1]-packet[1]["e"][1])/DELTA
    T = (packet[0]["q"][2]+packet[1]["q"][2])/2
    H = (packet[0]["e"][2]-packet[1]["e"][2])/DELTA
    if D.lo <= 0:
        raise ArithmeticError("nonpositive matching slope enclosure")
    sources = {}
    for source in ("sstar", "bv"):
        jq, jqp, jep = [], [], []
        for row in packet:
            q, qp, _ = row["q"]
            e, ep, _ = row["e"]
            s, sp, _ = row[source]
            qs, qsp, _ = row["q"+source]
            _, esp, _ = row["e"+source]
            jq.append(qs-q*s)
            jqp.append(qsp-qp*s-q*sp)
            jep.append(esp-ep*s-e*sp)
        jQ, jQp = sum(jq)/2, sum(jqp)/2
        jYp = (jep[0]-jep[1])/DELTA
        terms = {
            "direct": jYp/D,
            "root_motion": -H*jQ/(D**2),
            "slope_source": -B*jQp/(D**2),
            "slope_root": B*T*jQ/(D**3),
        }
        sources[source] = {"reduced_V": sum(terms.values()), "terms": terms, "root_tangent": -jQ/D}
    return {"reduced_U": B/D, "D": D, "sources": sources}


def git_value(spec):
    return subprocess.check_output(["git", "rev-parse", spec], cwd=ROOT, text=True).strip()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--counts-dir", type=Path, help="read saved Sstar/Bvac profiles without enumeration")
    args = parser.parse_args()
    started = time.perf_counter()
    start_utc = datetime.now(timezone.utc).isoformat()
    contract = json.loads(CONTRACT.read_text())
    if contract["N"] != N or contract["geometries"] != [[5, 0], [4, 3]]:
        raise ValueError("this scorer is restricted to the frozen N25 pair")
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    receipts, compile_command, compiler_version = [], None, None
    if args.counts_dir:
        paths = [args.counts_dir.resolve()/name for name in ("axis.csv", "tilted.csv")]
    else:
        paths = [out/name for name in ("axis.csv", "tilted.csv")]
        build = Path(tempfile.mkdtemp(prefix="p337-closed-source-n25-"))
        binary = build/"enumerate"
        compiler_version = subprocess.check_output(["/usr/bin/clang++", "--version"], text=True).splitlines()[0]
        compile_command = ["/usr/bin/clang++", "-O3", "-std=c++17", str(CPP), "-o", str(binary)]
        subprocess.run(compile_command, check=True)

        def enumerate_one(design):
            (a, b), path = design
            command = [str(binary), str(a), str(b), str(path)]
            completed = subprocess.run(command, check=True, capture_output=True, text=True)
            receipt = json.loads(completed.stdout)
            receipt.update(geometry=[a, b], command=command, binary_sha256=sha(binary),
                           exit_code=completed.returncode)
            print(json.dumps(receipt), flush=True)
            return receipt

        with ThreadPoolExecutor(max_workers=contract["workers"]) as pool:
            receipts = list(pool.map(enumerate_one, zip(contract["geometries"], paths)))
        (out/"enumeration.json").write_text(json.dumps(receipts, indent=2)+"\n")
    pair = [profiles(path) for path in paths]
    root = root_interval(pair, contract["root_bisection_steps"])
    scored = score(pair, root)
    with localcontext() as context:
        context.prec = 70
        A = Decimal(N)**(Decimal(13)/Decimal(8))/2
        decimation_factor = Decimal(2)**(Decimal(13)/Decimal(8))
        values = {s: A*decimal_fraction(middle(v["reduced_V"]))
                  for s, v in scored["sources"].items()}
        approximate = {
            "root": float(middle(root)),
            "U25": float(A*decimal_fraction(middle(scored["reduced_U"]))),
            "V25_Sstar": float(values["sstar"]), "V25_Bvac": float(values["bv"]),
            "V50_endpoint_F4_from_Bvac": float(decimation_factor*values["bv"]),
        }
    excludes_zero = interval_json(scored["sources"]["sstar"]["reduced_V"])["excludes_zero"]
    result = {
        "schema": "matching-one.p337-closed-source-n25.score.v1",
        "status": "completed_exact_finite_enumeration", "contract": contract,
        "contract_freeze_commit": "5598812612d176e30c0e9ee50d2fd78f382db632",
        "code_commit": git_value("HEAD"),
        "parent_provenance": {
            "enumerator_commit": git_value("2bfe9b90"),
            "enumerator_blob": git_value("2bfe9b90:scripts/exact_decimation_plaquette_u.cpp"),
            "scorer_commit": git_value("c76b038b"),
            "scorer_blob": git_value("c76b038b:scripts/analyze_decimation_plaquette_u.py"),
            "changes": "identical quotient/traversal; leaf statistics changed to Sstar/Bvac; corrected rational interval scorer reused",
        },
        "input_counts": [{"name": p.name, "sha256": sha(p), "configurations": sum(data["count"])}
                         for p, data in zip(paths, pair)],
        "root_enclosure": interval_json(root),
        "matching_slope_enclosure": interval_json(scored["D"]),
        "U25_over_A_enclosure": interval_json(scored["reduced_U"]),
        "source_enclosures": {
            s: {"V_over_A": interval_json(v["reduced_V"]),
                "root_tangent": interval_json(v["root_tangent"]),
                "terms_over_A": {k: interval_json(w) for k, w in v["terms"].items()}}
            for s, v in scored["sources"].items()},
        "numerical_values": approximate,
        "decision": "closed_Sstar_common_thermal_alias_rejected" if excludes_zero else
                    "not_resolved_by_frozen_exact_root_enclosure",
        "interval_scope": "exact rational enclosures conditional on graph/counting definitions; positive irrational area factor displayed numerically",
        "endpoint_scope": "V50_endpoint_F4=2^(13/8)*V25_Bvac under the stated decimation dictionary; no parent simulation",
        "sampling_error": "none", "new_random_samples": 0, "old_C_or_F4_rescored": False,
    }
    (out/"latest.json").write_text(json.dumps(result, indent=2)+"\n")
    report = f"""# Exact closed source Sstar reaches the N25 original global U

**V_Sstar={approximate['V25_Sstar']:.14g}; V_Bvac={approximate['V25_Bvac']:.14g}.**
The exact rational enclosure of V_Sstar/A {'excludes' if excludes_zero else 'contains'} zero.
Decision: `{result['decision']}`.

## Fixed mechanism, not a source scan

Sstar=Ctot+F4+Bvac is selected by T(Ctot)=Ctot+F4, T(F4)=Bvac, T(Bvac)=0.
It is fixed under this map, and bare Ctot reaches it after two steps.
The only companion source is Bvac, the next endpoint of F4. Both are bulk
sources exp(t*S); no amplitude or density normalization was fitted.

One full pass per geometry enumerated all 2^25 configurations of (5,0) and
(4,3). The inherited geometry, component updates and traversal were unchanged.
At each leaf Bvac=50−4K+occupied_NN_edges and Sstar=Ctot+F4+Bvac.
Only their sufficient statistics plus q/E were saved; old C/F4 source
responses were not recomputed. There were no random samples or cloud jobs.

## Root-complete finite response

| Quantity | Numerical evaluation of exact coefficients |
|---|---:|
| Fresh exact pooled root | {approximate['root']:.16g} |
| Native U25 | {approximate['U25']:.14g} |
| V25 for Sstar | {approximate['V25_Sstar']:.14g} |
| V25 for Bvac | {approximate['V25_Bvac']:.14g} |
| N50 F4 endpoint predicted from Bvac | {approximate['V50_endpoint_F4_from_Bvac']:.14g} |

Let Q=mean(q), Y=P4(E), D=Q_p, A=25^(13/8)/2 and j_O=Cov(O,S).
Every source uses the exact original-U derivative

`V_S/A = jY_p/D − Y_pp*jQ/D² − Y_p*jQ_p/D² + Y_p*Q_pp*jQ/D³`.

Root motion, slope motion and per-geometry covariance centering are all included.
The exact integer counts contain their configuration multiplicities; no extra
binomial factor is applied. The root uses 128 rational bisections in [11/20,13/20].
`latest.json` preserves outward rational enclosures for the root, V/A, and each
of its four terms. These are computational bounds, not confidence intervals.

## What changes

A nonzero Sstar response excludes common-thermal invisibility of this uniquely
closed source **for this finite global observable**. Closure of the source
under decimation does not force it to disappear after root normalization.
The companion predicts the next F4 endpoint through
`V50_endpoint_F4=2^(13/8)*V25_Bvac`; this is a transported derivative, not
an independently simulated parent endpoint.

The axis Z5xZ5 and tilted Z25 quotients have different Smith classes. Nothing
here identifies a continuum field or an asymptotic exponent. The calculation
does not change the independent larger-N F4 experiment or revive a lag-one
source. If the reported enclosure contains zero, the finite alias is unresolved.

## Provenance and reproduction

Contract freeze: `5598812612d176e30c0e9ee50d2fd78f382db632`.
Enumerator parent: `2bfe9b90:scripts/exact_decimation_plaquette_u.cpp`.
Corrected scorer parent: `c76b038b:scripts/analyze_decimation_plaquette_u.py`.
Full commit/blob IDs, source hashes and enumeration receipts are saved alongside
the exact profiles. Run once with
`python scripts/p337_closed_source_score.py --output-dir NEW_DIRECTORY`.
For future inspection only, `--counts-dir results/p337-closed-source-n25`
consumes these saved new-source profiles without enumerating again.
"""
    (out/"REPORT.md").write_text(report)
    receipt = {
        "started_utc": start_utc, "finished_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": time.perf_counter()-started, "command": sys.argv,
        "python": sys.version, "machine": platform.machine(), "compiler": compiler_version,
        "compile_command": compile_command, "enumeration_receipts": receipts,
        "hashes": {"contract": sha(CONTRACT), "cpp": sha(CPP), "python": sha(__file__),
                   "result": sha(out/"latest.json")},
        "new_random_samples": 0, "tests_run": 0, "cloud_jobs": 0,
    }
    (out/"run.json").write_text(json.dumps(receipt, indent=2)+"\n")
    print(json.dumps({"decision": result["decision"], "values": approximate,
                      "elapsed_seconds": receipt["elapsed_seconds"]}), flush=True)


if __name__ == "__main__":
    main()

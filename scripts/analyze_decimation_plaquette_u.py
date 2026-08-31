#!/usr/bin/env python3
"""Exact finite F4 -> original-U transmission on the fixed Gaussian N25 pair.

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
CONTRACT = ROOT / "analysis/decimation_plaquette_u_contract.json"
CPP = ROOT / "scripts/exact_decimation_plaquette_u.cpp"
N = 25
DELTA = F(1152, 625)
FIELDS = ("count", "q", "e", "c", "f4", "qc", "ec", "qf4", "ef4")


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
        faces = N*math.comb(N-4, k-4) if k >= 4 else 0
        if data["f4"][-1] != faces:
            raise ValueError("unit-face count is inconsistent with the fixed graph")
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
    for source in ("c", "f4"):
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--counts-dir", type=Path, help="reuse the two exact integer CSVs, without enumeration")
    args = parser.parse_args()
    started = time.perf_counter()
    start_utc = datetime.now(timezone.utc).isoformat()
    contract = json.loads(CONTRACT.read_text())
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    receipts = []
    compile_command = None
    compiler_version = None
    if args.counts_dir:
        paths = [args.counts_dir.resolve()/name for name in ("axis.csv", "tilted.csv")]
    else:
        paths = [out/name for name in ("axis.csv", "tilted.csv")]
        build = Path(tempfile.mkdtemp(prefix="decimation-plaquette-u-"))
        binary = build/"enumerate"
        compiler_version = subprocess.check_output(["/usr/bin/clang++", "--version"], text=True).splitlines()[0]
        compile_command = ["/usr/bin/clang++", "-O3", "-std=c++17", str(CPP), "-o", str(binary)]
        subprocess.run(compile_command, check=True)

        def enumerate_one(design):
            (a, b), path = design
            command = [str(binary), str(a), str(b), str(path)]
            completed = subprocess.run(command, check=True, capture_output=True, text=True)
            receipt = json.loads(completed.stdout)
            receipt.update(geometry=[a, b], command=command, binary_sha256=sha(binary), exit_code=completed.returncode)
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
        factor = Decimal(2)**(Decimal(13)/Decimal(8))
        u = A*decimal_fraction(middle(scored["reduced_U"]))
        values = {s: A*decimal_fraction(middle(v["reduced_V"])) for s, v in scored["sources"].items()}
        correction, bare = factor*values["f4"], factor*values["c"]
        complete = bare+correction
        approximate = {"root": float(middle(root)), "U25": float(u), "V25_cluster": float(values["c"]),
                       "V25_F4": float(values["f4"]), "V50_endpoint_cluster_bare": float(bare),
                       "V50_endpoint_cluster_complete": float(complete), "V50_endpoint_forced_correction": float(correction)}
    primary = scored["sources"]["f4"]["reduced_V"]
    excludes_zero = interval_json(primary)["excludes_zero"]
    result = {
        "schema": "matching-one.decimation-plaquette-u.v1",
        "status": "completed_exact_finite_enumeration",
        "contract": contract,
        "code_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "input_counts": [{"name": p.name, "sha256": sha(p), "configurations": sum(data["count"])} for p, data in zip(paths, pair)],
        "root_enclosure": interval_json(root),
        "matching_slope_enclosure": interval_json(scored["D"]),
        "U25_over_A_enclosure": interval_json(scored["reduced_U"]),
        "source_enclosures": {s: {"V_over_A": interval_json(v["reduced_V"]),
                                  "root_tangent": interval_json(v["root_tangent"]),
                                  "terms_over_A": {k: interval_json(w) for k, w in v["terms"].items()}}
                              for s, v in scored["sources"].items()},
        "numerical_values": approximate,
        "decision": "F4_thermal_alias_and_bare_cluster_endpoint_U_transport_rejected" if excludes_zero else "not_resolved_by_current_exact_root_enclosure",
        "interval_scope": "exact rational propagation from exhaustive integer coefficients and the physical root bracket; irrational positive area prefactors shown numerically, not Monte Carlo confidence intervals",
        "uncertainty": "no_sampling_error; computational conclusion conditional on stated graph/counting definitions",
        "new_random_samples": 0,
    }
    (out/"latest.json").write_text(json.dumps(result, indent=2)+"\n")
    report = f"""# The decimation-forced plaquette source reaches the original global U

## Finite mechanism result

The fixed N25 Gaussian pair gives **V_F4={approximate['V25_F4']:.12g}** in bulk source units.
The exact rational enclosure of V_F4/A {'excludes' if excludes_zero else 'contains'} zero.
Decision: `{result['decision']}`. This compares the forced decimation operator with
a thermal-only alias in the same original U, not a new fit or a lag1 rescue.

## Source and observer held fixed

All 2^25 configurations of each quotient (5,0) and (4,3) were enumerated: 67,108,864
configurations in total. Sources are Ctot=CB_NN+CW_matching and the number F4 of
fully occupied unit faces, with normalized weights exp(t*S) and no density factor.
q=CB-CW-(K-T_NN+F4), E=q^2 follow digital Alexander on these honest unit-cell tori.
P4 uses the exact direction difference 1152/625. The root is the new pair's pooled
matching root, not an old production calibration; A=25^(13/8)/2.

## Exact coefficients give a visible endpoint correction

| Quantity | Value (numerical evaluation of exact coefficients) |
|---|---:|
| Pooled root | {approximate['root']:.15g} |
| Native U25 | {approximate['U25']:.12g} |
| V25 for bulk cluster source | {approximate['V25_cluster']:.12g} |
| V25 for bulk F4 source | {approximate['V25_F4']:.12g} |
| Bare-cluster prediction for N50 endpoint | {approximate['V50_endpoint_cluster_bare']:.12g} |
| Complete forced-source prediction for N50 endpoint | {approximate['V50_endpoint_cluster_complete']:.12g} |
| Missing endpoint correction | {approximate['V50_endpoint_forced_correction']:.12g} |

The endpoint dictionary is C_parent=C_child+F4. Consequently the difference between
complete and bare predictions is exactly 2^(13/8)*V25_F4, including root and slope
motion. Parent generators are (5,5),(1,7); the complement sign and period rotation
cancel in U. These are theorem-transported endpoint derivatives, not parent
endpoint simulations. A nonzero F4 term makes the configuration-level failure of
bare source closure visible to the specified global observable.

## Calculation includes normalized covariance and the moving root

For each geometry j_q=Cov(q,S), j_E=Cov(E,S). Define Q=mean(q), Y=P4(E),
D=Q_p, r=Y_p/D. The evaluated response is

`V_S/A = jY_p/D - Y_pp*jQ/D^2 - Y_p*jQ_p/D^2 + Y_p*Q_pp*jQ/D^3`.

The four terms, source root shifts and full rational bounds are retained in
`latest.json`. Per-K integer sums already include binomial multiplicity.
All thermal derivatives include the derivatives of covariance centering.
The physical root is bracketed by exact rational bisection; rational interval
arithmetic then encloses the reduced response. Positive irrational area factors
are applied only for numerical presentation and do not change zero exclusion.

## Scope and uncertainty

No sampling errors or confidence levels apply to this exhaustive calculation.
The rational bounds propagate the root interval conditional on the supplied exact
graph counts. The two quotients have different Smith classes (Z5xZ5 and Z25);
the result is a finite-pair mechanism counterexample, not an asymptotic H4 law,
the N65/N85 production family, a continuum field identity or interior saturation
curvature. The prior P154/P334 stop decisions remain unchanged.

## Subsequent exact source completion and next question

The accompanying `notes/decimation-closed-source-and-global-u.md` now completes
the dictionary: F_parent=T_child-4K_child+2M and T_parent=4M-4K_child. Thus
S_hat=C+F+T-4K+2N is exactly unchanged by this endpoint decimation. This forced
finite source family closes without a fitted correction or another descriptor.
The note contains the proof; this numerical calculation establishes its first
otherwise-missing F4 contribution to the specified global U.

An interior transmission law for this same closed source remains open. The
endpoint identity alone is not that law, and no V_T value or interior curve is
claimed here. Repeating this F4 calculation or reopening a failed P154/P334
parameterization is not the next target. No new production block is launched.

## Source and reproduction

The dictionary is pinned at execution commit 207436518db46dd13ef0ec91168cb1c99d52eaea,
`notes/p337-checkerboard-decimation-global-u.md`; the topology proof is
56838d5f068f6f0ba7795926dc9343229bdd28ce, `notes/square-checkerboard-endpoint-homology.md`.
The contract and both scripts are pinned by `code_commit` in `latest.json` and
hashes in `run.json`. The stored score predates the subsequent source-completion
note; that explanatory addition does not change its counts or numerical result.
The table is used for exact value lookup rather than an inferred trend.

Run `python scripts/analyze_decimation_plaquette_u.py --output-dir NEW_DIRECTORY`.
Add `--counts-dir results/decimation-plaquette-u` to consume the saved integer
profiles without enumerating configurations again. Existing outputs are not overwritten.
"""
    (out/"REPORT.md").write_text(report)
    receipt = {"started_utc": start_utc, "finished_utc": datetime.now(timezone.utc).isoformat(),
               "elapsed_seconds": time.perf_counter()-started, "command": sys.argv,
               "python": sys.version, "machine": platform.machine(), "compiler": compiler_version,
               "compile_command": compile_command, "enumeration_receipts": receipts,
               "hashes": {"contract": sha(CONTRACT), "cpp": sha(CPP), "python": sha(__file__), "result": sha(out/"latest.json")},
               "new_random_samples": 0, "tests_run": 0, "cloud_jobs": 0}
    (out/"run.json").write_text(json.dumps(receipt, indent=2)+"\n")
    print(json.dumps({"decision": result["decision"], "values": approximate, "elapsed_seconds": receipt["elapsed_seconds"]}), flush=True)


if __name__ == "__main__":
    main()

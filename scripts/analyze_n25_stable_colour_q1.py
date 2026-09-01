#!/usr/bin/env python3
"""The fixed N25 stable-colour Q1 transmission and its full log-Q jet.

Uses saved seam counts and a saved Q1 root. Exact rational intervals carry
the root enclosure; first-order dual numbers differentiate log Q. Thermal
jets through order3 retain mixed root/slope motion without another solve.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import datetime, timezone
from decimal import Decimal, localcontext
from fractions import Fraction as F
import hashlib
import io
import json
import math
from pathlib import Path
import platform
import subprocess
import sys
import time

from analyze_decimation_plaquette_u import Interval as I, interval_json, middle

ROOT = Path(__file__).resolve().parents[1]
N, DELTA = 25, F(1152, 625)


@dataclass(frozen=True)
class Dual:
    value: I
    dot: I

    @staticmethod
    def of(x):
        return x if isinstance(x, Dual) else Dual(I.of(x), I.of(0))

    def __add__(self, x):
        x = self.of(x)
        return Dual(self.value + x.value, self.dot + x.dot)

    __radd__ = __add__

    def __neg__(self):
        return Dual(-self.value, -self.dot)

    def __sub__(self, x):
        return self + -self.of(x)

    def __rsub__(self, x):
        return self.of(x) + -self

    def __mul__(self, x):
        x = self.of(x)
        return Dual(self.value * x.value, self.dot * x.value + self.value * x.dot)

    __rmul__ = __mul__

    def __truediv__(self, x):
        x = self.of(x)
        return Dual(self.value / x.value,
                    (self.dot * x.value - self.value * x.dot) / (x.value**2))

    def __pow__(self, n):
        out = Dual.of(1)
        for _ in range(n):
            out = out * self
        return out


def raw_jet(coefficients, h):
    result = []
    for d in range(4):
        value = Dual.of(0)
        for k in range(N, d - 1, -1):
            value = value * h + coefficients[k] * math.comb(k, d) * math.factorial(d)
        result.append(value)
    return result


def normalized_jet(raw, z):
    result = []
    for n, value in enumerate(raw):
        for k in range(1, n + 1):
            value = value - math.comb(n, k) * z[k] * result[n-k]
        result.append(value / z[0])
    return result


def geometry(rows, h):
    names = ("z", "q", "e", "f12", "f21", "explicit")
    coeffs = {key: [Dual.of(0) for _ in range(N+1)] for key in names}
    accepted_rank1 = {(0, 0), (1, 1), (1, 2), (0, 1)}
    for row in rows:
        k, g, q, bad2, n3, count = (int(row[key]) for key in ("k", "g", "q", "bad2", "n_bad3", "count"))
        if q == 0 and (bad2, n3) not in accepted_rank1:
            raise ValueError("count packet is outside the declared N25 packing completion")
        weight = Dual(I.of(count), I.of(F(-g*count, 2)))
        for key, multiplier in (("z", 1), ("q", q), ("e", q*q)):
            coeffs[key][k] = coeffs[key][k] + weight * multiplier
        i12 = q == 0 and (bad2, n3) == (1, 2)
        i21 = q == 0 and (bad2, n3) == (0, 1)
        for key, selected, beta_dot in (("f12", i12, F(3, 2)), ("f21", i21, F(1, 2))):
            if selected:
                coeffs[key][k] = coeffs[key][k] + weight * Dual(I.of(-1), I.of(beta_dot))
                coeffs["explicit"][k] = coeffs["explicit"][k] + weight * beta_dot
    z = raw_jet(coeffs["z"], h)
    packet = {name: normalized_jet(raw_jet(coeffs[name], h), z) for name in names if name != "z"}
    packet["f"] = [a+b for a, b in zip(packet["f12"], packet["f21"])]
    return packet, coeffs


def comoving(packet, root_logQ):
    return {key: [Dual(jet[n].value, jet[n].dot + jet[n+1].value * root_logQ)
                  for n in range(3)] for key, jet in packet.items()}


def response(pair, field):
    q1, q2 = [p["q"] for p in pair]
    e1, e2 = [p["e"] for p in pair]
    f1, f2 = [p[field] for p in pair]
    D = (q1[1]+q2[1])/2
    Mhh = (q1[2]+q2[2])/2
    Yh, Yhh = (e1[1]-e2[1])/DELTA, (e1[2]-e2[2])/DELTA
    deltaM = -(q1[0]*f1[0]+q2[0]*f2[0])/2
    deltaMh = -(q1[1]*f1[0]+q1[0]*f1[1]+q2[1]*f2[0]+q2[0]*f2[1])/2
    deltaYh = -(e1[1]*f1[0]+e1[0]*f1[1]-e2[1]*f2[0]-e2[0]*f2[1])/DELTA
    terms = {"direct": deltaYh/D, "root_motion": -Yhh*deltaM/(D**2),
             "slope_source": -Yh*deltaMh/(D**2), "slope_root": Yh*Mhh*deltaM/(D**3)}
    return sum(terms.values()), terms


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, default=ROOT/"analysis/n25_stable_colour_q1_contract.json")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    contract_bytes = args.contract.read_bytes()
    contract = json.loads(contract_bytes)
    source_files = []

    def source(commit, path):
        data = subprocess.check_output(["git", "show", commit+":"+path], cwd=ROOT)
        source_files.append({"commit": commit, "path": path, "sha256": hashlib.sha256(data).hexdigest()})
        return data

    saved = json.loads(source(contract["root"]["commit"], contract["root"]["path"]))
    rp = saved["root_enclosure"]
    p = I(F(rp["lower_fraction"]), F(rp["upper_fraction"]))
    h = p/(1-p)
    packet, coefficient_packets = [], []
    for name in ("axis", "tilted"):
        data = source(contract["counts"]["commit"], contract["counts"]["directory"]+"/"+name+".csv")
        rows = list(csv.DictReader(io.StringIO(data.decode())))
        one, coeffs = geometry(rows, h)
        packet.append(one)
        coefficient_packets.append({"geometry": name, "coefficients_h": {
            key: [{"Q1": str(x.value.lo), "logQ_jet": str(x.dot.lo)} for x in values]
            for key, values in coeffs.items()}})
    D = (packet[0]["q"][1].value+packet[1]["q"][1].value)/2
    Mdot = (packet[0]["q"][0].dot+packet[1]["q"][0].dot)/2
    root_logQ = -Mdot/D
    moved = [comoving(p, root_logQ) for p in packet]
    B, terms = response(moved, "f")
    B12, _ = response(moved, "f12")
    B21, _ = response(moved, "f21")
    explicit, _ = response(moved, "explicit")
    quantities = {"B1": B.value, "B1_logQ": B.dot,
                  "B1_I12": B12.value, "B1_I21": B21.value,
                  "B1_logQ_explicit_beta": explicit.value,
                  "B1_logQ_measure_root_slope": B.dot-explicit.value}
    with localcontext() as ctx:
        ctx.prec = 65
        A = Decimal(N)**(Decimal(13)/8)/2
        values = {key: float(A*Decimal(middle(x).numerator)/Decimal(middle(x).denominator))
                  for key, x in quantities.items()}
    bounds = {key: interval_json(x) for key, x in quantities.items()}
    decisions = {key: "zero_rejected" if bounds[key]["excludes_zero"] else "not_resolved_by_saved_root_enclosure"
                 for key in ("B1", "B1_logQ")}
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    result = {"schema": "matching-one.n25-stable-colour-q1.result.v1", "status": "completed_saved_count_exact_reduction",
              "definition_commit": commit, "contract": contract, "decisions": decisions,
              "values": values, "enclosures_over_A": bounds,
              "four_terms_over_A": {key: {"B1": interval_json(x.value), "B1_logQ": interval_json(x.dot)} for key, x in terms.items()},
              "saved_p_root": rp, "h_root": interval_json(h), "h_root_logQ": interval_json(root_logQ),
              "geometry_trace": [{key: {"fraction_Q1": interval_json(p[key][0].value),
                                         "fixed_h_logQ": interval_json(p[key][0].dot),
                                         "fixed_Q_h_derivative": interval_json(p[key][1].value)}
                                  for key in ("f", "f12", "f21")} for p in packet],
              "source_files": source_files, "coefficient_packets": coefficient_packets,
              "new_lattice_enumerations": 0, "new_random_samples": 0, "new_root_searches": 0,
              "Q4_score_rerun": False, "tests_run": 0, "cloud_jobs": 0,
              "elapsed_seconds": time.perf_counter()-started,
              "boundary": contract["boundary"]}
    report = ["# The stable N25 colour trace at Q1", "",
              "One specified generic-Q full-central completion, scored at the saved original Q1 pooled root.",
              "The annular-packing proof makes the existing seam counts sufficient; no new occupation pass was needed.", "",
              "| Complete original-U quantity | Value |", "|---|---:|"]
    for key, value in values.items():
        report.append(f"| {key} | {value:+.16g} |")
    report += ["", "Primary B1 is the epsilon response, not the total logQ derivative of U.",
               "B1_logQ differentiates that response along the original Q family, including the explicit beta(Q), measure, root and thermal slope.",
               f"Zero decisions: B1={decisions['B1']}; B1_logQ={decisions['B1_logQ']}.",
               "Exact rational enclosures, four-term decompositions, full coefficient jets and input hashes are in [latest.json](latest.json).", "",
               "This finite colour-sector completion is not a unique local four-leg field or a continuum exponent measurement.",
               "The count population is shared with the completed Q4 calculation; it adds no independent statistical vote.",
               "Proof: [N25 stable completion](../../notes/n25-stable-colour-completion.md).", ""]
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    content = {"latest.json": json.dumps(result, indent=2)+"\n", "REPORT.md": "\n".join(report)}
    for name, data in content.items():
        (out/name).write_text(data)
    receipt = {"schema": "matching-one.n25-stable-colour-q1.run.v1", "definition_commit": commit,
               "command": sys.argv, "python_executable": sys.executable, "python": sys.version,
               "machine": platform.machine(), "created_utc": datetime.now(timezone.utc).isoformat(),
               "elapsed_seconds": result["elapsed_seconds"], "source_files": source_files,
               "contract_sha256": hashlib.sha256(contract_bytes).hexdigest(),
               "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               "interval_helper_sha256": hashlib.sha256((ROOT/"scripts/analyze_decimation_plaquette_u.py").read_bytes()).hexdigest(),
               "output_sha256": {name: hashlib.sha256(data.encode()).hexdigest() for name, data in content.items()},
               "new_enumerations": 0, "new_random_samples": 0, "root_searches": 0, "tests_run": 0, "cloud_jobs": 0}
    (out/"run.json").write_text(json.dumps(receipt, indent=2)+"\n")
    print(json.dumps({"decisions": decisions, "values": values, "elapsed_seconds": result["elapsed_seconds"]}))


if __name__ == "__main__":
    main()

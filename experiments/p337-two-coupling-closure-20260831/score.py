#!/usr/bin/env python3
"""Fixed exact four-profile tangent test; reads saved counts, never enumerates.

Run only after CONTRACT.md, inputs and this implementation are committed.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from fractions import Fraction as F
import hashlib
from itertools import combinations
import json
import math
from pathlib import Path
import platform
import re
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
N = 25
FIELDS = ("count", "q", "e", "sstar", "qsstar", "esstar")
ROWS = ("q_first", "E_first", "q_second", "E_second")


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_table(path):
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) != N + 1:
        raise ValueError(f"incomplete K profile: {path}")
    data = {field: [] for field in FIELDS}
    for k, row in enumerate(rows):
        if int(row["k"]) != k or int(row["count"]) != math.comb(N, k):
            raise ValueError(f"K/count mismatch: {path}, K={k}")
        for field in FIELDS:
            data[field].append(int(row[field if field == "count" else "sum_" + field]))
    return data


def parent_endpoint(child):
    # Complement child occupancy k -> 25-k; q and q*S change sign.
    return {field: [(-1 if field in ("q", "qsstar") else 1) * value
                    for value in reversed(values)] for field, values in child.items()}


def polynomial(coefficients, p):
    # Coefficients are population sums, already including binomial multiplicity.
    return sum(value * p**k * (1-p)**(N-k) for k, value in enumerate(coefficients))


def moments(data, p):
    pp = [p**k for k in range(N + 1)]
    zz = [(1-p)**k for k in range(N + 1)]
    result = {field: [Interval.of(0), Interval.of(0)] for field in FIELDS}
    for k in range(N + 1):
        m = N-k
        w = pp[k] * zz[m]
        wp = (k*pp[k-1]*zz[m] if k else 0) - (m*pp[k]*zz[m-1] if m else 0)
        for field in FIELDS:
            result[field][0] += data[field][k] * w
            result[field][1] += data[field][k] * wp
    # The verified count polynomial is identically (p+(1-p))**25 = 1.
    # Source normalization is therefore the per-geometry covariance below.
    return result


def det3(matrix):
    a, b, c = matrix
    return (a[0]*(b[1]*c[2]-b[2]*c[1])
            - a[1]*(b[0]*c[2]-b[2]*c[0])
            + a[2]*(b[0]*c[1]-b[1]*c[0]))


def square(x):
    # Generic interval multiplication loses the nonnegative lower bound at zero.
    return Interval(F(0) if x.lo <= 0 <= x.hi else min(x.lo*x.lo, x.hi*x.hi),
                    max(x.lo*x.lo, x.hi*x.hi))


def positive_witness(x):
    """A short exact positive lower bound even if the display grid rounds to zero."""
    if x.lo <= 0:
        return None
    value = F(1)
    while value > x.lo:
        value /= 10
    return str(value)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--freeze-commit", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if not re.fullmatch(r"[0-9a-f]{40}", args.freeze_commit):
        raise ValueError("provide the full committed contract/code SHA")
    started = time.perf_counter()
    started_utc = datetime.now(timezone.utc).isoformat()
    sources = json.loads((ROOT / "SOURCES.json").read_text())
    for item in sources["inputs"]:
        if sha(ROOT / item["local_path"]) != item["sha256"]:
            raise ValueError(f"input hash mismatch: {item['local_path']}")
    # Also establish that this exact contract, scorer and inputs occur at F0.
    repository = Path(subprocess.check_output(
        ["git", "rev-parse", "--show-toplevel"], cwd=ROOT, text=True).strip())
    frozen_paths = ["CONTRACT.md", "SOURCES.json", "score.py"] + [
        item["local_path"] for item in sources["inputs"]]
    for name in frozen_paths:
        path = ROOT / name
        frozen = subprocess.check_output(
            ["git", "show", f"{args.freeze_commit}:{path.relative_to(repository).as_posix()}"],
            cwd=repository)
        if hashlib.sha256(frozen).hexdigest() != sha(path):
            raise ValueError(f"file differs from freeze commit: {name}")
    # Import only the pinned arithmetic class/serializer, never the old main/score.
    global Interval, interval_json
    sys.path.insert(0, str(ROOT / "vendor"))
    from p337_closed_source_score import Interval, interval_json

    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    child = [read_table(ROOT / f"inputs/endpoint-{name}.csv") for name in ("axis", "tilted")]
    endpoint = [parent_endpoint(data) for data in child]
    defect = [read_table(ROOT / f"inputs/defect-{name}.csv") for name in ("first", "second")]
    saved = json.loads((ROOT / "inputs/endpoint-root.json").read_text())["root_enclosure"]
    root = Interval(1-F(saved["upper_fraction"]), 1-F(saved["lower_fraction"]))
    if not 0 < root.lo < root.hi < 1:
        raise ValueError("invalid complemented parent root bracket")
    qmean_coefficients = [(a+b)/F(2) for a, b in zip(endpoint[0]["q"], endpoint[1]["q"])]
    root_signs = [polynomial(qmean_coefficients, p) for p in (root.lo, root.hi)]
    if not root_signs[0] < 0 < root_signs[1]:
        raise ValueError("saved parent root does not have exact opposite endpoint signs")
    intact = [moments(data, root) for data in endpoint]
    one_hole = [moments(data, root) for data in defect]
    matrix = []
    for base, hole in zip(intact, one_hole):
        for field in ("q", "e"):
            T = base[field][1]
            C = base[field+"sstar"][0] - base[field][0]*base["sstar"][0]
            H = N*(1-root)*(hole[field][0]-base[field][0])
            matrix.append([T, C, H])
    slope = (matrix[0][0]+matrix[2][0])/2
    if slope.lo <= 0:
        raise ValueError("pooled root slope is not bounded strictly above zero")
    minors = [(indices, det3([matrix[i] for i in indices])) for indices in combinations(range(4), 3)]
    D3 = sum((square(value) for _, value in minors), Interval.of(0))
    minors2 = [(indices, matrix[indices[0]][0]*matrix[indices[1]][1]
                - matrix[indices[0]][1]*matrix[indices[1]][0])
               for indices in combinations(range(4), 2)]
    D2 = sum((square(value) for _, value in minors2), Interval.of(0))
    alpha_C = ((matrix[0][1]+matrix[2][1])/2)/slope
    alpha_H = ((matrix[0][2]+matrix[2][2])/2)/slope
    comoving = [[T, C-T*alpha_C, H-T*alpha_H] for T, C, H in matrix]
    rejected = D3.lo > 0
    result = {
        "schema": "matching-one.p337-two-coupling-closure.v1",
        "freeze_commit": args.freeze_commit, "sources": sources,
        "decision": "finite_root_two_coupling_profile_closure_rejected" if rejected else
                    "not_excluded_by_fixed_root_enclosure",
        "parent_N": 50, "generators": [[5, 5], [1, 7]], "rows": ROWS,
        "columns": ["T=partial_p_f", "C=Cov(f,S)", "H=partial_epsilon_f"],
        "root": interval_json(root), "pooled_slope": interval_json(slope),
        "root_endpoint_checks": {"exact_signs": [-1, 1], "qmean": [
            interval_json(Interval.of(value)) for value in root_signs]},
        "count_checks": {"each_K": "binomial(25,K)", "table_totals": [
            sum(data["count"]) for data in child+defect], "normalizer_identically_one": True},
        "parent_mapping": "child reversed k; q and qS negated; E and S unchanged",
        "matrix": [[interval_json(value) for value in row] for row in matrix],
        "minors3": [{"rows": [ROWS[i] for i in indices], "determinant": interval_json(value)}
                    for indices, value in minors],
        "D3": interval_json(D3), "D3_positive_rational_witness": positive_witness(D3),
        "auxiliary_endpoint_rank": {"minors2": [
            {"rows": [ROWS[i] for i in indices], "determinant": interval_json(value)}
            for indices, value in minors2], "sum_squared_minors": interval_json(D2),
            "rank_two_certified": D2.lo > 0},
        "root_comoving": {"interpretation": "same matrix under column operations; no extra test",
            "columns": ["T", "Cbar", "Hbar"], "alpha_C": interval_json(alpha_C),
            "alpha_H": interval_json(alpha_H),
            "matrix": [[interval_json(value) for value in row] for row in comoving]},
        "scope": "exact finite counts; computational enclosures, not confidence intervals; "
                 "no inference about scalar U alone, continuum fields or independent confirmation",
        "new_samples": 0, "new_enumerations": 0, "cloud_jobs": 0,
    }
    (out / "latest.json").write_text(json.dumps(result, indent=2)+"\n")
    number = lambda value: f"{float((value.lo+value.hi)/2):.14g}"
    table = "\n".join(f"| {name} | " + " | ".join(map(number, row)) + " |"
                      for name, row in zip(ROWS, matrix))
    conclusion = ("The positive rational lower bound rejects common thermal-plus-S profile closure "
                  "at this finite root, including p-dependent coordinate changes. At least one "
                  "additional finite response direction is needed beyond the two endpoint tangents."
                  if rejected else "The fixed root enclosure does not exclude the candidate. "
                  "A zero obstruction at one root would not prove closure on a neighborhood.")
    (out / "REPORT.md").write_text(f"""# Fixed two-coupling endpoint profile test

**{result['decision']}**. {conclusion}

| Profile | T | C | H |
|---|---:|---:|---:|
{table}

All four predeclared 3x3 minors enter D3; D3 ≈ {number(D3)}.
Its outward rational enclosure is [{result['D3']['lower_fraction']}, {result['D3']['upper_fraction']}].
The positive lower-bound witness, if any, is `{result['D3_positive_rational_witness']}`.
The auxiliary sum of squared (T,C) minors is ≈ {number(D2)};
endpoint rank two certified: {D2.lo > 0}. Full bounds and the equivalent root-comoving columns
are in `latest.json`. D3 magnitude depends on units; its zero criterion does not.

The four pinned tables pass every K/count check. Parent endpoint coefficients use the
fixed child complement/reversal. Exact pooled-root endpoint signs are negative/positive
and the slope is strictly positive. All arithmetic before serialization uses Fraction
intervals; reported bounds are rounded outward on the 1e-40 rational grid. The primary
decision uses the unrounded rational lower bound, with a short exact positive witness.

This tests one common map for all four profiles, not scalar U alone or separate maps
per profile. An effective source shift can survive an earlier scalar-gain rejection:
U(epsilon,t)=u(t+c*epsilon) gives U*U_t_epsilon-U_epsilon*U_t=c*(u*u''-u'^2),
which need not vanish. Thus that earlier rejection does not settle this tangent test.

The result is conditional on the published finite graph counts and their different
quotient Smith classes. It is not a continuum field count, an independent confirmation,
or a revision of previous P154/P334/F4 stop decisions. No samples, enumeration or cloud
jobs were added. Freeze commit: `{args.freeze_commit}`.
""")
    receipt = {"started_utc": started_utc, "finished_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": time.perf_counter()-started, "command": sys.argv,
        "python": sys.version, "machine": platform.machine(), "freeze_commit": args.freeze_commit,
        "hashes": {name: sha(ROOT/name) for name in frozen_paths},
        "output_hashes": {name: sha(out/name) for name in ("latest.json", "REPORT.md")},
        "new_samples": 0, "new_enumerations": 0, "cloud_jobs": 0, "exit_code": 0}
    (out / "run.json").write_text(json.dumps(receipt, indent=2)+"\n")
    print(json.dumps({"decision": result["decision"], "D3": result["D3"],
                      "elapsed_seconds": receipt["elapsed_seconds"]}))


if __name__ == "__main__":
    main()

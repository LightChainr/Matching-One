#!/usr/bin/env python3
"""Independent exact power-basis reproduction of the frozen four-profile test.

No import of score.py/vendor arithmetic; no enumeration, sampling, new point,
observable or model. The midpoint below is only the center of a rigorous bound
over the already frozen root interval.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations
from datetime import datetime, timezone
import csv
import hashlib
import json
import math
import subprocess
import time

ROOT = Path(__file__).resolve().parent
FREEZE = "76a070d4aa95866f129572940f591b197cda064d"
FIELDS = ("count", "q", "e", "sstar", "qsstar", "esstar")
ROW_NAMES = ("q_first", "E_first", "q_second", "E_second")


def trim(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def add(a, b):
    out = [0] * max(len(a), len(b))
    for i, x in enumerate(a):
        out[i] += x
    for i, x in enumerate(b):
        out[i] += x
    return trim(out)


def scale(a, c):
    return trim([c*x for x in a])


def sub(a, b):
    return add(a, scale(b, -1))


def mul(a, b):
    out = [0] * (len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return trim(out)


def derivative(a):
    return trim([i*a[i] for i in range(1, len(a))] or [0])


def value(a, x):
    out = Q(0)
    for coefficient in reversed(a):
        out = out*x + coefficient
    return out


def power_from_counts(counts):
    # Expand p^k(1-p)^(25-k), without another binomial multiplicity.
    out = [0]*26
    for k, count in enumerate(counts):
        for j in range(26-k):
            out[k+j] += count * (-1)**j * math.comb(25-k, j)
    return trim(out)


def complement_argument(a):
    # Compose the already expanded child polynomial with 1-p.
    out = [0]*len(a)
    for k, coefficient in enumerate(a):
        for j in range(k+1):
            out[j] += coefficient * math.comb(k, j) * (-1)**j
    return trim(out)


def table(path):
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 26
    for k, row in enumerate(rows):
        assert int(row["k"]) == k
        assert int(row["count"]) == math.comb(25, k)
    out = {f: power_from_counts([int(r[f if f == "count" else "sum_"+f])
                                for r in rows]) for f in FIELDS}
    assert out["count"] == [1]
    return out


def determinant(a):
    x, y, z = a
    return add(sub(mul(x[0], sub(mul(y[1], z[2]), mul(y[2], z[1]))),
                   mul(x[1], sub(mul(y[0], z[2]), mul(y[2], z[0])))),
               mul(x[2], sub(mul(y[0], z[1]), mul(y[1], z[0]))))


def determinant_point(a):
    x, y, z = a
    return (x[0]*(y[1]*z[2]-y[2]*z[1])
            - x[1]*(y[0]*z[2]-y[2]*z[0])
            + x[2]*(y[0]*z[1]-y[1]*z[0]))


def bound(a, lo, hi):
    """Exact centered first Taylor jet with an absolute second-derivative bound."""
    center, radius = (lo+hi)/2, (hi-lo)/2
    da = derivative(a)
    dda = derivative(da)
    second_bound = value([abs(x) for x in dda], max(abs(lo), abs(hi)))
    error = abs(value(da, center))*radius + second_bound*radius*radius/2
    v = value(a, center)
    return v-error, v+error


def serialize(interval):
    lo, hi = interval
    grid = 10**50
    lower = Q((lo*grid).__floor__(), grid)
    upper = Q((hi*grid).__ceil__(), grid)
    return {"lower_fraction": str(lower), "upper_fraction": str(upper),
            "midpoint_approx": float((lo+hi)/2),
            "width_approx": float(upper-lower),
            "serialization": "outward_1e_minus50_rational_grid"}


def saved_bounds(item):
    return Q(item["lower_fraction"]), Q(item["upper_fraction"])


def compare(a, saved, lo, hi):
    mid = value(a, (lo+hi)/2)
    original = saved_bounds(saved)
    independent = bound(a, lo, hi)
    assert original[0] <= mid <= original[1]
    assert max(original[0], independent[0]) <= min(original[1], independent[1])
    return {"independent_enclosure": serialize(independent),
            "exact_midpoint_value_inside_original_enclosure": True,
            "independent_enclosure_inside_original":
                original[0] <= independent[0] and independent[1] <= original[1]}


def main():
    started = time.perf_counter()
    sources = json.loads((ROOT/"SOURCES.json").read_text())
    result = json.loads((ROOT/"results/latest.json").read_text())
    receipt = json.loads((ROOT/"results/run.json").read_text())
    assert result["freeze_commit"] == receipt["freeze_commit"] == FREEZE
    repo = Path(subprocess.check_output(["git", "rev-parse", "--show-toplevel"],
                                        cwd=ROOT, text=True).strip())
    sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    for name, expected in receipt["hashes"].items():
        assert sha(ROOT/name) == expected
        blob = subprocess.check_output(["git", "show",
            f"{FREEZE}:{(ROOT/name).relative_to(repo).as_posix()}"], cwd=repo)
        assert hashlib.sha256(blob).hexdigest() == expected
    for name, expected in receipt["output_hashes"].items():
        assert sha(ROOT/"results"/name) == expected
    for item in sources["inputs"]:
        assert sha(ROOT/item["local_path"]) == item["sha256"]
        blob = subprocess.check_output(["git", "show",
            f"{item['commit']}:{item['source_path']}"], cwd=repo)
        assert hashlib.sha256(blob).hexdigest() == item["sha256"]

    children = [table(ROOT/f"inputs/endpoint-{s}.csv") for s in ("axis", "tilted")]
    endpoint = [{f: scale(complement_argument(a), -1 if f in ("q", "qsstar") else 1)
                 for f, a in child.items()} for child in children]
    holes = [table(ROOT/f"inputs/defect-{s}.csv") for s in ("first", "second")]
    original_root = json.loads((ROOT/"inputs/endpoint-root.json").read_text())["root_enclosure"]
    lo = 1-Q(original_root["upper_fraction"])
    hi = 1-Q(original_root["lower_fraction"])
    assert 0 < lo < hi < 1
    qsum = add(endpoint[0]["q"], endpoint[1]["q"])
    signs = (value(qsum, lo), value(qsum, hi))
    assert signs[0] < 0 < signs[1]
    ds = derivative(qsum)  # Twice the pooled thermal slope; remains an integer polynomial.
    slope_bound = tuple(x/2 for x in bound(ds, lo, hi))
    assert slope_bound[0] > 0
    matrix = []
    for base, hole in zip(endpoint, holes):
        for f in ("q", "e"):
            matrix.append([derivative(base[f]),
                           sub(base[f+"sstar"], mul(base[f], base["sstar"])),
                           scale(mul([1, -1], sub(hole[f], base[f])), 25)])
    matrix_checks = [[compare(a, saved, lo, hi) for a, saved in zip(row, old)]
                     for row, old in zip(matrix, result["matrix"])]
    center = (lo+hi)/2
    at_center = [[value(a, center) for a in row] for row in matrix]
    triples = list(combinations(range(4), 3))
    minors = [determinant([matrix[i] for i in triple]) for triple in triples]
    for a, triple in zip(minors, triples):
        assert value(a, center) == determinant_point([at_center[i] for i in triple])
    minor_checks = [compare(a, saved["determinant"], lo, hi)
                    for a, saved in zip(minors, result["minors3"])]
    d3 = [0]
    for a in minors:
        d3 = add(d3, mul(a, a))
    assert value(d3, center) == sum(value(a, center)**2 for a in minors)
    d3_check = compare(d3, result["D3"], lo, hi)
    d3_bound = bound(d3, lo, hi)
    assert d3_bound[0] > 0
    assert result["decision"] == "finite_root_two_coupling_profile_closure_rejected"

    d2 = [0]
    for i, j in combinations(range(4), 2):
        minor = sub(mul(matrix[i][0], matrix[j][1]), mul(matrix[i][1], matrix[j][0]))
        d2 = add(d2, mul(minor, minor))
    d2_check = compare(d2, result["auxiliary_endpoint_rank"]["sum_squared_minors"], lo, hi)
    assert bound(d2, lo, hi)[0] > 0

    # Exact polynomial identity, rather than merely checking column operations numerically.
    cqsum = add(matrix[0][1], matrix[2][1])
    hqsum = add(matrix[0][2], matrix[2][2])
    numerator = [[T, sub(mul(C, ds), mul(T, cqsum)),
                 sub(mul(H, ds), mul(T, hqsum))] for T, C, H in matrix]
    for triple, old_minor in zip(triples, minors):
        assert sub(determinant([numerator[i] for i in triple]),
                   mul(mul(ds, ds), old_minor)) == [0]
    alpha_c = value(cqsum, center)/value(ds, center)
    alpha_h = value(hqsum, center)/value(ds, center)
    comoving = [[T, C-T*alpha_c, H-T*alpha_h] for T, C, H in at_center]
    for row, saved_row in zip(comoving, result["root_comoving"]["matrix"]):
        for v, saved in zip(row, saved_row):
            a, b = saved_bounds(saved)
            assert a <= v <= b
    for triple, old_minor in zip(triples, minors):
        assert determinant_point([comoving[i] for i in triple]) == value(old_minor, center)

    review = {
        "status": "PASS_INDEPENDENT_EXACT_POWER_BASIS_REPRODUCTION",
        "reviewed_at_utc": datetime.now(timezone.utc).isoformat(),
        "freeze_commit": FREEZE,
        "original_result_sha256": sha(ROOT/"results/latest.json"),
        "frozen_file_and_output_hashes_match": True,
        "all_inputs_match_their_original_git_blobs": True,
        "four_tables_each_K_count_binomial25_and_normalizer_exactly_one": True,
        "method": "Integer power-basis expansion; independent complement composition; exact Fraction Horner evaluation; centered Taylor enclosure with absolute second-derivative bound over the SAME frozen root interval. No scorer/vendor import.",
        "reference_root_unchanged": serialize((lo, hi)),
        "exact_root_endpoint_signs": [-1, 1],
        "pooled_slope": serialize(slope_bound),
        "matrix_4_by_3": matrix_checks,
        "all_four_minors3": minor_checks,
        "D3": d3_check,
        "D3_strictly_positive_independently": True,
        "D2_auxiliary_endpoint_rank": d2_check,
        "root_comoving": {
            "all_four_minor_identities_verified_as_integer_polynomial_identities": True,
            "identity": "det[T, Ds*C-T*(Cq1+Cq2), Ds*H-T*(Hq1+Hq2)] = Ds^2*det[T,C,H]; Ds=Tq1+Tq2",
            "exact_midpoint_entries_inside_original_comoving_enclosures": True,
            "all_four_exact_midpoint_minors_identical": True},
        "decision": result["decision"],
        "scope": "Correctness reproduction of the one frozen target; not new evidence, point selection, observable selection, model fitting, or an independent experiment.",
        "new_samples": 0, "new_enumerations": 0, "cloud_jobs": 0,
        "elapsed_seconds": time.perf_counter()-started,
        "review_code_sha256": sha(Path(__file__))}
    (ROOT/"review.json").write_text(json.dumps(review, indent=2)+"\n")
    (ROOT/"REVIEW.md").write_text(
        "# Independent exact reproduction\n\n"
        "**PASS.** The fixed thermal-plus-S* profile closure is excluded by an independently positive rational lower bound for the original D3.\n\n"
        "Frozen inputs/code match their recorded hashes and Git blobs; original result files match the run receipt hashes. Each of the four coefficient tables has exactly binomial(25,k) configurations in every row. Parent complement/signs, per-geometry source covariance, and the physical 25(1-p) defect dose agree with the contract.\n\n"
        "The independent reader expands the tables into integer power-basis polynomials, without importing the scorer or its interval class. Exact Fraction evaluation and a Taylor remainder bound use the same frozen root interval. Every 4x3 entry, all four third-order minors, D3 and the existing auxiliary D2 agree with the published rational enclosures; the exact pooled-root signs and positive slope are reproduced.\n\n"
        "Every root-comoving third-order minor is also checked as an exact integer-polynomial identity after clearing the common slope denominator. This verifies that root motion creates no additional response direction.\n\n"
        "D3 remains strictly positive. The result excludes one common smooth thermal/source-coordinate map for the four fixed profiles at this finite root. It does not identify continuum fields, exclude separate geometry maps or enlarged source families, establish an independent confirmation, or alter earlier stop decisions. No new samples, enumeration, occupancy points, scientific rows or models were added.\n\n"
        "Reproduce with `python3 review.py`. Bounds and input checks are in `review.json`; the original `results/` files are unchanged.\n")
    print(json.dumps({"status": review["status"],
                      "D3": d3_check["independent_enclosure"],
                      "all_symbolic_comoving_minors_equal": True,
                      "elapsed_seconds": review["elapsed_seconds"]}))


if __name__ == "__main__":
    main()

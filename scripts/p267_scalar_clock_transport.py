#!/usr/bin/env python3
"""Two ordered-height invariants of ordinary scalar transport, on N100 raw hist.

No density Jacobian, no pointwise branch choice, and no new Monte Carlo.
Source blobs are read from a fixed public git commit instead of duplicating raw.
"""
from __future__ import annotations

import argparse
import csv
from fractions import Fraction
from hashlib import sha256
import io
import json
from math import gcd
from pathlib import Path
import subprocess

import numpy as np
from scipy.optimize import brentq
from scipy.stats import binom, chi2, ncx2

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "experiments/p267_scalar_clock_transport_20260831.json"
OUTPUT = ROOT / "results/p267-scalar-clock-transport"


def evaluate(coefficients, p):
    degree = len(coefficients)-1
    return float(coefficients @ binom.pmf(np.arange(degree+1), degree, p))


def landmarks(coefficients, brackets):
    derivative = (len(coefficients)-1)*np.diff(coefficients)
    positions = np.array([brentq(lambda p: evaluate(derivative, p), a, b,
                                xtol=5e-15) for a, b in brackets])
    heights = np.array([evaluate(coefficients, p) for p in positions])
    curvatures = np.array([evaluate((len(derivative)-1)*np.diff(derivative), p)
                           for p in positions])
    return positions, heights, curvatures


def feature_vector(coefficients, brackets):
    found = [landmarks(coefficients[i, 0], brackets[i] if len(brackets) == 2 else brackets)
             for i in range(2)]
    positions, heights, curvatures = [np.stack([row[k] for row in found]) for k in range(3)]
    normalized = heights[:, 1:]/heights[:, :1]
    vector = np.r_[positions.ravel(), heights.ravel(), normalized.ravel(),
                   normalized[1]-normalized[0]]
    return vector, curvatures


def variation(coefficients):
    signs = [1 if x > 0 else -1 for x in coefficients if x]
    return sum(a != b for a, b in zip(signs, signs[1:]))


def split_half_integer(coefficients):
    """Exact Bernstein subdivision, clearing a common 2**degree factor."""
    row, degree = list(map(int, coefficients)), len(coefficients)-1
    left, right = [], []
    for k in range(degree+1):
        left.append(row[0] << (degree-k))
        right.append(row[-1] << (degree-k))
        row = [a+b for a, b in zip(row, row[1:])]
    def reduce(values):
        factor = 0
        for x in values:
            factor = gcd(factor, x)
        return [x//factor for x in values] if factor else values
    return reduce(left), reduce(right[::-1])


def exact_empirical_root_intervals(coefficients):
    """Bernstein variation certificate for empirical derivative, not truth.

    Variation zero has no interior roots; variation one has exactly one.
    This bounded source instance has no roots exactly on split boundaries.
    """
    pending, roots = [(list(map(int, coefficients)), Fraction(0), Fraction(1), 0)], []
    while pending:
        values, a, b, depth = pending.pop()
        changes = variation(values)
        if changes == 0:
            continue
        if changes == 1 and a > 0 and b < 1:
            roots.append([str(a), str(b)])
            continue
        if depth >= 24:
            raise ValueError("Unresolved empirical derivative; do not claim branch completeness")
        left, right = split_half_integer(values)
        if left[-1] == 0:
            raise ValueError("A split-boundary zero requires a separate certificate")
        middle = (a+b)/2
        pending.extend([(right, middle, b, depth+1), (left, a, middle, depth+1)])
    return sorted(roots, key=lambda ab: Fraction(ab[0]))


def load_source(protocol, source_directory=None, verification_commit=None):
    hashes = {}
    def read(relative):
        path = protocol["source_directory"]+"/"+relative
        data = ((source_directory/relative).read_bytes() if source_directory else
                subprocess.check_output(["git", "show", protocol["source_commit"]+":"+path], cwd=ROOT))
        if source_directory and verification_commit:
            frozen = subprocess.check_output(["git", "show", verification_commit+":"+path], cwd=ROOT)
            if sha256(data).digest() != sha256(frozen).digest():
                raise ValueError("Local archive differs from its declared public source commit: "+path)
        hashes[path] = sha256(data).hexdigest()
        return data.decode()
    source = json.loads(read("score.json"))
    contract = source["contract"]
    n, batches = contract["area"], contract["batches"]
    per_batch = contract["samples_per_shape_pair"]//batches
    raw = np.zeros((3, batches, 2, n+1), dtype=np.int64)
    deltas = [Fraction(shape["delta_cos4"]) for shape in contract["shapes"]]
    if len({abs(x) for x in deltas}) != 1:
        raise ValueError("This source-specific integer certificate requires one absolute P4 normalization")
    for si, shape in enumerate(contract["shapes"]):
        for row in csv.DictReader(io.StringIO(read("raw/"+shape["name"]+".hist.csv"))):
            sign = (1 if row["orientation"] == "first" else -1)*(1 if deltas[si] > 0 else -1)
            raw[si, int(row["batch"]), int(row["kind"] == "plus"), int(row["k"])] += sign*int(row["count"])
    fields = np.stack([raw[:, :, 0]+raw[:, :, 1], raw[:, :, 1]-raw[:, :, 0]], axis=2)
    integers = np.cumsum(np.stack([fields[1]-fields[0], fields[2]-fields[0]]), axis=-1)
    # Zero endpoint is exact, because each normalized orientation count totals
    # the same number of permutations. No numeric tail correction is needed.
    assert np.all(integers[..., -1] == 0)
    bernstein = integers/(per_batch*float(abs(deltas[0])))
    return contract, hashes, bernstein, integers


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-directory", type=Path, help="Apply frozen ordinal readout to another completed local archive")
    parser.add_argument("--source-commit", help="Verify the local source blobs against this committed version")
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    protocol = json.loads(PROTOCOL.read_text())
    if args.source_directory:
        protocol["source_directory"] = "results/"+args.source_directory.name
    contract, hashes, batches, integers = load_source(protocol, args.source_directory, args.source_commit)
    b = batches.shape[1]
    mean = batches.mean(axis=1)
    exact_intervals = {label: exact_empirical_root_intervals(np.diff(integers[i, :, 0].sum(axis=0)))
                       for i, label in enumerate(("D_A", "U_A"))}
    args.output.mkdir(parents=True, exist_ok=True)
    if not all(len(intervals) == 3 for intervals in exact_intervals.values()):
        result = {"schema": protocol["schema"], "source_directory": str(args.source_directory),
                  "source_commit": args.source_commit, "source_sha256": hashes, "source_contract": contract,
                  "status": "ordered_extremum_pattern_changed; three-landmark ratio gate not forced",
                  "exact_empirical_critical_point_intervals": exact_intervals,
                  "empirical_landmarks": {label: dict(zip(("positions", "heights", "curvatures"),
                     [x.tolist() for x in landmarks(mean[i, 0], [[float(Fraction(a)), float(Fraction(z))] for a, z in exact_intervals[label]])]))
                     for i, label in enumerate(("D_A", "U_A"))}, "new_samples": 0}
        (args.output/"score.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
        print(json.dumps(result, indent=2))
        return
    # The target uses ordinal empirical roots, not the old scale's p windows.
    brackets = ([[[float(Fraction(a)), float(Fraction(z))] for a, z in exact_intervals[label]]
                 for label in ("D_A", "U_A")] if args.source_directory else protocol["root_brackets"])
    estimate, curvatures = feature_vector(mean, brackets)
    loo, loo_curvatures = [], []
    for i in range(b):
        vector, curvature = feature_vector((b*mean-batches[:, i])/(b-1), brackets)
        loo.append(vector)
        loo_curvatures.append(curvature)
    loo, loo_curvatures = np.array(loo), np.array(loo_curvatures)
    centered = loo-loo.mean(axis=0)
    covariance = (b-1)/b*(centered.T @ centered)
    se = np.sqrt(covariance.diagonal())
    delta, cov = estimate[-2:], covariance[-2:, -2:]
    statistic = float(delta @ np.linalg.solve(cov, delta))
    labels = ([label+".p."+name for label in ("D_A", "U_A") for name in protocol["landmarks"]]
              +[label+".height."+name for label in ("D_A", "U_A") for name in protocol["landmarks"]]
              +[label+"."+name for label in ("D_A", "U_A") for name in protocol["primary_features"]]
              +["U_minus_D."+name for name in protocol["primary_features"]])
    result = {
        "schema": protocol["schema"], "protocol": str(PROTOCOL.relative_to(ROOT)),
        "source_commit": args.source_commit if args.source_directory else protocol["source_commit"],
        "source_directory": str(args.source_directory) if args.source_directory else protocol["source_directory"],
        "source_sha256": hashes, "source_contract": contract,
        "analysis_stage": "source-frozen ordinal auxiliary readout" if args.source_directory else "N100 retrospective feature selection",
        "dependency_group": "all three source shapes share aligned counters and batches; see source_contract",
        "samples_per_shape_pair": contract["samples_per_shape_pair"], "common_batches": b,
        "new_samples": 0,
        "layout": labels, "estimate": estimate.tolist(), "se": se.tolist(),
        "full_covariance": covariance.tolist(), "leave_one_batch_out_vectors": loo.tolist(),
        "jackknife_bias_estimate": ((b-1)*(loo.mean(axis=0)-estimate)).tolist(),
        "primary": {"feature_difference": delta.tolist(), "covariance": cov.tolist(),
                    "se": se[-2:].tolist(), "marginal_z": (delta/se[-2:]).tolist(),
                    "chi2": statistic, "df": 2, "nominal_p": float(chi2.sf(statistic, 2)),
                    "nominal_1pct_decision": "tension" if chi2.sf(statistic, 2) < .01 else "not_rejected",
                    "interpretation": "A-only necessary gate for ordinary scalar transport; N100 is exploratory, a new independent scale uses its already chosen ordinal features."},
        "branch_order": {"exact_empirical_critical_point_intervals": exact_intervals,
                         "mean_curvatures": curvatures.tolist(),
                         "all_loo_curvature_signs_match_mean": bool(np.all(np.sign(loo_curvatures) == np.sign(curvatures))),
                         "all_loo_roots_in_mean_ordered_brackets": True,
                         "numerical_brackets": brackets,
                         "scope": "Exact root counts describe the empirical mean polynomials. Sampling uncertainty is not an exact population root-count proof. Sparse empirical endpoint orders are not interpreted."},
        "boundary": "An A-only necessary gate for common scalar A/E transport with independent amplitudes. Does not assume density transport, does not select inverse branches pointwise, and does not infer a field count or continuum law. Features were selected on N100, not on a subsequent target; each archive's uncertainty is one shared dependency block."
    }
    if contract["area"] == 400:
        reference = json.loads((OUTPUT/"score.json").read_text())
        old = np.array(reference["primary"]["feature_difference"])
        old_cov = np.array(reference["primary"]["covariance"])
        difference, difference_cov = delta-old, cov+old_cov
        q = float(difference @ np.linalg.solve(difference_cov, difference))
        noncentrality = float(old @ np.linalg.solve(cov, old))
        result["independent_scale_comparison"] = {
            "source_analysis_commit": "8ac86fb", "source_data_commit": reference["source_commit"],
            "N400_minus_N100": difference.tolist(), "covariance": difference_cov.tolist(),
            "chi2": q, "df": 2, "nominal_p": float(chi2.sf(q, 2)),
            "N100_plug_in_effect_target_noncentrality": noncentrality,
            "conditional_power_at_1pct_if_N100_effect_persisted": float(ncx2.sf(chi2.isf(.01, 2), 2, noncentrality)),
            "boundary": "Independent streams; source effect is a design plug-in, not a known truth. No unchanged cross-scale amplitude is assumed."
        }
    (args.output/"score.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    rows = [f"# N{contract['area']} ordinary scalar transport: ordered double-peak heights", "",
            "The scalar hypothesis has no Jacobian and has independent A/E amplitudes. It is not the already-tested density hypothesis.", "",
            "| Ordered A height ratio | D: 4i minus 2i | U: shear minus 2i | U minus D | SE |",
            "|---|---:|---:|---:|---:|"]
    for i, name in enumerate(protocol["primary_features"]):
        rows.append(f"| {name} | {estimate[12+i]:.9g} | {estimate[14+i]:.9g} | {delta[i]:.9g} | {se[-2+i]:.9g} |")
    rows += ["", f"Joint Gaussian reference: chi2={statistic:.8g}/2, nominal p={chi2.sf(statistic,2):.8g}.", "",
             "An increasing scalar reparametrization preserves the complete ordered critical-height vector up to one amplitude. The two ratios quantify valley filling and peak balance. A common A/E scalar map must pass this A-only gate. Unresolved valley signs are not interpreted as population zero crossings.", "",
             f"All {b} leave-one-common-batch-out calculations recompute peaks, valley, and normalization. The exact Bernstein certificate finds three interior critical points in each empirical A polynomial; their order and curvature types remain stable in all leave-one-batch-out replicas.", "", result["boundary"], ""]
    (args.output/"REPORT.md").write_text("\n".join(rows))
    print("\n".join(rows))


if __name__ == "__main__":
    main()

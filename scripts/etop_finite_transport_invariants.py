#!/usr/bin/env python3
"""Test common finite profile transport, beyond any velocity truncation.

Independent A/E amplitudes are allowed. Relative cumulative moments remove
the entire unknown monotone coordinate map, without dividing by U_E's area.
"""
from __future__ import annotations

import csv
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.special import betainc, roots_legendre
from scipy.stats import binom

from etop_thermal_transport import ROOT, SOURCE, load_fields, moment_kernel, describe

OUT = ROOT / "results/etop-finite-transport-invariants"


def odd_sign_certificate(contract, source_dir=SOURCE):
    """Integer Bernstein proof for these empirical curves, not their population."""
    nums = []
    for shape in contract["shapes"]:
        a = [0]*(contract["area"]+1)
        # This source's three magnitudes are all 1152/625; only signs differ.
        sign = -1 if shape["delta_cos4"].startswith("-") else 1
        with (source_dir/"raw"/(shape["name"]+".hist.csv")).open() as stream:
            for row in csv.DictReader(stream):
                a[int(row["k"])] += sign*(1 if row["orientation"] == "first" else -1)*int(row["count"])
        nums.append(a)
    certificates = {}
    for name, threshold in (("D_A", [a-b for a, b in zip(nums[1], nums[0])]),
                            ("minus_U_A", [a-b for a, b in zip(nums[0], nums[2])])):
        bernstein, total = [], 0
        for count in threshold:
            total += count
            bernstein.append(total)
        n = len(bernstein)-1
        work = bernstein.copy()
        left, right = [work[0] << n], [work[-1] << n]
        for depth in range(1, n+1):
            work = [a+b for a, b in zip(work, work[1:])]
            left.append(work[0] << (n-depth))
            right.append(work[-1] << (n-depth))
        right.reverse()
        certificates[name] = {"bernstein_integer_numerators": bernstein,
            "half_interval_integer_numerators": [[str(x) for x in left], [str(x) for x in right]],
            "half_interval_extra_denominator": str(1 << n),
            "strictly_positive_on_open_unit_interval": bool(
                min(left) >= 0 and min(right) >= 0 and max(left) > 0 and max(right) > 0 and left[-1] > 0)}
    return {"common_positive_scale": f"625/(1152*{contract['samples_per_shape_pair']})",
            "method": "One exact de Casteljau subdivision at p=1/2; nonnegative coefficients on both halves and a positive midpoint prove strict positivity inside (0,1).",
            "curves": certificates,
            "consequence": "If both certificates hold, these empirical A curves admit a unique interior-regular cumulative-quantile map after area normalization. Endpoint derivatives may degenerate. This proves neither population-wide positivity nor a finite physical coordinate identification."}


def main(source_dir=SOURCE, output_dir=OUT, source_commit="7b30648"):
    source, fields = load_fields(source_dir)
    n = source["contract"]["area"]
    b = len(fields)
    max_m = 6
    # Degree of F_A^m D_E is at most m(N+1)+N. This rule integrates the
    # finite histogram polynomials exactly in exact arithmetic.
    nodes = (max_m*(n+1)+n+2)//2
    x, weights = roots_legendre(nodes)
    p, weights = (x+1)/2, weights/2
    k = np.arange(n+1, dtype=float)
    tails = binom.sf(k[:, None]-1, n, p[None, :])
    integrals = p[None, :]*tails-(k[:, None]/(n+1))*betainc(
        k[:, None]+1, n+1-k[:, None], p[None, :])
    area_kernel = moment_kernel(n, 0., 1., 0)[:, 0]
    mean = fields.mean(axis=0)
    loo = (b*mean-fields)/(b-1)

    def value(averaged, full=False):
        d, u = averaged[1]-averaged[0], averaged[2]-averaged[0]
        area_d, area_u = d@area_kernel, u@area_kernel
        r = area_u/area_d
        cd, cu = (d@integrals)[0]/area_d[0], (u@integrals)[0]/area_u[0]
        ed, eu = (d@tails)[1], (u@tails)[1]
        jd = np.array([np.dot(weights*cd**m, ed) for m in range(1, max_m+1)])
        ju = np.array([np.dot(weights*cu**m, eu) for m in range(1, max_m+1)])
        residual = ju-r[1]*jd
        # Omega=2 integral F_A D_E - A_0 E_0; this difference permits both
        # independent amplitudes and is exactly 2 U_A0 times residual[0].
        omega_d = area_d[0]*(2*jd[0]-area_d[1])
        omega_u = area_u[0]*(2*ju[0]-area_u[1])
        omega_residual = omega_u-r[0]*r[1]*omega_d
        result = {"residual": residual, "r": r, "jd": jd, "ju": ju,
                  "area_d": area_d, "area_u": area_u,
                  "omega": np.array([omega_residual]),
                  "omega_D_U": [float(omega_d), float(omega_u)]}
        if full:
            result["A_normalized_cumulative_range_D_U"] = [float(cd.min()), float(cd.max()), float(cu.min()), float(cu.max())]
        return result

    point = value(mean, True)
    deleted = [value(a) for a in loo]
    score = describe(point["residual"], np.array([a["residual"] for a in deleted]), with_score=True)
    omega = describe(point["omega"], np.array([a["omega"] for a in deleted]), with_score=True)
    r_score = describe(point["r"], np.array([a["r"] for a in deleted]))
    result = {"status": "Finite common-transport readout of the identified source block; no additional samples generated by this scorer",
        "source_commit": source_commit, "source_path": str((source_dir/"score.json").relative_to(ROOT)),
        "input_prediction_freeze_commit": source.get("prediction_freeze_commit"),
        "source_sha256": hashlib.sha256((source_dir/"score.json").read_bytes()).hexdigest(),
        "dependency": source["contract"]["sampling"],
        "candidate": "U_j(p)=r_j phi'(p) D_j(phi(p)), j=A,E, one common increasing endpoint-fixing phi and independent signed nonzero r_A,r_E",
        "definitions": ["D=Y(4i)-Y(2i), U=Y(1/2+i)-Y(2i)",
                        "F_A(p)=integral_0^p D_A(t) dt; A0=integral_0^1 D_A",
                        "J_m(D)=integral_0^1 [F_A(p)/A0]^m D_E(p) dp",
                        "Common finite transport implies J_m(U)=r_E J_m(D) for all m; no velocity truncation",
                        "Omega(D)=integral(F_A D_E-F_E D_A); Omega(U)=r_A*r_E Omega(D)"],
        "areas_D_A_E": point["area_d"].tolist(), "areas_U_A_E": point["area_u"].tolist(),
        "area_amplitudes_A_E": r_score, "moment_orders": list(range(1, max_m+1)),
        "source_moments_D": point["jd"].tolist(), "target_moments_U": point["ju"].tolist(),
        "finite_transport_remainders": score, "oriented_area_remainder": omega,
        "Omega_D_U": point["omega_D_U"],
        "empirical_odd_sign_certificate": odd_sign_certificate(source["contract"], source_dir),
        "quadrature": {"nodes": nodes, "max_integrand_degree": max_m*(n+1)+n,
                       "exact_arithmetic_scope": "Gauss-Legendre integrates these histogram polynomials exactly; displayed computation is double precision"},
        "scope": ["Allows arbitrary finite increasing coordinate map, not just translation/dilation or polynomial velocity.",
                  "A and E have separate area amplitudes; no equality of clock and lifetime transfer is imposed.",
                  "No division by the weak U_E area is performed; r_E divides only by D_E area.",
                  "Common transport is stronger than giving each observer its own map. A single positive normalized curve alone always has a quantile map.",
                  "This is signed-profile/Jacobian transport. Ordinary scalar-observable coordinate relabelling without a Jacobian is a different class.",
                  "One paired block, all source/target moments and fitted area amplitudes are delete-one jackknifed. Nominal Gaussian-reference scores are exploratory, not independent validation.",
                  "Neither a finite-transport failure nor its extra shape coordinates count continuum fields."]}
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir/"invariants.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    sign_proved = all(c["strictly_positive_on_open_unit_interval"] for c in result["empirical_odd_sign_certificate"]["curves"].values())
    sign_text = ("For the actual empirical histograms, exact integer Bernstein subdivision proves D_A>0 and U_A<0 for every0<p<1. After dividing by their signed areas, both are positive normalized profiles. Thus the odd data alone admit a unique cumulative-quantile map. The joint readout tests whether that same map also transports E. This existence statement allows degenerate endpoint derivatives and does not prove positivity of the underlying population curves."
                 if sign_proved else "One exact Bernstein subdivision does not certify the required signs for both empirical A curves at this scale. This is not a proof of a sign change. The cumulative-moment null remains a necessary common-transport condition without this positivity certificate; a quantile-map existence claim is not made for this target.")
    lines = ["# Does one finite thermal coordinate explain both odd and even profiles?", "",
        "A single sign-definite profile can always be matched to another equal-area profile by a cumulative-quantile map. The substantive question is whether the SAME map explains A and E, even after allowing each its own amplitude. The cumulative cross-moments remove that map exactly.", "",
        "Candidate: U_j(p)=r_j phi'(p)D_j(phi(p)), j=A,E, with phi increasing and fixing0/1. Define J_m(D)=integral(F_A/A0)^m D_E. Then J_m(U)-r_E J_m(D)=0 for every m, without a small-warp or polynomial approximation.", "",
        f"Measured area amplitudes: r_A={point['r'][0]:.9g}, r_E={point['r'][1]:.9g}. The even target area is not used as a denominator.", "",
        "| cumulative power m | target-source remainder | paired jackknife SE |", "|---:|---:|---:|"]
    for i in range(max_m):
        lines.append(f"| {i+1} | {point['residual'][i]:.9g} | {score['se'][i]:.6g} |")
    lines += ["", f"Joint: chi-square={score['chi_square']:.8g}/{score['df']}, nominal p={score['nominal_p_value']:.6g}.", "",
        f"The single oriented-area remainder Omega(U)-r_A*r_E Omega(D) is {point['omega'][0]:.9g} +/- {omega['se'][0]:.6g}; chi-square={omega['chi_square']:.7g}/1, nominal p={omega['nominal_p_value']:.6g}.", "",
        f"The six finite-polynomial integrals use {nodes}-node Gauss-Legendre (degree at most {max_m*(n+1)+n}), with all source and target uncertainty propagated by deleting each of {b} aligned batches. This scorer performs no additional sampling.", "",
        "## Why the odd curve alone is not enough", "",
        sign_text, "",
        "## Interpretation boundary", "",
        "These are necessary conditions for one common signed-profile transport, not for every imaginable redefinition of two observables. A separate map per observer, ordinary observable relabelling without its Jacobian, or a genuine extra response component remain different mechanisms. The oriented area and six moments reuse one data block and must not be counted as separate confirmations.", ""]
    (output_dir/"REPORT.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--output", type=Path, default=OUT)
    parser.add_argument("--source-commit", default="7b30648")
    args = parser.parse_args()
    main(args.source.resolve(), args.output.resolve(), args.source_commit)

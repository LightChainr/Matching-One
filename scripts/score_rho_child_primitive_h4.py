#!/usr/bin/env python3
"""Score the frozen non-A_top complex response on the three rho children."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from math import gcd
from pathlib import Path

import mpmath as mp

from derive_hexagonal_degree2_hecke import normalized_eisenstein
from pinson_arguin_primitive import engine_to_paper, primitive_probability_direct
from rho_child_primitive_h4_mc import CHILD_DESIGNS, physical_phase


def tau_from_matrix(matrix) -> mp.mpc:
    (a, b), (c, d) = matrix
    return mp.mpc(b, d) / mp.mpc(a, c)


def canonical_lines(cutoff: int):
    for u in range(-cutoff, cutoff + 1):
        for v in range(-cutoff, cutoff + 1):
            if (u, v) == (0, 0) or gcd(abs(u), abs(v)) != 1:
                continue
            if u < 0 or (u == 0 and v < 0):
                continue
            yield u, v


def continuum_harmonic(matrix, *, cutoff: int, dps: int) -> mp.mpc:
    tau = tau_from_matrix(matrix)
    total = mp.mpc(0)
    for line in canonical_lines(cutoff):
        sector = engine_to_paper(line)
        probability = primitive_probability_direct(*sector, tau, dps=dps)
        total += probability * physical_phase(matrix, line)
    return +total


def real_vector(values):
    return mp.matrix([coordinate for value in values for coordinate in (mp.re(value), mp.im(value))])


def complex_design(vectors):
    rows = []
    for child in range(len(vectors[0])):
        real_row = []
        imag_row = []
        for vector in vectors:
            value = vector[child]
            real_row.extend([mp.re(value), -mp.im(value)])
            imag_row.extend([mp.im(value), mp.re(value)])
        rows.extend([real_row, imag_row])
    return mp.matrix(rows)


def gls(y, covariance, design):
    precision = covariance ** -1
    gram = design.T * precision * design
    beta = mp.lu_solve(gram, design.T * precision * y)
    residual = y - design * beta
    chi2 = (residual.T * precision * residual)[0]
    dof = len(y) - design.cols
    p = mp.gammainc(mp.mpf(dof) / 2, chi2 / 2, mp.inf) / mp.gamma(mp.mpf(dof) / 2) if dof else mp.mpf(1)
    return {
        "chi_square": mp.nstr(chi2, 20),
        "dof": dof,
        "survival_p": mp.nstr(p, 20),
        "coefficients": [mp.nstr(value, 20) for value in beta],
        "residual": [mp.nstr(value, 20) for value in residual],
    }


def score(run: dict, manifest: dict, *, cutoff: int = 9, dps: int = 60) -> dict:
    if not manifest.get("production_authorized"):
        raise ValueError("manifest is not authorized")
    if run.get("manifest_runner_commit") != manifest.get("runner_commit"):
        raise ValueError("runner commit mismatch")
    acquisition = manifest["acquisition"]
    for key in ("samples", "batches", "workers", "seed"):
        if run[key] != acquisition[key]:
            raise ValueError(f"acquisition mismatch: {key}")
    summary = run["summary"]
    if not summary["all_invariant_failures_zero"]:
        raise ValueError("homology invariant failures are nonzero")
    point = [mp.mpf(str(value)) for value in summary["lattice_H4_point_re_im"]]
    covariance = mp.matrix([
        [mp.mpf(str(value)) for value in row]
        for row in summary["full_common_field_covariance_6x6"]
    ])
    baselines = []
    convergence = []
    shapes = {"E4_r1": [], "E6_r0": [], "E4sq_r2": []}
    for _, matrix in CHILD_DESIGNS:
        baseline = continuum_harmonic(matrix, cutoff=cutoff, dps=dps)
        previous = continuum_harmonic(matrix, cutoff=cutoff - 2, dps=dps)
        baselines.append(baseline)
        convergence.append(abs(baseline - previous))
        tau = tau_from_matrix(matrix)
        e4, _, _ = normalized_eisenstein(4, tau, dps=dps)
        e6, _, _ = normalized_eisenstein(6, tau, dps=dps)
        shapes["E4_r1"].append(e4)
        shapes["E6_r0"].append(e6)
        shapes["E4sq_r2"].append(e4 * e4)
    y = mp.matrix(point) - real_vector(baselines)
    zero_chi2 = (y.T * (covariance ** -1) * y)[0]
    zero_p = mp.gammainc(3, zero_chi2 / 2, mp.inf) / mp.gamma(3)
    pure_scores = {
        name: gls(y, covariance, complex_design([vector]))
        for name, vector in shapes.items()
    }
    pairs = (("E4_r1", "E6_r0"), ("E4_r1", "E4sq_r2"), ("E6_r0", "E4sq_r2"))
    mixture_scores = {
        "+".join(pair): gls(y, covariance, complex_design([shapes[pair[0]], shapes[pair[1]]]))
        for pair in pairs
    }
    zeta = mp.exp(2 * mp.pi * mp.j / 3)
    complex_y = [mp.mpc(y[2 * j], y[2 * j + 1]) for j in range(3)]
    dft = [sum(complex_y[j] * zeta ** (-r * j) for j in range(3)) / 3 for r in range(3)]
    dft_transform_rows = []
    for r in range(3):
        real_row = []
        imag_row = []
        for j in range(3):
            coefficient = zeta ** (-r * j) / 3
            real_row.extend([mp.re(coefficient), -mp.im(coefficient)])
            imag_row.extend([mp.im(coefficient), mp.re(coefficient)])
        dft_transform_rows.extend([real_row, imag_row])
    dft_transform = mp.matrix(dft_transform_rows)
    dft_covariance = dft_transform * covariance * dft_transform.T
    dft_components = []
    for r, value in enumerate(dft):
        component_covariance = mp.matrix([
            [dft_covariance[2*r+i, 2*r+j] for j in range(2)] for i in range(2)
        ])
        component = mp.matrix([mp.re(value), mp.im(value)])
        component_chi2 = (component.T * (component_covariance ** -1) * component)[0]
        component_p = mp.exp(-component_chi2 / 2)
        dft_components.append({
            "r": r,
            "value_re_im": [mp.nstr(mp.re(value), 20), mp.nstr(mp.im(value), 20)],
            "covariance_2x2": [
                [mp.nstr(component_covariance[i, j], 20) for j in range(2)]
                for i in range(2)
            ],
            "zero_chi_square": mp.nstr(component_chi2, 20),
            "zero_dof": 2,
            "zero_survival_p": mp.nstr(component_p, 20),
        })
    return {
        "schema": "matching-one/rho-child-primitive-h4-score/v1",
        "status": "frozen_reveal",
        "issues": [156, 205, 226, 250, 275],
        "observer": "delta_H4_primitive = E[1_rank1 u(ell)^4] - Pinson/Arguin H4(tau)",
        "observer_sector": "topology_typed_rank1_polarization_non_A_top",
        "child_order": [name for name, _ in CHILD_DESIGNS],
        "complex_frame": "u(ell)=P*ell/|P*ell| in the common laboratory square-lattice frame",
        "continuum_baseline": [[mp.nstr(mp.re(v), 30), mp.nstr(mp.im(v), 30)] for v in baselines],
        "baseline_cutoff": cutoff,
        "cutoff_minus_2_change_abs": [mp.nstr(v, 8) for v in convergence],
        "delta_H4_re_im": [mp.nstr(value, 20) for value in y],
        "full_covariance_6x6": [[mp.nstr(covariance[i, j], 20) for j in range(6)] for i in range(6)],
        "ideal_C3_DFT_components": dft_components,
        "ideal_C3_DFT_full_covariance_6x6": [
            [mp.nstr(dft_covariance[i, j], 20) for j in range(6)] for i in range(6)
        ],
        "null": {"chi_square": mp.nstr(zero_chi2, 20), "dof": 6, "survival_p": mp.nstr(zero_p, 20)},
        "pell_transported_pure_models": pure_scores,
        "two_character_mixture_opponents": mixture_scores,
        "saturated_mixture": {"parameters": 6, "dof": 0, "chi_square": "0"},
        "claim_boundary": [
            "The result classifies one topology-typed complex child vector; it does not reopen ordinary H4/H8/H12 voting.",
            "Pure labels use area-normalized E4, E6 and E4 squared at the actual Pell-child moduli, not the exact-rho limit alone.",
            "A mixture win means this observable is not one ordinary modular-ring character at N112; it is not a continuum no-go theorem.",
        ],
    }


def render(report: dict) -> str:
    lines = [
        "# Three-rho-child non-A_top complex C3 pilot", "",
        f"Observer: `{report['observer']}`.", "",
        "## Frozen scores", "",
        "| model | chi-square/df | p |", "|---|---:|---:|",
        f"| null | {report['null']['chi_square']}/{report['null']['dof']} | {report['null']['survival_p']} |",
    ]
    for name, value in report["pell_transported_pure_models"].items():
        lines.append(f"| {name} | {value['chi_square']}/{value['dof']} | {value['survival_p']} |")
    for name, value in report["two_character_mixture_opponents"].items():
        lines.append(f"| mixture {name} | {value['chi_square']}/{value['dof']} | {value['survival_p']} |")
    lines += ["", "## Complex response", ""]
    for index, name in enumerate(report["child_order"]):
        lines.append(f"- {name}: `{report['delta_H4_re_im'][2*index]} + i {report['delta_H4_re_im'][2*index+1]}`")
    lines += ["", "## Ideal C3 DFT diagnostics", ""]
    for component in report["ideal_C3_DFT_components"]:
        lines.append(
            f"- r={component['r']}: `{component['value_re_im']}`, "
            f"zero chi-square `{component['zero_chi_square']}/2`, "
            f"p `{component['zero_survival_p']}`"
        )
    lines += ["", "## Boundary", ""] + [f"- {line}" for line in report["claim_boundary"]] + [""]
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    parser.add_argument("--cutoff", type=int, default=9)
    parser.add_argument("--dps", type=int, default=60)
    args = parser.parse_args()
    report = score(json.loads(args.run.read_text()), json.loads(args.manifest.read_text()), cutoff=args.cutoff, dps=args.dps)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    args.markdown.write_text(render(report))
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

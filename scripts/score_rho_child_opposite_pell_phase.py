#!/usr/bin/env python3
"""Score the frozen opposite-Pell normalization-free complex phase test."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mpmath as mp

from derive_hexagonal_degree2_hecke import normalized_eisenstein
from rho_child_primitive_h4_mc import N60_CHILD_DESIGNS
from score_rho_child_primitive_h4 import continuum_harmonic, tau_from_matrix


def dot(left, right):
    return sum(a * b for a, b in zip(left, right))


def cross(left, right):
    return [
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    ]


def block_diagonal(first, second):
    output = mp.zeros(first.rows + second.rows, first.cols + second.cols)
    for i in range(first.rows):
        for j in range(first.cols):
            output[i, j] = first[i, j]
    for i in range(second.rows):
        for j in range(second.cols):
            output[first.rows + i, first.cols + j] = second[i, j]
    return output


def combine_runs(runs, manifest):
    order = [name for name, _ in N60_CHILD_DESIGNS]
    by_child = {run.get("selected_child"): run for run in runs}
    if sorted(by_child) != sorted(order):
        raise ValueError("exactly three N60 child runs are required")
    point = []
    covariance = mp.zeros(6, 6)
    for index, child in enumerate(order):
        run = by_child[child]
        expected = manifest["runs"][child]
        for key in ("family", "child", "samples", "batches", "workers", "seed", "replica_offset"):
            observed = run["selected_family"] if key == "family" else run["selected_child"] if key == "child" else run[key]
            if observed != expected[key]:
                raise ValueError(f"{child}: frozen acquisition mismatch for {key}")
        if run.get("manifest_runner_commit") != manifest["runner_commit"]:
            raise ValueError(f"{child}: runner commit mismatch")
        if run["child_gate"].get("family") != "N60_Dminus2":
            raise ValueError(f"{child}: wrong Pell family")
        if not run["child_gate"]["direction_alias_gate"]["all_rank_two"]:
            raise ValueError(f"{child}: direction alias gate failed")
        summary = run["summary"]
        if summary["primary_order"] != [f"{child}_re", f"{child}_im"]:
            raise ValueError(f"{child}: coordinate order mismatch")
        if not summary["all_invariant_failures_zero"]:
            raise ValueError(f"{child}: homology invariant failure")
        point.extend(mp.mpf(str(value)) for value in summary["lattice_H4_point_re_im"])
        block = summary["full_common_field_covariance_6x6"]
        if len(block) != 2 or any(len(row) != 2 for row in block):
            raise ValueError(f"{child}: covariance is not 2x2")
        for i in range(2):
            for j in range(2):
                covariance[2*index+i, 2*index+j] = mp.mpf(str(block[i][j]))
    return mp.matrix(point), covariance


def annihilator_coefficient(delta, covariance, designs, *, dps):
    e4, e6 = [], []
    for _, matrix in designs:
        tau = tau_from_matrix(matrix)
        e4.append(normalized_eisenstein(4, tau, dps=dps)[0])
        e6.append(normalized_eisenstein(6, tau, dps=dps)[0])
    w = cross(e4, e6)
    denominator = dot(w, [value * value for value in e4])
    if abs(dot(w, e4)) > mp.mpf("1e-40") or abs(dot(w, e6)) > mp.mpf("1e-40"):
        raise ArithmeticError("E4/E6 annihilator did not close")
    coefficients = [value / denominator for value in w]
    transform = mp.matrix([
        [coordinate for value in coefficients for coordinate in (mp.re(value), -mp.im(value))],
        [coordinate for value in coefficients for coordinate in (mp.im(value), mp.re(value))],
    ])
    coefficient = transform * delta
    coefficient_covariance = transform * covariance * transform.T
    return coefficient, coefficient_covariance, w, denominator


def _ray_objective(reference, new, covariance, sign, log_scale):
    scale = mp.exp(log_scale)
    design = mp.matrix([
        [1, 0], [0, 1], [sign * scale, 0], [0, sign * scale]
    ])
    values = mp.matrix([reference[0], reference[1], new[0], new[1]])
    precision = covariance**-1
    beta = mp.lu_solve(design.T * precision * design, design.T * precision * values)
    residual = values - design * beta
    return (residual.T * precision * residual)[0], scale, beta, residual


def ray_score(reference, reference_covariance, new, new_covariance, sign):
    covariance = block_diagonal(reference_covariance, new_covariance)
    lower, upper = mp.mpf("-13.815510557964274"), mp.mpf("13.815510557964274")
    grid = [lower + (upper - lower) * index / 240 for index in range(241)]
    values = [_ray_objective(reference, new, covariance, sign, value)[0] for value in grid]
    best = min(range(len(grid)), key=values.__getitem__)
    left = grid[max(0, best - 1)]
    right = grid[min(len(grid) - 1, best + 1)]
    golden = (mp.sqrt(5) - 1) / 2
    x1 = right - golden * (right - left)
    x2 = left + golden * (right - left)
    f1 = _ray_objective(reference, new, covariance, sign, x1)[0]
    f2 = _ray_objective(reference, new, covariance, sign, x2)[0]
    for _ in range(100):
        if f1 > f2:
            left, x1, f1 = x1, x2, f2
            x2 = left + golden * (right - left)
            f2 = _ray_objective(reference, new, covariance, sign, x2)[0]
        else:
            right, x2, f2 = x2, x1, f1
            x1 = right - golden * (right - left)
            f1 = _ray_objective(reference, new, covariance, sign, x1)[0]
    result = _ray_objective(reference, new, covariance, sign, (left + right) / 2)
    chi2, scale, beta, residual = result
    p = mp.erfc(mp.sqrt(chi2 / 2))
    return {
        "ray_sign": sign,
        "best_positive_scale": mp.nstr(scale, 20),
        "chi_square": mp.nstr(chi2, 20),
        "dof": 1,
        "survival_p": mp.nstr(p, 20),
        "latent_reference_re_im": [mp.nstr(value, 20) for value in beta],
        "joint_residual": [mp.nstr(value, 20) for value in residual],
    }


def score(runs, manifest, reference_score, *, cutoff=9, dps=60):
    mp.mp.dps = dps
    if not manifest.get("production_authorized"):
        raise ValueError("manifest is not authorized")
    point, covariance = combine_runs(runs, manifest)
    baselines, convergence = [], []
    for _, matrix in N60_CHILD_DESIGNS:
        baseline = continuum_harmonic(matrix, cutoff=cutoff, dps=dps)
        previous = continuum_harmonic(matrix, cutoff=cutoff-2, dps=dps)
        baselines.append(baseline)
        convergence.append(abs(baseline - previous))
    baseline_vector = mp.matrix([coordinate for value in baselines for coordinate in (mp.re(value), mp.im(value))])
    delta60 = point - baseline_vector
    c60, cov60, w60, denominator60 = annihilator_coefficient(delta60, covariance, N60_CHILD_DESIGNS, dps=dps)

    delta112 = mp.matrix([mp.mpf(value) for value in reference_score["delta_H4_re_im"]])
    cov112_full = mp.matrix([[mp.mpf(value) for value in row] for row in reference_score["full_covariance_6x6"]])
    from rho_child_primitive_h4_mc import CHILD_DESIGNS
    c112, cov112, w112, denominator112 = annihilator_coefficient(delta112, cov112_full, CHILD_DESIGNS, dps=dps)

    preserve = ray_score(c112, cov112, c60, cov60, +1)
    flip = ray_score(c112, cov112, c60, cov60, -1)
    alpha = mp.mpf(str(manifest["decision_alpha"]))
    preserve_pass = mp.mpf(preserve["survival_p"]) >= alpha
    flip_pass = mp.mpf(flip["survival_p"]) >= alpha
    decision = (
        "E4_squared_phase_preserved" if preserve_pass and not flip_pass else
        "signed_Pell_phase_flip" if flip_pass and not preserve_pass else
        "phase_incoherent_both_rejected" if not preserve_pass and not flip_pass else
        "phase_unresolved_both_survive"
    )
    ratio = mp.mpc(c60[0], c60[1]) / mp.mpc(c112[0], c112[1])
    null_chi2 = (c60.T * cov60**-1 * c60)[0]
    null_p = mp.exp(-null_chi2 / 2)
    return {
        "schema": "matching-one/rho-child-opposite-pell-phase-score/v1",
        "status": "frozen_opposite_pell_reveal",
        "decision": decision,
        "decision_alpha": manifest["decision_alpha"],
        "primary": "positive-real ray c60=a*c112 with a>0 (E4-squared phase preservation)",
        "opponent": "negative-real ray c60=-a*c112 with a>0 (signed-Pell phase flip)",
        "c112_re_im": [mp.nstr(value, 20) for value in c112],
        "c60_re_im": [mp.nstr(value, 20) for value in c60],
        "c112_covariance_2x2": [[mp.nstr(cov112[i,j], 20) for j in range(2)] for i in range(2)],
        "c60_covariance_2x2": [[mp.nstr(cov60[i,j], 20) for j in range(2)] for i in range(2)],
        "raw_complex_ratio_re_im": [mp.nstr(mp.re(ratio), 20), mp.nstr(mp.im(ratio), 20)],
        "raw_phase_difference_rad": mp.nstr(mp.arg(ratio), 20),
        "phase_preserving_ray": preserve,
        "phase_flipping_ray": flip,
        "N60_zero": {"chi_square": mp.nstr(null_chi2, 20), "dof": 2, "survival_p": mp.nstr(null_p, 20)},
        "N60_delta_H4_re_im": [mp.nstr(value, 20) for value in delta60],
        "N60_full_covariance_6x6": [[mp.nstr(covariance[i,j], 20) for j in range(6)] for i in range(6)],
        "N60_continuum_baseline": [[mp.nstr(mp.re(v), 30), mp.nstr(mp.im(v), 30)] for v in baselines],
        "baseline_cutoff": cutoff,
        "cutoff_minus_2_change_abs": [mp.nstr(value, 8) for value in convergence],
        "annihilator_checks": {
            "N112_w_dot_E4sq_re_im": [mp.nstr(mp.re(denominator112), 20), mp.nstr(mp.im(denominator112), 20)],
            "N60_w_dot_E4sq_re_im": [mp.nstr(mp.re(denominator60), 20), mp.nstr(mp.im(denominator60), 20)],
            "N112_w": [[mp.nstr(mp.re(v), 20), mp.nstr(mp.im(v), 20)] for v in w112],
            "N60_w": [[mp.nstr(mp.re(v), 20), mp.nstr(mp.im(v), 20)] for v in w60],
        },
        "claim_boundary": "This tests the completion phase of one typed primitive-line row; it is not a new pure-harmonic vote.",
    }


def render(result):
    return "\n".join([
        "# Opposite-Pell N60 normalization-free phase reveal", "",
        f"Decision: **{result['decision']}** at alpha `{result['decision_alpha']}`.", "",
        f"- c112: `{result['c112_re_im']}`", f"- c60: `{result['c60_re_im']}`",
        f"- raw c60/c112: `{result['raw_complex_ratio_re_im']}`; phase `{result['raw_phase_difference_rad']}` rad", "",
        f"- phase-preserving E4-squared ray: `{result['phase_preserving_ray']['chi_square']}/1`, p `{result['phase_preserving_ray']['survival_p']}`",
        f"- signed-Pell flipping ray: `{result['phase_flipping_ray']['chi_square']}/1`, p `{result['phase_flipping_ray']['survival_p']}`",
        f"- N60 coefficient zero: `{result['N60_zero']['chi_square']}/2`, p `{result['N60_zero']['survival_p']}`", "",
        result["claim_boundary"], "",
    ])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path, nargs=3)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("reference_score", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args()
    result = score(
        [json.loads(path.read_text()) for path in args.run],
        json.loads(args.manifest.read_text()),
        json.loads(args.reference_score.read_text()),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    args.markdown.write_text(render(result))
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

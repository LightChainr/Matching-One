#!/usr/bin/env python3
"""One counterfactual: delete stationary current, keep width-eight sources."""
from __future__ import annotations

import os
for _key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_key] = "1"

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import linalg, sparse
from scipy.optimize import brentq
from threadpoolctl import threadpool_limits

from p398_width8_projected_memory import decomposition, exponential_row, series
from p398_width8_source_spectrum import complex_display

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "analysis/p398_width8_reversible_current_control.json"
OUT = ROOT / "results/p398-width8-reversible-current-control/latest.json"


def maximum(array):
    return float(np.max(np.abs(array), initial=0))


def propagation(matrix, source, lags, reversible, tolerance):
    source = source / linalg.norm(source)
    if reversible:
        values, vectors = linalg.eigh(matrix)
        residues = np.abs(vectors.conj().T @ source) ** 2
    else:
        values, residues = exponential_row(matrix, source.conj(), source)
    visible = np.abs(residues) > tolerance * np.max(np.abs(residues))
    first = np.flatnonzero(visible)[np.argmin(values[visible].real)]
    omega = complex(np.vdot(source, matrix @ source))
    second = complex(np.vdot(source, matrix @ matrix @ source))
    complement = linalg.null_space(source.conj()[None, :])
    back = source.conj() @ matrix @ complement
    force = complement.conj().T @ matrix @ source
    hidden = complement.conj().T @ matrix @ complement
    integral = complex(back @ linalg.solve(hidden, force))
    first_moment = complex(back @ linalg.solve(hidden, linalg.solve(hidden, force)))
    kernel_values, kernel_residues = exponential_row(hidden, back, force)
    u = lambda t: series(values, residues, t)
    row = {
        "lowest_source_visible_mass": float(values[first].real),
        "lowest_source_visible_mass_imaginary": float(np.imag(values[first])),
        "lowest_source_visible_residue_re_im": complex_display(residues[first]),
        "initial_decay_re_im": complex_display(omega),
        "initial_second_derivative_re_im": complex_display(second),
        "initial_log_curvature": float((second - omega**2).real),
        "projected_memory_integral": float(integral.real),
        "projected_memory_first_moment": float(first_moment.real),
        "projected_memory_mean_time": float((first_moment / integral).real),
        "samples": [{"s": t, "u": u(t),
                     "instantaneous_decay": -series(values, residues, t, 1) / u(t),
                     "projected_memory": series(kernel_values, kernel_residues, t)}
                    for t in lags],
        "mass_spectrum_re_im": complex_display(values),
        "normalized_residues_re_im": complex_display(residues),
        "visible_mode_count_at_frozen_tolerance": int(np.sum(visible)),
        "maximum_imaginary_correlation_on_fixed_grid": max(
            maximum(np.sum(residues * np.exp(-values*t)).imag) for t in lags),
    }
    return row, u


def build_result():
    protocol = json.loads(PROTOCOL.read_text())
    lags = protocol["distances"]
    tolerance = protocol["visibility_relative_numerical_tolerance"]
    states, mass, f, _, q, h, sources, pi, weight, dual_basis, phase = decomposition()

    # Backward row generators. pi is a column but stationary on the left.
    g = -mass
    g_star = sparse.diags(1/pi) @ g.T @ sparse.diags(pi)
    s = (g + g_star) / 2
    pi_s = sparse.diags(pi) @ s
    offdiag = s - sparse.diags(s.diagonal())
    checks = {
        "states": len(states), "minimum_stationary_probability": float(pi.min()),
        "stationary_sum": float(pi.sum()),
        "original_rowsum_max": maximum(g @ np.ones(len(pi))),
        "reversible_rowsum_max": maximum(s @ np.ones(len(pi))),
        "reversible_stationarity_max": maximum(s.T @ pi),
        "detailed_balance_max": maximum((pi_s - pi_s.T).data),
        "unchanged_diagonal_max": maximum(s.diagonal() - g.diagonal()),
        "minimum_reversible_offdiagonal": float(offdiag.data.min()),
        "original_stationary_current_max": maximum((sparse.diags(pi) @ g - g.T @ sparse.diags(pi)).data),
    }
    assert checks["minimum_stationary_probability"] > 0
    assert checks["minimum_reversible_offdiagonal"] >= 0
    assert checks["reversible_rowsum_max"] < 1e-9
    assert checks["detailed_balance_max"] < 1e-12
    assert checks["unchanged_diagonal_max"] == 0

    rows = []
    functions = {"original": [], "reversible": []}
    for sign, ray_name in ((-1, "minus"), (1, "plus")):
        v = dual_basis[:, :93] if sign < 0 else dual_basis[:, 93:]
        ray_mass = v.conj().T @ h @ v
        ray_weight = v.conj().T @ weight @ v
        upper = linalg.cholesky(ray_weight, lower=False)
        whitened = upper @ ray_mass @ linalg.inv(upper)
        symmetric = (whitened + whitened.conj().T) / 2
        antisymmetric = (whitened - whitened.conj().T) / 2
        physical_source = (f[:, 0] + sign * phase * f[:, 1]) / np.sqrt(2)
        source = upper @ v.conj().T @ sources @ np.array([1, sign * phase]) / np.sqrt(2)
        variance = float(np.vdot(source, source).real)
        normalized = source / np.sqrt(variance)
        original, original_u = propagation(whitened, source, lags, False, tolerance)
        reversible, reversible_u = propagation(symmetric, source, lags, True, tolerance)
        functions["original"].append(original_u)
        functions["reversible"].append(reversible_u)
        current_force = float(linalg.norm(antisymmetric @ normalized)**2)
        current_mean = complex(np.vdot(normalized, antisymmetric @ normalized))
        curvature_difference = reversible["initial_log_curvature"] - original["initial_log_curvature"]
        predicted_difference = current_force - abs(current_mean)**2
        direct_second = float(np.vdot(physical_source, pi * (s @ (s @ physical_source))).real / variance)
        ray_checks = {
            "variance_matches_full_state": abs(variance - np.vdot(physical_source, pi * physical_source).real),
            "time_reverse_second_derivative_matches_full_state": abs(direct_second - reversible["initial_second_derivative_re_im"][0]),
            "real_initial_decay_difference": reversible["initial_decay_re_im"][0] - original["initial_decay_re_im"][0],
            "current_force_variance": current_force,
            "current_mean_re_im": complex_display(current_mean),
            "observed_initial_log_curvature_increase": curvature_difference,
            "predicted_initial_log_curvature_increase": predicted_difference,
            "curvature_identity_absolute_error": abs(curvature_difference - predicted_difference),
        }
        assert ray_checks["curvature_identity_absolute_error"] < 1e-10
        assert ray_checks["time_reverse_second_derivative_matches_full_state"] < 1e-9
        rows.append({"ray": ray_name, "sign": sign, "variance": variance,
                     "original": original, "reversible": reversible, "checks": ray_checks})

    comparisons = {}
    bracket = protocol["crossing_bracket"]
    for name, (minus_u, plus_u) in functions.items():
        difference = lambda t: minus_u(t) - plus_u(t)
        crossing = float(brentq(difference, *bracket)) if difference(bracket[0]) * difference(bracket[1]) < 0 else None
        comparisons[name] = {
            "crossing_in_frozen_bracket": crossing,
            "difference_at_bracket_ends": [difference(t) for t in bracket],
            "plus_to_minus_at_fixed_lags": [{"s": t, "u_plus_over_u_minus": plus_u(t)/minus_u(t)} for t in lags],
            "plus_starts_faster": rows[1][name]["initial_decay_re_im"][0] > rows[0][name]["initial_decay_re_im"][0],
            "plus_has_lower_visible_mass": rows[1][name]["lowest_source_visible_mass"] < rows[0][name]["lowest_source_visible_mass"],
        }
    return {
        "schema": "matching-one/p398-width8-reversible-current-control/v1",
        "definition_commit": "5e47bdb6", "protocol": str(PROTOCOL.relative_to(ROOT)),
        "generator_checks": checks, "ray_rows": rows, "comparison": comparisons,
        "interpretation_boundary": "A single finite-process current-deletion counterfactual; not a local square-bond row word, continuum field count, or Jordan test.",
        "arithmetic": "float64 deterministic finite matrices; one BLAS thread; no fit, MC, or parameter scan",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path, default=OUT)
    args = parser.parse_args()
    with threadpool_limits(limits=1):
        value = build_result()
    inputs = (PROTOCOL, Path(__file__), ROOT/"scripts/p398_width8_projected_memory.py",
              ROOT/"scripts/p398_width8_source_spectrum.py")
    value["input_sha256"] = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(value, indent=2) + "\n")
    for row in value["ray_rows"]:
        print(row["ray"], "variance", row["variance"])
        for name in ("original", "reversible"):
            print(name, {k: row[name][k] for k in ("lowest_source_visible_mass", "initial_decay_re_im", "initial_log_curvature", "projected_memory_integral")})
        print("curvature", row["checks"])
    print("comparison", json.dumps(value["comparison"], indent=2))


if __name__ == "__main__":
    main()

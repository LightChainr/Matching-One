#!/usr/bin/env python3
"""Frozen spectral-energy/variance readout, not unconditional high-pass means."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.stats import norm
from p437_high_pass_mc import COEFFICIENTS, NAMES


def score(directory):
    run = json.loads((directory / "run.json").read_text())
    path = directory / "batches.csv"
    if hashlib.sha256(path.read_bytes()).hexdigest() != run["batch_sha256"]:
        raise ValueError("batch hash mismatch")
    with path.open() as stream:
        rows = list(csv.DictReader(stream))
    x = np.array([[float(row[name]) / int(row["samples"]) for name in NAMES] for row in rows])
    mean = x.mean(axis=0)
    covariance = np.cov(x, rowvar=False, ddof=1) / len(rows)
    coefficients = np.array(COEFFICIENTS, dtype=float)
    transform = np.zeros((8, len(NAMES)))
    for j in range(6):
        transform[j, j:36:6] = coefficients
    transform[6, -2] = 1
    transform[7, -1] = 1
    point, cov = transform @ mean, transform @ covariance @ transform.T
    names = ("raw_hp_F_re_zero", "raw_hp_F_im_zero", "spectral_energy", "energy_im_zero",
             "euler_energy_zero", "degree5_energy", "unfiltered_variance", "euler_unfiltered_variance")
    summaries = {}
    for j, name in enumerate(names):
        target = 9765 / 32768 if name == "degree5_energy" else 0
        se = float(np.sqrt(cov[j, j]))
        z = float((point[j] - target) / se) if se else None
        summaries[name] = {"mean": float(point[j]), "se": se, "exact_control_target": target if j != 2 and j < 6 else None,
                           "z_vs_target": z, "p_two_sided": float(2 * norm.sf(abs(z))) if z is not None else None}
    energy, variance = point[2], point[6]
    hp_snr2 = float(energy ** 2 / cov[2, 2])
    raw_snr2 = float(variance ** 2 / cov[6, 6])
    cpu = run["cpu_seconds"]
    raw_cpu = run["comparator_classification_cpu_seconds"]
    gradient = np.array([1 / variance, -energy / variance ** 2])
    fraction_se = float(np.sqrt(gradient @ cov[np.ix_([2, 6], [2, 6])] @ gradient))
    count = run["run"]["samples"]
    universal_energy_bound = 1 / 9  # F is 0 or a sixth root of unity divided by 3.
    optimistic_maximum_snr2 = float(universal_energy_bound ** 2 / cov[2, 2])
    return {"schema": "matching-one/p437-high-pass-score/v1", "samples": run["run"]["samples"],
            "coordinate_order": NAMES, "level_point": mean.tolist(), "full_covariance_of_mean": covariance.tolist(),
            "readout_order": names, "readout_covariance": cov.tolist(), "readouts": summaries,
            "high_degree_weighted_energy_fraction": {"point": float(energy / variance), "delta_se": fraction_se,
                "population_range": [0, 1], "finite_pilot_estimate_not_clipped": True},
            "variance_inflation_vs_unfiltered": float(cov[2, 2] / cov[6, 6]),
            "primary_p_one_sided_positive_energy": float(norm.sf(summaries["spectral_energy"]["z_vs_target"])),
            "pilot_variance_power_envelope": {
                "exact_universal_energy_upper_bound": "1/9",
                "reason": "|F|<=1/3 and 0<=H<=I imply 0<=<F,HF><=VarF<=1/9",
                "five_sigma_samples_even_at_universal_maximum": float(count * 25 / optimistic_maximum_snr2),
                "five_sigma_samples_if_fraction_of_observed_VarF": {
                    str(fraction): float(count * 25 * cov[2, 2] / (fraction * variance) ** 2)
                    for fraction in (1, 0.1, 0.01)},
                "best_case_relative_efficiency_at_universal_maximum": float((optimistic_maximum_snr2 / cpu) / (raw_snr2 / raw_cpu)),
                "degree5_known_signal_five_sigma_sample_projection": float(count * 25 * cov[5, 5] / (9765 / 32768) ** 2),
                "boundary": "Fixed estimator and pilot variance extrapolation, not a universal impossibility bound on high-pass acquisition; no production approval."},
            "efficiency": {"highpass_snr_squared": hp_snr2, "unfiltered_snr_squared": raw_snr2,
                "highpass_cpu_seconds": cpu, "comparator_classification_cpu_seconds": raw_cpu,
                "highpass_snr_squared_per_cpu_second": hp_snr2 / cpu,
                "unfiltered_snr_squared_per_classification_cpu_second": raw_snr2 / raw_cpu,
                "point_relative_efficiency": (hp_snr2 / cpu) / (raw_snr2 / raw_cpu),
                "cpu_boundary": "baseline excludes RNG/control/output overhead; highpass uses all worker CPU, conservative relative cost"},
            "five_sigma_sample_projection_from_noisy_point": float(run["run"]["samples"] * 25 / hp_snr2) if energy > 0 else None,
            "decision": "resolved_high_order_energy_at_1pct" if summaries["spectral_energy"]["z_vs_target"] > norm.isf(.01) else "underpowered_high_order_energy",
            "stop": "fixed20k_complete_no_automatic_extension",
            "boundary": "Self-source spectral energy is not an independent source/observer bridge; estimator variance is not filtered population energy. Raw HP mean is exact zero. No field name or state count."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(score(args.directory), indent=2, sort_keys=True, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()

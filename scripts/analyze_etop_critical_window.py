#!/usr/bin/env python3
"""Locate the N100 clock-quotient redistribution using archived histograms.

This is a post-reveal finite-N decomposition, not a new hypothesis test or
Monte Carlo run. It consumes immutable PR484 blobs; fetch that branch first
if its commit is absent. No branch checkout or merge is required.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import io
import json
from math import comb
from pathlib import Path
import platform
import subprocess
import sys
import time

import numpy as np
import scipy
from scipy.special import betainc


ROOT = Path(__file__).resolve().parents[1]
SOURCE_COMMIT = "894b3d800c5aeaad3dd8b0f893b6f17d85d234c6"
SOURCE_ROOT = "results/etop-n100-three-modulus"


def archived(path, commit, inputs):
    raw = subprocess.check_output(["git", "show", f"{commit}:{path}"], cwd=ROOT)
    inputs.append({"commit": commit, "path": path, "bytes": len(raw),
                   "sha256": hashlib.sha256(raw).hexdigest()})
    return raw.decode()


def kernel(n, pc, scale, lo, hi, power):
    """Integral z**power B_{l,n}(p) dp, by incomplete beta identities.

    The formula is analytic; evaluation uses double-precision betainc,
    not an interval-arithmetic exact numeric certificate or quadrature.
    """
    l = np.arange(n + 1, dtype=float)
    output = np.zeros(n + 1)
    factor = np.full(n + 1, 1.0 / (n + 1))
    for m in range(power + 1):
        if m:
            factor *= (l + m) / (n + 1 + m)
        a, b = l + m + 1, n - l + 1
        # Use the smaller-tail difference to limit cancellation near 1.
        interval = np.where(lo + hi > 1,
                            betainc(b, a, 1-lo) - betainc(b, a, 1-hi),
                            betainc(a, b, hi) - betainc(a, b, lo))
        output += comb(power, m) * (-pc)**(power-m) * factor * interval
    return scale**power * output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-commit", default=SOURCE_COMMIT)
    parser.add_argument("--output-dir", type=Path,
                        default=ROOT / "results/etop-critical-window-n100")
    args = parser.parse_args()
    out = args.output_dir
    if any((out / name).exists() for name in ("latest.json", "REPORT.md")):
        raise SystemExit("Use a fresh output directory; existing results are not overwritten.")
    start = time.perf_counter()
    inputs = []
    source = json.loads(archived(f"{SOURCE_ROOT}/score.json", args.source_commit, inputs))
    contract = source["contract"]
    n, batches = contract["area"], contract["batches"]
    per_batch = contract["samples_per_shape_pair"] // batches
    counts = np.zeros((3, batches, 2, n + 1))
    for s, shape in enumerate(contract["shapes"]):
        raw = archived(f"{SOURCE_ROOT}/raw/{shape['name']}.hist.csv",
                       args.source_commit, inputs)
        delta = float(Fraction(shape["delta_cos4"]))
        totals = np.zeros((batches, 2, 2), dtype=np.int64)
        for row in csv.DictReader(io.StringIO(raw)):
            batch, k, value = int(row["batch"]), int(row["k"]), int(row["count"])
            if row["orientation"] not in ("first", "second") or row["kind"] not in ("minus", "plus"):
                raise ValueError("Unknown archived orientation/threshold convention")
            orientation = int(row["orientation"] == "second")
            event = int(row["kind"] == "plus")
            if int(row["n"]) != n or int(row["samples"]) != per_batch:
                raise ValueError("Unexpected archived area/batch size")
            totals[batch, orientation, event] += value
            counts[s, batch, event, k] += (1 - 2*orientation) * value / (per_batch*delta)
        if not np.all(totals == per_batch):
            raise ValueError("Incomplete marginal histogram")

    # A=F1+F2-1 and E=1+F2-F1; constants cancel in each orientation contrast.
    # CDF-tail coefficients become Bernstein coefficients by cumulative sum.
    fields = np.stack((counts[:, :, 0]+counts[:, :, 1],
                       counts[:, :, 1]-counts[:, :, 0]), axis=2)
    coeff = np.cumsum(fields, axis=-1)
    clocks = fields[:, :, 0] @ (np.arange(n+1)/(2*(n+1)))
    d, u = coeff[1]-coeff[0], coeff[2]-coeff[0]
    dc, uc = clocks[1]-clocks[0], clocks[2]-clocks[0]
    ratio = uc.mean()/dc.mean()
    ratio_loo = (uc.sum()-uc)/(dc.sum()-dc)
    residual = u.mean(axis=0)-ratio*d.mean(axis=0)
    loo = ((u.sum(axis=0)-u)-ratio_loo[:, None, None]*(d.sum(axis=0)-d))/(batches-1)

    pc, scale = contract["fixed_p"], n**(3/8)
    regions = [{"id": "full", "lo": 0.0, "hi": 1.0}]
    for width in (0.5, 1.0, 1.5):
        lo, hi = max(0.0, pc-width/scale), min(1.0, pc+width/scale)
        for name, a, b in (("lower", 0.0, lo), ("core", lo, hi), ("upper", hi, 1.0)):
            regions.append({"id": f"w{width:g}_{name}", "half_width_z": width,
                            "region": name, "lo": a, "hi": b})
    kernels = np.stack([kernel(n, pc, scale, row["lo"], row["hi"], j)
                        for row in regions for j in range(3)], axis=-1)
    # F1=(A-E)/2, F2=(A+E)/2. These redundant coordinates are kept together,
    # not counted as independent evidence.
    transform = np.array([[1, 0], [0, 1], [.5, -.5], [.5, .5]])
    names = ("A", "E", "F1", "F2")
    values = (transform @ residual @ kernels).reshape(-1)
    leave_one = np.einsum("fe,bek,kl->bfl", transform, loo, kernels).reshape(batches, -1)
    labels = [f"{field}:{row['id']}:j{j}" for field in names for row in regions for j in range(3)]
    centered = leave_one-leave_one.mean(axis=0)
    covariance = (batches-1)/batches * centered.T @ centered
    se = np.sqrt(np.maximum(0, np.diag(covariance)))
    ratio_se = np.sqrt((batches-1)/batches*np.sum((ratio_loo-ratio_loo.mean())**2))
    # Algebraic additivity is recorded as numerical closure, not another test
    # campaign. The exact odd full area is zero by definition of the clock.
    view = values.reshape(4, len(regions), 3)
    closure = {f"w{w:g}": float(np.max(np.abs(view[:, 1+3*i:4+3*i].sum(axis=1)-view[:, 0])))
               for i, w in enumerate((.5, 1., 1.5))}
    dipole_indices = [labels.index(f"A:{name}:j1")
                      for name in ("full", "w1_lower", "w1_core", "w1_upper")]
    result = {
        "schema": "matching-one.etop-critical-window.v1",
        "status": "completed_post_reveal_existing_data_analysis",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_pr": 484, "source_commit": args.source_commit, "inputs": inputs,
        "dependency_group": "N100_three_modulus_seed20260831125401_offset267100000000",
        "new_samples": 0, "batches": batches,
        "definitions": {"contrast": "D=Y(4i)-Y(2i), U=Y(1/2+i)-Y(2i), R=U-r_C D",
            "clock": "C=(E[K1]+E[K2])/(2(N+1)), r_C=U_C/D_C",
            "coordinate": "z=N^(3/8)*(p-p_ref); finite-N convention, no exponent fit",
            "moment": "integral_lo^hi z^j R_field(p) dp, j=0,1,2; integration is dp, not dz",
            "births": "K1=K_minus, K2=K_plus; R_F1=(R_A-R_E)/2, R_F2=(R_A+R_E)/2",
            "primary_descriptive_window": "|z|<=1; half-widths .5,1.5 are sensitivity descriptions, not extra tests",
            "uncertainty": "delete the same batch across all shapes/orientations, recomputing r_C in every deletion; full delete-one covariance"},
        "area": n, "p_ref": pc, "scale": scale, "regions": regions,
        "clock_ratio": float(ratio), "clock_ratio_se": float(ratio_se),
        "clock_ratio_leave_one": ratio_loo.tolist(),
        "labels": labels, "mean": values.tolist(), "se": se.tolist(),
        "covariance": covariance.tolist(), "delete_one": leave_one.tolist(),
        "primary_odd_dipole": {"labels": [labels[i] for i in dipole_indices],
            "mean": values[dipole_indices].tolist(), "se": se[dipole_indices].tolist(),
            "covariance": covariance[np.ix_(dipole_indices, dipole_indices)].tolist()},
        "identities": ["full integral R_A = 0 by clock calibration",
                       "lower + core + upper = full at each width and moment",
                       "F1+F2=A and F2-F1=E; birth rows are dependent readouts"],
        "numerical_closure": {"odd_full_area": float(view[0, 0, 0]), "window_additivity_max_abs": closure},
        "boundaries": ["A single finite N localizes response; it does not determine a scaling limit, thin-geometry cause or continuum operator.",
            "All windows/moments/birth splits reuse PR484's same random block; no p-values or independent evidence counts are combined.",
            "Signed dipoles can cancel. They are not nonnegative mass fractions or variance shares.",
            "The analysis follows disclosure of the full-p curve and is retrospective; widths are not claimed prospectively frozen.",
            "Incomplete-beta formulas are evaluated in float64, not certified exact arithmetic."],
        "execution": {"command": " ".join(sys.argv), "python": platform.python_version(),
            "platform": platform.platform(), "machine": platform.machine(),
            "numpy": np.__version__, "scipy": scipy.__version__,
            "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "elapsed_seconds": time.perf_counter()-start}}
    lines = ["# N100: locate the thermal redistribution before assigning its mechanism", "",
             "This consumes PR484's archived 200 aligned batches. No new sampling, curve fitting, model vote or transfer-engine run is added.", "",
             "## Primary descriptive window and odd dipole", "",
             f"The center window is p in [{regions[5]['lo']:.9f}, {regions[5]['hi']:.9f}], or |z|<=1 with z=N^(3/8)(p-p_ref). The exponent is a coordinate convention, not a new estimate.", "",
             "| region | integral z R_A dp | delete-one SE |", "|---|---:|---:|"]
    for i in dipole_indices:
        lines.append(f"| {labels[i].split(':')[1]} | {values[i]:.10g} | {se[i]:.6g} |")
    lines += ["", "These are additive **signed** contributions, not shares of positive signal mass. They retain all cross-region covariance and the uncertainty in the same-stream clock ratio.", "",
              "## First/second activation in the center window", "",
              "| readout | integral R dp | delete-one SE |", "|---|---:|---:|"]
    for field in names:
        i = labels.index(f"{field}:w1_core:j0")
        lines.append(f"| {field} | {values[i]:.10g} | {se[i]:.6g} |")
    lines += ["", "F1=(A-E)/2 and F2=(A+E)/2; the birth rows are a change of coordinates, not two independent replications.", "",
              "## Width sensitivity of the odd dipole", "",
              "| z half-width | lower | core | upper |", "|---|---:|---:|---:|"]
    for width in (.5, 1., 1.5):
        items = [labels.index(f"A:w{width:g}_{name}:j1") for name in ("lower", "core", "upper")]
        lines.append(f"| {width:g} | " + " | ".join(f"{values[i]:.8g} ± {se[i]:.3g}" for i in items) + " |")
    lines += ["", "## Interpretation and next discriminant", "",
              "This is a localization of the existing finite-N response. A nonzero center-window response does not identify a critical field; a larger outer contribution does not establish a thin-geometry mechanism. The next scale comparison should preserve the three homothetic shapes and these named regions, then ask whether the distribution contracts in p or remains at an off-critical location. PR484's N400 design is an existing acquisition option, not a new approval gate.", "",
              "No transport test is repeated: PR484 already contains the failed joint finite-Jacobian transport comparison. A-only quantile transport is not a discriminating model.", "",
              "## Reproduction and dependence", "",
              f"Source: open PR484, `{args.source_commit}`; all three shape pairs share the same 2,000,000 permutation counters, seed 20260831125401, offset 267100000000. Histograms and their SHA256 values are listed in latest.json.", "",
              f"Clock ratio {ratio:.10g} ± {ratio_se:.6g}. Each delete-one removes one aligned batch from all three shape pairs and refits this ratio. JSON stores all 200 vectors, the complete covariance, and window bounds; singular directions and overlapping windows are preserved.", "",
              "The Bernstein integrals use incomplete-beta identities evaluated in float64, not quadrature and not a rigorous numeric certificate. This follows disclosure of the full-p curve; all window claims are retrospective.", "",
              "```bash", "git fetch origin analysis/etop-modulus-survivors-20260831",
              "python3 scripts/analyze_etop_critical_window.py --output-dir /path/to/fresh-output", "```", ""]
    out.mkdir(parents=True, exist_ok=True)
    (out / "latest.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    (out / "REPORT.md").write_text("\n".join(lines))
    print(json.dumps({"primary_odd_dipole": result["primary_odd_dipole"],
                      "clock_ratio": ratio, "closure": result["numerical_closure"],
                      "seconds": result["execution"]["elapsed_seconds"],
                      "output": str(out)}, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Current orientation in the single fixed pair: ray source and pair density."""
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
from scipy import linalg
from threadpoolctl import threadpool_limits

from p398_width8_projected_memory import decomposition
from p398_width8_source_spectrum import complex_display

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "analysis/p398_width8_pair_density_current.json"
OUT = ROOT / "results/p398-width8-pair-density-current/latest.json"


def pair_rows(mass, e, partner, lags):
    pair = np.column_stack((e, partner))
    values, vectors = linalg.eig(mass)
    left = pair.conj().T @ vectors
    right = linalg.solve(vectors, pair)
    rows = []
    for s in lags:
        correlation = (left * np.exp(-s*values)[None, :]) @ right
        anti = correlation - correlation.conj().T
        rows.append({"s": s, "cross_e_partner_re_im": complex_display(correlation[0,1]),
                     "cross_partner_e_re_im": complex_display(correlation[1,0]),
                     "anti_cross_re_im": complex_display(anti[0,1]),
                     "anti_cross_absolute": float(abs(anti[0,1]))})
    current = -(mass - mass.conj().T) / 2
    derivative = 2 * np.vdot(e, current @ partner)
    return {"anti_derivative_at_zero_re_im": complex_display(derivative),
            "anti_derivative_at_zero_absolute": float(abs(derivative)),
            "samples": rows}


def build_result():
    protocol = json.loads(PROTOCOL.read_text())
    _, _, _, t2, q, h, sources, _, weight, dual_basis, phase = decomposition()
    reduced_t2 = q.conj().T @ t2
    rows = []
    for sign, name in ((-1, "minus"), (1, "plus")):
        v = dual_basis[:, :93] if sign < 0 else dual_basis[:, 93:]
        upper = linalg.cholesky(v.conj().T @ weight @ v, lower=False)
        mass = upper @ (v.conj().T @ h @ v) @ linalg.inv(upper)
        psi = upper @ v.conj().T @ sources @ np.array([1, sign*phase]) / np.sqrt(2)
        t = upper @ v.conj().T @ reduced_t2
        source_variance = float(np.vdot(psi, psi).real)
        t_variance = float(np.vdot(t, t).real)
        e = psi / np.sqrt(source_variance)
        coefficient = np.vdot(e, t)
        orthogonal = t - e*coefficient
        orthogonal_variance = float(np.vdot(orthogonal, orthogonal).real)
        eta = orthogonal / np.sqrt(orthogonal_variance)
        current = -(mass - mass.conj().T) / 2
        j_source = current @ e
        norm_j = float(np.vdot(j_source, j_source).real)
        projection_j = np.vdot(eta, j_source)
        raw = pair_rows(mass, e, t/np.sqrt(t_variance), protocol["distances"])
        orth = pair_rows(mass, e, eta, protocol["distances"])
        rows.append({
            "ray": name, "sign": sign, "source_variance": source_variance,
            "projected_T2_variance": t_variance,
            "source_T2_inner_product_re_im": complex_display(coefficient),
            "orthogonal_T2_variance": orthogonal_variance,
            "orthogonal_fraction_of_T2_variance": orthogonal_variance/t_variance,
            "current_source_squared_norm": norm_j,
            "current_source_along_orthogonal_T2_re_im": complex_display(projection_j),
            "current_source_squared_norm_captured": float(abs(projection_j)**2),
            "current_source_norm_fraction_captured": float(abs(projection_j)**2/norm_j),
            "raw_normalized_pair": raw, "orthonormal_pair": orth,
            "time_reversed_control": "negative of every anti-cross and derivative, by adjoint identity",
            "reversible_control": "identically zero, by self-adjointness; not numerically refitted",
        })
    return {"schema": protocol["schema"], "protocol": str(PROTOCOL.relative_to(ROOT)),
            "parent": protocol["parent"], "ray_rows": rows,
            "arithmetic": "same inherited integer process and numerical stationary law, one BLAS thread, no new samples",
            "direction_gauge": "physical T2 site origin and fixed psi phases; eta is only orthogonalized/divided by a positive real norm"}


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
        print(row["ray"], "captured", row["current_source_norm_fraction_captured"],
              "derivative", row["orthonormal_pair"]["anti_derivative_at_zero_re_im"])
        for sample in row["orthonormal_pair"]["samples"]:
            print(sample["s"], sample["anti_cross_re_im"])


if __name__ == "__main__":
    main()

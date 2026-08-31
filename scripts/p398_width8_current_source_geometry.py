#!/usr/bin/env python3
"""Locate the unique current-source direction in existing named geometry."""
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
from scipy.optimize import brentq
from threadpoolctl import threadpool_limits

from p398_width8_geometric_compression import features, response
from p398_width8_pair_density_current import pair_rows
from p398_width8_projected_memory import decomposition, series
from p398_width8_source_spectrum import complex_display, kreweras

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "analysis/p398_width8_current_source_geometry.json"
GEOMETRY = ROOT / "results/p398-width8-t4-post-reveal/latest.json"
FULL = ROOT / "results/p398-width8-reversible-current-control/latest.json"
OUT = ROOT / "results/p398-width8-current-source-geometry/latest.json"


def build_result():
    protocol = json.loads(PROTOCOL.read_text())
    geometry = json.loads(GEOMETRY.read_text())
    full = json.loads(FULL.read_text())
    states, _, f, t2, q, h, source, _, weight, _, phase = decomposition()
    index = {state: j for j, state in enumerate(states)}
    complement = [index[kreweras(state)] for state in states]
    named = features(states, f, t2)
    upper = linalg.cholesky(weight, lower=False)
    inverse = linalg.inv(upper)
    mass = upper @ h @ inverse
    current = -(mass-mass.conj().T)/2
    k = upper @ (q.conj().T @ q[complement, :]).toarray() @ inverse
    involution = k/phase
    _, dual_basis = linalg.eigh((involution+involution.conj().T)/2)
    rows = []
    projected_functions = {"original": [], "reversible": []}
    for ray_index, sign in enumerate((-1, 1)):
        psi = upper @ source @ np.array([1, sign*phase])/np.sqrt(2)
        e = psi/linalg.norm(psi)
        je = current @ e
        source_overlap = np.vdot(e, je)
        orthogonal_je = je-e*source_overlap
        norm = float(linalg.norm(orthogonal_je))
        r = orthogonal_je/norm
        bases, captures = {}, []
        for dimension, stage in zip((7, 8), geometry["stages"]):
            names = stage["exact_named_columns"]["selected_columns_in_declaration_order"]
            columns = [named[name[2:]][complement] if name.startswith("K_") else named[name] for name in names]
            z = upper @ np.asarray(q.conj().T @ np.column_stack(columns))
            z = (z+sign*k @ z/phase)/2
            basis = linalg.svd(z, full_matrices=False)[0][:, :dimension]
            bases[dimension] = basis
            capture = float(linalg.norm(basis.conj().T @ r)**2)
            captures.append({"dimension": dimension, "fraction_of_current_direction_squared_norm": capture})
        t4 = upper @ np.asarray(q.conj().T @ named["T4"])
        t4 = (t4+sign*k @ t4/phase)/2
        residual = t4-bases[7] @ (bases[7].conj().T @ t4)
        t4_direction = residual/linalg.norm(residual)
        increment = float(abs(np.vdot(t4_direction, r))**2)

        # Same fixed two columns for both dynamics; r is not redefined under S.
        pair = np.column_stack((e, r))
        projected = pair.conj().T @ mass @ pair
        two_dimensional = {}
        for process, matrix in (("original", projected), ("reversible", (projected+projected.conj().T)/2)):
            score, values, residues = response(matrix, np.array([1., 0.]))
            projected_functions[process].append((values, residues))
            reference = full["ray_rows"][ray_index][process]
            score["mass_relative_error"] = score["mass"]/reference["lowest_source_visible_mass"]-1
            score["samples"] = [{"s": row["s"], "u": series(values, residues, row["s"]),
                                  "reference_u": row["u"], "relative_error": series(values, residues, row["s"])/row["u"]-1}
                                 for row in reference["samples"]]
            score["orthonormal_projected_mass_re_im"] = complex_display(matrix)
            score["initial_memory_k0"] = float((matrix[0,1]*matrix[1,0]).real)
            two_dimensional[process] = score
        v = dual_basis[:, :93] if sign < 0 else dual_basis[:, 93:]
        cross = pair_rows(v.conj().T @ mass @ v, v.conj().T @ e, v.conj().T @ r,
                          [0.]+protocol["distances"])
        rows.append({"ray": full["ray_rows"][ray_index]["ray"], "sign": sign,
                     "current_source_squared_norm": norm**2,
                     "source_current_overlap_re_im": complex_display(source_overlap),
                     "geometry_captures": captures,
                     "T4_residual_current_capture_increment": increment,
                     "captured_increment_from_projector_difference": captures[1]["fraction_of_current_direction_squared_norm"]-captures[0]["fraction_of_current_direction_squared_norm"],
                     "direct_current_pair": cross,
                     "maximal_initial_anti_cross_magnitude_by_definition": 2*norm,
                     "two_dimensional_projection": two_dimensional})
    comparisons = {}
    for process, curves in projected_functions.items():
        difference = lambda s: series(*curves[0],s)-series(*curves[1],s)
        lo, hi = protocol["crossing_bracket"]
        comparisons[process] = {"crossing_in_fixed_bracket": float(brentq(difference,lo,hi)) if difference(lo)*difference(hi)<0 else None,
                                "difference_at_fixed_bracket_ends": [difference(lo), difference(hi)]}
    return {"schema": protocol["schema"], "protocol": str(PROTOCOL.relative_to(ROOT)),
            "ray_rows": rows, "two_dimensional_ray_comparison": comparisons,
            "scope": "same finite process, unique generator-defined observer, existing geometry only; no new samples or state-count claim"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path, default=OUT)
    args = parser.parse_args()
    with threadpool_limits(limits=1):
        value = build_result()
    inputs = (PROTOCOL, Path(__file__), GEOMETRY, FULL,
              ROOT/"scripts/p398_width8_geometric_compression.py", ROOT/"scripts/p398_width8_pair_density_current.py")
    value["input_sha256"] = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(value, indent=2)+"\n")
    for row in value["ray_rows"]:
        print(row["ray"], row["geometry_captures"], "T4 increment", row["T4_residual_current_capture_increment"])
        print("derivative", row["direct_current_pair"]["anti_derivative_at_zero_re_im"])
        for process in ("original", "reversible"):
            score = row["two_dimensional_projection"][process]
            print(process, "mass", score["mass"], "k0", score["initial_memory_k0"], "t4 error", score["samples"][-1]["relative_error"])
    print(value["two_dimensional_ray_comparison"])


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Fixed named geometry x stationary-current deletion, with no basis search."""
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

from p398_width8_geometric_compression import features
from p398_width8_projected_memory import decomposition
from p398_width8_source_spectrum import complex_display, kreweras

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = ROOT / "analysis/p398_width8_geometry_current_interaction.json"
GEOMETRY = ROOT / "results/p398-width8-geometric-compression/latest.json"
T4 = ROOT / "results/p398-width8-t4-post-reveal/latest.json"
FULL = ROOT / "results/p398-width8-reversible-current-control/latest.json"
OUT = ROOT / "results/p398-width8-geometry-current-interaction/latest.json"


def build_result():
    protocol = json.loads(PROTOCOL.read_text())
    prior = json.loads(GEOMETRY.read_text())
    extension = json.loads(T4.read_text())
    full = json.loads(FULL.read_text())
    stages = {row["stage"]: row for row in prior["stages"]+extension["stages"]}
    states, _, f, t2, q, h, source, _, weight, _, phase = decomposition()
    index = {state: j for j, state in enumerate(states)}
    complement = [index[kreweras(state)] for state in states]
    named = features(states, f, t2)
    upper = linalg.cholesky(weight, lower=False)
    inverse = linalg.inv(upper)
    mass = upper @ h @ inverse
    k = upper @ (q.conj().T @ q[complement, :]).toarray() @ inverse
    rows = []
    full_crossing = {name: full["comparison"][name]["crossing_in_frozen_bracket"]
                     for name in ("original", "reversible")}
    for stage_name, dimension in zip(protocol["fixed_spans"], protocol["exact_dimensions_per_ray"]):
        archived = stages[stage_name]
        names = archived["exact_named_columns"]["selected_columns_in_declaration_order"]
        columns = [named[name[2:]][complement] if name.startswith("K_") else named[name] for name in names]
        z = upper @ np.asarray(q.conj().T @ np.column_stack(columns))
        stage = {"stage": stage_name, "dimension_per_ray": dimension,
                 "archived_exact_selected_columns": names, "rays": []}
        correlations = []
        for ray_index, sign in enumerate((-1, 1)):
            projected = (z + sign*k @ z/phase) / 2
            vectors, _, _ = linalg.svd(projected, full_matrices=False)
            basis = vectors[:, :dimension]
            b = basis.conj().T @ mass @ basis
            symmetric = (b + b.conj().T) / 2
            psi = upper @ source @ np.array([1, sign*phase]) / np.sqrt(2)
            a = basis.conj().T @ psi
            values, eigenvectors = linalg.eigh(symmetric)
            residues = np.abs(eigenvectors.conj().T @ a)**2 / np.vdot(a, a).real
            visible = np.flatnonzero(residues > 1e-10 * residues.max())
            lowest = int(visible[0])
            u = lambda s, ev=values, r=residues: float(np.sum(r*np.exp(-ev*s)))
            correlations.append(u)
            reference = full["ray_rows"][ray_index]
            old = archived["rays"][ray_index]
            original = {"mass": old["mass"],
                        "mass_relative_error": old["mass"]/reference["original"]["lowest_source_visible_mass"]-1,
                        "samples": [{"s": row["t"], "u": row["u_compressed"],
                                     "reference_u": row["u_reference"], "relative_error": row["relative_error"]}
                                    for row in old["samples"]]}
            reversible = {
                "mass": float(values[lowest]),
                "mass_relative_error": float(values[lowest]/reference["reversible"]["lowest_source_visible_mass"]-1),
                "samples": [{"s": row["s"], "u": u(row["s"]), "reference_u": row["u"],
                             "relative_error": u(row["s"])/row["u"]-1}
                            for row in reference["reversible"]["samples"]],
                "all_masses": values.tolist(), "normalized_residues": residues.tolist()}
            stage["rays"].append({"ray": reference["ray"], "sign": sign,
                "original": original, "reversible": reversible,
                "orthonormal_projected_original_mass_re_im": complex_display(b),
                "orthonormal_projected_source_re_im": complex_display(a)})
        lo, hi = protocol["crossing_bracket"]
        difference = lambda s: correlations[0](s)-correlations[1](s)
        crossing = float(brentq(difference, lo, hi)) if difference(lo)*difference(hi) < 0 else None
        stage["crossing"] = {
            "original": {"value": archived["crossing_in_frozen_bracket"],
                         "relative_error": archived["crossing_in_frozen_bracket"]/full_crossing["original"]-1},
            "reversible": {"value": crossing,
                           "relative_error": crossing/full_crossing["reversible"]-1 if crossing is not None else None}}
        rows.append(stage)

    def metric_map(stage, process):
        output = {"crossing": stage["crossing"][process]["relative_error"]}
        for ray in stage["rays"]:
            output[ray["ray"]+"_mass"] = ray[process]["mass_relative_error"]
            for sample in ray[process]["samples"]:
                if sample["s"] in protocol["primary_tail_distances"]:
                    output[ray["ray"]+"_u_t"+str(int(sample["s"]))] = sample["relative_error"]
        return output

    four_cells = {process: {str(row["dimension_per_ray"]): metric_map(row, process) for row in rows[-2:]}
                  for process in ("original", "reversible")}
    interaction = []
    for metric in four_cells["original"]["7"]:
        item = {"metric": metric}
        for process in ("original", "reversible"):
            e7 = four_cells[process]["7"][metric]
            e8 = four_cells[process]["8"][metric]
            item[process] = {"relative_error_7": e7, "relative_error_8": e8,
                             "signed_error_change_8_minus_7": e8-e7,
                             "absolute_error_reduction": abs(e7)-abs(e8)}
        item["signed_error_difference_in_differences_G_minus_S"] = item["original"]["signed_error_change_8_minus_7"] - item["reversible"]["signed_error_change_8_minus_7"]
        item["absolute_improvement_difference_G_minus_S"] = item["original"]["absolute_error_reduction"] - item["reversible"]["absolute_error_reduction"]
        interaction.append(item)
    return {"schema": protocol["schema"], "protocol": str(PROTOCOL.relative_to(ROOT)),
            "full_reference_crossings": full_crossing,
            "full_reference_masses": [{"ray": row["ray"], **{p: row[p]["lowest_source_visible_mass"] for p in ("original", "reversible")}} for row in full["ray_rows"]],
            "stages": rows, "geometry_T4_by_current_two_by_two": interaction,
            "original_process_outputs": "Reused from archived scores; not refitted/recomputed",
            "counterfactual": "Only take Hermitian part of each unchanged orthonormal projected mass",
            "scope": "fixed observable Galerkin spaces, not Markov state chains; descriptive finite-process interaction, no new samples"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path, default=OUT)
    args = parser.parse_args()
    with threadpool_limits(limits=1):
        value = build_result()
    inputs = (PROTOCOL, Path(__file__), GEOMETRY, T4, FULL,
              ROOT/"scripts/p398_width8_geometric_compression.py", ROOT/"scripts/p398_width8_projected_memory.py")
    value["input_sha256"] = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(value, indent=2)+"\n")
    for row in value["stages"]:
        print("dimension", row["dimension_per_ray"], "crossing", row["crossing"])
        for ray in row["rays"]:
            print(ray["ray"], {p: {"mass": ray[p]["mass"], "mass_error": ray[p]["mass_relative_error"],
                                   "tail_errors": [r["relative_error"] for r in ray[p]["samples"][-2:]]}
                               for p in ("original", "reversible")})
    print("interaction", json.dumps(value["geometry_T4_by_current_two_by_two"], indent=2))


if __name__ == "__main__":
    main()

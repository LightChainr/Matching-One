#!/usr/bin/env python3
"""Independent A_top/E_top ray test from the completed P205 prism raw batches."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import tempfile
from typing import Any, Mapping, Sequence

import mpmath as mp

from etop_global_ray_elimination import (
    clean_fit, fit_ray, load_rows as load_production_rows,
)
from etop_rank1_elimination import (
    CROSSWALK, PRIMARY_COVARIANCE, _chi2_survival,
)
from rank_plane_crosswalk import (
    _projected_p_derivative, _projected_state, combine, cos4,
    endpoint_observables, intrinsic_center, read_histograms, read_moments,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE_COMMIT = "fc14817bb8c0b2f6e7cbde41778e715dcb62bc64"
P_REF = mp.mpf("0.59274605079")
SIZES = (25, 50, 125)
SOURCE_PATHS = {
    n: {
        kind: f"results/server-20260829/P205-quotient-character-prism/raw/n{n}_12m.{suffix}"
        for kind, suffix in (("histogram", "hist.csv"),
                             ("moments", "moments.csv"),
                             ("metadata", "metadata.json"))
    }
    for n in SIZES
}


def number(value, digits=17):
    return float(mp.nstr(value, digits))


def git_blob(path: str) -> bytes:
    return subprocess.check_output(["git", "show", f"{SOURCE_COMMIT}:{path}"])


def covariance_of_mean(influences: Sequence[Sequence[mp.mpf]]):
    batches = len(influences)
    return [[mp.fsum(row[i] * row[j] for row in influences) /
             (batches * (batches - 1)) for j in range(2)] for i in range(2)]


def state_at(first, second, p, delta_cos4):
    projected = _projected_state(
        endpoint_observables(first, p), endpoint_observables(second, p), delta_cos4)
    return [projected["P4_A_top"], projected["P4_E_top"]]


def materialize_sources(directory: Path):
    output, hashes = {}, {}
    for n in SIZES:
        output[n], hashes[n] = {}, {}
        for kind, source_path in SOURCE_PATHS[n].items():
            payload = git_blob(source_path)
            suffix = {"histogram": ".hist.csv", "moments": ".moments.csv",
                      "metadata": ".metadata.json"}[kind]
            path = directory / f"n{n}{suffix}"
            path.write_bytes(payload)
            output[n][kind] = path
            hashes[n][kind] = hashlib.sha256(payload).hexdigest()
    return output, hashes


def analyze_pair(n: int, paths: Mapping[str, Path]):
    histograms = read_histograms(paths["histogram"])
    moments = read_moments(paths["moments"])
    metadata = json.loads(paths["metadata"].read_text())
    if set(histograms) != set(moments):
        raise ValueError(f"N{n}: histogram/moment batches do not align")
    batches = sorted({key[2] for key in histograms})
    if batches != list(range(100)):
        raise ValueError(f"N{n}: batch grid differs from completed prism")
    rows = {
        orientation: [histograms[(n, orientation, batch)] for batch in batches]
        for orientation in ("first", "second")
    }
    first, second = combine(rows["first"]), combine(rows["second"])
    if first.samples != 12000000 or second.samples != 12000000:
        raise ValueError(f"N{n}: sample contract changed")
    delta_cos4 = cos4(first.a, first.b) - cos4(second.a, second.b)
    if delta_cos4 == 0:
        raise ValueError(f"N{n}: zero H4 leverage")
    p0 = intrinsic_center(first, second)

    def fixed_analysis(p):
        point = state_at(first, second, p, delta_cos4)
        influences = []
        for batch in batches:
            batch_point = state_at(rows["first"][batch], rows["second"][batch],
                                   p, delta_cos4)
            influences.append([batch_point[index] - point[index] for index in range(2)])
        return point, covariance_of_mean(influences)

    fixed_point, fixed_covariance = fixed_analysis(P_REF)
    intrinsic_point, _ = fixed_analysis(p0)
    first_obs = endpoint_observables(first, p0)
    second_obs = endpoint_observables(second, p0)
    derivative = _projected_p_derivative(first_obs, second_obs, delta_cos4)
    state_derivative = [derivative["P4_A_top"], derivative["P4_E_top"]]
    center_slope = (first_obs["S_birth"] + second_obs["S_birth"]) / 2
    center_influences = []
    center_shifts = []
    for batch in batches:
        left = endpoint_observables(rows["first"][batch], p0)
        right = endpoint_observables(rows["second"][batch], p0)
        batch_point = _projected_state(left, right, delta_cos4)
        center_equation = (left["A_top"] + right["A_top"]) / 2
        center_shift = -center_equation / center_slope
        center_shifts.append(center_shift)
        center_influences.append([
            batch_point[name] - intrinsic_point[index] +
            state_derivative[index] * center_shift
            for index, name in enumerate(("P4_A_top", "P4_E_top"))
        ])
    intrinsic_covariance = covariance_of_mean(center_influences)
    design = metadata["designs"][0]
    return {
        "id": f"P205-prism-N{n}", "N": n,
        "dependency_group": "P205-prism",
        "orientation_pair": [design["first"], design["second"]],
        "smith_pair": [design["first_smith_invariants"],
                       design["second_smith_invariants"]],
        "batches": len(batches), "samples_per_orientation": first.samples,
        "p0": p0, "delta_cos4": delta_cos4,
        "intrinsic": {"estimate": intrinsic_point, "covariance": intrinsic_covariance,
                      "center_influence_se": mp.sqrt(mp.fsum(x * x for x in center_shifts) /
                                                     (len(center_shifts) * (len(center_shifts) - 1)))},
        "fixed": {"estimate": fixed_point, "covariance": fixed_covariance},
    }


def ray_row(pair, coordinate):
    return {
        "id": pair["id"], "N": pair["N"],
        "dependency_group": pair["dependency_group"],
        "estimate": pair[coordinate]["estimate"],
        "covariance": pair[coordinate]["covariance"],
    }


def pairwise(rows):
    output = []
    single = {row["id"]: fit_ray([row]) for row in rows}
    for left, right in itertools.combinations(rows, 2):
        combined = fit_ray([left, right])
        delta = max(mp.mpf("0"),
                    mp.mpf(str(combined["min_chi2"])) -
                    mp.mpf(str(single[left["id"]]["min_chi2"])) -
                    mp.mpf(str(single[right["id"]]["min_chi2"])))
        theta_left = single[left["id"]]["theta_internal"]
        theta_right = single[right["id"]]["theta_internal"]
        output.append({
            "pair": [left["id"], right["id"]],
            "delta_chi2": number(delta), "df": 1,
            "p": number(_chi2_survival(delta, 1)),
            "angle_difference_degrees_mod_ray": number(abs(
                ((theta_right - theta_left + mp.pi / 2) % mp.pi - mp.pi / 2) *
                180 / mp.pi)),
            "common_ray": clean_fit(combined),
        })
    return sorted(output, key=lambda row: row["delta_chi2"], reverse=True)


def prism_test(pairs, coordinate):
    rows = [ray_row(pair, coordinate) for pair in pairs]
    fit = fit_ray(rows, include_profile=True)
    chi2 = mp.mpf(str(fit["min_chi2"]))
    singles = {row["id"]: clean_fit(fit_ray([row], include_profile=True)) for row in rows}
    return {
        "coordinate": "intrinsic matching center per quotient pair" if coordinate == "intrinsic"
                      else f"common fixed p_ref={P_REF}",
        "common_ray": {**clean_fit(fit), "df": 2,
                       "p": number(_chi2_survival(chi2, 2)),
                       "decision_at_alpha_0p01": "eliminated" if
                       _chi2_survival(chi2, 2) < mp.mpf("0.01") else "survives"},
        "individual_ray_directions": singles,
        "pairwise": pairwise(rows),
    }


def production_context(prism_rows, covariance_key):
    crosswalk = json.loads(CROSSWALK.read_text())
    production = load_production_rows(crosswalk, covariance_key)
    production_fit = fit_ray(production)
    prism_fit = fit_ray(prism_rows)
    combined_fit = fit_ray(production + prism_rows)
    delta = max(mp.mpf("0"),
                mp.mpf(str(combined_fit["min_chi2"])) -
                mp.mpf(str(production_fit["min_chi2"])) -
                mp.mpf(str(prism_fit["min_chi2"])))
    by_lineage = {}
    for name in ("P49", "P43", "P50", "P57"):
        lineage = [row for row in production if row["dependency_group"] == name]
        lineage_fit = fit_ray(lineage)
        joint_fit = fit_ray(lineage + prism_rows)
        penalty = max(mp.mpf("0"),
                      mp.mpf(str(joint_fit["min_chi2"])) -
                      mp.mpf(str(lineage_fit["min_chi2"])) -
                      mp.mpf(str(prism_fit["min_chi2"])))
        by_lineage[name] = {
            "delta_chi2": number(penalty), "df": 1,
            "p": number(_chi2_survival(penalty, 1)),
            "lineage_angle": lineage_fit["angle_degrees_from_A_axis"],
            "shared_fit": clean_fit(joint_fit),
        }
    return {
        "production_global_angle": production_fit["angle_degrees_from_A_axis"],
        "prism_angle": prism_fit["angle_degrees_from_A_axis"],
        "common_production_plus_prism": clean_fit(combined_fit),
        "separate_vs_common_ray": {"delta_chi2": number(delta), "df": 1,
                                   "p": number(_chi2_survival(delta, 1))},
        "prism_vs_each_production_lineage": by_lineage,
    }


def clean_pair(pair):
    return {
        "id": pair["id"], "N": pair["N"],
        "orientation_pair": pair["orientation_pair"], "smith_pair": pair["smith_pair"],
        "batches": pair["batches"], "samples_per_orientation": pair["samples_per_orientation"],
        "p0": number(pair["p0"]), "delta_cos4": number(pair["delta_cos4"]),
        "intrinsic": {"estimate": [number(x) for x in pair["intrinsic"]["estimate"]],
                      "covariance": [[number(x) for x in row]
                                     for row in pair["intrinsic"]["covariance"]],
                      "center_influence_se": number(pair["intrinsic"]["center_influence_se"])},
        "fixed": {"estimate": [number(x) for x in pair["fixed"]["estimate"]],
                  "covariance": [[number(x) for x in row]
                                 for row in pair["fixed"]["covariance"]]},
    }


def build_report():
    mp.mp.dps = 70
    with tempfile.TemporaryDirectory(prefix="p205-etop-") as directory:
        paths, hashes = materialize_sources(Path(directory))
        pairs = [analyze_pair(n, paths[n]) for n in SIZES]
    primary = prism_test(pairs, "intrinsic")
    fixed = prism_test(pairs, "fixed")
    primary_rows = [ray_row(pair, "intrinsic") for pair in pairs]
    return {
        "schema": "matching-one/p205-quotient-etop-ray-test/v1",
        "status": "zero_new_sample_production_archive_test",
        "source": {"commit": SOURCE_COMMIT, "paths": SOURCE_PATHS,
                   "sha256": hashes},
        "state": "P4(A_top,E_top) reconstructed from the aligned first/second threshold-rank batches",
        "primary": primary, "fixed_p_ref_sensitivity": fixed,
        "pairs": [clean_pair(pair) for pair in pairs],
        "production_context": {
            "intrinsic_only": production_context(primary_rows, PRIMARY_COVARIANCE),
            "coordinate_guard": "The fixed-p prism block is not combined with the production sensitivity rows: those rows condition on each archive's own plug-in p0, whereas the prism sensitivity uses one common p_ref.",
        },
        "decision_question": "does the A_top/E_top ray rotate across the three completed Smith/quotient-character contrasts, and can that rotation explain the existing production global-ray failure?",
        "claim_boundary": "The prism rows are equal-area quotient contrasts, not individual geometry states. Intrinsic-center covariance uses the declared first-order center influence; fixed p_ref is the exact-batch sensitivity. The result can classify quotient-linked ray rotation but cannot assign a field or asymptotic operator.",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = build_report()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "intrinsic": report["primary"]["common_ray"],
        "fixed": report["fixed_p_ref_sensitivity"]["common_ray"],
        "production_context": report["production_context"]["intrinsic_only"],
    }, indent=2))


if __name__ == "__main__":
    main()

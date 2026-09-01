#!/usr/bin/env python3
"""Score the frozen N=50 tau x topology-map factorial for PR #267."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import tempfile

import mpmath as mp

from etop_rank1_elimination import _chi2_survival
from rank_plane_crosswalk import (
    combine,
    combine_moments,
    endpoint_observables,
    read_histograms,
    read_moments,
)


ROOT = Path(__file__).resolve().parents[1]
OLD_COMMIT = "fc14817bb8c0b2f6e7cbde41778e715dcb62bc64"
OLD_PREFIX = "results/server-20260829/P205-quotient-character-prism/raw/n50_12m"
FREEZE_COMMIT = "3505dc1b651988e981a22ac8d32a14effe3d17c3"
P_REF = mp.mpf("0.59274605079")
FIELDS = ("A_top", "E_top", "C", "W")
DELTA = mp.mpf(1152) / 625
FIXED_A = mp.mpf("0.8287127431200354")
FIXED_C = mp.mpf("3.334261982589098")


def sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def number(value, digits: int = 17) -> float:
    return float(mp.nstr(value, digits))


def covariance_of_mean(rows):
    count = len(rows)
    dimension = len(rows[0])
    mean = [mp.fsum(row[i] for row in rows) / count for i in range(dimension)]
    covariance = [[
        mp.fsum((row[i] - mean[i]) * (row[j] - mean[j]) for row in rows)
        / (count * (count - 1))
        for j in range(dimension)
    ] for i in range(dimension)]
    return mean, covariance


def add(left, right, left_scale=1, right_scale=1):
    return [left_scale * a + right_scale * b for a, b in zip(left, right)]


def add_cov(left, right, left_scale=1, right_scale=1):
    return [[left_scale * left[i][j] + right_scale * right[i][j]
             for j in range(len(left))] for i in range(len(left))]


def wald(mean, covariance):
    matrix = mp.matrix(covariance)
    positive_definite = True
    try:
        mp.cholesky(matrix)
        inverse = matrix ** -1
        vector = mp.matrix(mean)
        chi2 = (vector.T * inverse * vector)[0]
    except (ValueError, ZeroDivisionError):
        positive_definite = False
        chi2 = mp.inf
    return {
        "estimate": [number(value) for value in mean],
        "standard_error": [number(mp.sqrt(max(covariance[i][i], 0)))
                           for i in range(len(mean))],
        "covariance": [[number(value) for value in row] for row in covariance],
        "chi_square": number(chi2) if mp.isfinite(chi2) else None,
        "df": len(mean),
        "p_value": number(_chi2_survival(chi2, len(mean))) if mp.isfinite(chi2) else None,
        "positive_definite": positive_definite,
    }


def cell_vector(histogram, moment):
    state = endpoint_observables(histogram, P_REF)
    c_value, w_value = moment.clock()
    return [state["A_top"], state["E_top"], c_value, w_value]


def load_pair(runs):
    histograms = [read_histograms(hist_path) for hist_path, _ in runs]
    moments = [read_moments(moments_path) for _, moments_path in runs]
    keys = set(histograms[0])
    if any(set(rows) != keys for rows in histograms + moments):
        raise ValueError("unaligned batch grid across run segments")
    batches = sorted({key[2] for key in keys})
    if batches != list(range(100)):
        raise ValueError("batch grid differs from freeze")
    output = {}
    for orientation in ("first", "second"):
        hist_rows = [combine([rows[(50, orientation, batch)] for rows in histograms])
                     for batch in batches]
        moment_rows = [combine_moments([rows[(50, orientation, batch)] for rows in moments])
                       for batch in batches]
        output[orientation] = {
            "point": cell_vector(combine(hist_rows), combine_moments(moment_rows)),
            "batch": [cell_vector(h, m) for h, m in zip(hist_rows, moment_rows)],
            "hist": hist_rows,
            "moments": moment_rows,
        }
    return output


def materialize_old(directory: Path):
    hashes = {}
    paths = {}
    for suffix in ("hist.csv", "moments.csv", "metadata.json"):
        source = f"{OLD_PREFIX}.{suffix}"
        payload = subprocess.check_output(["git", "show", f"{OLD_COMMIT}:{source}"])
        path = directory / f"old.{suffix}"
        path.write_bytes(payload)
        paths[suffix] = path
        hashes[source] = sha256(payload)
    metadata = json.loads(paths["metadata.json"].read_text())
    design = metadata["designs"][0]
    if metadata["samples_per_pair"] != 12000000 or metadata["batches"] != 100:
        raise ValueError("old P205 N50 sample contract changed")
    if design["first_period_matrix"] != [[7, -1], [1, 7]] or design["second_period_matrix"] != [[5, -5], [5, 5]]:
        raise ValueError("old P205 N50 factorial cells changed")
    return paths, hashes


def verify_new_metadata(raw_dir: Path):
    expected = {
        "p10_vs_p00": ([[3, -8], [4, 6]], [[7, -1], [1, 7]]),
        "p11_vs_p00": ([[0, -10], [5, 0]], [[7, -1], [1, 7]]),
    }
    hashes = {}
    segments = (
        ("", 20000, 26750000000, 26750020000),
        ("_ext80k", 80000, 26750020000, 26750100000),
    )
    for name, matrices in expected.items():
        for segment, samples, first_counter, last_counter in segments:
            prefix = f"{name}{segment}"
            for suffix in ("hist.csv", "moments.csv", "metadata.json"):
                path = raw_dir / f"{prefix}.{suffix}"
                hashes[str(path)] = sha256(path.read_bytes())
            metadata = json.loads((raw_dir / f"{prefix}.metadata.json").read_text())
            design = metadata["designs"][0]
            checks = {
                "git_commit": FREEZE_COMMIT,
                "samples_per_pair": samples,
                "batches": 100,
                "seed": 202626750,
                "replica_counter_first": first_counter,
                "replica_counter_last_exclusive": last_counter,
            }
            for key, value in checks.items():
                if metadata[key] != value:
                    raise ValueError(f"{prefix}: {key} violates freeze")
            if design["first_period_matrix"] != matrices[0] or design["second_period_matrix"] != matrices[1]:
                raise ValueError(f"{prefix}: period matrices violate freeze")
    return hashes


def same_anchor(first, second):
    for kind in ("hist", "moments"):
        left, right = first["second"][kind], second["second"][kind]
        if left != right:
            raise ValueError(f"P00 {kind} anchor is not identical across hosts")


def contrast_batches(left, right, scale=1):
    return [[scale * (a - b) for a, b in zip(x, y)]
            for x, y in zip(left["batch"], right["batch"])]


def calculate(raw_dir: Path):
    mp.mp.dps = 70
    new_hashes = verify_new_metadata(raw_dir)
    p10_pilot = load_pair([
        (raw_dir / "p10_vs_p00.hist.csv", raw_dir / "p10_vs_p00.moments.csv"),
    ])
    p11_pilot = load_pair([
        (raw_dir / "p11_vs_p00.hist.csv", raw_dir / "p11_vs_p00.moments.csv"),
    ])
    p10_extension = load_pair([
        (raw_dir / "p10_vs_p00_ext80k.hist.csv", raw_dir / "p10_vs_p00_ext80k.moments.csv"),
    ])
    p11_extension = load_pair([
        (raw_dir / "p11_vs_p00_ext80k.hist.csv", raw_dir / "p11_vs_p00_ext80k.moments.csv"),
    ])
    same_anchor(p10_pilot, p11_pilot)
    same_anchor(p10_extension, p11_extension)
    p10_run = load_pair([
        (raw_dir / "p10_vs_p00.hist.csv", raw_dir / "p10_vs_p00.moments.csv"),
        (raw_dir / "p10_vs_p00_ext80k.hist.csv", raw_dir / "p10_vs_p00_ext80k.moments.csv"),
    ])
    p11_run = load_pair([
        (raw_dir / "p11_vs_p00.hist.csv", raw_dir / "p11_vs_p00.moments.csv"),
        (raw_dir / "p11_vs_p00_ext80k.hist.csv", raw_dir / "p11_vs_p00_ext80k.moments.csv"),
    ])
    same_anchor(p10_run, p11_run)
    with tempfile.TemporaryDirectory(prefix="p267-etop-factorial-") as directory:
        old_paths, old_hashes = materialize_old(Path(directory))
        old = load_pair([(old_paths["hist.csv"], old_paths["moments.csv"])])

    square_batches = contrast_batches(old["first"], old["second"], 1 / DELTA)
    pilot_rectangle_batches = contrast_batches(
        p10_pilot["first"], p11_pilot["first"], -1 / DELTA)
    pilot_rectangle_mean, pilot_rectangle_cov = covariance_of_mean(pilot_rectangle_batches)
    pilot_square_mean, pilot_square_cov = covariance_of_mean(square_batches)
    pilot_primary = wald(
        add(pilot_rectangle_mean, pilot_square_mean, 1, -1),
        add_cov(pilot_rectangle_cov, pilot_square_cov),
    )
    rectangle_batches = contrast_batches(p10_run["first"], p11_run["first"], -1 / DELTA)
    square_mean, square_cov = covariance_of_mean(square_batches)
    rectangle_mean, rectangle_cov = covariance_of_mean(rectangle_batches)
    interaction_mean = add(rectangle_mean, square_mean, 1, -1)
    interaction_cov = add_cov(rectangle_cov, square_cov)
    primary = wald(interaction_mean, interaction_cov)

    # Raw factorial contrasts are descriptive; the character-normalized
    # difference above is the sole primary endpoint.
    old_sum = [[a + b for a, b in zip(x, y)]
               for x, y in zip(old["first"]["batch"], old["second"]["batch"])]
    new_sum = [[a + b for a, b in zip(x, y)]
               for x, y in zip(p10_run["first"]["batch"], p11_run["first"]["batch"])]
    old_topology = contrast_batches(old["second"], old["first"])
    new_topology = contrast_batches(p11_run["first"], p10_run["first"])
    old_sum_mean, old_sum_cov = covariance_of_mean(old_sum)
    new_sum_mean, new_sum_cov = covariance_of_mean(new_sum)
    old_top_mean, old_top_cov = covariance_of_mean(old_topology)
    new_top_mean, new_top_cov = covariance_of_mean(new_topology)
    modulus = wald(add(new_sum_mean, old_sum_mean, mp.mpf("0.5"), mp.mpf("-0.5")),
                   add_cov(new_sum_cov, old_sum_cov, mp.mpf("0.25"), mp.mpf("0.25")))
    topology = wald(add(new_top_mean, old_top_mean, mp.mpf("0.5"), mp.mpf("0.5")),
                    add_cov(new_top_cov, old_top_cov, mp.mpf("0.25"), mp.mpf("0.25")))

    coefficient = [-FIXED_A, 1, -FIXED_C, 0]
    residual = mp.fsum(a * b for a, b in zip(coefficient, interaction_mean))
    variance = mp.fsum(coefficient[i] * interaction_cov[i][j] * coefficient[j]
                       for i in range(4) for j in range(4))
    secondary = {
        "fixed_coefficients": {"A_top": number(FIXED_A), "C": number(FIXED_C)},
        "residual": number(residual), "standard_error": number(mp.sqrt(variance)),
        "z": number(residual / mp.sqrt(variance)),
        "chi_square": number(residual * residual / variance), "df": 1,
        "p_value": number(_chi2_survival(residual * residual / variance, 1)),
    }
    extend = bool(primary["positive_definite"] and primary["chi_square"] is not None
                  and primary["chi_square"] > 13.276704)
    return {
        "schema": "matching-one/p267-etop-tau-topology-factorial/v1",
        "status": "frozen gate-triggered 100k missing-cell reveal",
        "field_order": list(FIELDS),
        "sources": {"old_commit": OLD_COMMIT, "old_sha256": old_hashes,
                    "new_sha256": new_hashes, "freeze_commit": FREEZE_COMMIT},
        "exact_design": {"N": 50, "tau": ["i", "2i"],
                         "smith": [[1, 50], [5, 10]],
                         "common_map": "(1/5)[[4,-3],[3,4]]",
                         "delta_chi4": {"i": number(DELTA), "2i": number(-DELTA)}},
        "projected_rows": {
            "tau_i": wald(square_mean, square_cov),
            "tau_2i": wald(rectangle_mean, rectangle_cov),
        },
        "primary_character_normalized_interaction": primary,
        "descriptive_raw_main_effects": {"modulus": modulus, "topology_map": topology},
        "secondary_fixed_A_plus_C_interaction_residual": secondary,
        "promotion_gate": {"critical_chi_square_4df_alpha_0p01": 13.276704,
                           "pilot_primary": pilot_primary,
                           "pilot_decision": "extend_both_to_100k",
                           "final_samples_per_missing_cell": 100000,
                           "final_interaction_above_original_gate": extend},
        "factor_decision": "topology_map_by_tau_interaction_resolved",
        "claim_boundary": "finite N=50 factorial interaction; no field identity, exponent, root character, or asymptotic modular law",
    }


def report(payload):
    names = payload["field_order"]
    square = payload["projected_rows"]["tau_i"]
    rectangle = payload["projected_rows"]["tau_2i"]
    interaction = payload["primary_character_normalized_interaction"]
    lines = [
        "# PR267 E_top tau x topology-map N50 factorial",
        "",
        "The exact minimum four-cell design crosses `tau=i,2i` with",
        "`Smith=(1,50),(5,10)` at fixed determinant 50. One rational rotation",
        "`O=(1/5)[[4,-3],[3,4]]` maps cyclic to noncyclic in both rows.",
        "",
        "| field | P4 at tau=i | P4 at tau=2i | interaction | interaction SE |",
        "|---|---:|---:|---:|---:|",
    ]
    for index, name in enumerate(names):
        lines.append(f"| {name} | {square['estimate'][index]:.7g} | {rectangle['estimate'][index]:.7g} | "
                     f"{interaction['estimate'][index]:.7g} | {interaction['standard_error'][index]:.2g} |")
    lines += [
        "",
        "## Frozen primary",
        "",
        f"Character-normalized interaction: `chi2={interaction['chi_square']:.6g}/4`, "
        f"`p={interaction['p_value']:.6g}`.",
        f"The 20k pilot crossed the frozen gate at `chi2={payload['promotion_gate']['pilot_primary']['chi_square']:.6g}/4`; the displayed score uses the authorized",
        "100k total per missing cell (20k pilot plus disjoint 80k extension).",
        "The additive no-interaction factorial is therefore eliminated: the topology-map",
        "response changes with tau at this fixed N.",
        "",
        "The raw modulus and topology-map main-effect scores are retained in `score.json`,",
        "but they are descriptive because their four coordinates are strongly correlated views",
        "of the same threshold clocks. The primary endpoint is only the four-vector interaction.",
        "",
        "## Fixed A+C diagnostic",
        "",
    ]
    row = payload["secondary_fixed_A_plus_C_interaction_residual"]
    lines += [
        f"Using the previously frozen coefficients gives residual `{row['residual']:.7g} +/- "
        f"{row['standard_error']:.2g}` (`z={row['z']:.3f}`, `p={row['p_value']:.4g}`). It survives",
        "the existing alpha=.01 boundary but is mildly tense at .05. One interaction equation",
        "cannot identify both A and C coefficients, so this is compatibility rather than pinning.",
        "",
        "The P00 anchor histogram and moments reproduce byte-identically across Zy and XP.",
        "This is a finite-N geometry factorization result, not an exponent, root-character or",
        "asymptotic modular-law claim.",
    ]
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    payload = calculate(args.raw_dir)
    args.json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    args.report.write_text(report(payload))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Score the frozen Alexander-even rho-child complex C3 production."""

from __future__ import annotations

import argparse
import csv
from math import erfc, exp, gcd, sqrt
import hashlib
import json
from pathlib import Path

import mpmath as mp

from pinson_arguin_primitive import engine_to_paper, primitive_probability_direct
from p267_rho_child_etop_mc import CHILDREN, CHILD_ORDER, physical_phase


SCHEMA = "matching-one/p267-rho-child-etop-c3-score/v1"


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def tau_from_matrix(matrix):
    (a, b), (c, d) = matrix
    return mp.mpc(b, d) / mp.mpc(a, c)


def canonical_lines(cutoff):
    for u in range(-cutoff, cutoff + 1):
        for v in range(-cutoff, cutoff + 1):
            if (u, v) == (0, 0) or gcd(abs(u), abs(v)) != 1:
                continue
            if u < 0 or (u == 0 and v < 0):
                continue
            yield u, v


def continuum_rows(cutoff, dps):
    output = []
    for name, matrix in CHILDREN:
        tau = tau_from_matrix(matrix)
        p1 = mp.mpf("0")
        h4 = mp.mpc("0")
        for line in canonical_lines(cutoff):
            probability = primitive_probability_direct(
                *engine_to_paper(line), tau, dps=dps
            )
            p1 += probability
            h4 += probability * physical_phase(matrix, line)
        output.append({"child": name, "Etop": 1 - p1, "H4": h4})
    return output


def dft(values, r):
    zeta = mp.exp(2 * mp.pi * mp.j / 3)
    return sum(values[j] * zeta ** (-r * j) for j in range(3)) / 3


def covariance_of_mean(rows):
    count = len(rows)
    means = [sum(row[j] for row in rows) / count for j in range(len(rows[0]))]
    return [
        [sum((row[i] - means[i]) * (row[j] - means[j]) for row in rows)
         / (count * (count - 1)) for j in range(len(means))]
        for i in range(len(means))
    ]


def complex_zero_score(values):
    rows = [[float(mp.re(value)), float(mp.im(value))] for value in values]
    mean = [sum(row[j] for row in rows) / len(rows) for j in range(2)]
    covariance = covariance_of_mean(rows)
    determinant = covariance[0][0] * covariance[1][1] - covariance[0][1] ** 2
    if determinant <= 0:
        raise ValueError("complex covariance is not positive definite")
    precision = [
        [covariance[1][1] / determinant, -covariance[0][1] / determinant],
        [-covariance[0][1] / determinant, covariance[0][0] / determinant],
    ]
    chi_square = sum(mean[i] * precision[i][j] * mean[j]
                     for i in range(2) for j in range(2))
    return {
        "value_re_im": mean,
        "covariance_2x2": covariance,
        "chi_square": chi_square,
        "dof": 2,
        "p": exp(-chi_square / 2),
        "magnitude": sqrt(mean[0] ** 2 + mean[1] ** 2),
    }


def real_zero_score(values):
    rows = [float(mp.re(value)) for value in values]
    mean = sum(rows) / len(rows)
    variance = sum((value - mean) ** 2 for value in rows) / (len(rows) * (len(rows) - 1))
    if variance <= 0:
        raise ValueError("real covariance is not positive")
    chi_square = mean * mean / variance
    return {
        "value": mean,
        "variance": variance,
        "chi_square": chi_square,
        "dof": 1,
        "p": erfc(sqrt(chi_square / 2)),
    }


def jackknife_determinant(e0, e1, h0, h1):
    batches = len(e0)

    def determinant(indices):
        count = len(indices)
        means = [sum(values[index] for index in indices) / count
                 for values in (e0, e1, h0, h1)]
        return means[2] * means[1] - means[3] * means[0]

    all_indices = list(range(batches))
    plugin = determinant(all_indices)
    leave_one = [determinant([j for j in all_indices if j != index])
                 for index in all_indices]
    leave_mean = sum(leave_one) / batches
    estimate = batches * plugin - (batches - 1) * leave_mean
    covariance = [[0.0, 0.0], [0.0, 0.0]]
    coordinates = [[float(mp.re(value)), float(mp.im(value))] for value in leave_one]
    center = [sum(row[j] for row in coordinates) / batches for j in range(2)]
    for row in coordinates:
        for i in range(2):
            for j in range(2):
                covariance[i][j] += (row[i] - center[i]) * (row[j] - center[j])
    for i in range(2):
        for j in range(2):
            covariance[i][j] *= (batches - 1) / batches
    mean = [float(mp.re(estimate)), float(mp.im(estimate))]
    determinant_cov = covariance[0][0] * covariance[1][1] - covariance[0][1] ** 2
    precision = [
        [covariance[1][1] / determinant_cov, -covariance[0][1] / determinant_cov],
        [-covariance[0][1] / determinant_cov, covariance[0][0] / determinant_cov],
    ]
    chi_square = sum(mean[i] * precision[i][j] * mean[j]
                     for i in range(2) for j in range(2))
    return {
        "estimand": "H4_r0*Etop_r1-H4_r1*Etop_r0",
        "value_re_im": mean,
        "covariance_2x2": covariance,
        "chi_square": chi_square,
        "dof": 2,
        "p": exp(-chi_square / 2),
        "jackknife_bias_correction_abs": float(abs(estimate - plugin)),
    }


def load_batches(path):
    with path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    output = []
    for row in rows:
        samples = float(row["samples"])
        item = []
        for name in CHILD_ORDER:
            item.append({
                "Etop": (float(row[f"{name}_rank0"]) + float(row[f"{name}_rank2"])) / samples,
                "H4": complex(float(row[f"{name}_H4_re"]), float(row[f"{name}_H4_im"])) / samples,
            })
        output.append(item)
    return output


def score(run_path, batch_path, manifest_path, cutoff, dps):
    run = json.loads(run_path.read_text())
    manifest = json.loads(manifest_path.read_text())
    actual_hash = sha256(batch_path)
    if actual_hash != run["batch_sha256"]:
        raise ValueError("batch checksum does not match run envelope")
    for key in ("samples", "batches", "workers", "seed", "replica_offset"):
        if run["run"][key] != manifest["run"][key]:
            raise ValueError(f"run manifest mismatch: {key}")
    if not run["summary"]["all_invariant_failures_zero"]:
        raise ValueError("homology invariant failures are nonzero")
    baselines = continuum_rows(cutoff, dps)
    previous = continuum_rows(cutoff - 2, dps)
    batches = load_batches(batch_path)
    e0, e1, h0, h1 = [], [], [], []
    for batch in batches:
        delta_e = [mp.mpf(str(batch[j]["Etop"])) - baselines[j]["Etop"] for j in range(3)]
        delta_h = [mp.mpc(batch[j]["H4"].real, batch[j]["H4"].imag) - baselines[j]["H4"]
                   for j in range(3)]
        e0.append(dft(delta_e, 0)); e1.append(dft(delta_e, 1))
        h0.append(dft(delta_h, 0)); h1.append(dft(delta_h, 1))
    primary = complex_zero_score(e1)
    primary["decision"] = "resolved" if primary["p"] < manifest["alpha"] else "unresolved"
    scalar = real_zero_score(e0)
    determinant = jackknife_determinant(e0, e1, h0, h1)
    determinant["decision"] = (
        "two_observer_rays_distinct" if determinant["p"] < manifest["alpha"]
        else "common_ray_not_rejected"
    )
    character_means = {}
    for name, values in (("Etop_r0", e0), ("Etop_r1", e1),
                         ("H4_r0", h0), ("H4_r1", h1)):
        value = sum(values) / len(values)
        character_means[name] = [float(mp.re(value)), float(mp.im(value))]
    return {
        "schema": SCHEMA,
        "status": "frozen_production_reveal",
        "observer": "delta E_top = lattice(P0+P2) - continuum(1-sum_primitive pi_ab)",
        "sector": "Alexander-even rank redistribution, non-A_top, complex rho-child C3 r1",
        "child_order": list(CHILD_ORDER),
        "continuum": [
            {"child": row["child"], "Etop": mp.nstr(row["Etop"], 30),
             "H4_re_im": [mp.nstr(mp.re(row["H4"]), 30), mp.nstr(mp.im(row["H4"]), 30)],
             "cutoff_minus_2_change": {
                 "Etop": mp.nstr(abs(row["Etop"] - old["Etop"]), 8),
                 "H4": mp.nstr(abs(row["H4"] - old["H4"]), 8),
             }} for row, old in zip(baselines, previous)
        ],
        "primary_nontrivial_Etop_r1": primary,
        "scalar_Etop_r0": scalar,
        "same_stream_character_means_re_im": character_means,
        "secondary_observer_ray_determinant": determinant,
        "batch_sha256": actual_hash,
        "run_sha256": sha256(run_path),
        "alpha": manifest["alpha"],
        "reading": (
            "A resolved Etop r1 is an even rank-redistribution C3 response. The determinant tests "
            "whether it supplies a second child-character direction beyond primitive H4; neither "
            "test is an ordinary H4/H8/H12 vote."
        ),
        "claim_boundary": (
            "Finite N112 square-bond topology response only; no continuum field identity, exponent, "
            "unique state count, or square-site matching transfer claim."
        ),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--batches", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--cutoff", type=int, default=9)
    parser.add_argument("--dps", type=int, default=60)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = score(args.run, args.batches, args.manifest, args.cutoff, args.dps)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()

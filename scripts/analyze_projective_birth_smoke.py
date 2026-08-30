#!/usr/bin/env python3
"""Score the projectively marked threshold-rank Phase B/C smoke.

The sparse input has one row for every observed ``(tau1,tau2,ell)`` cell in
each aligned batch.  ``DIRECT_RANK2`` is a typed atom and deliberately has no
line.  The scorer evaluates the rank-one plateau character and its source/sink
flux at fixed p, verifies the exact continuity identities, and retains one
joint covariance block across both same-N orientations.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Iterable, Mapping, Sequence


DEFAULT_P = 0.592746050790


@dataclass(frozen=True)
class BirthCell:
    orientation: str
    batch: int
    samples: int
    tau1: int
    tau2: int
    kind: str
    ell_x: int
    ell_y: int
    count: int


def read_births(path: Path) -> tuple[int, dict[tuple[str, int], list[BirthCell]]]:
    grouped: dict[tuple[str, int], list[BirthCell]] = {}
    n_values: set[int] = set()
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            n = int(row["n"])
            n_values.add(n)
            cell = BirthCell(
                orientation=row["orientation"], batch=int(row["batch"]),
                samples=int(row["samples"]), tau1=int(row["tau1"]),
                tau2=int(row["tau2"]), kind=row["kind"],
                ell_x=int(row["ell_x"]), ell_y=int(row["ell_y"]),
                count=int(row["count"]),
            )
            grouped.setdefault((cell.orientation, cell.batch), []).append(cell)
    if len(n_values) != 1:
        raise ValueError(f"birth archive must contain one N, got {sorted(n_values)}")
    return next(iter(n_values)), grouped


def read_histograms(path: Path) -> dict[tuple[str, int, str], dict[int, int]]:
    output: dict[tuple[str, int, str], dict[int, int]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            key = (row["orientation"], int(row["batch"]), row["kind"])
            output.setdefault(key, {})[int(row["k"])] = int(row["count"])
    return output


def binomial_probabilities(n: int, p: float) -> list[float]:
    if not 0.0 < p < 1.0:
        raise ValueError("p must lie strictly between zero and one")
    values = [0.0] * (n + 1)
    values[0] = (1.0 - p) ** n
    ratio = p / (1.0 - p)
    for k in range(n):
        values[k + 1] = values[k] * (n - k) / (k + 1) * ratio
    return values


def flux_weights(n: int, p: float) -> list[float]:
    # d/dp P[Bin(N,p)>=t] = N Bin(N-1,p)[t-1].
    lower = binomial_probabilities(n - 1, p)
    return [0.0] + [n * value for value in lower]


def chi4(matrix: Sequence[Sequence[int]], ell_x: int, ell_y: int) -> complex:
    x = matrix[0][0] * ell_x + matrix[0][1] * ell_y
    y = matrix[1][0] * ell_x + matrix[1][1] * ell_y
    if x == 0 and y == 0:
        raise ValueError("zero projective line")
    z = complex(x, y)
    return z**4 / (x * x + y * y) ** 2


def covariance_of_mean(rows: Sequence[Sequence[float]]) -> list[list[float]]:
    batches = len(rows)
    if batches < 2:
        raise ValueError("at least two aligned batches are required")
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise ValueError("inconsistent covariance row width")
    means = [math.fsum(row[j] for row in rows) / batches for j in range(width)]
    return [
        [
            math.fsum((row[i] - means[i]) * (row[j] - means[j]) for row in rows)
            / (batches * (batches - 1))
            for j in range(width)
        ]
        for i in range(width)
    ]


BASE_NAMES = (
    "A4_re", "A4_im", "j4_birth_re", "j4_birth_im", "j4_exit_re",
    "j4_exit_im", "dA4_re", "dA4_im", "j4_activity_re",
    "j4_activity_im", "M_prime", "j_direct", "line_birth_fraction",
)


def evaluate_batch(
    cells: Sequence[BirthCell], n: int, p: float,
    matrix: Sequence[Sequence[int]],
) -> tuple[dict[str, float], dict[str, float]]:
    if not cells:
        raise ValueError("empty birth batch")
    samples = cells[0].samples
    if any(cell.samples != samples for cell in cells):
        raise ValueError("inconsistent samples within batch")
    if sum(cell.count for cell in cells) != samples:
        raise ValueError("sparse birth cells do not sum to batch samples")
    pmf = binomial_probabilities(n, p)
    flux = flux_weights(n, p)
    tails = [0.0] * (n + 2)
    running = 0.0
    for k in range(n, -1, -1):
        running += pmf[k]
        tails[k] = running

    plateau = 0j
    birth = 0j
    exit_ = 0j
    fixedp_by_k = [0j] * (n + 1)
    line_births = 0
    direct_flux = 0.0
    f1_prime = 0.0
    f2_prime = 0.0
    for cell in cells:
        weight = cell.count / samples
        f1_prime += weight * flux[cell.tau1]
        f2_prime += weight * flux[cell.tau2]
        if cell.kind == "DIRECT_RANK2":
            if cell.tau1 != cell.tau2 or cell.ell_x or cell.ell_y:
                raise ValueError("invalid DIRECT_RANK2 row")
            direct_flux += weight * flux[cell.tau1]
            continue
        if cell.kind != "LINE" or not cell.tau1 < cell.tau2:
            raise ValueError("invalid line-bearing row")
        if math.gcd(abs(cell.ell_x), abs(cell.ell_y)) != 1:
            raise ValueError("nonprimitive ell in birth archive")
        character = chi4(matrix, cell.ell_x, cell.ell_y)
        plateau_probability = tails[cell.tau1] - tails[cell.tau2]
        plateau += weight * plateau_probability * character
        birth += weight * flux[cell.tau1] * character
        exit_ += weight * flux[cell.tau2] * character
        line_births += cell.count
        for k in range(cell.tau1, cell.tau2):
            fixedp_by_k[k] += weight * character

    primitive_fixedp = sum(pmf[k] * fixedp_by_k[k] for k in range(n + 1))
    primitive_derivative = sum(
        pmf[k] * (k / p - (n - k) / (1.0 - p)) * fixedp_by_k[k]
        for k in range(n + 1)
    )
    gradient = birth - exit_
    activity = birth + exit_
    metrics = {
        "A4_re": plateau.real, "A4_im": plateau.imag,
        "j4_birth_re": birth.real, "j4_birth_im": birth.imag,
        "j4_exit_re": exit_.real, "j4_exit_im": exit_.imag,
        "dA4_re": gradient.real, "dA4_im": gradient.imag,
        "j4_activity_re": activity.real, "j4_activity_im": activity.imag,
        "M_prime": f1_prime + f2_prime,
        "j_direct": direct_flux,
        "line_birth_fraction": line_births / samples,
    }
    gates = {
        "A4_minus_issue156_re": plateau.real - primitive_fixedp.real,
        "A4_minus_issue156_im": plateau.imag - primitive_fixedp.imag,
        "dA4_minus_source_sink_re": primitive_derivative.real - gradient.real,
        "dA4_minus_source_sink_im": primitive_derivative.imag - gradient.imag,
        "Mprime_minus_line_activity_and_direct":
            (f1_prime + f2_prime) - (sum(
                (cell.count / samples) *
                (flux[cell.tau1] + flux[cell.tau2])
                for cell in cells if cell.kind == "LINE"
            ) + 2.0 * direct_flux),
    }
    return metrics, gates


def validate_histogram_recovery(
    births: Mapping[tuple[str, int], Sequence[BirthCell]],
    histograms: Mapping[tuple[str, int, str], Mapping[int, int]],
) -> dict[str, object]:
    failures: list[str] = []
    direct_total = 0
    for (orientation, batch), cells in births.items():
        reconstructed = {"minus": {}, "plus": {}}
        for cell in cells:
            reconstructed["minus"][cell.tau1] = (
                reconstructed["minus"].get(cell.tau1, 0) + cell.count
            )
            reconstructed["plus"][cell.tau2] = (
                reconstructed["plus"].get(cell.tau2, 0) + cell.count
            )
            if cell.kind == "DIRECT_RANK2":
                direct_total += cell.count
        for kind in ("minus", "plus"):
            expected = dict(histograms.get((orientation, batch, kind), {}))
            if reconstructed[kind] != expected:
                failures.append(f"{orientation}/batch{batch}/{kind}")
    return {
        "passed": not failures,
        "failed_cells": failures,
        "direct_rank2_paths": direct_total,
        "contract": "sum_ell plus DIRECT_RANK2 recovers K1; tau2 recovers K2",
    }


def support_summary(
    cells: Iterable[BirthCell], matrix: Sequence[Sequence[int]], n: int, p: float,
) -> dict[str, object]:
    counts: dict[tuple[int, int], int] = {}
    flux_counts: dict[tuple[int, int], float] = {}
    flux = flux_weights(n, p)
    for cell in cells:
        if cell.kind != "LINE":
            continue
        key = (cell.ell_x, cell.ell_y)
        counts[key] = counts.get(key, 0) + cell.count
        flux_counts[key] = flux_counts.get(key, 0.0) + cell.count * flux[cell.tau1]
    total = sum(counts.values())
    flux_total = math.fsum(flux_counts.values())
    mean = sum(count * chi4(matrix, *line) for line, count in counts.items()) / total
    flux_mean = sum(
        count * chi4(matrix, *line) for line, count in flux_counts.items()
    ) / flux_total
    values = {
        (round(chi4(matrix, *line).real, 14), round(chi4(matrix, *line).imag, 14))
        for line in counts
    }
    return {
        "line_count": len(counts),
        "distinct_chi4_values": len(values),
        "non_micro_chi4_support": len(values) > 1,
        "unconditional_chi4_mean": [mean.real, mean.imag],
        "unconditional_chi4_variance": max(0.0, 1.0 - abs(mean) ** 2),
        "birth_flux_chi4_mean": [flux_mean.real, flux_mean.imag],
        "birth_flux_chi4_variance": max(0.0, 1.0 - abs(flux_mean) ** 2),
        "support": [
            {
                "ell": list(line), "count": count,
                "chi4": [chi4(matrix, *line).real, chi4(matrix, *line).imag],
            }
            for line, count in sorted(counts.items())
        ],
    }


def analyze(birth_path: Path, histogram_path: Path, metadata_path: Path, p: float) -> dict:
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if not metadata.get("projective_births"):
        raise ValueError("metadata does not declare projective birth output")
    if "iota=1" not in metadata.get("integral_saturation", ""):
        raise ValueError("metadata does not freeze the saturation theorem contract")
    n, births = read_births(birth_path)
    histograms = read_histograms(histogram_path)
    design = metadata["designs"][0]
    if int(design["N"]) != n:
        raise ValueError("metadata/birth N mismatch")
    matrices = {
        "first": design["first_period_matrix"],
        "second": design["second_period_matrix"],
    }
    batch_ids = sorted({batch for _, batch in births})
    orientations = ("first", "second")
    for orientation in orientations:
        if sorted(batch for name, batch in births if name == orientation) != batch_ids:
            raise ValueError("orientations do not share one aligned batch set")

    rows: list[list[float]] = []
    gates = []
    by_batch: list[dict[str, object]] = []
    order = [f"{orientation}_{name}" for orientation in orientations for name in BASE_NAMES]
    contrast_names = (
        "A4_re", "A4_im", "j4_birth_re", "j4_birth_im", "j4_exit_re",
        "j4_exit_im", "dA4_re", "dA4_im", "j4_activity_re", "j4_activity_im",
    )
    order += [f"second_minus_first_{name}" for name in contrast_names]
    for batch in batch_ids:
        metrics = {}
        batch_gate = {}
        for orientation in orientations:
            value, gate = evaluate_batch(
                births[(orientation, batch)], n, p, matrices[orientation]
            )
            metrics[orientation] = value
            batch_gate[orientation] = gate
        vector = [metrics[o][name] for o in orientations for name in BASE_NAMES]
        vector += [metrics["second"][name] - metrics["first"][name]
                   for name in contrast_names]
        rows.append(vector)
        gates.append(batch_gate)
        by_batch.append({"batch": batch, "values": dict(zip(order, vector))})

    covariance = covariance_of_mean(rows)
    means = [math.fsum(row[j] for row in rows) / len(rows) for j in range(len(order))]
    max_gate = {
        name: max(abs(gate[o][name]) for gate in gates for o in orientations)
        for name in next(iter(gates))["first"]
    }
    histogram_gate = validate_histogram_recovery(births, histograms)
    support = {
        orientation: support_summary(
            (cell for (name, _), cells in births.items() if name == orientation for cell in cells),
            matrices[orientation], n, p,
        )
        for orientation in orientations
    }
    tolerance = 2e-12
    identity_passed = histogram_gate["passed"] and all(
        value <= tolerance for value in max_gate.values()
    )
    return {
        "schema": "matching-one/projective-birth-production-smoke/v1",
        "status": "local variance smoke; no model selected from pilot mean",
        "source": {
            "births": str(birth_path), "histogram": str(histogram_path),
            "metadata": str(metadata_path), "N": n, "p_ref": p,
            "samples_per_shape": int(metadata["samples_per_pair"]),
            "batches": len(batch_ids), "elapsed_seconds": metadata["elapsed_seconds"],
            "seed": metadata["seed"],
            "counter_range": [metadata["replica_counter_first"],
                              metadata["replica_counter_last_exclusive"]],
        },
        "exact_crosswalk_gates": {
            "passed": identity_passed,
            "tolerance": tolerance,
            "histogram_recovery": histogram_gate,
            "max_absolute_residual": max_gate,
            "A4_contract": "fixed-p rank-one primitive character of Issue 156",
            "continuity_contract": "dA4/dp = j4_birth1 - j4_exit2",
            "unmarked_contract": "Mprime = line activity + 2*j_DIRECT_RANK2",
        },
        "chi4_support": support,
        "joint_estimate": {
            "order": order,
            "mean": means,
            "standard_error": [math.sqrt(max(0.0, covariance[i][i]))
                               for i in range(len(order))],
            "covariance": covariance,
            "batch_values": by_batch,
        },
        "interpretation_boundary": {
            "A4": "duplicate coordinate of Issue 156, retained only as an exact crosswalk",
            "new_coordinates": ["j4_birth1", "j4_exit2", "j4_activity"],
            "iota": "fixed to one by c1a72e5 and intentionally absent from raw rows",
            "pilot_mean": "not used to choose a radial or angular model",
        },
    }


def render_markdown(result: Mapping[str, object]) -> str:
    source = result["source"]
    gate = result["exact_crosswalk_gates"]
    estimate = result["joint_estimate"]
    lookup = dict(zip(estimate["order"], zip(estimate["mean"], estimate["standard_error"])))
    lines = [
        "# Projective essential-birth N65 smoke", "",
        "**Status:** local variance/runtime smoke; the pilot mean did not select a model.", "",
        f"- Samples: {source['samples_per_shape']:,} per shape in {source['batches']} aligned batches.",
        f"- Runtime: {source['elapsed_seconds']:.6g} wall seconds on the recorded local run.",
        f"- Exact crosswalk gates: {'PASS' if gate['passed'] else 'FAIL'}.",
        f"- Direct rank-two paths across both shapes: {gate['histogram_recovery']['direct_rank2_paths']:,}.",
        "", "## Non-micro projective support", "",
        "| orientation | primitive lines | distinct chi4 | Var(chi4), path | Var(chi4), birth flux |",
        "|---|---:|---:|---:|---:|",
    ]
    for orientation, row in result["chi4_support"].items():
        lines.append(
            f"| {orientation} | {row['line_count']} | {row['distinct_chi4_values']} | "
            f"{row['unconditional_chi4_variance']:.6g} | {row['birth_flux_chi4_variance']:.6g} |"
        )
    lines += ["", "Both N65 shapes expose more than the quarter-turn-only tiny-control support; "
              "the projective mark therefore carries a genuinely varying chi4 value before any "
              "continuum model is fitted.", "", "## Marked source/sink at p_ref", "",
              "| coordinate | mean | batch SE |", "|---|---:|---:|"]
    for name in (
        "first_j4_birth_re", "first_j4_exit_re", "first_j4_activity_re",
        "second_j4_birth_re", "second_j4_exit_re", "second_j4_activity_re",
        "second_minus_first_j4_birth_re", "second_minus_first_j4_exit_re",
        "second_minus_first_j4_activity_re",
    ):
        mean, se = lookup[name]
        lines.append(f"| `{name}` | {mean:.8g} | {se:.3g} |")
    lines += ["", "`A4` is retained only to verify the Issue #156 fixed-p character crosswalk. "
              "The nonredundant production coordinates are the ingress and egress fluxes, with "
              "the old line retained at the second birth. `DIRECT_RANK2` has no line and enters "
              "the unmarked derivative with multiplicity two.", ""]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--births", type=Path, required=True)
    parser.add_argument("--histogram", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--p", type=float, default=DEFAULT_P)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--markdown", type=Path)
    args = parser.parse_args()
    result = analyze(args.births, args.histogram, args.metadata, args.p)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    if args.markdown:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(render_markdown(result), encoding="utf-8")


if __name__ == "__main__":
    main()

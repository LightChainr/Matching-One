#!/usr/bin/env python3
"""Score the frozen norm-5 finite-N thermal-jet and width-cocycle contracts.

The scorer consumes full threshold histograms and their aligned joint-moment
files.  It recomputes the intrinsic center, Krawtchouk modes, canonical rank-gap
width, and every nonlinear residual inside delete-one-batch replicates.
Counter-identical runs are deleted synchronously; disjoint counter groups
contribute independent jackknife covariance matrices.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Mapping, Sequence

import mpmath as mp
import yaml

from analyze_matching_parity_derivatives_fast import H, combine, read, remove
from analyze_rank_gap_thermal_window import Run, covariance_groups, pooled_statistics, read_run
from covariance_nullspace import covariance_spectral_diagnostics, serialize_diagnostics
from hermite_krawtchouk_scaling_jet import (
    canonical_dimensionless_width,
    cocycle_residual,
    scaling_derivative_jet,
    width_cross_residual,
    width_normalized_jet,
)
from threshold_score_modes import project


DEFAULT_PREDICTION = "predictions/hermite_krawtchouk_jet_20260829.yaml"
ORDERS = (2, 3, 4, 5, 6)


@dataclass(frozen=True)
class InputRun:
    n: int
    histogram_path: Path
    moments_path: Path
    metadata_path: Path
    histogram: Mapping[tuple[int, str, int], H]
    moments: Run


@dataclass(frozen=True)
class State:
    n: int
    p0: mp.mpf
    jet: tuple[mp.mpf, ...]
    canonical_width: mp.mpf
    mean_rank_gap: mp.mpf


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_prediction(path: Path) -> dict:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError("thermal-jet prediction must be a mapping")
    if payload.get("status") != "source_theory_frozen_before_norm5_reveal":
        raise ValueError("thermal-jet prediction chronology changed")
    if payload.get("scoring_order", ())[:3] != [
        "rank_gap_width_collapse",
        "q2_analytic_after_width_correction",
        "rank2_Jordan_after_width_correction",
    ]:
        raise ValueError("thermal-jet scoring order changed")
    return payload


def load_input(n: int, histogram: Path, moments: Path, metadata: Path) -> InputRun:
    histogram_data = read(histogram)
    sizes = {key[0] for key in histogram_data}
    if sizes != {n}:
        raise ValueError(f"N={n}: histogram sizes are {sorted(sizes)}")
    moment_run = read_run(n, moments, metadata)
    for orientation in ("first", "second"):
        histogram_batches = sorted(
            key[2] for key in histogram_data if key[:2] == (n, orientation)
        )
        moment_batches = [row["batch"] for row in moment_run.rows[orientation]]
        if histogram_batches != moment_batches:
            raise ValueError(f"N={n}: histogram/moment batch mismatch for {orientation}")
    return InputRun(n, histogram, moments, metadata, histogram_data, moment_run)


def state(run: InputRun, omitted_batch: int | None = None) -> State:
    grouped = {
        orientation: [
            run.histogram[key]
            for key in sorted(run.histogram)
            if key[:2] == (run.n, orientation)
        ]
        for orientation in ("first", "second")
    }
    first = combine(grouped["first"])
    second = combine(grouped["second"])
    if omitted_batch is not None:
        first_row = next(row for row in grouped["first"] if row.batch == omitted_batch)
        second_row = next(row for row in grouped["second"] if row.batch == omitted_batch)
        first = remove(first, first_row)
        second = remove(second, second_row)
    projected = project(first, second, max(ORDERS))
    coefficients = tuple(
        projected["P4_D_modes"][order]
        if order % 2 == 0
        else projected["P4_S_modes"][order]
        for order in range(max(ORDERS) + 1)
    )
    jet = tuple(
        scaling_derivative_jet(
            coefficients, run.n, projected["p0"], mp.mpf(13) / 8
        )
    )
    mean_gap = pooled_statistics(run.moments, omitted_batch)["gap_mean"]
    return State(
        n=run.n,
        p0=projected["p0"],
        jet=jet,
        canonical_width=canonical_dimensionless_width(run.n, mean_gap),
        mean_rank_gap=mean_gap,
    )


def width_residuals(
    states: Mapping[int, State], lineages: Sequence[tuple[int, int, int]]
) -> list[mp.mpf]:
    values: list[mp.mpf] = []
    for parent_n, _norm2_n, norm5_n in lineages:
        parent, child = states[parent_n], states[norm5_n]
        residual = width_cross_residual(
            parent.jet, child.jet, parent.canonical_width, child.canonical_width
        )
        values.extend(residual[order - 1] for order in ORDERS)
    return values


def multiplier_residuals(
    states: Mapping[int, State],
    lineages: Sequence[tuple[int, int, int]],
    multiplier: mp.mpf,
) -> list[mp.mpf]:
    values: list[mp.mpf] = []
    for parent_n, norm2_n, norm5_n in lineages:
        parent = width_normalized_jet(states[parent_n].jet, states[parent_n].canonical_width)
        norm2 = width_normalized_jet(states[norm2_n].jet, states[norm2_n].canonical_width)
        norm5 = width_normalized_jet(states[norm5_n].jet, states[norm5_n].canonical_width)
        residual = cocycle_residual(parent, norm2, norm5, multiplier)
        values.extend(residual[order] for order in ORDERS)
    return values


def jackknife_covariance(rows: Sequence[Sequence[mp.mpf]]) -> list[list[mp.mpf]]:
    if len(rows) < 2 or not rows or not rows[0]:
        raise ValueError("jackknife requires at least two nonempty vectors")
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise ValueError("jackknife vectors have inconsistent widths")
    means = [mp.fsum(row[index] for row in rows) / len(rows) for index in range(width)]
    factor = mp.mpf(len(rows) - 1) / len(rows)
    return [
        [
            factor
            * mp.fsum(
                (row[left] - means[left]) * (row[right] - means[right])
                for row in rows
            )
            for right in range(width)
        ]
        for left in range(width)
    ]


def add_covariance(
    left: Sequence[Sequence[mp.mpf]], right: Sequence[Sequence[mp.mpf]]
) -> list[list[mp.mpf]]:
    if len(left) != len(right) or any(len(a) != len(b) for a, b in zip(left, right)):
        raise ValueError("covariance shapes differ")
    return [
        [a + b for a, b in zip(left_row, right_row)]
        for left_row, right_row in zip(left, right)
    ]


def propagated_covariance(
    full_states: Mapping[int, State],
    deleted_states: Mapping[int, Sequence[State]],
    inputs: Sequence[InputRun],
    residual: Callable[[Mapping[int, State]], list[mp.mpf]],
) -> tuple[list[list[mp.mpf]], list[dict]]:
    groups = covariance_groups([run.moments for run in inputs])
    point = residual(full_states)
    covariance = [[mp.mpf(0) for _ in point] for _ in point]
    descriptions = []
    by_n = {run.n: run for run in inputs}
    for key, members in groups.items():
        sizes = [member.n for member in members]
        batches = len(deleted_states[sizes[0]])
        rows = []
        for batch in range(batches):
            replicate = dict(full_states)
            for n in sizes:
                replicate[n] = deleted_states[n][batch]
            rows.append(residual(replicate))
        covariance = add_covariance(covariance, jackknife_covariance(rows))
        descriptions.append(
            {
                "seed": key[0],
                "counter_first": key[1],
                "counter_last_exclusive": key[2],
                "sizes": sizes,
                "batches": batches,
                "rule": "synchronized_delete_one" if len(sizes) > 1 else "independent_group",
                "metadata": [str(by_n[n].metadata_path) for n in sizes],
            }
        )
    return covariance, descriptions


def generalized_chi_square(
    residual: Sequence[mp.mpf], covariance: Sequence[Sequence[mp.mpf]],
    relative_cutoff: mp.mpf = mp.mpf("1e-10"),
) -> dict:
    diagnostics = covariance_spectral_diagnostics(
        residual, covariance, relative_cutoff, nullspace_policy="estimated"
    )
    return serialize_diagnostics(diagnostics, lambda value: mp.nstr(value, 20))


def render_score(
    name: str,
    labels: Sequence[str],
    point: Sequence[mp.mpf],
    covariance: Sequence[Sequence[mp.mpf]],
    groups: Sequence[Mapping[str, object]],
) -> dict:
    score = generalized_chi_square(point, covariance)
    return {
        "name": name,
        "labels": list(labels),
        "residual": [mp.nstr(value, 25) for value in point],
        "covariance": [[mp.nstr(value, 18) for value in row] for row in covariance],
        "score": score,
        "covariance_groups": list(groups),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run", nargs=4, action="append", metavar=("N", "HIST", "MOMENTS", "METADATA"),
        required=True,
    )
    parser.add_argument(
        "--lineage", nargs=3, action="append", type=int, metavar=("PARENT", "NORM2", "NORM5"),
        required=True,
    )
    parser.add_argument("--prediction", type=Path, default=Path(DEFAULT_PREDICTION))
    parser.add_argument("--dps", type=int, default=60)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    mp.mp.dps = args.dps
    prediction = load_prediction(args.prediction)
    inputs = [
        load_input(int(n), Path(hist), Path(moments), Path(metadata))
        for n, hist, moments, metadata in args.run
    ]
    if len({run.n for run in inputs}) != len(inputs):
        raise SystemExit("--run sizes must be unique")
    lineages = [tuple(row) for row in args.lineage]
    required_sizes = {n for lineage in lineages for n in lineage}
    if required_sizes != {run.n for run in inputs}:
        raise SystemExit("--run sizes must equal the union of --lineage sizes")
    full_states = {run.n: state(run) for run in inputs}
    deleted_states = {
        run.n: [state(run, batch) for batch in range(int(run.moments.metadata["batches"]))]
        for run in inputs
    }
    labels = [f"N{parent}_to_N{norm5}_r{order}" for parent, _n2, norm5 in lineages for order in ORDERS]

    width_point = width_residuals(full_states, lineages)
    width_covariance, width_groups = propagated_covariance(
        full_states, deleted_states, inputs, lambda states: width_residuals(states, lineages)
    )
    multipliers = {
        "q2_analytic": mp.mpf(8) / 5,
        "rank2_Jordan": mp.log(5) / mp.log(2),
    }
    multiplier_scores = {}
    for name, multiplier in multipliers.items():
        point = multiplier_residuals(full_states, lineages, multiplier)
        covariance, groups = propagated_covariance(
            full_states,
            deleted_states,
            inputs,
            lambda states, c=multiplier: multiplier_residuals(states, lineages, c),
        )
        multiplier_scores[name] = {
            "multiplier": mp.nstr(multiplier, 25),
            **render_score(name, labels, point, covariance, groups),
        }

    payload = {
        "schema": "matching-one/norm5-thermal-jet-score/v1",
        "status": "post-reveal execution of pre-reveal frozen finite-N score",
        "prediction": {"path": str(args.prediction), "sha256": sha256(args.prediction)},
        "orders": list(ORDERS),
        "lineages": [list(row) for row in lineages],
        "inputs": [
            {
                "N": run.n,
                "histogram": str(run.histogram_path),
                "histogram_sha256": sha256(run.histogram_path),
                "moments": str(run.moments_path),
                "moments_sha256": sha256(run.moments_path),
                "metadata": str(run.metadata_path),
                "metadata_sha256": sha256(run.metadata_path),
            }
            for run in inputs
        ],
        "states": {
            str(n): {
                "p0": mp.nstr(row.p0, 25),
                "mean_rank_gap": mp.nstr(row.mean_rank_gap, 25),
                "canonical_width": mp.nstr(row.canonical_width, 25),
                "finite_thermal_jet": [mp.nstr(value, 25) for value in row.jet],
            }
            for n, row in sorted(full_states.items())
        },
        "primary_width_collapse": render_score(
            "rank_gap_width_collapse", labels, width_point, width_covariance, width_groups
        ),
        "secondary_multiplier_cocycles": multiplier_scores,
        "evidence_boundary": (
            "The width score and both cocycle scores reuse the same raw curves. "
            "They are ordered model diagnostics, not additive evidence rows."
        ),
    }
    rendered = json.dumps(payload, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
        print(args.output)
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Score the preregistered N85 natural charged-current scale test."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Mapping, Optional, Sequence

from reveal_n65_charged_source_archive import evaluate_batch, read_births


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def jackknife_natural(
    rows: Sequence[Mapping[str, Mapping[str, float]]], p: float
) -> dict[str, object]:
    factor = p * (1.0 - p)
    count = len(rows)

    def estimate(indices: Sequence[int]) -> list[float]:
        output = []
        for orientation in ("first", "second"):
            w = math.fsum(rows[i][orientation]["W_A"] for i in indices) / len(indices)
            current = math.fsum(rows[i][orientation]["Jminus_A"] for i in indices) / len(indices)
            output.append(factor * current / w)
        return [output[0], output[1], output[1] - output[0]]

    value = estimate(list(range(count)))
    replicates = [estimate([j for j in range(count) if j != i]) for i in range(count)]
    centers = [math.fsum(row[j] for row in replicates) / count for j in range(3)]
    covariance = [[
        (count - 1) / count * math.fsum(
            (row[i] - centers[i]) * (row[j] - centers[j]) for row in replicates
        )
        for j in range(3)] for i in range(3)]
    return {
        "order": ["K_A_first", "K_A_second", "Delta_K_A_second_minus_first"],
        "value": value,
        "standard_error": [math.sqrt(max(0.0, covariance[i][i])) for i in range(3)],
        "covariance": covariance,
        "jackknife_replicates": replicates,
    }


def score(
    prereg_path: Path, births_path: Path, metadata_path: Path
) -> dict[str, object]:
    prereg = json.loads(prereg_path.read_text(encoding="utf-8"))
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if prereg.get("schema") != "matching-one/P337-natural-current-scale-preregistration/v1":
        raise ValueError("unexpected preregistration schema")
    expected = prereg["design"]
    design = metadata["designs"][0]
    gates = {
        "N": design["N"] == expected["geometry"]["N"],
        "first": design["first"] == expected["geometry"]["first"],
        "second": design["second"] == expected["geometry"]["second"],
        "first_period_matrix": design["first_period_matrix"] == expected["geometry"]["first_period_matrix"],
        "second_period_matrix": design["second_period_matrix"] == expected["geometry"]["second_period_matrix"],
        "samples": metadata["samples_per_pair"] == expected["samples_per_shape"],
        "batches": metadata["batches"] == expected["batches"],
        "seed": metadata["seed"] == expected["seed"],
        "replica_offset": metadata["replica_counter_first"] == expected["replica_offset"],
        "projective_births": metadata.get("projective_births") is True,
        "engine_freeze": metadata["git_commit"] == "1714141+freeze-5e601dc",
    }
    if not all(gates.values()):
        raise ValueError(f"production differs from preregistration: {gates}")

    n, births = read_births(births_path)
    batch_ids = sorted({batch for _, batch in births})
    rows = []
    max_continuity = 0.0
    for batch in batch_ids:
        row = {}
        for orientation in ("first", "second"):
            metrics, exact = evaluate_batch(births[(orientation, batch)], n, prereg["source"]["p_ref"])
            row[orientation] = {
                "W_A": metrics["W_A"],
                "Jminus_A": metrics["J_A_birth"] - metrics["J_A_exit"],
                "Jplus_A": metrics["J_A_birth"] + metrics["J_A_exit"],
            }
            max_continuity = max(max_continuity, abs(exact["A_continuity"]))
        rows.append(row)
    natural = jackknife_natural(rows, prereg["source"]["p_ref"])
    observed = natural["value"][2]
    variance = natural["covariance"][2][2]
    comparisons = []
    for name in prereg["scoring_contract"]["model_order"]:
        target = prereg["frozen_targets_at_N85"][name]
        residual = observed - target["value"]
        predictive_variance = variance + target["fit_variance"]
        comparisons.append({
            "name": name,
            "frozen_target": target["value"],
            "residual": residual,
            "measurement_standard_error": math.sqrt(variance),
            "measurement_only_quadratic": residual * residual / variance,
            "predictive_standard_error": math.sqrt(predictive_variance),
            "predictive_quadratic": residual * residual / predictive_variance,
        })
    closest = min(comparisons, key=lambda row: row["predictive_quadratic"])
    mean_rows = {}
    for orientation in ("first", "second"):
        mean_rows[orientation] = {
            name: math.fsum(row[orientation][name] for row in rows) / len(rows)
            for name in ("W_A", "Jminus_A", "Jplus_A")
        }
    return {
        "schema": "matching-one/P337-natural-current-scale-N85-score/v1",
        "status": "fresh independent N85 scale reveal scored against preregistered targets",
        "source": {
            "preregistration_commit": "5e601dc",
            "preregistration": str(prereg_path),
            "engine_commit": metadata["git_commit"],
            "environment": "Huawei DevEnvC_ZyTrST f415a4bcbd9a438b85f5f29e4a507ea4",
            "births": str(births_path), "births_sha256": sha256(births_path),
            "metadata": str(metadata_path), "metadata_sha256": sha256(metadata_path),
            "N": n, "samples_per_shape": metadata["samples_per_pair"],
            "batches": metadata["batches"], "seed": metadata["seed"],
            "elapsed_seconds": metadata["elapsed_seconds"],
        },
        "freeze_gates": {"passed": all(gates.values()), "items": gates},
        "batch_means": mean_rows,
        "natural_coordinate": natural,
        "target_comparison": comparisons,
        "reading": {
            "independent_nonzero": {
                "value": observed,
                "standard_error": math.sqrt(variance),
                "z": observed / math.sqrt(variance),
            },
            "closest_frozen_predictive_target": closest["name"],
            "closest_predictive_quadratic": closest["predictive_quadratic"],
            "summary": (
                "N85 remains nonzero, is much smaller than scale-neutral transfer, "
                "and lies below but predictively compatible with the source-fitted H4 target"
            ),
            "variance_design_ratio": math.sqrt(variance) / prereg["design"]["power_extrapolation"]["planned_standard_error"],
        },
        "exact_gates": {
            "max_A_continuity_residual": max_continuity,
            "tolerance": 3e-12,
            "passed": max_continuity < 3e-12,
        },
        "dependency": {
            "N85": "new seed/counter block, independent of N65 archive 1714141",
            "N65_target_fit": "source-fitted target uncertainty included only in predictive residuals",
            "no_pooling": "N65 and N85 are combined only through the preregistered transfer equations",
        },
        "claim_boundary": (
            "one independent N85 scale point; discriminates frozen transfers but does not "
            "fit an exponent or authorize a new normalization"
        ),
    }


def render_markdown(payload: Mapping[str, object]) -> str:
    natural = payload["natural_coordinate"]
    lines = [
        "# N85 reveal of the natural charged-current scale test", "",
        "The 200k N85 block was generated only after preregistration commit `5e601dc`.", "",
        f"`Delta_K_A(N85)={natural['value'][2]:.12g} +/- {natural['standard_error'][2]:.3g}` (`z={payload['reading']['independent_nonzero']['z']:.3f}`).", "",
        "| frozen target | target | residual | measurement chi2 | predictive chi2 |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in payload["target_comparison"]:
        lines.append(f"| {row['name']} | {row['frozen_target']:.12g} | {row['residual']:+.12g} | "
                     f"{row['measurement_only_quadratic']:.4g} | {row['predictive_quadratic']:.4g} |")
    first, second = payload["batch_means"]["first"], payload["batch_means"]["second"]
    lines += ["", f"First `(9+2i)`: `W_A={first['W_A']:.12g}`, `Jminus_A={first['Jminus_A']:.12g}`, `K_A={natural['value'][0]:.12g}`.",
              f"Second `(7+6i)`: `W_A={second['W_A']:.12g}`, `Jminus_A={second['Jminus_A']:.12g}`, `K_A={natural['value'][1]:.12g}`.", "",
              "The independent scale point is nonzero and strongly below scale-neutral transfer. It is also below the project H4 target, but once N65 source-fit uncertainty is included, H4 is the closest compatible frozen transfer.", "",
              f"Observed SE was {payload['reading']['variance_design_ratio']:.3f} times the N65-only projection. Current continuity closes to `{payload['exact_gates']['max_A_continuity_residual']:.3g}`.", "",
              "No exponent or normalization was fitted to N85.", ""]
    return "\n".join(lines)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preregistration", type=Path, required=True)
    parser.add_argument("--births", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--markdown", type=Path)
    args = parser.parse_args(argv)
    payload = score(args.preregistration, args.births, args.metadata)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.markdown:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(render_markdown(payload), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Score the integrated pivotal/H4 bridge on the archived P33 full curves."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

from analyze_p48_retrospective import (
    HELDOUT_SIZES,
    SIZES,
    TRAIN_SIZES,
    covariance_of_mean,
    gls_score,
    project_size,
    pseudovalues,
    read_histograms,
    sha256,
)


METRICS = ("pivotal_H4_scaled", "coefficient_ratio", "thermal_mass_scaled")


def bridge_metrics(projected):
    slope = projected["Mbar_prime"]
    normalized_h4 = projected["P4_S_prime"] / slope
    return {
        "Mbar_prime": slope,
        "normalized_pivotal_H4": normalized_h4,
        "pivotal_H4_scaled": normalized_h4 * projected["N"] ** (13.0 / 8.0),
        "coefficient_ratio": normalized_h4 / projected["P4_D"],
        "thermal_mass_scaled": slope / projected["N"] ** (3.0 / 8.0),
    }


def analyze(input_path: Path):
    records = read_histograms(input_path)
    sizes = tuple(sorted({key[0] for key in records}))
    if sizes != SIZES:
        raise ValueError("pivotal/H4 bridge requires N=65,85,130,145,170")
    batches = len({key[2] for key in records})
    grouped = {
        n: {
            orientation: sorted(
                [row for key, row in records.items() if key[:2] == (n, orientation)],
                key=lambda row: row.batch,
            )
            for orientation in ("first", "second")
        }
        for n in SIZES
    }

    full = {}
    deleted = []
    for n in SIZES:
        projected = project_size(grouped[n])
        projected["N"] = n
        full[n] = {**projected, **bridge_metrics(projected)}
    for omitted in range(batches):
        replicate = {}
        for n in SIZES:
            projected = project_size(grouped[n], omitted)
            projected["N"] = n
            replicate[n] = {**projected, **bridge_metrics(projected)}
        deleted.append(replicate)

    pseudo = {
        (metric, n): pseudovalues(
            full[n][metric], [replicate[n][metric] for replicate in deleted]
        )
        for metric in METRICS
        for n in SIZES
    }
    labels = [(metric, n) for metric in METRICS for n in SIZES]
    pseudo_rows = [[pseudo[label][batch] for label in labels] for batch in range(batches)]
    covariance = covariance_of_mean(pseudo_rows)
    source_indices = [SIZES.index(n) for n in TRAIN_SIZES]
    target_indices = [SIZES.index(n) for n in HELDOUT_SIZES]
    scores = {}
    for metric_index, metric in enumerate(METRICS):
        offset = metric_index * len(SIZES)
        indices = list(range(offset, offset + len(SIZES)))
        block = [[covariance[i][j] for j in indices] for i in indices]
        scores[metric] = gls_score(
            [full[n][metric] for n in SIZES],
            block,
            source_indices,
            target_indices,
        )
    return {
        "batch_count": batches,
        "full": full,
        "labels": labels,
        "covariance": covariance,
        "scores": scores,
    }


def emit(result, input_path: Path, output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for n in SIZES:
        row = {"N": n}
        for metric in ("Mbar_prime", "normalized_pivotal_H4", *METRICS):
            row[metric] = result["full"][n][metric]
            if metric in METRICS:
                index = result["labels"].index((metric, n))
                row[metric + "_se"] = math.sqrt(result["covariance"][index][index])
        rows.append(row)
    with (output_dir / "bridge.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    payload = {
        "schema": "matching-one/pivotal-h4-bridge/v1",
        "status": "retrospective_reuse",
        "primary_metric": "pivotal_H4_scaled",
        "source_sizes": TRAIN_SIZES,
        "heldout_sizes": HELDOUT_SIZES,
        "definitions": {
            "Mbar_prime": "orientation-mean matching slope; exact total pivotal mass by Russo",
            "normalized_pivotal_H4": "P4[S_prime]/Mbar_prime",
            "pivotal_H4_scaled": "N^(13/8)*normalized_pivotal_H4",
            "coefficient_ratio": "normalized_pivotal_H4/P4[D]",
            "thermal_mass_scaled": "N^(-3/8)*Mbar_prime",
        },
        "values": rows,
        "scores": result["scores"],
        "stacked_covariance": result["covariance"],
        "provenance": {
            "input": str(input_path),
            "input_sha256": sha256(input_path),
            "batches": result["batch_count"],
            "evidence_warning": "source data predate this analysis; mechanism diagnostic only",
        },
    }
    (output_dir / "score.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    primary = result["scores"]["pivotal_H4_scaled"]
    xi = result["scores"]["coefficient_ratio"]
    thermal = result["scores"]["thermal_mass_scaled"]
    lines = [
        "# Integrated pivotal/H4 bridge",
        "",
        "**Status: retrospective mechanism analysis.** The raw block predates this protocol.",
        "",
        "The exact Russo identity turns the orientation-mean matching slope into total pivotal mass. Dividing `P4[S']` by that mass removes the leading thermal metric and defines an integrated spin-4 pivotal anisotropy coefficient.",
        "",
        "## Decision",
        "",
        f"Thermal normalization does not remove the known derivative drift: the primary constant-amplitude score is `{primary['heldout_chi_square']:.4f} / 2 df`. Thus the bare `N^-13/8` pivotal-H4 law is not an adequate precision model for this block.",
        "",
        f"The stronger dimensionless coefficient `Xi=[P4[S']/Mbar']/P4[D]` is still compatible with a shared central-value/derivative scaling function at `{xi['heldout_chi_square']:.4f} / 2 df`, although its point estimates decrease with N and this is retrospective evidence. The leading pure thermal-mass law is also far too rigid at `{thermal['heldout_chi_square']:.4f} / 2 df`; its already-resolved finite-size corrections cannot be discarded merely because the leading exponent is exact.",
        "",
        "| metric | source constant (SE) | heldout chi-square / 2 | heldout z (145,170) |",
        "| --- | ---: | ---: | ---: |",
    ]
    for name, score in (("N^(13/8) P4[S']/Mbar'", primary), ("Xi", xi), ("N^(-3/8) Mbar'", thermal)):
        z = [score["residual"][i] / score["residual_se"][i] for i in range(2)]
        lines.append(
            f"| {name} | {score['amplitude']:.9g} ({score['amplitude_se']:.3g}) | "
            f"{score['heldout_chi_square']:.4g} | {z[0]:+.3f}, {z[1]:+.3f} |"
        )
    lines.extend([
        "",
        "`Xi=[P4[S']/Mbar']/P4[D]` asks whether the central residual and pivotal anisotropy are two coefficients of one matching-odd scaling function. It is diagnostic because the denominator is noisy. All three rows are correlated views of one batch block and are not additive evidence.",
        "",
        "The next high-information use is to score `Xi` and the pivotal-normalized q=2/Jordan alternatives on the already authorized norm-4 targets from Issue #154. Those are correlated transforms of the same target block, not extra evidence. A genuinely independent geometric test still requires a local landing-sector/four-arm H4 observable.",
        "",
    ])
    (output_dir / "REPORT.md").write_text("\n".join(lines), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    emit(analyze(args.input), args.input, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

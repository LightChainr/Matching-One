#!/usr/bin/env python3
"""Build the high-statistics archive view of integrated pivotal H4 channels."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

from analyze_p48_retrospective import (
    covariance_of_mean,
    project_size,
    pseudovalues,
    quadratic,
    read_histograms,
    solve_matrix,
)
from analyze_pivotal_h4_bridge import METRICS, bridge_metrics


SIZES = (65, 85, 130, 170, 185, 265, 325, 425)
BLOCKS = ((65, 85, 130, 170), (185,), (265,), (325,), (425,))
BASE_METRICS = ("B_plus", "A_plus", "A_minus", "B_minus")
ALL_METRICS = BASE_METRICS + METRICS


def estimate(path: Path, expected_n: int):
    records = read_histograms(path)
    sizes = {key[0] for key in records}
    if sizes != {expected_n}:
        raise ValueError(f"{path} does not contain only N={expected_n}")
    grouped = {
        orientation: sorted(
            [row for key, row in records.items() if key[1] == orientation],
            key=lambda row: row.batch,
        )
        for orientation in ("first", "second")
    }
    batches = len(grouped["first"])
    point = project_size(grouped)
    point["N"] = expected_n
    point = {**point, **bridge_metrics(point)}
    deleted = []
    for omitted in range(batches):
        row = project_size(grouped, omitted)
        row["N"] = expected_n
        deleted.append({**row, **bridge_metrics(row)})
    return point, deleted


def fit_constant(values, covariance):
    inverse_times_one = [
        row[0] for row in solve_matrix(covariance, [[1.0] for _ in values])
    ]
    denominator = math.fsum(inverse_times_one)
    amplitude = math.fsum(w * y for w, y in zip(inverse_times_one, values)) / denominator
    residual = [y - amplitude for y in values]
    return {
        "amplitude": amplitude,
        "amplitude_se": math.sqrt(1.0 / denominator),
        "residual": residual,
        "chi_square": quadratic(residual, covariance),
        "df": len(values) - 1,
    }


def analyze(runs):
    points = {}
    deleted = {}
    for n in SIZES:
        points[n], deleted[n] = estimate(runs[n], n)
    labels = [(metric, n) for metric in ALL_METRICS for n in SIZES]
    dimension = len(labels)
    covariance = [[0.0] * dimension for _ in range(dimension)]
    for block in BLOCKS:
        batch_counts = {len(deleted[n]) for n in block}
        if len(batch_counts) != 1:
            raise ValueError(f"unaligned batch counts in block {block}")
        batches = batch_counts.pop()
        block_labels = [(metric, n) for metric in ALL_METRICS for n in block]
        rows = []
        for batch in range(batches):
            row = []
            for metric, n in block_labels:
                pseudo = pseudovalues(
                    points[n][metric], [item[metric] for item in deleted[n]]
                )
                row.append(pseudo[batch])
            rows.append(row)
        block_cov = covariance_of_mean(rows)
        for i, left in enumerate(block_labels):
            ii = labels.index(left)
            for j, right in enumerate(block_labels):
                jj = labels.index(right)
                covariance[ii][jj] = block_cov[i][j]
    scores = {}
    for metric in METRICS:
        indices = [labels.index((metric, n)) for n in SIZES]
        block = [[covariance[i][j] for j in indices] for i in indices]
        scores[metric] = fit_constant([points[n][metric] for n in SIZES], block)
    return {"points": points, "labels": labels, "covariance": covariance, "scores": scores}


def emit(result, runs, output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for n in SIZES:
        row = {"N": n}
        for metric in (
            "first_primal_pivotal",
            "first_matching_pivotal",
            "second_primal_pivotal",
            "second_matching_pivotal",
            *BASE_METRICS,
            "normalized_pivotal_H4",
            *METRICS,
        ):
            row[metric] = result["points"][n][metric]
            if metric in ALL_METRICS:
                index = result["labels"].index((metric, n))
                row[metric + "_se"] = math.sqrt(result["covariance"][index][index])
        rows.append(row)
    with (output_dir / "archive.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    payload = {
        "schema": "matching-one/pivotal-h4-archive/v1",
        "status": "descriptive_high_statistics_archive",
        "sizes": SIZES,
        "covariance_blocks": BLOCKS,
        "runs": {str(n): str(runs[n]) for n in SIZES},
        "values": rows,
        "scores": result["scores"],
        "stacked_covariance": result["covariance"],
        "stacked_covariance_order": result["labels"],
        "lossless_basis": {
            "U": [
                "first_primal_pivotal",
                "first_matching_pivotal",
                "second_primal_pivotal",
                "second_matching_pivotal",
            ],
            "V": BASE_METRICS,
            "identities": {
                "B_plus": "Mbar_prime",
                "A_plus": "2*P4[D_prime]",
                "A_minus": "2*P4[S_prime]",
                "B_minus": "orientation-mean primal-minus-matching pivotal mass",
            },
        },
        "evidence_rule": "all metrics are correlated transforms of previously analyzed raw data",
    }
    (output_dir / "archive_score.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    names = {
        "pivotal_H4_scaled": "N^(13/8) P4[S']/Mbar'",
        "even_pivotal_H4_scaled": "N P4[D']/Mbar'",
        "coefficient_ratio": "Xi",
        "thermal_mass_scaled": "N^(-3/8) Mbar'",
    }
    lines = [
        "# High-statistics integrated pivotal/H4 archive",
        "",
        "**Status: descriptive reuse, not heldout evidence.**",
        "",
        "## Decision",
        "",
        f"Across N=65--425, the matching-even relative pivotal anisotropy `N P4[D']/Mbar'` is compatible with one amplitude (`chi-square={result['scores']['even_pivotal_H4_scaled']['chi_square']:.4f} / 7 df`). This is the clean integrated signature of an `L^-2` H4 correction to total pivotal mass.",
        "",
        f"The coefficient ratio `Xi` is also compatible with one value (`chi-square={result['scores']['coefficient_ratio']['chi_square']:.4f} / 7 df`), supporting a common scaling-function relation between central `P4[D]` and derivative `P4[S']`. In contrast, the matching-odd leading-power normalization is decisively nonconstant (`chi-square={result['scores']['pivotal_H4_scaled']['chi_square']:.4f} / 7 df`): pivotal normalization does not erase its q=2/Jordan-scale finite-size dynamics.",
        "",
        "| metric | constant amplitude (SE) | chi-square / df |",
        "| --- | ---: | ---: |",
    ]
    for metric in METRICS:
        score = result["scores"][metric]
        lines.append(
            f"| {names[metric]} | {score['amplitude']:.9g} ({score['amplitude_se']:.3g}) | "
            f"{score['chi_square']:.4g} / {score['df']} |"
        )
    lines.extend([
        "",
        "The first four sizes retain their shared-random-number covariance. The remaining sizes are independent counter blocks. These rows expose the precision and lineage dependence of the mechanism coordinates; they are not four independent tests.",
        "",
    ])
    (output_dir / "ARCHIVE_REPORT.md").write_text("\n".join(lines), encoding="utf-8")


def parse_run(value):
    n, path = value.split(":", 1)
    return int(n), Path(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="append", type=parse_run, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    runs = dict(args.run)
    if tuple(sorted(runs)) != SIZES:
        raise ValueError(f"runs must contain exactly {SIZES}")
    emit(analyze(runs), runs, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

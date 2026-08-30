#!/usr/bin/env python3
"""Orthogonally split frozen charged currents into net and activity modes."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Mapping, Optional, Sequence

from reveal_n65_charged_source_archive import covariance_of_mean, quadratic_n


ORIENTATIONS = ("first", "second")
CHANNELS = ("A", "D")
TRIPLET = ("W", "J_minus", "J_plus")


def correlation(covariance: Sequence[Sequence[float]]) -> list[list[float]]:
    scale = [math.sqrt(value) for value in (covariance[i][i] for i in range(len(covariance)))]
    return [[covariance[i][j] / (scale[i] * scale[j])
             for j in range(len(covariance))] for i in range(len(covariance))]


def symmetric_eigen(matrix: Sequence[Sequence[float]]) -> tuple[list[float], list[list[float]]]:
    """Jacobi eigenpairs for a small real symmetric matrix, ascending in value."""
    n = len(matrix)
    work = [list(row) for row in matrix]
    vectors = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    for _ in range(100):
        p, q = max(((i, j) for i in range(n) for j in range(i + 1, n)),
                   key=lambda pair: abs(work[pair[0]][pair[1]]))
        if abs(work[p][q]) < 1e-15:
            break
        angle = 0.5 * math.atan2(2.0 * work[p][q], work[q][q] - work[p][p])
        cosine, sine = math.cos(angle), math.sin(angle)
        app, aqq, apq = work[p][p], work[q][q], work[p][q]
        for k in range(n):
            if k in (p, q):
                continue
            akp, akq = work[k][p], work[k][q]
            work[k][p] = work[p][k] = cosine * akp - sine * akq
            work[k][q] = work[q][k] = sine * akp + cosine * akq
        work[p][p] = cosine * cosine * app - 2 * sine * cosine * apq + sine * sine * aqq
        work[q][q] = sine * sine * app + 2 * sine * cosine * apq + cosine * cosine * aqq
        work[p][q] = work[q][p] = 0.0
        for k in range(n):
            vkp, vkq = vectors[k][p], vectors[k][q]
            vectors[k][p] = cosine * vkp - sine * vkq
            vectors[k][q] = sine * vkp + cosine * vkq
    pairs = sorted((work[i][i], [vectors[row][i] for row in range(n)]) for i in range(n))
    values, columns = [], []
    for value, vector in pairs:
        pivot = max(range(n), key=lambda i: abs(vector[i]))
        if vector[pivot] < 0:
            vector = [-entry for entry in vector]
        values.append(value)
        columns.append(vector)
    return values, columns


def block_score(
    vector: Sequence[float], covariance: Sequence[Sequence[float]]
) -> dict[str, object]:
    standard_error = [math.sqrt(max(0.0, covariance[i][i])) for i in range(3)]
    z = [vector[i] / standard_error[i] for i in range(3)]
    subsets = {
        "W_only": (0,), "net_only": (1,), "activity_only": (2,),
        "W_net": (0, 1), "W_activity": (0, 2),
        "net_activity": (1, 2), "full": (0, 1, 2),
    }
    quadratic = {}
    for name, indices in subsets.items():
        subvector = [vector[i] for i in indices]
        subcovariance = [[covariance[i][j] for j in indices] for i in indices]
        quadratic[name] = quadratic_n(subvector, subcovariance)
    quadratic["increment_net_given_W_activity"] = quadratic["full"] - quadratic["W_activity"]
    quadratic["increment_activity_given_W_net"] = quadratic["full"] - quadratic["W_net"]
    quadratic["increment_W_given_net_activity"] = quadratic["full"] - quadratic["net_activity"]

    corr = correlation(covariance)
    eigenvalues, eigenvectors = symmetric_eigen(corr)
    modes = []
    for eigenvalue, eigenvector in zip(eigenvalues, eigenvectors):
        projection = math.fsum(a * b for a, b in zip(z, eigenvector))
        modes.append({
            "eigenvalue": eigenvalue,
            "standardized_eigenvector": eigenvector,
            "signal_projection": projection,
            "quadratic_contribution": projection * projection / eigenvalue,
        })
    return {
        "order": list(TRIPLET), "value": list(vector),
        "standard_error": standard_error, "marginal_z": z,
        "covariance": [list(row) for row in covariance],
        "correlation": corr, "quadratic": quadratic,
        "correlation_eigenmodes": modes,
    }


def score(reveal_path: Path) -> dict[str, object]:
    parent = json.loads(reveal_path.read_text(encoding="utf-8"))
    if parent.get("schema") != "matching-one/N65-F3-charged-source-archive-reveal/v1":
        raise ValueError("unexpected charged-source reveal schema")
    parent_order = parent["joint_estimate"]["order"]
    batches = parent["joint_estimate"]["batch_values"]
    transformed_order = [f"{orientation}_{channel}_{name}"
                         for orientation in ORIENTATIONS for channel in CHANNELS
                         for name in TRIPLET]
    transformed_order += [f"second_minus_first_{channel}_{name}"
                          for channel in CHANNELS for name in TRIPLET]
    rows = []
    transformed_batches = []
    for row in batches:
        parent_values = row["values"]
        values: dict[str, float] = {}
        for orientation in ORIENTATIONS:
            for channel in CHANNELS:
                birth = parent_values[f"{orientation}_J_{channel}_birth"]
                exit_ = parent_values[f"{orientation}_J_{channel}_exit"]
                values[f"{orientation}_{channel}_W"] = parent_values[f"{orientation}_W_{channel}"]
                values[f"{orientation}_{channel}_J_minus"] = birth - exit_
                values[f"{orientation}_{channel}_J_plus"] = birth + exit_
        for channel in CHANNELS:
            for name in TRIPLET:
                values[f"second_minus_first_{channel}_{name}"] = (
                    values[f"second_{channel}_{name}"] - values[f"first_{channel}_{name}"]
                )
        vector = [values[name] for name in transformed_order]
        rows.append(vector)
        transformed_batches.append({"batch": row["batch"], "values": values})

    covariance = covariance_of_mean(rows)
    means = [math.fsum(row[i] for row in rows) / len(rows)
             for i in range(len(transformed_order))]
    errors = [math.sqrt(max(0.0, covariance[i][i])) for i in range(len(transformed_order))]
    by_channel = {}
    for channel in CHANNELS:
        names = [f"second_minus_first_{channel}_{name}" for name in TRIPLET]
        indices = [transformed_order.index(name) for name in names]
        vector = [means[i] for i in indices]
        block = [[covariance[i][j] for j in indices] for i in indices]
        by_channel[channel] = block_score(vector, block)

    a = by_channel["A"]["quadratic"]
    d = by_channel["D"]["quadratic"]
    return {
        "schema": "matching-one/N65-F3-charged-activity-net/v1",
        "status": "same-archive orthogonal activity/net decomposition; no new simulation",
        "source": {
            "parent_commit": "83207c7",
            "parent_certificate": str(reveal_path),
            "archive_commit": parent["source"]["archive_commit"],
            "dependency_group": "N65 projective-birth 20k shared-counter block at 1714141",
            "batches": parent["source"]["batches"],
            "samples_per_shape": parent["source"]["samples_per_shape"],
        },
        "frozen_transform": {
            "W": "rank-one charged plateau susceptibility",
            "J_minus": "J_birth-J_exit=dW/dp, net derivative",
            "J_plus": "J_birth+J_exit, common activity",
            "matrix_from_W_birth_exit": [[1, 0, 0], [0, 1, -1], [0, 1, 1]],
            "selection_changed": False,
        },
        "orientation_contrast": by_channel,
        "mechanism_reading": {
            "A": (
                f"net derivative carries the contrast: marginal z={by_channel['A']['marginal_z'][1]:.6g}, "
                f"net-only quadratic={a['net_only']:.6g}, activity-only={a['activity_only']:.6g}; "
                f"adding activity after W+net contributes only {a['increment_activity_given_W_net']:.6g}"
            ),
            "D": (
                f"no resolved orientation contrast: full quadratic={d['full']:.6g}/3; "
                f"net-only={d['net_only']:.6g}, activity-only={d['activity_only']:.6g}"
            ),
            "classification": "A net-timing counterflow; no D mode; common activity cancels between orientations",
        },
        "p334_crosswalk": {
            "common_activity": "J_plus is the positive birth-plus-exit traffic used by the #334 orbit-flux phase diagram",
            "counterflow": "J_minus is the small source-sink residual; here its A/B1 orientation contrast survives while J_plus does not",
            "sector_difference": "this reveal is intra-axis A/B1 and intra-diagonal D/B2; #334 77aa3fe is axis-versus-diagonal H4 orbit composition",
            "same_archive_coordinates": "381984d conditional line sorting and this score reuse 1714141 and are not independent evidence",
            "different_dependency_group": "77aa3fe exact N13/N17 orbit tables are a qualitative mechanism crosswalk only; no quadratic or significance is pooled",
        },
        "joint_estimate": {
            "order": transformed_order, "mean": means, "standard_error": errors,
            "covariance": covariance, "batch_values": transformed_batches,
        },
        "claim_boundary": (
            "fixed linear decomposition of one exploratory 20k dependency block; "
            "mechanism localization, not additive evidence or large-N scaling"
        ),
    }


def render_markdown(payload: Mapping[str, object]) -> str:
    lines = [
        "# Charged activity/net decomposition of the existing N65 archive", "",
        "The frozen transform is `J_minus=birth-exit` and `J_plus=birth+exit`; no source or data block changed.", "",
        "| channel | coordinate | contrast | SE | z | marginal quadratic |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for channel in CHANNELS:
        block = payload["orientation_contrast"][channel]
        for i, name in enumerate(TRIPLET):
            lines.append(f"| {channel} | {name} | {block['value'][i]:.12g} | "
                         f"{block['standard_error'][i]:.3g} | {block['marginal_z'][i]:.3f} | "
                         f"{block['quadratic'][{'W':'W_only','J_minus':'net_only','J_plus':'activity_only'}[name]]:.4g} |")
    a = payload["orientation_contrast"]["A"]["quadratic"]
    d = payload["orientation_contrast"]["D"]["quadratic"]
    lines += ["", f"A full triplet remains `{a['full']:.4g} / 3 df`; net alone carries `{a['net_only']:.4g}`, while activity alone carries `{a['activity_only']:.4g}`.",
              f"Conditioned on `(W,net)`, activity adds only `{a['increment_activity_given_W_net']:.4g}`. A is therefore a net-timing counterflow, not a common-activity amplitude.", "",
              f"D gives `{d['full']:.4g} / 3 df`; neither net nor activity resolves an orientation response.", "",
              "## Covariance eigenmodes", ""]
    for channel in CHANNELS:
        lines.append(f"{channel}, standardized `(W,J_minus,J_plus)` correlation modes:")
        for mode in payload["orientation_contrast"][channel]["correlation_eigenmodes"]:
            vector = ", ".join(f"{value:+.3f}" for value in mode["standardized_eigenvector"])
            lines.append(f"- lambda={mode['eigenvalue']:.4g}, vector=({vector}), quadratic contribution={mode['quadratic_contribution']:.4g}.")
        lines.append("")
    lines += ["## #334 crosswalk", "",
              "`J_plus` is the common positive activity and `J_minus` the source-sink residual. The A result matches the counterflow morphology: the large common activity cancels between orientations and information remains in the smaller net timing current.", "",
              "The sectors and dependency groups remain distinct. This A/B1 coordinate is intra-axis; the exact #334 N13/N17 result is axis-versus-diagonal H4 orbit composition. The exact tables provide a mechanism crosswalk only. The conditional N65 line-sorting score reuses this same 1714141 block and is not additive evidence.", ""]
    return "\n".join(lines)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reveal", type=Path, required=True)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--markdown", type=Path)
    args = parser.parse_args(argv)
    payload = score(args.reveal)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.markdown:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(render_markdown(payload), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

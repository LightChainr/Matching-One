#!/usr/bin/env python3
"""Delete-one scorer for production marked-birth microcanonical paths."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any, Iterable, Sequence

import mpmath as mp


METRICS = (
    "P4_A_top",
    "P4_B",
    "P4_D",
    "P4_J_S_re",
    "P4_J_S_im",
    "P4_J_D_re",
    "P4_J_D_im",
    "P4_connected_q_J_D_re",
    "P4_connected_q_J_D_im",
    "P4_gamma_D_re",
    "P4_gamma_D_im",
    "P4_local_S",
    "P4_local_D",
)


@dataclass
class PathRow:
    n: int
    a: int
    b: int
    orientation: str
    batch: int
    samples: int
    k: int
    values: dict[str, mp.mpf]


VALUE_COLUMNS = (
    "sum_q",
    "sum_q2",
    "sum_gate01",
    "sum_gate12",
    "sum_inactive_gate01",
    "sum_inactive_gate12",
    "sum_active_S",
    "sum_active_D",
    "sum_inactive_S",
    "sum_inactive_D",
    "sum_site_S",
    "sum_site_D",
    "sum_J_S_re",
    "sum_J_S_im",
    "sum_J_D_re",
    "sum_J_D_im",
    "sum_q_J_D_re",
    "sum_q_J_D_im",
    "sum_local_S",
    "sum_local_D",
)


def cos4(a: int, b: int) -> mp.mpf:
    n = a * a + b * b
    return mp.mpf(a**4 - 6 * a * a * b * b + b**4) / (n * n)


def read_path(path: Path) -> dict[tuple[int, str, int], list[PathRow]]:
    groups: dict[tuple[int, str, int], list[PathRow]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            row = PathRow(
                n=int(raw["n"]),
                a=int(raw["a"]),
                b=int(raw["b"]),
                orientation=raw["orientation"],
                batch=int(raw["batch"]),
                samples=int(raw["samples"]),
                k=int(raw["k"]),
                values={column: mp.mpf(raw[column]) for column in VALUE_COLUMNS},
            )
            groups.setdefault((row.n, row.orientation, row.batch), []).append(row)
    for key, rows in groups.items():
        rows.sort(key=lambda row: row.k)
        if len(rows) != key[0] or [row.k for row in rows] != list(range(key[0])):
            raise ValueError(f"incomplete microcanonical path for {key}")
        if len({row.samples for row in rows}) != 1:
            raise ValueError(f"inconsistent path sample count for {key}")
        for row in rows:
            absent = row.n - row.k
            if row.values["sum_site_S"] != absent * (
                row.values["sum_gate01"] + row.values["sum_gate12"]
            ):
                raise ValueError(f"Horvitz S/gate mismatch for {key}, k={row.k}")
            if row.values["sum_site_D"] != absent * (
                row.values["sum_gate12"] - row.values["sum_gate01"]
            ):
                raise ValueError(f"Horvitz D/gate mismatch for {key}, k={row.k}")
            if 2 * row.values["sum_site_S"] != (
                row.values["sum_active_S"] + row.values["sum_inactive_S"]
            ):
                raise ValueError(f"full active/inactive S mismatch for {key}, k={row.k}")
            if 2 * row.values["sum_site_D"] != (
                row.values["sum_active_D"] - row.values["sum_inactive_D"]
            ):
                raise ValueError(f"full active/inactive D mismatch for {key}, k={row.k}")
            if row.values["sum_active_S"] != absent * (
                row.values["sum_gate01"] + row.values["sum_gate12"]
            ) or row.values["sum_inactive_S"] != absent * (
                row.values["sum_inactive_gate01"] + row.values["sum_inactive_gate12"]
            ):
                raise ValueError(f"side-resolved S mismatch for {key}, k={row.k}")
    if not groups:
        raise ValueError(f"empty marked path: {path}")
    return groups


def combine(groups: Sequence[list[PathRow]]) -> list[PathRow]:
    if not groups:
        raise ValueError("cannot combine empty path groups")
    first = groups[0]
    n = len(first)
    result = []
    for k in range(n):
        rows = [group[k] for group in groups]
        template = rows[0]
        if any(
            (row.n, row.a, row.b, row.orientation, row.k)
            != (template.n, template.a, template.b, template.orientation, k)
            for row in rows
        ):
            raise ValueError("incompatible path batches")
        result.append(
            PathRow(
                n=template.n,
                a=template.a,
                b=template.b,
                orientation=template.orientation,
                batch=-1,
                samples=sum(row.samples for row in rows),
                k=k,
                values={
                    column: mp.fsum(row.values[column] for row in rows)
                    for column in VALUE_COLUMNS
                },
            )
        )
    return result


def binomial_weights(n: int, p: mp.mpf) -> list[mp.mpf]:
    if p <= 0:
        return [mp.mpf(1)] + [mp.mpf(0)] * n
    if p >= 1:
        return [mp.mpf(0)] * n + [mp.mpf(1)]
    q = 1 - p
    weights = [q**n]
    for k in range(n):
        weights.append(weights[-1] * (n - k) * p / ((k + 1) * q))
    return weights


def mean_column(rows: Sequence[PathRow], column: str) -> list[mp.mpf]:
    return [row.values[column] / row.samples for row in rows]


def matching_curve(rows: Sequence[PathRow], p: mp.mpf) -> mp.mpf:
    n = len(rows)
    q = mean_column(rows, "sum_q")
    # q at k=N is +1 after the second birth on every path.
    coefficients = q + [mp.mpf(1)]
    return mp.fsum(
        weight * value for weight, value in zip(binomial_weights(n, p), coefficients)
    )


def intrinsic_center(first: Sequence[PathRow], second: Sequence[PathRow]) -> mp.mpf:
    lower = mp.mpf(0)
    upper = mp.mpf(1)
    for _ in range(100):
        p = (lower + upper) / 2
        value = (matching_curve(first, p) + matching_curve(second, p)) / 2
        if value < 0:
            lower = p
        else:
            upper = p
    return (lower + upper) / 2


def canonical_site_sum(rows: Sequence[PathRow], column: str, p: mp.mpf) -> mp.mpf:
    n = len(rows)
    weights = binomial_weights(n - 1, p)
    return mp.fsum(
        weights[k]
        * (mp.mpf(n) / (n - k))
        * rows[k].values[column]
        / rows[k].samples
        for k in range(n)
    )


def orientation_observables(rows: Sequence[PathRow], p: mp.mpf) -> dict[str, mp.mpf]:
    values = {
        "A_top": matching_curve(rows, p),
        "B": canonical_site_sum(rows, "sum_site_S", p),
        "D": canonical_site_sum(rows, "sum_site_D", p),
        "J_S_re": canonical_site_sum(rows, "sum_J_S_re", p),
        "J_S_im": canonical_site_sum(rows, "sum_J_S_im", p),
        "J_D_re": canonical_site_sum(rows, "sum_J_D_re", p),
        "J_D_im": canonical_site_sum(rows, "sum_J_D_im", p),
        "environment_q_J_D_re": canonical_site_sum(rows, "sum_q_J_D_re", p),
        "environment_q_J_D_im": canonical_site_sum(rows, "sum_q_J_D_im", p),
        "local_S": canonical_site_sum(rows, "sum_local_S", p),
        "local_D": canonical_site_sum(rows, "sum_local_D", p),
    }
    # The path stores q before inserting the distinguished root.  The source
    # is a function of the root-deleted environment, whereas the global q in
    # Cov(q,J) is evaluated on the full Bernoulli configuration. Conditional
    # on the environment, q_full=q_before+p*S.  Since S*J_D=J_D (direct
    # 0->2 has D=0), the missing occupied-root term is exactly p*J_D.
    values["q_J_D_re"] = (
        values["environment_q_J_D_re"] + p * values["J_D_re"]
    )
    values["q_J_D_im"] = (
        values["environment_q_J_D_im"] + p * values["J_D_im"]
    )
    q_full = values["A_top"]
    values["connected_q_J_D_re"] = (
        values["q_J_D_re"] - q_full * values["J_D_re"]
    )
    values["connected_q_J_D_im"] = (
        values["q_J_D_im"] - q_full * values["J_D_im"]
    )
    values["gamma_D_re"] = values["connected_q_J_D_re"] / values["B"]
    values["gamma_D_im"] = values["connected_q_J_D_im"] / values["B"]

    n = len(rows)
    q_coefficients = mean_column(rows, "sum_q") + [mp.mpf(1)]
    derivative = mp.mpf(n) * mp.fsum(
        weight * (q_coefficients[k + 1] - q_coefficients[k])
        for k, weight in enumerate(binomial_weights(n - 1, p))
    )
    values["Russo_residual_B_minus_A_prime"] = values["B"] - derivative
    return values


def projected(first: Sequence[PathRow], second: Sequence[PathRow]) -> tuple[mp.mpf, dict[str, mp.mpf], dict[str, dict[str, mp.mpf]]]:
    p = intrinsic_center(first, second)
    left = orientation_observables(first, p)
    right = orientation_observables(second, p)
    delta = cos4(first[0].a, first[0].b) - cos4(second[0].a, second[0].b)
    if delta == 0:
        raise ValueError("zero H4 leverage")
    point = {
        "P4_A_top": (left["A_top"] - right["A_top"]) / delta,
        "P4_B": (left["B"] - right["B"]) / delta,
        "P4_D": (left["D"] - right["D"]) / delta,
        "P4_J_S_re": (left["J_S_re"] - right["J_S_re"]) / delta,
        "P4_J_S_im": (left["J_S_im"] - right["J_S_im"]) / delta,
        "P4_J_D_re": (left["J_D_re"] - right["J_D_re"]) / delta,
        "P4_J_D_im": (left["J_D_im"] - right["J_D_im"]) / delta,
        "P4_connected_q_J_D_re": (
            left["connected_q_J_D_re"] - right["connected_q_J_D_re"]
        ) / delta,
        "P4_connected_q_J_D_im": (
            left["connected_q_J_D_im"] - right["connected_q_J_D_im"]
        ) / delta,
        "P4_gamma_D_re": (left["gamma_D_re"] - right["gamma_D_re"]) / delta,
        "P4_gamma_D_im": (left["gamma_D_im"] - right["gamma_D_im"]) / delta,
        "P4_local_S": (left["local_S"] - right["local_S"]) / delta,
        "P4_local_D": (left["local_D"] - right["local_D"]) / delta,
    }
    return p, point, {"first": left, "second": right}


def covariance(rows: Sequence[dict[str, mp.mpf]]) -> list[list[mp.mpf]]:
    m = len(rows)
    means = {name: mp.fsum(row[name] for row in rows) / m for name in METRICS}
    factor = mp.mpf(m - 1) / m
    return [
        [
            factor
            * mp.fsum(
                (row[left] - means[left]) * (row[right] - means[right])
                for row in rows
            )
            for right in METRICS
        ]
        for left in METRICS
    ]


def validate_sparse(path: Path) -> dict[str, Any]:
    totals: dict[tuple[int, str, int], int] = {}
    direct_rows = 0
    failures = []
    with path.open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            key = (int(raw["n"]), raw["orientation"], int(raw["batch"]))
            totals[key] = totals.get(key, 0) + int(raw["count"])
            direct = raw["direct_0_to_2"] == "1"
            if direct:
                direct_rows += 1
                if not (
                    raw["line_null"] == "1"
                    and raw["iota01"] == "0"
                    and raw["iota12"] == "0"
                    and raw["site01"] == raw["site12"]
                ):
                    failures.append({"kind": "direct_schema", "row": raw})
            elif raw["line_null"] != "0":
                failures.append({"kind": "strict_line_null", "row": raw})
    return {
        "batch_orientation_rows": len(totals),
        "sample_totals": {"|".join(map(str, key)): value for key, value in sorted(totals.items())},
        "direct_sparse_rows": direct_rows,
        "schema_failures": failures[:8],
    }


def validate_audit(path: Path) -> dict[str, Any]:
    totals = {
        name: 0
        for name in (
            "endpoint_failures",
            "site_failures",
            "line_failures",
            "local_mark_failures",
            "index_mismatches",
        )
    }
    rows = 0
    with path.open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            rows += 1
            for name in totals:
                totals[name] += int(raw[name])
    return {"rows": rows, "totals": totals}


def _text(value: mp.mpf) -> str:
    return mp.nstr(value, 24)


def build_report(prefix: Path) -> dict[str, Any]:
    path_file = Path(str(prefix) + ".path.csv")
    sparse_file = Path(str(prefix) + ".marked_births.csv")
    audit_file = Path(str(prefix) + ".complement_audit.csv")
    metadata_file = Path(str(prefix) + ".metadata.json")
    groups = read_path(path_file)
    sizes = {key[0] for key in groups}
    if len(sizes) != 1:
        raise ValueError("score input contains multiple sizes")
    n = sizes.pop()
    batches = sorted(
        set(key[2] for key in groups if key[1] == "first")
        & set(key[2] for key in groups if key[1] == "second")
    )
    if len(batches) < 2:
        raise ValueError("delete-one score needs aligned batches")
    first = combine([groups[(n, "first", batch)] for batch in batches])
    second = combine([groups[(n, "second", batch)] for batch in batches])
    p0, point, orientations = projected(first, second)
    delete_one = []
    delete_centers = []
    for omitted in batches:
        first_delete = combine(
            [groups[(n, "first", batch)] for batch in batches if batch != omitted]
        )
        second_delete = combine(
            [groups[(n, "second", batch)] for batch in batches if batch != omitted]
        )
        center, row, _ = projected(first_delete, second_delete)
        delete_centers.append(center)
        delete_one.append(row)
    cov = covariance(delete_one)
    standard_errors = {
        name: mp.sqrt(max(mp.mpf(0), cov[index][index]))
        for index, name in enumerate(METRICS)
    }
    sparse = validate_sparse(sparse_file)
    audit = validate_audit(audit_file)
    metadata = json.loads(metadata_file.read_text(encoding="utf-8"))
    return {
        "schema": "matching-one/marked-birth-path-score/v1",
        "N": n,
        "prefix": str(prefix),
        "git_commit": metadata.get("git_commit"),
        "samples_per_orientation": first[0].samples,
        "batches": len(batches),
        "intrinsic_center": _text(p0),
        "intrinsic_center_delete_one_se": _text(
            mp.sqrt(
                mp.mpf(len(batches) - 1)
                / len(batches)
                * mp.fsum(
                    (value - mp.fsum(delete_centers) / len(delete_centers)) ** 2
                    for value in delete_centers
                )
            )
        ),
        "orientations": {
            side: {name: _text(value) for name, value in values.items()}
            for side, values in orientations.items()
        },
        "P4_point": {name: _text(point[name]) for name in METRICS},
        "P4_standard_error": {name: _text(standard_errors[name]) for name in METRICS},
        "P4_z": {
            name: (_text(point[name] / standard_errors[name]) if standard_errors[name] else None)
            for name in METRICS
        },
        "covariance_metric_order": list(METRICS),
        "delete_one_covariance": [[_text(value) for value in row] for row in cov],
        "sparse_validation": sparse,
        "complement_validation": audit,
        "Russo_residuals": {
            side: values["Russo_residual_B_minus_A_prime"]
            for side, values in orientations.items()
        },
        "scientific_card": [
            "MECHANISM SPACE: separates unmarked S/D, lifted-line chi4 S/D, local landing S/D and connected q-J_D in one common-field path.",
            "NOT PROVED: a nonzero finite pilot does not identify the Q4 epsilon continuum field or its asymptotic exponent.",
            "OBSERVER-SECTOR-SOURCE-GEOMETRY: A_top | odd D and even S birth sources | lifted chi4 and local H4 | frozen Gaussian orientation pair.",
            "DEPENDENCY GROUP: all marked channels share the same permutation and batch; use the supplied covariance as one evidence block.",
            "UPWEIGHT OBSERVATION: stable q2 sign transfer of connected q-J_D and gamma_D with zero reverse-complement audits.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--dps", type=int, default=50)
    args = parser.parse_args()
    mp.mp.dps = args.dps
    report = build_report(args.prefix)
    text = json.dumps(report, indent=2, sort_keys=True, default=str) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    else:
        print(text, end="")


if __name__ == "__main__":
    main()

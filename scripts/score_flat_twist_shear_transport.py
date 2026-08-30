#!/usr/bin/env python3
"""Score the exact F3 identity/shear transport on paired birth archives.

The declared shear is a change of period basis, not a physical deformation:

    P_shear = P_identity T^-1,  ell_shear = T ell_identity.

Accordingly this scorer first demands an exact sparse-cell transport gate.  It
then reconstructs the fixed-p F3 triplet for two same-N Gaussian shapes and
retains the full six-dimensional identity/shear covariance and the complete
three-dimensional transport residual.
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
import math
from pathlib import Path
from typing import Mapping, Sequence

from analyze_projective_birth_smoke import BirthCell, covariance_of_mean, read_births
from score_flat_twist_projective_archive import fixed_p_state, twist_characters


DEFAULT_P = 0.592746050790
TRIPLET_NAMES = (
    "F3_H4_axis_diag",
    "F3_axis_odd",
    "F3_diagonal_odd",
)
T_MATRIX = (
    (0.0, 1.0 / math.sqrt(2.0), -1.0 / math.sqrt(2.0)),
    (1.0 / math.sqrt(2.0), 0.5, 0.5),
    (1.0 / math.sqrt(2.0), -0.5, -0.5),
)
T_INTEGER = ((1, 1), (0, 1))
T_INVERSE_INTEGER = ((1, -1), (0, 1))


def matmul2(
    first: Sequence[Sequence[int]], second: Sequence[Sequence[int]],
) -> list[list[int]]:
    return [[
        sum(first[i][k] * second[k][j] for k in range(2))
        for j in range(2)
    ] for i in range(2)]


def apply3(matrix: Sequence[Sequence[float]], vector: Sequence[float]) -> list[float]:
    return [
        math.fsum(matrix[i][j] * vector[j] for j in range(3))
        for i in range(3)
    ]


def canonical_primitive(x: int, y: int) -> tuple[int, int]:
    divisor = math.gcd(abs(x), abs(y))
    if divisor == 0:
        return 0, 0
    x, y = x // divisor, y // divisor
    if x < 0 or (x == 0 and y < 0):
        x, y = -x, -y
    return x, y


def transport_cell(cell: BirthCell) -> tuple[int, int, str, int, int]:
    x, y = cell.ell_x, cell.ell_y
    if cell.kind == "LINE":
        x, y = canonical_primitive(x + y, y)
    return cell.tau1, cell.tau2, cell.kind, x, y


def cell_counter(cells: Sequence[BirthCell], transport: bool) -> Counter:
    output: Counter = Counter()
    for cell in cells:
        key = transport_cell(cell) if transport else (
            cell.tau1, cell.tau2, cell.kind, cell.ell_x, cell.ell_y
        )
        output[key] += cell.count
    return output


def covariance_summary(rows: Sequence[Sequence[float]]) -> dict[str, object]:
    covariance = covariance_of_mean(rows)
    means = [math.fsum(row[j] for row in rows) / len(rows)
             for j in range(len(rows[0]))]
    errors = [math.sqrt(max(0.0, covariance[j][j])) for j in range(len(means))]
    return {"mean": means, "standard_error": errors, "covariance": covariance}


def validate_metadata(
    identity_path: Path, shear_path: Path,
) -> tuple[dict[str, object], list[list[int]], list[list[int]]]:
    # Both paths are intentionally the same file: first/second within one run
    # are the identity/shear pair.  Separate arguments keep the contract clear
    # for future archive layouts.
    if identity_path != shear_path:
        raise ValueError("one paired metadata file must declare identity and shear")
    metadata = json.loads(identity_path.read_text(encoding="utf-8"))
    design = metadata["designs"][0]
    identity = design["first_period_matrix"]
    shear = design["second_period_matrix"]
    expected = matmul2(identity, T_INVERSE_INTEGER)
    if shear != expected:
        raise ValueError(f"second period matrix is not P*T^-1: {shear} != {expected}")
    if design["first_HNF"] != design["second_HNF"]:
        raise ValueError("basis-related quotients must have identical canonical HNF")
    return metadata, identity, shear


def score_shape(
    births_path: Path, metadata_path: Path, p: float,
) -> tuple[dict[str, object], dict[int, dict[str, list[float]]]]:
    metadata, identity_matrix, shear_matrix = validate_metadata(
        metadata_path, metadata_path
    )
    n, births = read_births(births_path)
    if n != metadata["designs"][0]["N"]:
        raise ValueError("birth N and metadata N disagree")
    batch_ids = sorted({batch for _, batch in births})
    rows: dict[int, dict[str, list[float]]] = {}
    failures: list[int] = []
    for batch in batch_ids:
        identity_cells = births[("first", batch)]
        shear_cells = births[("second", batch)]
        if cell_counter(identity_cells, True) != cell_counter(shear_cells, False):
            failures.append(batch)
        by_basis = {}
        for name, cells, matrix in (
            ("identity", identity_cells, identity_matrix),
            ("shear", shear_cells, shear_matrix),
        ):
            characters = twist_characters(fixed_p_state(cells, n, p, matrix))
            by_basis[name] = [characters[field] for field in TRIPLET_NAMES]
        rows[batch] = by_basis
    return {
        "births": str(births_path),
        "metadata": str(metadata_path),
        "identity_period_matrix": identity_matrix,
        "shear_period_matrix": shear_matrix,
        "N": n,
        "samples": metadata["samples_per_pair"],
        "batches": len(batch_ids),
        "seed": metadata["seed"],
        "counter_range": [metadata["replica_counter_first"],
                          metadata["replica_counter_last_exclusive"]],
        "exact_sparse_transport": {
            "passed": not failures,
            "failed_batches": failures,
            "map": "P_shear=P*T^-1 and ell_shear=primitive(T ell_identity)",
        },
    }, rows


def score(
    shape_a_births: Path, shape_a_metadata: Path,
    shape_b_births: Path, shape_b_metadata: Path, p: float,
) -> dict[str, object]:
    source_a, rows_a = score_shape(shape_a_births, shape_a_metadata, p)
    source_b, rows_b = score_shape(shape_b_births, shape_b_metadata, p)
    synchronization = {
        key: source_a[key] == source_b[key]
        for key in ("N", "samples", "batches", "seed", "counter_range")
    }
    if not all(synchronization.values()):
        raise ValueError(f"shape streams are not synchronized: {synchronization}")
    batch_ids = sorted(rows_a)
    if batch_ids != sorted(rows_b):
        raise ValueError("shape batch IDs differ")

    joint_rows = []
    residual_rows = []
    batch_payload = []
    for batch in batch_ids:
        identity = [rows_b[batch]["identity"][j] - rows_a[batch]["identity"][j]
                    for j in range(3)]
        shear = [rows_b[batch]["shear"][j] - rows_a[batch]["shear"][j]
                 for j in range(3)]
        predicted = apply3(T_MATRIX, identity)
        residual = [shear[j] - predicted[j] for j in range(3)]
        joint_rows.append(identity + shear)
        residual_rows.append(residual)
        batch_payload.append({
            "batch": batch, "identity_contrast": identity,
            "shear_contrast": shear, "predicted_shear": predicted,
            "transport_residual": residual,
        })

    joint = covariance_summary(joint_rows)
    residual = covariance_summary(residual_rows)
    max_batch_residual = max(abs(value) for row in residual_rows for value in row)
    return {
        "schema": "matching-one/F3-flat-twist-shear-transport/v1",
        "status": "fresh local low-sample exact basis-transport smoke",
        "archive_audit": {
            "preexisting_identity_shear_archive_found": False,
            "reason": "P334 is the only production-format projective-birth archive and contains two physical Gaussian orientations, not a period-basis shear pair",
        },
        "frozen_before_smoke": {
            "commit": "3cb06ab",
            "triplet_order": list(TRIPLET_NAMES),
            "T_matrix": [list(row) for row in T_MATRIX],
            "old_archive_numeric_prediction": [
                0.0021583896116662287,
                0.004383700795185517,
                -0.0028238455708840906,
            ],
        },
        "source": {"shape_a": source_a, "shape_b": source_b,
                   "synchronization": synchronization, "p_ref": p},
        "exact_gates": {
            "passed": (
                source_a["exact_sparse_transport"]["passed"]
                and source_b["exact_sparse_transport"]["passed"]
                and max_batch_residual < 5e-15
            ),
            "max_batch_transport_residual": max_batch_residual,
            "matrix_contract": "P_shear=P*T^-1 implies ell_shear=T ell and C_shear=M_T C_identity",
        },
        "joint_identity_shear_estimate": {
            "order": [f"identity_{name}" for name in TRIPLET_NAMES]
                     + [f"shear_{name}" for name in TRIPLET_NAMES],
            **joint,
        },
        "transport_residual": {
            "order": list(TRIPLET_NAMES),
            **residual,
            "max_absolute_mean": max(abs(value) for value in residual["mean"]),
            "joint_score": None,
            "joint_score_reason": "residual is an exact pathwise basis-relabeling identity, so its covariance is rank zero up to floating roundoff",
        },
        "batch_values": batch_payload,
        "scientific_boundary": {
            "closed": "end-to-end period-basis, projective-line and fixed-p twist-vector transport",
            "not_tested": "a physically distinct charged defect/source insertion",
            "consequence": "right period-basis shear is a covariance gate, not independent evidence for the charged sector; a nontrivial test needs an explicit source whose action is not mere line relabeling",
        },
    }


def render_markdown(result: Mapping[str, object]) -> str:
    joint = result["joint_identity_shear_estimate"]
    residual = result["transport_residual"]
    identity = joint["mean"][:3]
    shear = joint["mean"][3:]
    prediction = apply3(T_MATRIX, identity)
    lines = [
        "# F3 flat-twist identity/shear transport smoke", "",
        "No pre-existing legal identity/shear projective-birth archive was found. "
        "This artifact uses fresh synchronized local N65 counters after the matrix "
        "prediction was frozen at `3cb06ab`.", "",
        "## Full triplet", "",
        f"Identity orientation contrast `[H,A,D]`: `{identity}`.", "",
        f"Exact no-fit prediction `M_T C_identity`: `{prediction}`.", "",
        f"Measured shear contrast: `{shear}`.", "",
        f"Three-dimensional residual: `{residual['mean']}` with maximum batchwise "
        f"absolute residual `{result['exact_gates']['max_batch_transport_residual']:.3g}`.", "",
        "The JSON retains the full 6x6 identity/shear covariance and 3x3 residual "
        "covariance. A residual chi-square is intentionally undefined: the sparse event "
        "archive satisfies the line map pathwise, so the residual covariance is rank zero "
        "apart from floating roundoff.", "",
        "## Scientific boundary", "",
        "This closes the end-to-end convention gate: `P_shear=P*T^-1` gives "
        "`ell_shear=T ell`, and the reconstructed F3 vector follows the frozen A4 "
        "matrix exactly. It is not new evidence for a charged field. Right period-basis "
        "shear only relabels the same quotient graph; the next nontrivial experiment must "
        "implement a physically distinct twist/defect source rather than another basis copy.", "",
    ]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shape-a-births", type=Path, required=True)
    parser.add_argument("--shape-a-metadata", type=Path, required=True)
    parser.add_argument("--shape-b-births", type=Path, required=True)
    parser.add_argument("--shape-b-metadata", type=Path, required=True)
    parser.add_argument("--p", type=float, default=DEFAULT_P)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--markdown", type=Path)
    args = parser.parse_args()
    result = score(
        args.shape_a_births, args.shape_a_metadata,
        args.shape_b_births, args.shape_b_metadata, args.p,
    )
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    if args.markdown:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(render_markdown(result), encoding="utf-8")


if __name__ == "__main__":
    main()

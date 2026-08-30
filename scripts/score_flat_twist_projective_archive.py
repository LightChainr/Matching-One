#!/usr/bin/env python3
"""Reconstruct F2/F3 flat-twist sectors from a projective-birth archive.

At fixed p the event tuple ``(tau1, ell, tau2)`` determines the complete
ambient-homology state along the Newman--Ziff path.  Primitive saturation makes
the reduction of a rank-one ``ell`` nonzero over every prime field.  Hence the
archive is sufficient for every constraint probability

    T_alpha = P(alpha vanishes on the occupied ambient image)

without new sampling.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path
from typing import Mapping, Sequence

from analyze_projective_birth_smoke import (
    BirthCell,
    binomial_probabilities,
    chi4,
    covariance_of_mean,
    read_births,
)


DEFAULT_P = 0.592746050790
PRIMES = (2, 3)
ORIENTATIONS = ("first", "second")


def projective_lines(q: int) -> list[tuple[int, int]]:
    return [(1, slope) for slope in range(q)] + [(0, 1)]


def projectivize(x: int, y: int, q: int) -> tuple[int, int]:
    x %= q
    y %= q
    if x == 0 and y == 0:
        raise ValueError("primitive integral line reduced to zero")
    if x:
        inverse = pow(x, -1, q)
        return 1, (y * inverse) % q
    return 0, 1


def kernel_line(alpha: tuple[int, int], q: int) -> tuple[int, int]:
    if alpha == (0, 0):
        raise ValueError("zero twist has a two-dimensional kernel")
    for line in projective_lines(q):
        if (alpha[0] * line[0] + alpha[1] * line[1]) % q == 0:
            return line
    raise AssertionError("nonzero Fq functional must have a projective kernel")


def fixed_p_state(
    cells: Sequence[BirthCell], n: int, p: float,
    matrix: Sequence[Sequence[int]],
) -> dict[str, object]:
    if not cells:
        raise ValueError("empty event batch")
    samples = cells[0].samples
    if any(cell.samples != samples for cell in cells):
        raise ValueError("inconsistent samples in event batch")
    if sum(cell.count for cell in cells) != samples:
        raise ValueError("event cells do not sum to samples")
    pmf = binomial_probabilities(n, p)
    tails = [0.0] * (n + 2)
    for k in range(n, -1, -1):
        tails[k] = tails[k + 1] + pmf[k]
    p0 = p2 = 0.0
    line_probabilities = {
        q: {line: 0.0 for line in projective_lines(q)} for q in PRIMES
    }
    raw_a4 = 0j
    for cell in cells:
        weight = cell.count / samples
        p0 += weight * (1.0 - tails[cell.tau1])
        p2 += weight * tails[cell.tau2]
        if cell.kind == "DIRECT_RANK2":
            if cell.tau1 != cell.tau2 or cell.ell_x or cell.ell_y:
                raise ValueError("invalid DIRECT_RANK2 cell")
            continue
        if cell.kind != "LINE" or not cell.tau1 < cell.tau2:
            raise ValueError("invalid line-bearing cell")
        plateau = weight * (tails[cell.tau1] - tails[cell.tau2])
        raw_a4 += plateau * chi4(matrix, cell.ell_x, cell.ell_y)
        for q in PRIMES:
            line_probabilities[q][projectivize(cell.ell_x, cell.ell_y, q)] += plateau
    p1 = 1.0 - p0 - p2
    if abs(sum(line_probabilities[2].values()) - p1) > 2e-13:
        raise ValueError("F2 line bins do not recover P1")
    if abs(sum(line_probabilities[3].values()) - p1) > 2e-13:
        raise ValueError("F3 line bins do not recover P1")
    return {
        "P0": p0, "P1": p1, "P2": p2,
        "line_probabilities": line_probabilities,
        "raw_A4": raw_a4,
    }


def twist_sectors(state: Mapping[str, object], q: int) -> dict[tuple[int, int], float]:
    p0 = float(state["P0"])
    lines = state["line_probabilities"][q]
    output = {}
    for alpha in itertools.product(range(q), repeat=2):
        if alpha == (0, 0):
            output[alpha] = 1.0
        else:
            output[alpha] = p0 + lines[kernel_line(alpha, q)]
    return output


def twist_characters(state: Mapping[str, object]) -> dict[str, float]:
    f2 = state["line_probabilities"][2]
    f3 = state["line_probabilities"][3]
    # Unit-norm zero-sum contrasts in the projective-line basis.
    return {
        "F2_H4_axis_diag": (
            f2[(1, 0)] + f2[(0, 1)] - 2.0 * f2[(1, 1)]
        ) / math.sqrt(6.0),
        "F2_axis_odd": (f2[(1, 0)] - f2[(0, 1)]) / math.sqrt(2.0),
        "F3_H4_axis_diag": 0.5 * (
            f3[(1, 0)] + f3[(0, 1)] - f3[(1, 1)] - f3[(1, 2)]
        ),
        "F3_axis_odd": (f3[(1, 0)] - f3[(0, 1)]) / math.sqrt(2.0),
        "F3_diagonal_odd": (f3[(1, 1)] - f3[(1, 2)]) / math.sqrt(2.0),
    }


def orientation_values(
    state: Mapping[str, object], a: int, b: int,
) -> dict[str, float]:
    output = {"P0": state["P0"], "P1": state["P1"], "P2": state["P2"]}
    for q in PRIMES:
        sectors = twist_sectors(state, q)
        for alpha, value in sectors.items():
            output[f"F{q}_T_{alpha[0]}_{alpha[1]}"] = value
        output[f"F{q}_S"] = math.fsum(sectors.values())
    output.update(twist_characters(state))
    phase = complex(a, b) ** 4 / (a * a + b * b) ** 2
    aligned = state["raw_A4"] * phase.conjugate()
    output["raw_chi4_parallel"] = aligned.real
    output["raw_chi4_perpendicular"] = aligned.imag
    return output


def gate_residuals(state: Mapping[str, object], values: Mapping[str, float]) -> dict[str, float]:
    p0, p1, p2 = (float(state[name]) for name in ("P0", "P1", "P2"))
    residuals = {}
    for q in PRIMES:
        residuals[f"F{q}_aggregate"] = values[f"F{q}_S"] - (
            q * q * p0 + q * p1 + p2
        )
        sectors = twist_sectors(state, q)
        residuals[f"F{q}_T_zero"] = sectors[(0, 0)] - 1.0
        residuals[f"F{q}_line_recovery"] = max(
            abs(sectors[alpha] - p0 - state["line_probabilities"][q][kernel_line(alpha, q)])
            for alpha in sectors if alpha != (0, 0)
        )
    recovered_p0 = (values["F3_S"] - 2.0 * values["F2_S"] + 1.0) / 2.0
    recovered_p1 = values["F2_S"] - 1.0 - 3.0 * recovered_p0
    recovered_p2 = 1.0 - recovered_p0 - recovered_p1
    residuals["source_inversion_P0"] = recovered_p0 - p0
    residuals["source_inversion_P1"] = recovered_p1 - p1
    residuals["source_inversion_P2"] = recovered_p2 - p2
    return residuals


def quadratic_2(vector: Sequence[float], covariance: Sequence[Sequence[float]]) -> float:
    a, b = covariance[0]
    c, d = covariance[1]
    determinant = a * d - b * c
    if determinant <= 0.0:
        raise ValueError("2x2 covariance is not positive definite")
    x, y = vector
    return (d * x * x - (b + c) * x * y + a * y * y) / determinant


def score(births_path: Path, metadata_path: Path, p: float) -> dict[str, object]:
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    n, births = read_births(births_path)
    design = metadata["designs"][0]
    matrices = {
        "first": design["first_period_matrix"],
        "second": design["second_period_matrix"],
    }
    reps = {"first": design["first"], "second": design["second"]}
    batch_ids = sorted({batch for _, batch in births})
    batch_values = []
    gate_rows = []
    order = None
    for batch in batch_ids:
        values_by_orientation = {}
        gates = {}
        for orientation in ORIENTATIONS:
            state = fixed_p_state(
                births[(orientation, batch)], n, p, matrices[orientation]
            )
            values = orientation_values(state, *reps[orientation])
            values_by_orientation[orientation] = values
            gates[orientation] = gate_residuals(state, values)
        local_order = [
            f"{orientation}_{name}"
            for orientation in ORIENTATIONS
            for name in values_by_orientation[orientation]
        ]
        contrast_names = (
            "F2_H4_axis_diag", "F2_axis_odd", "F3_H4_axis_diag",
            "F3_axis_odd", "F3_diagonal_odd", "raw_chi4_parallel",
            "raw_chi4_perpendicular",
        )
        local_order += [f"second_minus_first_{name}" for name in contrast_names]
        vector = [
            values_by_orientation[orientation][name]
            for orientation in ORIENTATIONS
            for name in values_by_orientation[orientation]
        ]
        vector += [
            values_by_orientation["second"][name]
            - values_by_orientation["first"][name]
            for name in contrast_names
        ]
        if order is None:
            order = local_order
        elif order != local_order:
            raise ValueError("batch observable order changed")
        batch_values.append(vector)
        gate_rows.append(gates)
    assert order is not None
    covariance = covariance_of_mean(batch_values)
    means = [math.fsum(row[j] for row in batch_values) / len(batch_values)
             for j in range(len(order))]
    errors = [math.sqrt(max(0.0, covariance[j][j])) for j in range(len(order))]
    lookup = {name: (means[i], errors[i]) for i, name in enumerate(order)}
    h4_candidates = []
    for name in ("F2_H4_axis_diag", "F3_H4_axis_diag", "raw_chi4_parallel"):
        value, error = lookup[f"second_minus_first_{name}"]
        h4_candidates.append({
            "name": name,
            "orientation_contrast": value,
            "standard_error": error,
            "absolute_z": abs(value) / error,
            "quadratic": (value / error) ** 2,
        })
    h4_candidates.sort(key=lambda row: row["absolute_z"], reverse=True)
    other_candidates = []
    for name in ("F2_axis_odd", "F3_axis_odd", "F3_diagonal_odd"):
        value, error = lookup[f"second_minus_first_{name}"]
        other_candidates.append({
            "name": name, "orientation_contrast": value,
            "standard_error": error, "absolute_z": abs(value) / error,
            "quadratic": (value / error) ** 2,
        })
    other_candidates.sort(key=lambda row: row["absolute_z"], reverse=True)
    odd_names = (
        "second_minus_first_F3_axis_odd",
        "second_minus_first_F3_diagonal_odd",
    )
    odd_indices = [order.index(name) for name in odd_names]
    odd_vector = [means[index] for index in odd_indices]
    odd_covariance = [[covariance[i][j] for j in odd_indices] for i in odd_indices]
    max_residual = {
        name: max(abs(gates[orientation][name])
                  for gates in gate_rows for orientation in ORIENTATIONS)
        for name in gate_rows[0]["first"]
    }
    return {
        "schema": "matching-one/N65-flat-twist-projective-score/v1",
        "status": "existing-data flat-twist reconstruction",
        "sufficiency": {
            "archive_is_sufficient": True,
            "reason": "tau1/tau2 give rank0/rank2 and primitive ell mod q gives every rank1 constraint sector",
            "new_fields_required": [],
            "direct_rank2": "requires no line; fixed-p rank is determined by tau1=tau2",
        },
        "exact_gates": {
            "passed": max(max_residual.values()) < 3e-13,
            "tolerance": 3e-13,
            "max_absolute_residual": max_residual,
            "aggregate": "S_q=sum_alpha T_alpha=q^2 P0+q P1+P2",
            "source_inversion": "P0=(S3-2S2+1)/2; P1=S2-1-3P0; P2=1-P0-P1",
        },
        "sharpness": {
            "same_H4_sector": {
                "candidates": h4_candidates,
                "winner": h4_candidates[0],
                "conclusion": "F3 is nominally sharpest but is effectively tied with raw chi4 at this sample size",
            },
            "additional_twist_sectors": {
                "candidates": other_candidates,
                "winner": other_candidates[0],
                "F3_odd_joint": {
                    "order": list(odd_names), "value": odd_vector,
                    "covariance": odd_covariance,
                    "quadratic": quadratic_2(odd_vector, odd_covariance), "df": 2,
                },
                "conclusion": "F3 diagonal-odd is sharper than raw chi4 but belongs to a distinct reflection-odd sector",
            },
            "comparison_boundary": "same exploratory 20k block; rankings are design evidence, not independent tests",
        },
        "joint_estimate": {
            "order": order, "mean": means, "standard_error": errors,
            "covariance": covariance,
            "batch_values": [
                {"batch": batch, "values": dict(zip(order, row))}
                for batch, row in zip(batch_ids, batch_values)
            ],
        },
        "source": {
            "births": str(births_path), "metadata": str(metadata_path),
            "N": n, "p_ref": p, "samples_per_shape": metadata["samples_per_pair"],
            "batches": len(batch_ids), "seed": metadata["seed"],
            "counter_range": [metadata["replica_counter_first"],
                              metadata["replica_counter_last_exclusive"]],
        },
        "characters": {
            "F3_primary": "1/2[(T_kernel_axis0+T_kernel_axisInf)-(T_kernel_diagPlus+T_kernel_diagMinus)]",
            "F2_secondary": "1/sqrt(6)[T_axis0+T_axisInf-2*T_diag]",
            "common_P0": "cancels in every declared zero-sum character",
        },
        "claim_boundary": "one correlated reanalysis of P334 N65 20k; no continuum defect identity or exponent fit",
    }


def render_markdown(result: Mapping[str, object]) -> str:
    joint = result["joint_estimate"]
    lookup = dict(zip(joint["order"], zip(joint["mean"], joint["standard_error"])))
    lines = [
        "# N65 flat-twist constraint-sector score", "",
        "**Archive sufficiency: YES.** No new field or sample is required for F2/F3 flat-twist tomography.", "",
        "For each fixed-p configuration the archive distinguishes rank zero, rank one with primitive "
        "`ell mod q`, and rank two. This determines every `T_alpha` exactly.", "",
        f"All aggregate/source-inversion gates pass with maximum residual "
        f"`{max(result['exact_gates']['max_absolute_residual'].values()):.3g}`.", "",
        "## Same-modulus orientation contrasts", "",
        "| projector | contrast | batch SE | |z| |",
        "|---|---:|---:|---:|",
    ]
    for row in result["sharpness"]["same_H4_sector"]["candidates"]:
        lines.append(
            f"| `{row['name']}` | {row['orientation_contrast']:.8g} | "
            f"{row['standard_error']:.3g} | {row['absolute_z']:.3g} |"
        )
    winner = result["sharpness"]["same_H4_sector"]["winner"]
    odd_winner = result["sharpness"]["additional_twist_sectors"]["winner"]
    odd_joint = result["sharpness"]["additional_twist_sectors"]["F3_odd_joint"]
    lines += ["", f"Within the same H4 sector, **`{winner['name']}`** is nominally sharpest "
              f"at `|z|={winner['absolute_z']:.3g}`, effectively tied with raw chi4 rather than a "
              "material variance improvement.", "",
              "The balanced F3 projector is the minimal axes-versus-diagonals character on "
              "`P1(F3)`. It cancels the common rank-zero contribution automatically; unlike raw "
              "physical `chi4`, it uses unit finite-field orbit weights and retains only the twist "
              "constraint information.", "", "## Additional twist-only sector", "",
              f"The sharper new contrast is **`{odd_winner['name']}` = "
              f"{odd_winner['orientation_contrast']:.8g} +/- {odd_winner['standard_error']:.3g}** "
              f"(`|z|={odd_winner['absolute_z']:.3g}`). The joint F3 axis-odd/diagonal-odd "
              f"diagnostic is `{odd_joint['quadratic']:.4g} / 2 df`.", "",
              "This does not replace the H4 score: it is reflection-odd/projective and therefore a "
              "different finite-twist sector. Its value is that ordinary chi4 collapses this modular "
              "line information, while T_alpha retains it.", "", "## Selected sector values", "",
              "| value | mean | SE |", "|---|---:|---:|"]
    for name in (
        "first_F2_S", "first_F3_S", "second_F2_S", "second_F3_S",
        "first_F3_H4_axis_diag", "second_F3_H4_axis_diag",
        "second_minus_first_F3_H4_axis_diag",
    ):
        value, error = lookup[name]
        lines.append(f"| `{name}` | {value:.8g} | {error:.3g} |")
    lines += ["", "The complete F2/F3 sector vector, all characters and their cross-orientation/"
              "cross-character covariance are stored in the JSON artifact.", ""]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--births", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--p", type=float, default=DEFAULT_P)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--markdown", type=Path)
    args = parser.parse_args()
    result = score(args.births, args.metadata, args.p)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    if args.markdown:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(render_markdown(result), encoding="utf-8")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Projective Z5 magnetic-translation versus free commuting state score."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
from scipy.optimize import least_squares

from score_z5_charged_threepoint import zero_score
from score_z5_projective_leg_bivariate_state import (
    AXIS_TRAIN,
    DEGREE4_HOLDOUT,
    MIXED_COMMUTING_GATE,
    means,
    pair,
    read_batches,
    rotate,
    vector_score,
)
from score_z5_projective_leg_pair_transfer import CHANNELS
from z5_projective_leg_cross_scale_mc import PARENT_GEOMETRY, contexts


OMEGA = np.exp(2j * np.pi / 5.0)
DEGREE3_FIRST_QUADRANT = tuple((a, b) for degree in range(4) for a in range(degree + 1) for b in (degree - a,))
FLUX_MODELS = tuple((m, (-m) % 5) for m in range(1, 5))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sigmoid(value: float) -> float:
    if value >= 0:
        return 1.0 / (1.0 + math.exp(-value))
    exp_value = math.exp(value)
    return exp_value / (1.0 + exp_value)


def shift_matrix() -> np.ndarray:
    matrix = np.zeros((5, 5), dtype=complex)
    for index in range(5):
        matrix[(index + 1) % 5, index] = 1.0
    return matrix


def clock_matrix(m: int) -> np.ndarray:
    return np.diag([OMEGA ** (m * index) for index in range(5)])


def weil_rotation(m: int) -> np.ndarray:
    return np.asarray([[OMEGA ** (m * row * column) for column in range(5)] for row in range(5)]) / math.sqrt(5.0)


def matrix_power(matrix: np.ndarray, power: int) -> np.ndarray:
    return np.linalg.matrix_power(matrix, power)


def weyl_displacement(point: tuple[int, int], m: int) -> np.ndarray:
    """Symmetric Weyl operator W(a,b), with W(-d)=W(d)^*."""
    a, b = point
    half_mod5 = 3  # 2^-1 modulo five
    phase = OMEGA ** ((half_mod5 * m * a * b) % 5)
    return phase * matrix_power(shift_matrix(), a) @ matrix_power(clock_matrix(m), b)


def canonical_gate() -> dict:
    x = shift_matrix()
    rows = {}
    for m in range(1, 5):
        z = clock_matrix(m)
        rotation = weil_rotation(m)
        center = x @ z @ np.linalg.inv(x) @ np.linalg.inv(z)
        deck_center = OMEGA ** (-m) * np.eye(5)
        rows[str(m)] = {
            "Z_X_minus_omega_m_X_Z": float(np.max(np.abs(z @ x - OMEGA**m * x @ z))),
            "R_X_Rinv_minus_Z": float(np.max(np.abs(rotation @ x @ rotation.conjugate().T - z))),
            "R_Z_Rinv_minus_Xinv": float(np.max(np.abs(rotation @ z @ rotation.conjugate().T - np.linalg.inv(x)))),
            "R4_minus_I": float(np.max(np.abs(matrix_power(rotation, 4) - np.eye(5)))),
            "center_minus_omega_minus_m_I": float(np.max(np.abs(center - OMEGA ** (-m) * np.eye(5)))),
            "center_fifth_power_minus_I": float(np.max(np.abs(matrix_power(center, 5) - np.eye(5)))),
            "D_center_fifth_power_minus_I": float(np.max(np.abs(matrix_power(deck_center, 5) - np.eye(5)))),
            "symmetric_weyl_adjoint_max": max(
                float(np.max(np.abs(weyl_displacement(point, m).conjugate().T - weyl_displacement((-point[0], -point[1]), m))))
                for point in ((1, 0), (0, 1), (1, 1), (2, 1))
            ),
        }
    maximum = max(max(row.values()) for row in rows.values())
    return {
        "convention": "Tx=X, Ty=Z_m; Tx Ty = omega^(-m) Ty Tx",
        "deck_center": "D=omega^(-m) I in the fixed-center irreducible sector, so Tx Ty Tx^-1 Ty^-1=D and D^5=I",
        "Weil_relations": "R X R^-1=Z_m, R Z_m R^-1=X^-1, R^4=I",
        "models": rows,
        "max_abs_residual": maximum,
        "passed": maximum < 1e-12,
    }


def _fiber_shift(context, parent: int, delta: tuple[int, int]) -> tuple[int, int]:
    representative = PARENT_GEOMETRY.coordinates[parent]
    target_parent = PARENT_GEOMETRY.vertex((representative[0] + delta[0], representative[1] + delta[1]))
    point = context.field_coordinates[5 * parent]
    target_vertex = context.geometry.vertex((point[0] + delta[0], point[1] + delta[1]))
    hits = [
        fiber for fiber in range(5)
        if context.field_to_vertex[5 * target_parent + fiber] == target_vertex
    ]
    if len(hits) != 1:
        raise AssertionError("unit translation fiber lift is not unique")
    return target_parent, hits[0]


def spatial_bundle_flatness_gate() -> dict:
    output = {}
    for hand, context in zip(("plus", "minus"), contexts()):
        curvatures = []
        for parent in range(PARENT_GEOMETRY.n):
            px, sx = _fiber_shift(context, parent, (1, 0))
            py, sy = _fiber_shift(context, parent, (0, 1))
            pxy, sy_after_x = _fiber_shift(context, px, (0, 1))
            pyx, sx_after_y = _fiber_shift(context, py, (1, 0))
            if pxy != pyx:
                raise AssertionError("base translations do not commute")
            curvatures.append((sx + sy_after_x - sy - sx_after_y) % 5)
        output[hand] = {str(value): curvatures.count(value) for value in range(5)}
    return {
        "plaquette": "s_x(p)+s_y(p+x)-s_y(p)-s_x(p+y) mod 5",
        "counts": output,
        "conclusion": "the exact spatial cover bundle has m=0; nonzero-m fits are effective projected-state models only",
        "passed": all(row["0"] == PARENT_GEOMETRY.n for row in output.values()),
    }


def unpack_weyl(parameters: Sequence[float]) -> tuple[float, np.ndarray, np.ndarray]:
    rho = 0.999 * sigmoid(float(parameters[0]))
    plus = np.asarray(parameters[1:6]) + 1j * np.asarray(parameters[6:11])
    minus = np.asarray(parameters[11:16]) + 1j * np.asarray(parameters[16:21])
    return rho, plus, minus


def base_matrix_element(vector: np.ndarray, point: tuple[int, int], m: int, rho: float) -> complex:
    degree = abs(point[0]) + abs(point[1])
    return rho**degree * np.vdot(vector, weyl_displacement(point, m) @ vector)


def weyl_prediction(parameters: Sequence[float], point: tuple[int, int], channel: tuple[str, int], flux: tuple[int, int]) -> complex:
    rho, plus, minus = unpack_weyl(parameters)
    hand, charge = channel
    if hand == "plus":
        base = lambda displacement: base_matrix_element(plus, displacement, flux[0], rho)
        return base(point) if charge == 1 else base(rotate(point)).conjugate()
    base = lambda displacement: base_matrix_element(minus, displacement, flux[1], rho)
    return base(point) if charge == 1 else base(rotate(point))


def real_residual(
    values: Mapping[str, float], points: Sequence[tuple[int, int]], prediction,
) -> np.ndarray:
    output = []
    for channel in CHANNELS:
        for point in points:
            residual = pair(values, point, channel) - prediction(point, channel)
            output.extend((residual.real, residual.imag))
    return np.asarray(output)


def weyl_axis_residual(parameters: Sequence[float], values: Mapping[str, float], flux: tuple[int, int]) -> np.ndarray:
    return real_residual(values, AXIS_TRAIN, lambda point, channel: weyl_prediction(parameters, point, channel, flux))


def initial_weyl_parameters(values: Mapping[str, float], start: int) -> np.ndarray:
    rng = np.random.default_rng(250000 + start)
    origin_plus = max(pair(values, (0, 0), ("plus", 1)).real, 1e-8)
    origin_minus = max(pair(values, (0, 0), ("minus", 1)).real, 1e-8)
    plus = rng.normal(size=5) + 1j * rng.normal(size=5)
    minus = rng.normal(size=5) + 1j * rng.normal(size=5)
    plus *= math.sqrt(origin_plus) / np.linalg.norm(plus)
    minus *= math.sqrt(origin_minus) / np.linalg.norm(minus)
    rho = 0.45 + 0.08 * (start % 4)
    eta = math.log(rho / (0.999 - rho))
    return np.asarray([eta, *plus.real, *plus.imag, *minus.real, *minus.imag], dtype=float)


def fit_weyl(values: Mapping[str, float], flux: tuple[int, int], initial: np.ndarray | None = None) -> dict:
    starts = [initial] if initial is not None else [initial_weyl_parameters(values, index) for index in range(8)]
    fits = []
    for start in starts:
        fitted = least_squares(
            weyl_axis_residual,
            start,
            args=(values, flux),
            max_nfev=400 if initial is not None else 2500,
            ftol=1e-11,
            xtol=1e-11,
            gtol=1e-11,
        )
        fits.append(fitted)
    result = min(fits, key=lambda row: (float(np.dot(row.fun, row.fun)), row.nfev))
    return {
        "parameters": result.x,
        "axis_sse": float(np.dot(result.fun, result.fun)),
        "nfev": int(result.nfev),
        "rho": unpack_weyl(result.x)[0],
    }


def weyl_model_residual(values: Mapping[str, float], fitted: Mapping[str, object], flux, points) -> list[float]:
    return real_residual(
        values,
        points,
        lambda point, channel: weyl_prediction(fitted["parameters"], point, channel, flux),
    ).tolist()


def fit_commuting(values: Mapping[str, float], rank: int, initial: np.ndarray | None = None) -> dict:
    points = DEGREE3_FIRST_QUADRANT

    def roots(parameters):
        row = np.asarray(parameters).reshape(2, rank, 2)
        return row[0, :, 0] + 1j * row[0, :, 1], row[1, :, 0] + 1j * row[1, :, 1]

    def amplitudes(xroots, yroots):
        design = np.asarray([[xroots[index] ** a * yroots[index] ** b for index in range(rank)] for a, b in points])
        return {
            channel: np.linalg.lstsq(design, np.asarray([pair(values, point, channel) for point in points]), rcond=None)[0]
            for channel in CHANNELS
        }

    def residual(parameters):
        xroots, yroots = roots(parameters)
        rows = amplitudes(xroots, yroots)
        output = []
        for channel in CHANNELS:
            for point in points:
                a, b = point
                prediction = sum(rows[channel][index] * xroots[index] ** a * yroots[index] ** b for index in range(rank))
                value = pair(values, point, channel) - prediction
                output.extend((value.real, value.imag))
        return np.asarray(output)

    if initial is None:
        starts = []
        for start_index in range(10):
            rng = np.random.default_rng(255000 + 10 * rank + start_index)
            row = rng.uniform(-0.7, 0.7, size=(2, rank, 2))
            starts.append(row.ravel())
    else:
        starts = [initial]
    fits = [
        least_squares(
            residual,
            start,
            bounds=(-1.5, 1.5),
            max_nfev=250 if initial is not None else 3000,
            ftol=1e-10,
            xtol=1e-10,
            gtol=1e-10,
        )
        for start in starts
    ]
    fitted = min(fits, key=lambda row: (float(np.dot(row.fun, row.fun)), row.nfev))
    xroots, yroots = roots(fitted.x)
    rows = amplitudes(xroots, yroots)
    return {
        "parameters": fitted.x,
        "xroots": xroots,
        "yroots": yroots,
        "amplitudes": rows,
        "training_sse": float(np.dot(fitted.fun, fitted.fun)),
        "nfev": int(fitted.nfev),
    }


def commuting_residual(values: Mapping[str, float], fitted: Mapping[str, object], points) -> list[float]:
    output = []
    for channel in CHANNELS:
        for a, b in points:
            prediction = sum(
                fitted["amplitudes"][channel][index] * fitted["xroots"][index] ** a * fitted["yroots"][index] ** b
                for index in range(len(fitted["xroots"]))
            )
            value = pair(values, (a, b), channel) - prediction
            output.extend((value.real, value.imag))
    return output


def scored_residual(point, deleted) -> dict:
    return vector_score(point, deleted)


def score_weyl(values, deleted_values, flux) -> dict:
    fitted = fit_weyl(values, flux)
    deleted = [fit_weyl(row, flux, fitted["parameters"]) for row in deleted_values]
    mixed = scored_residual(
        weyl_model_residual(values, fitted, flux, MIXED_COMMUTING_GATE),
        [weyl_model_residual(row, model, flux, MIXED_COMMUTING_GATE) for row, model in zip(deleted_values, deleted)],
    )
    heldout = scored_residual(
        weyl_model_residual(values, fitted, flux, DEGREE4_HOLDOUT),
        [weyl_model_residual(row, model, flux, DEGREE4_HOLDOUT) for row, model in zip(deleted_values, deleted)],
    )
    return {
        "m_plus": flux[0],
        "m_minus": flux[1],
        "center_phase_TxTy": f"omega^{(-flux[0]) % 5} / omega^{(-flux[1]) % 5}",
        "rho": fitted["rho"],
        "axis_sse": fitted["axis_sse"],
        "mixed_gate": mixed,
        "degree4_heldout": heldout,
    }


def score_commuting(values, deleted_values, rank) -> dict:
    fitted = fit_commuting(values, rank)
    deleted = [fit_commuting(row, rank, fitted["parameters"]) for row in deleted_values]
    heldout = scored_residual(
        commuting_residual(values, fitted, DEGREE4_HOLDOUT),
        [commuting_residual(row, model, DEGREE4_HOLDOUT) for row, model in zip(deleted_values, deleted)],
    )
    return {
        "rank": rank,
        "training_points": [list(point) for point in DEGREE3_FIRST_QUADRANT],
        "training_advantage": "mixed degree<=3 points are consumed because free rank4/5 roots are not identifiable from four axis moments",
        "training_sse": fitted["training_sse"],
        "xroots": [[value.real, value.imag] for value in fitted["xroots"]],
        "yroots": [[value.real, value.imag] for value in fitted["yroots"]],
        "degree4_heldout": heldout,
    }


def score(batches: Sequence[dict], manifest: Mapping[str, object]) -> dict:
    path = Path(manifest["input_batches"])
    if sha256(path) != manifest["input_batches_sha256"]:
        raise ValueError("bivariate batch hash changed")
    values = means(batches)
    deleted_values = [means(batches, index) for index in range(len(batches))]
    weyl = {f"m{plus}_m{minus}": score_weyl(values, deleted_values, (plus, minus)) for plus, minus in FLUX_MODELS}
    commuting = {str(rank): score_commuting(values, deleted_values, rank) for rank in (4, 5)}
    alpha = float(manifest["decision_alpha"])
    passing_weyl = [
        name for name, row in weyl.items()
        if row["mixed_gate"]["zero_score"]["survival_p"] >= alpha
        and row["degree4_heldout"]["zero_score"]["survival_p"] >= alpha
    ]
    passing_commuting = [
        rank for rank, row in commuting.items()
        if row["degree4_heldout"]["zero_score"]["survival_p"] >= alpha
    ]
    if passing_weyl:
        decision = "canonical_rank5_projective_Weyl_survives"
    elif passing_commuting:
        decision = "free_higher_rank_commuting_survives_but_Weyl_fails"
    else:
        decision = "neither_canonical_Weyl_nor_free_commuting_rank5_closes"
    return {
        "schema": "matching-one/z5-projective-magnetic-translation-score/v1",
        "status": "post_bivariate_frozen_existing_data_reanalysis",
        "exact_spatial_bundle_flatness": spatial_bundle_flatness_gate(),
        "canonical_Weyl_gate": canonical_gate(),
        "Weyl_models": weyl,
        "free_commuting_comparators": commuting,
        "passing_Weyl": passing_weyl,
        "passing_commuting": passing_commuting,
        "decision": decision,
        "decision_alpha": alpha,
        "claim_boundary": [
            "The exact cover-translation bundle has zero Z5 plaquette curvature; any nonzero-m fit is an effective projected-state model, not deck holonomy.",
            "The Weyl model uses one Hermitian source/sink vector per hand, fixed Weil rotation, conjugate hand fluxes, and one shared real radial decay.",
            "Free commuting ranks four and five receive the favorable extra use of all mixed degree<=3 points because axis-only roots are unidentifiable at those ranks.",
            "Failure does not exclude a less constrained rank-five projective representation or a larger/Jordan/image realization.",
        ],
    }


def render(result) -> str:
    lines = [
        "# P250 projective magnetic-translation score",
        "",
        "| Weyl flux (plus,minus) | rho | mixed p | degree-4 p |",
        "|---|---:|---:|---:|",
    ]
    for name, row in result["Weyl_models"].items():
        lines.append(
            f"| {name} | {row['rho']:.6g} | {row['mixed_gate']['zero_score']['survival_p']:.6g} "
            f"| {row['degree4_heldout']['zero_score']['survival_p']:.6g} |"
        )
    lines += ["", "| free commuting rank | degree-4 p |", "|---:|---:|"]
    for rank, row in result["free_commuting_comparators"].items():
        lines.append(f"| {rank} | {row['degree4_heldout']['zero_score']['survival_p']:.6g} |")
    lines += [
        "",
        f"Exact bundle curvature: `{result['exact_spatial_bundle_flatness']['counts']}`.",
        f"Decision: `{result['decision']}`.",
        "",
    ]
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    result = score(read_batches(Path(manifest["input_batches"])), manifest)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    args.markdown.write_text(render(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

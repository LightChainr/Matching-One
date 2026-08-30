#!/usr/bin/env python3
"""Cross-fitted existing-data audit of P250 pre-outcome state information.

The archived P250 N505 production retained a sampled nonzero parent residue,
the ordered response Rminus, and the four support descriptions D0, J0,
J_after_D, and D_after_J.  It did not retain a terminal full-field spectrum,
so Delta Pminus cannot be reconstructed.  This scorer therefore compares:

* M_spec: a small, training-only selected residue Fourier model, fit separately
  in each child hand; and
* M_state: the identical M_spec design plus an antisymmetric pre-outcome
  support-increment block Sminus.

Outer folds are whole production batches.  The same replica and its two child
views consequently remain in one fold.  No L or R value enters Sminus.
"""

from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
from scipy.stats import t as student_t

from p250_adaptive_spatial_joint_mc import PARENT_GEOMETRY, parent_residue


SCHEMA = "matching-one/p250-adaptive-spatial-state-reanalysis/v1"
HANDS = ("plus", "minus")
GROUP_ORDER = 101
BASIS_CATEGORIES = (
    "empty",
    "axis_x",
    "axis_y",
    "diagonal_positive",
    "diagonal_negative",
    "rank_two",
)
REQUIRED_COLUMNS = (
    "batch",
    "replica",
    "hand",
    "residue",
    "translation",
    "fiber",
    "orientation",
    "R_minus",
    "L_D",
    "L_J",
    "L_DJ",
    "L_JD",
    *(
        f"{support}_{field}"
        for support in ("D0", "J0", "J_after_D", "D_after_J")
        for field in (
            "site_field",
            "source_id",
            "target_id",
            "source_rank",
            "target_rank",
            "source_basis",
            "target_basis",
        )
    ),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_manifest(path: Path) -> dict:
    manifest = json.loads(path.read_text(encoding="utf-8"))
    if manifest.get("schema") != (
        "matching-one/p250-adaptive-spatial-state-reanalysis-freeze/v1"
    ):
        raise ValueError("wrong or missing reanalysis freeze schema")
    return manifest


def load_rows(path: Path, manifest: dict) -> list[dict]:
    actual_hash = sha256(path)
    expected_hash = manifest["input"]["defined_sha256"]
    if actual_hash != expected_hash:
        raise ValueError(f"defined.csv hash changed: {actual_hash} != {expected_hash}")
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        missing = sorted(set(REQUIRED_COLUMNS) - set(reader.fieldnames or ()))
        if missing:
            raise ValueError(f"defined.csv is missing frozen columns: {missing}")
        rows = list(reader)
    if len(rows) != manifest["input"]["defined_rows"]:
        raise ValueError("defined row count changed")
    for row in rows:
        for field in (
            "batch",
            "replica",
            "residue",
            "translation",
            "fiber",
            "orientation",
            "R_minus",
            "L_D",
            "L_J",
            "L_DJ",
            "L_JD",
        ):
            row[field] = int(row[field])
        for support in ("D0", "J0", "J_after_D", "D_after_J"):
            for field in (
                "site_field",
                "source_id",
                "target_id",
                "source_rank",
                "target_rank",
            ):
                row[f"{support}_{field}"] = int(row[f"{support}_{field}"])
    batches = {row["batch"] for row in rows}
    if batches != set(range(manifest["cross_fitting"]["batches"])):
        raise ValueError("batch ids are not the frozen consecutive 0..399 set")
    replica_batches: dict[int, set[int]] = {}
    for row in rows:
        replica_batches.setdefault(row["replica"], set()).add(row["batch"])
    if any(len(value) != 1 for value in replica_batches.values()):
        raise ValueError("a replica crosses production batches")
    if {row["hand"] for row in rows} != set(HANDS):
        raise ValueError("child hands changed")
    if any(not 1 <= row["residue"] < GROUP_ORDER for row in rows):
        raise ValueError("the archived nonzero-residue contract changed")
    return rows


def parse_basis(value: str) -> tuple[tuple[int, int], ...]:
    parsed = ast.literal_eval(value)
    if not isinstance(parsed, tuple):
        raise ValueError(f"basis is not a tuple: {value}")
    output = tuple(tuple(int(item) for item in vector) for vector in parsed)
    if any(len(vector) != 2 for vector in output):
        raise ValueError(f"basis vector is not two-dimensional: {value}")
    return output


def basis_category(value: str) -> str:
    basis = parse_basis(value)
    if not basis:
        return "empty"
    if len(basis) >= 2:
        return "rank_two"
    x, y = basis[0]
    if y == 0:
        return "axis_x"
    if x == 0:
        return "axis_y"
    if x * y > 0:
        return "diagonal_positive"
    if x * y < 0:
        return "diagonal_negative"
    return "other_rank_one"


PARENT_INDEX_TO_RESIDUE = np.asarray(
    [
        parent_residue(tuple(PARENT_GEOMETRY.coordinates[index]))
        for index in range(GROUP_ORDER)
    ],
    dtype=np.int64,
)
if set(PARENT_INDEX_TO_RESIDUE.tolist()) != set(range(GROUP_ORDER)):
    raise AssertionError("parent CRT index map is no longer bijective")


def phase_features(row: dict, support: str) -> tuple[float, float, float, float]:
    site_field = row[f"{support}_site_field"]
    parent_index, site_fiber = divmod(site_field, 5)
    relative_residue = int(
        (PARENT_INDEX_TO_RESIDUE[parent_index]
         - PARENT_INDEX_TO_RESIDUE[row["translation"]])
        % GROUP_ORDER
    )
    parent_angle = 2.0 * math.pi * relative_residue / GROUP_ORDER
    fiber_angle = 2.0 * math.pi * ((site_fiber - row["fiber"]) % 5) / 5.0
    return (
        math.cos(parent_angle),
        math.sin(parent_angle),
        math.cos(fiber_angle),
        math.sin(fiber_angle),
    )


STATE_FEATURE_NAMES = (
    "Sminus_source_rank_increment",
    "Sminus_target_rank_increment",
    "Sminus_site_parent_cos1_increment",
    "Sminus_site_parent_sin1_increment",
    "Sminus_site_fiber_cos1_increment",
    "Sminus_site_fiber_sin1_increment",
    "Sminus_site_changed",
    "Sminus_source_component_changed",
    "Sminus_target_component_changed",
    *(f"Sminus_source_basis_{category}" for category in BASIS_CATEGORIES),
    *(f"Sminus_target_basis_{category}" for category in BASIS_CATEGORIES),
)


def state_features(row: dict) -> np.ndarray:
    """Return [J-after-D - J0] - [D-after-J - D0] support increments."""

    j_after, j_base = "J_after_D", "J0"
    d_after, d_base = "D_after_J", "D0"

    def increment(after: str, base: str, field: str) -> float:
        return float(row[f"{after}_{field}"] - row[f"{base}_{field}"])

    values = [
        increment(j_after, j_base, "source_rank")
        - increment(d_after, d_base, "source_rank"),
        increment(j_after, j_base, "target_rank")
        - increment(d_after, d_base, "target_rank"),
    ]
    j_after_phase = phase_features(row, j_after)
    j_base_phase = phase_features(row, j_base)
    d_after_phase = phase_features(row, d_after)
    d_base_phase = phase_features(row, d_base)
    values.extend(
        (j_after_phase[index] - j_base_phase[index])
        - (d_after_phase[index] - d_base_phase[index])
        for index in range(4)
    )
    values.extend(
        [
            float(row[f"{j_after}_site_field"] != row[f"{j_base}_site_field"])
            - float(row[f"{d_after}_site_field"] != row[f"{d_base}_site_field"]),
            float(row[f"{j_after}_source_id"] != row[f"{j_base}_source_id"])
            - float(row[f"{d_after}_source_id"] != row[f"{d_base}_source_id"]),
            float(row[f"{j_after}_target_id"] != row[f"{j_base}_target_id"])
            - float(row[f"{d_after}_target_id"] != row[f"{d_base}_target_id"]),
        ]
    )
    for field in ("source_basis", "target_basis"):
        categories = {
            support: basis_category(row[f"{support}_{field}"])
            for support in (j_after, j_base, d_after, d_base)
        }
        for category in BASIS_CATEGORIES:
            values.append(
                float(categories[j_after] == category)
                - float(categories[j_base] == category)
                - float(categories[d_after] == category)
                + float(categories[d_base] == category)
            )
    output = np.asarray(values, dtype=np.float64)
    if output.shape != (len(STATE_FEATURE_NAMES),):
        raise AssertionError("Sminus feature contract changed")
    return output


def nuisance_matrix(rows: Sequence[dict]) -> tuple[np.ndarray, list[str]]:
    values = []
    for row in rows:
        values.append(
            [
                1.0,
                *(float(row["orientation"] == value) for value in (1, 2, 3)),
                *(float(row["fiber"] == value) for value in (1, 2, 3, 4)),
            ]
        )
    return np.asarray(values, dtype=np.float64), [
        "intercept",
        "orientation_1",
        "orientation_2",
        "orientation_3",
        "fiber_1",
        "fiber_2",
        "fiber_3",
        "fiber_4",
    ]


def fourier_matrix(rows: Sequence[dict], frequencies: Sequence[int]) -> tuple[np.ndarray, list[str]]:
    values = []
    for row in rows:
        angle = 2.0 * math.pi * row["residue"] / GROUP_ORDER
        values.append(
            [
                coordinate
                for frequency in frequencies
                for coordinate in (
                    math.cos(frequency * angle),
                    math.sin(frequency * angle),
                )
            ]
        )
    names = [
        name
        for frequency in frequencies
        for name in (f"residue_cos_{frequency}", f"residue_sin_{frequency}")
    ]
    return np.asarray(values, dtype=np.float64), names


def select_frequencies(
    rows: Sequence[dict],
    y: np.ndarray,
    candidates: Iterable[int],
    count: int,
) -> tuple[int, ...]:
    nuisance, _ = nuisance_matrix(rows)
    nuisance_projection = np.linalg.pinv(nuisance, rcond=1e-12)
    residual_y = y - nuisance @ (nuisance_projection @ y)
    scores = []
    for frequency in candidates:
        pair, _ = fourier_matrix(rows, (frequency,))
        residual_pair = pair - nuisance @ (nuisance_projection @ pair)
        gram = residual_pair.T @ residual_pair
        cross = residual_pair.T @ residual_y
        reduction = float(cross @ np.linalg.pinv(gram, rcond=1e-12) @ cross)
        scores.append((reduction, int(frequency)))
    scores.sort(key=lambda item: (-item[0], item[1]))
    return tuple(sorted(frequency for _, frequency in scores[:count]))


def design_matrix(
    rows: Sequence[dict], frequencies: Sequence[int], include_state: bool
) -> tuple[np.ndarray, list[str]]:
    nuisance, nuisance_names = nuisance_matrix(rows)
    fourier, fourier_names = fourier_matrix(rows, frequencies)
    blocks = [nuisance, fourier]
    names = [*nuisance_names, *fourier_names]
    if include_state:
        blocks.append(np.vstack([state_features(row) for row in rows]))
        names.extend(STATE_FEATURE_NAMES)
    return np.column_stack(blocks), names


def fit_standardized_ridge(x: np.ndarray, y: np.ndarray, ridge: float) -> dict:
    if x.ndim != 2 or x.shape[0] != len(y) or not np.allclose(x[:, 0], 1.0):
        raise ValueError("ridge design must start with an intercept")
    means = x[:, 1:].mean(axis=0)
    scales = x[:, 1:].std(axis=0)
    keep = scales > max(float(scales.max(initial=0.0)) * 1e-12, 1e-14)
    z = np.column_stack(
        [np.ones(len(x)), (x[:, 1:][:, keep] - means[keep]) / scales[keep]]
    )
    penalty = np.eye(z.shape[1], dtype=np.float64) * ridge
    penalty[0, 0] = 0.0
    try:
        beta = np.linalg.solve(z.T @ z + penalty, z.T @ y)
    except np.linalg.LinAlgError:
        beta = np.linalg.pinv(z.T @ z + penalty, rcond=1e-12) @ (z.T @ y)
    return {"means": means, "scales": scales, "keep": keep, "beta": beta}


def predict_standardized(model: dict, x: np.ndarray) -> np.ndarray:
    keep = model["keep"]
    z = np.column_stack(
        [
            np.ones(len(x)),
            (x[:, 1:][:, keep] - model["means"][keep])
            / model["scales"][keep],
        ]
    )
    return z @ model["beta"]


def cluster_summary(values: np.ndarray) -> dict:
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    mean = float(values.mean())
    standard_error = float(values.std(ddof=1) / math.sqrt(len(values)))
    statistic = mean / standard_error if standard_error else math.copysign(math.inf, mean)
    return {
        "clusters": int(len(values)),
        "mean": mean,
        "standard_error": standard_error,
        "t": statistic,
        "degrees_of_freedom": int(len(values) - 1),
        "one_sided_p_positive": (
            float(student_t.sf(statistic, len(values) - 1))
            if standard_error
            else (0.0 if mean > 0 else 1.0)
        ),
        "ci95": [mean - 1.96 * standard_error, mean + 1.96 * standard_error],
    }


def batch_loss_summary(
    rows: Sequence[dict],
    y: np.ndarray,
    pred_spec: np.ndarray,
    pred_state: np.ndarray,
    batches: int,
    hand: str | None = None,
) -> dict:
    mask = np.asarray(
        [hand is None or row["hand"] == hand for row in rows], dtype=bool
    )
    spec_error = np.square(y - pred_spec)
    state_error = np.square(y - pred_state)
    batch_spec = np.full(batches, np.nan)
    batch_state = np.full(batches, np.nan)
    for batch in range(batches):
        selected = mask & np.asarray([row["batch"] == batch for row in rows])
        if selected.any():
            batch_spec[batch] = float(spec_error[selected].mean())
            batch_state[batch] = float(state_error[selected].mean())
    difference = batch_spec - batch_state
    valid = np.isfinite(difference)
    return {
        "rows": int(mask.sum()),
        "held_out_mse_spec": float(spec_error[mask].mean()),
        "held_out_mse_state": float(state_error[mask].mean()),
        "row_weighted_relative_gain": float(
            1.0 - state_error[mask].mean() / spec_error[mask].mean()
        ),
        "equal_batch_loss_reduction": cluster_summary(difference[valid]),
        "batch_loss_spec": batch_spec[valid].tolist(),
        "batch_loss_state": batch_state[valid].tolist(),
        "batch_loss_reduction": difference[valid].tolist(),
    }


def residual_fourier(rows: Sequence[dict], residual: np.ndarray, hand: str) -> dict:
    indices = np.asarray([row["hand"] == hand for row in rows], dtype=bool)
    selected_rows = [row for row, keep in zip(rows, indices) if keep]
    selected_residual = residual[indices]
    amplitudes = {}
    for frequency in range(1, 51):
        pair, _ = fourier_matrix(selected_rows, (frequency,))
        design = np.column_stack([np.ones(len(pair)), pair])
        beta = np.linalg.pinv(design, rcond=1e-12) @ selected_residual
        amplitudes[str(frequency)] = float(math.hypot(beta[1], beta[2]))
    peak = max(amplitudes, key=lambda value: amplitudes[value])
    return {
        "mean": float(selected_residual.mean()),
        "rmse": float(np.sqrt(np.mean(np.square(selected_residual)))),
        "peak_frequency": int(peak),
        "peak_amplitude": amplitudes[peak],
        "amplitude_k1": amplitudes["1"],
        "amplitude_k10": amplitudes["10"],
        "all_frequency_amplitudes": amplitudes,
    }


def markdown_report(payload: dict) -> str:
    joint = payload["held_out"]["joint"]
    reduction = joint["equal_batch_loss_reduction"]
    decision = payload["decision"]
    lines = [
        "# P250 N505 existing-data state increment",
        "",
        "This analysis adds no samples. It cross-fits the archived 200k/400-batch",
        "defined-event table, with whole batches (and therefore shared replicas and",
        "both child views) kept together.",
        "",
        "## Result",
        "",
        f"- M_spec held-out MSE: `{joint['held_out_mse_spec']:.8g}`.",
        f"- M_state held-out MSE: `{joint['held_out_mse_state']:.8g}`.",
        f"- Row-weighted relative gain: `{100.0 * joint['row_weighted_relative_gain']:.3f}%`.",
        "- Equal-batch loss reduction (M_spec - M_state): "
        f"`{reduction['mean']:.8g} +/- {reduction['standard_error']:.3g}`; "
        f"one-sided cluster t p=`{reduction['one_sided_p_positive']:.4g}`.",
        f"- Frozen decision: **{decision['label']}**.",
        "",
        "| hand | rows | M_spec MSE | M_state MSE | relative gain | batch p |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for hand in HANDS:
        row = payload["held_out"][hand]
        lines.append(
            f"| {hand} | {row['rows']} | {row['held_out_mse_spec']:.8g} | "
            f"{row['held_out_mse_state']:.8g} | "
            f"{100.0 * row['row_weighted_relative_gain']:.3f}% | "
            f"{row['equal_batch_loss_reduction']['one_sided_p_positive']:.4g} |"
        )
    lines.extend(
        [
            "",
            "## What the two models contain",
            "",
            "M_spec is fit separately in each hand. In every outer fold it selects",
            "two residue Fourier frequencies using training rows only, then fits those",
            "four sine/cosine coordinates together with orientation and fibre controls.",
            "M_state uses the identical selected frequencies and adds only Sminus:",
            "the antisymmetric change of rank, canonical basis class, support site phase,",
            "and component/site-change indicators between the two ordered intermediate",
            "supports. No L or R coordinate is a predictor.",
            "",
            "Selected training-only frequencies by fold:",
            "",
        ]
    )
    for hand in HANDS:
        lines.append(
            f"- {hand}: "
            + ", ".join(
                f"k={frequency} ({count}/10 folds)"
                for frequency, count in payload["selected_frequency_counts"][hand].items()
            )
        )
    lines.extend(
        [
            "",
            "## Boundary",
            "",
            payload["spectral_boundary"],
            "",
            payload["claim_boundary"],
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--defined", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-md", type=Path, required=True)
    args = parser.parse_args()

    manifest = load_manifest(args.manifest)
    if sha256(Path(__file__)) != manifest["software"]["scorer_sha256"]:
        raise ValueError("scorer hash differs from the frozen implementation")
    rows = load_rows(args.defined, manifest)
    y = np.asarray([row["R_minus"] for row in rows], dtype=np.float64)
    batches = manifest["cross_fitting"]["batches"]
    folds = manifest["cross_fitting"]["outer_folds"]
    candidates = range(1, manifest["models"]["frequency_candidates_max"] + 1)
    frequency_count = manifest["models"]["training_selected_frequencies_per_hand"]
    ridge = manifest["models"]["ridge_lambda_standardized"]

    pred_spec = np.full(len(rows), np.nan, dtype=np.float64)
    pred_state = np.full(len(rows), np.nan, dtype=np.float64)
    selected_by_fold: dict[str, list[dict]] = {hand: [] for hand in HANDS}
    for fold in range(folds):
        test_fold = np.asarray([row["batch"] % folds == fold for row in rows])
        for hand in HANDS:
            hand_mask = np.asarray([row["hand"] == hand for row in rows])
            train = hand_mask & ~test_fold
            test = hand_mask & test_fold
            train_rows = [row for row, keep in zip(rows, train) if keep]
            test_rows = [row for row, keep in zip(rows, test) if keep]
            frequencies = select_frequencies(
                train_rows, y[train], candidates, frequency_count
            )
            selected_by_fold[hand].append(
                {"fold": fold, "frequencies": list(frequencies)}
            )
            x_train_spec, _ = design_matrix(train_rows, frequencies, False)
            x_test_spec, _ = design_matrix(test_rows, frequencies, False)
            x_train_state, _ = design_matrix(train_rows, frequencies, True)
            x_test_state, _ = design_matrix(test_rows, frequencies, True)
            model_spec = fit_standardized_ridge(x_train_spec, y[train], ridge)
            model_state = fit_standardized_ridge(x_train_state, y[train], ridge)
            pred_spec[test] = predict_standardized(model_spec, x_test_spec)
            pred_state[test] = predict_standardized(model_state, x_test_state)
    if not np.isfinite(pred_spec).all() or not np.isfinite(pred_state).all():
        raise AssertionError("cross-fitting did not predict every row exactly once")

    held_out = {
        "joint": batch_loss_summary(
            rows, y, pred_spec, pred_state, batches, hand=None
        )
    }
    for hand in HANDS:
        held_out[hand] = batch_loss_summary(
            rows, y, pred_spec, pred_state, batches, hand=hand
        )
    counts = {
        hand: {
            str(frequency): count
            for frequency, count in sorted(
                Counter(
                    frequency
                    for fold in selected_by_fold[hand]
                    for frequency in fold["frequencies"]
                ).items()
            )
        }
        for hand in HANDS
    }
    residuals = {
        "M_spec": {
            hand: residual_fourier(rows, y - pred_spec, hand) for hand in HANDS
        },
        "M_state": {
            hand: residual_fourier(rows, y - pred_state, hand) for hand in HANDS
        },
    }
    alpha = manifest["decision"]["alpha"]
    joint_increment = held_out["joint"]["equal_batch_loss_reduction"]
    per_hand_positive = all(
        held_out[hand]["equal_batch_loss_reduction"]["mean"] > 0
        for hand in HANDS
    )
    pass_gate = (
        joint_increment["mean"] > 0
        and joint_increment["one_sided_p_positive"] < alpha
        and per_hand_positive
    )
    payload = {
        "schema": SCHEMA,
        "status": "completed_existing_data_cross_fit",
        "freeze_commit": manifest["freeze_commit"],
        "input": {
            "path": str(args.defined),
            "sha256": sha256(args.defined),
            "rows": len(rows),
            "batches": batches,
            "replica_clusters": len({row["replica"] for row in rows}),
            "two_hand_replica_clusters": sum(
                len({row["hand"] for row in rows if row["replica"] == replica}) == 2
                for replica in {row["replica"] for row in rows}
            ),
        },
        "models": manifest["models"],
        "state_feature_order": list(STATE_FEATURE_NAMES),
        "selected_frequencies_by_fold": selected_by_fold,
        "selected_frequency_counts": counts,
        "held_out": held_out,
        "residual_diagnostics": residuals,
        "decision": {
            "alpha": alpha,
            "requires": (
                "positive equal-batch loss reduction jointly, one-sided cluster "
                "t p<alpha, and positive point reduction in each child"
            ),
            "passed": pass_gate,
            "label": (
                "pre-outcome antisymmetric state block adds held-out information"
                if pass_gate
                else "no certified held-out increment from the frozen state block"
            ),
        },
        "spectral_boundary": manifest["spectral_boundary"],
        "claim_boundary": manifest["claim_boundary"],
        "software": {
            "scorer": str(Path(__file__)),
            "scorer_sha256": sha256(Path(__file__)),
        },
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_md.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    args.output_md.write_text(markdown_report(payload), encoding="utf-8")


if __name__ == "__main__":
    main()

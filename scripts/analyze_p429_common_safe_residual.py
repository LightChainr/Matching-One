#!/usr/bin/env python3
"""Cross-fit the P429 common-safe clone-dependence remainder.

The analysis unit is one pre-common checkpoint row whose common update was
safe, together with its paired clone responses.  The code never expands the
clones into independent rows, and every held-out fold contains complete
batches from both orientations and both sizes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import scipy
from scipy.optimize import minimize
from scipy.special import expit


ROOT = Path(__file__).resolve().parents[1]
ENVIRONMENTS = ("N325_first", "N325_second", "N425_first", "N425_second")
SIZES = ("N325", "N425")
ORIENTATIONS = ("first", "second")
BINARY_FEATURES: set[str] = set()
REQUIRED_COLUMNS = {
    "n", "a", "b", "orientation", "batch", "replica", "k0", "age_steps",
    "ell_u", "ell_v", "essential_size", "essential_carriers",
    "occupied_frontier", "vacant_frontier", "boundary_cut_edges",
    "boundary_multicontact_sites", "boundary_contact_pairs", "core_vertices",
    "core_edges", "articulation_vertices", "bridges",
    "boundary_axis_imbalance", "boundary_corner_balance", "frontier_components",
    "largest_frontier_component", "frontier_component_sumsq", "H2",
    "H2_theta", "H2_figure8", "H2_separate", "H2_direction_positive",
    "H2_direction_negative", "H2_direction_mixed", "next_exit",
    "checkpoint_b1_safe_count", "branch_common_safe",
    "branch_clone1_survives", "branch_clone2_survives", "branch_both_survive",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def derive_features(frame: pd.DataFrame) -> pd.DataFrame:
    """Return the fixed, outcome-blind transforms declared in the contract."""
    out = frame.copy()
    n = out["n"].astype(float)
    length = np.sqrt(n)
    remaining = (out["n"] - out["k0"]).astype(float)
    h2 = out["H2"].astype(float)
    frontier = out["vacant_frontier"].astype(float).clip(lower=1.0)
    pairs = out["boundary_contact_pairs"].astype(float).clip(lower=1.0)
    cut = out["boundary_cut_edges"].astype(float).clip(lower=1.0)

    out["h2_rate"] = h2 / remaining
    out["h2_direction_balance"] = (
        out["H2_direction_positive"] - out["H2_direction_negative"]
    ) / h2.clip(lower=1.0)
    out["multicontact_fraction"] = out["boundary_multicontact_sites"] / frontier
    out["contact_pair_intensity"] = (
        out["boundary_contact_pairs"] /
        out["boundary_multicontact_sites"].astype(float).clip(lower=1.0)
    )
    out["corner_balance_ratio"] = out["boundary_corner_balance"] / pairs

    out["essential_fraction"] = out["essential_size"] / n
    out["vacant_frontier_L"] = out["vacant_frontier"] / length
    out["core_cycle_density"] = (out["core_edges"] - out["core_vertices"]) / n
    out["articulation_L"] = out["articulation_vertices"] / length
    out["axis_anisotropy"] = out["boundary_axis_imbalance"].abs() / cut
    out["frontier_components_L"] = out["frontier_components"] / length
    out["largest_frontier_L"] = out["largest_frontier_component"] / length
    out["frontier_concentration"] = out["frontier_component_sumsq"] / (frontier ** 2)
    out["age_fraction"] = out["age_steps"] / n

    # P has columns equal to the lifted period vectors.  For these four fixed
    # designs P=((a,b),(0,1)), so P*ell=(a*ell_u+b*ell_v, ell_v).
    physical_x = out["a"] * out["ell_u"] + out["b"] * out["ell_v"]
    physical_y = out["ell_v"]
    radius2 = physical_x.astype(float) ** 2 + physical_y.astype(float) ** 2
    _require(bool((radius2 > 0).all()), "rank-one physical line must be nonzero")
    out["line_chi4_re"] = (
        physical_x.astype(float) ** 4
        - 6.0 * physical_x.astype(float) ** 2 * physical_y.astype(float) ** 2
        + physical_y.astype(float) ** 4
    ) / (radius2 ** 2)
    out["line_chi4_im"] = (
        4.0 * physical_x.astype(float) * physical_y.astype(float)
        * (physical_x.astype(float) ** 2 - physical_y.astype(float) ** 2)
    ) / (radius2 ** 2)
    return out


def load_inputs(contract: dict, root: Path) -> tuple[pd.DataFrame, dict]:
    raw_lock_path = root / contract["inputs"]["raw_lock"]
    raw_lock = json.loads(raw_lock_path.read_text())
    parent_score_path = root / raw_lock["score"]["path"]
    parent_score_hash = sha256(parent_score_path)
    _require(parent_score_hash == raw_lock["score"]["sha256"],
             "parent production score hash drifted")
    parent_score = json.loads(parent_score_path.read_text())
    frames: list[pd.DataFrame] = []
    audit: dict = {
        "raw_lock": {"path": str(raw_lock_path.relative_to(root)),
                     "sha256": sha256(raw_lock_path)},
        "parent_score": {
            "path": str(parent_score_path.relative_to(root)),
            "sha256": parent_score_hash,
            "secondary_vector_order": parent_score[
                "secondary_successor_heterogeneity"]["vector_order"],
            "secondary_vector": parent_score[
                "secondary_successor_heterogeneity"]["vector"],
            "secondary_size_common_gap": parent_score[
                "secondary_successor_heterogeneity"]["size_common_gap"],
        },
        "sizes": {},
    }

    for size in SIZES:
        path = root / contract["inputs"][f"{size}_csv"]
        expected_hash = raw_lock["production"][size]["csv_sha256"]
        observed_hash = sha256(path)
        _require(observed_hash == expected_hash, f"{size} production CSV hash drifted")
        frame = pd.read_csv(path)
        missing = sorted(REQUIRED_COLUMNS - set(frame.columns))
        _require(not missing, f"{size} missing required columns: {missing}")
        _require(len(frame) == raw_lock["production"][size]["at_risk_rows"],
                 f"{size} at-risk row count drifted")
        _require(set(frame["orientation"]) == set(ORIENTATIONS),
                 f"{size} orientation support drifted")
        _require(set(frame["batch"]) == set(range(100)),
                 f"{size} batch support must be 0..99")
        _require(bool((frame["n"] == int(size[1:])).all()), f"{size} n drifted")
        _require(bool((frame["a"] == frame["n"]).all()),
                 f"{size} period-matrix a must equal n")
        for column in ("branch_common_safe", "branch_clone1_survives",
                       "branch_clone2_survives", "branch_both_survive"):
            _require(set(frame[column].unique()).issubset({0, 1}),
                     f"{size} {column} is not binary")
        _require(bool((frame["branch_both_survive"] ==
                       frame["branch_clone1_survives"] *
                       frame["branch_clone2_survives"]).all()),
                 f"{size} paired clone product drifted")
        _require(bool((frame["checkpoint_b1_safe_count"] ==
                       frame["n"] - frame["k0"] - frame["H2"]).all()),
                 f"{size} b1/H2 identity drifted")
        _require(bool((frame["H2"] == frame["H2_theta"] +
                       frame["H2_figure8"] + frame["H2_separate"]).all()),
                 f"{size} H2 type decomposition drifted")
        _require(bool((frame["H2"] == frame["H2_direction_positive"] +
                       frame["H2_direction_negative"] +
                       frame["H2_direction_mixed"]).all()),
                 f"{size} H2 direction decomposition drifted")

        safe = frame.loc[frame["branch_common_safe"] == 1].copy()
        _require(bool((safe["next_exit"] == 0).all()),
                 f"{size} common-safe rows contain a common exit")
        _require(bool((safe["H2_theta"] == safe["H2"]).all()),
                 f"{size} H2_theta no longer equals H2")
        _require(int(safe["H2_figure8"].abs().sum()) == 0,
                 f"{size} H2_figure8 no longer degenerate")
        _require(int(safe["H2_separate"].abs().sum()) == 0,
                 f"{size} H2_separate no longer degenerate")
        safe["size"] = size
        safe["environment"] = size + "_" + safe["orientation"].astype(str)
        safe["fold"] = safe["batch"] % 5
        safe["u"] = (safe["branch_clone1_survives"] +
                     safe["branch_clone2_survives"]) / 2.0
        safe = derive_features(safe)
        frames.append(safe)

        first_safe = set(safe.loc[safe["orientation"] == "first", "replica"])
        second_safe = set(safe.loc[safe["orientation"] == "second", "replica"])
        audit["sizes"][size] = {
            "path": str(path.relative_to(root)),
            "sha256": observed_hash,
            "at_risk_rows": int(len(frame)),
            "common_safe_rows": int(len(safe)),
            "common_safe_rows_by_orientation": {
                orientation: int((safe["orientation"] == orientation).sum())
                for orientation in ORIENTATIONS
            },
            "joint_orientation_common_safe_replicas": int(len(first_safe & second_safe)),
            "exact_line_support_by_orientation": {
                orientation: int(safe.loc[safe["orientation"] == orientation,
                                          ["ell_u", "ell_v"]].drop_duplicates().shape[0])
                for orientation in ORIENTATIONS
            },
            "top4_exact_line_fraction_by_orientation": {
                orientation: float(
                    safe.loc[safe["orientation"] == orientation]
                    .groupby(["ell_u", "ell_v"]).size().nlargest(4).sum()
                    / max(int((safe["orientation"] == orientation).sum()), 1)
                )
                for orientation in ORIENTATIONS
            },
            "rare_support_audit": {
                "H2_direction_mixed_nonzero_fraction": float(
                    (safe["H2_direction_mixed"] != 0).mean()),
                "H2_direction_mixed_max": int(safe["H2_direction_mixed"].max()),
                "multiple_essential_carriers_fraction": float(
                    (safe["essential_carriers"] > 1).mean()),
            },
            "batches": int(safe["batch"].nunique()),
        }

    data = pd.concat(frames, ignore_index=True)
    _require(set(data["environment"]) == set(ENVIRONMENTS),
             "environment support drifted")
    return data, audit


def _logit(probability: float) -> float:
    probability = min(max(probability, 1e-6), 1.0 - 1e-6)
    return math.log(probability / (1.0 - probability))


def fit_fractional_logit(
    design: np.ndarray,
    target: np.ndarray,
    weights: np.ndarray,
    penalty_mask: np.ndarray,
    ridge_lambda: float,
    initial: np.ndarray,
) -> tuple[np.ndarray, dict]:
    """Fit a fixed ridge-logit nuisance model to pair-average responses."""
    weights = np.asarray(weights, dtype=float)
    weights = weights / weights.sum()

    def objective(beta: np.ndarray) -> tuple[float, np.ndarray]:
        eta = design @ beta
        loss = float(np.sum(weights * (np.logaddexp(0.0, eta) - target * eta)))
        loss += 0.5 * ridge_lambda * float(np.sum((beta * penalty_mask) ** 2))
        gradient = design.T @ (weights * (expit(eta) - target))
        gradient += ridge_lambda * penalty_mask * beta
        return loss, gradient

    result = minimize(
        lambda beta: objective(beta)[0],
        initial,
        method="L-BFGS-B",
        jac=lambda beta: objective(beta)[1],
        options={"maxiter": 600, "ftol": 1e-12, "gtol": 1e-8, "maxls": 40},
    )
    _require(bool(result.success), f"fractional logit failed: {result.message}")
    return np.asarray(result.x), {
        "converged": bool(result.success),
        "iterations": int(result.nit),
        "objective": float(result.fun),
        "gradient_max_abs": float(np.max(np.abs(result.jac))),
    }


def make_design(
    train: pd.DataFrame,
    test: pd.DataFrame,
    features: list[str],
) -> tuple[np.ndarray, np.ndarray, list[str], dict]:
    train_columns = [(train["environment"] == env).to_numpy(dtype=float)
                     for env in ENVIRONMENTS]
    test_columns = [(test["environment"] == env).to_numpy(dtype=float)
                    for env in ENVIRONMENTS]
    names = [f"intercept:{env}" for env in ENVIRONMENTS]
    scaling: dict = {}
    for feature in features:
        train_values = np.zeros(len(train), dtype=float)
        test_values = np.zeros(len(test), dtype=float)
        scaling[feature] = {}
        for env in ENVIRONMENTS:
            train_mask = (train["environment"] == env).to_numpy()
            test_mask = (test["environment"] == env).to_numpy()
            raw_train = train.loc[train_mask, feature].to_numpy(dtype=float)
            raw_test = test.loc[test_mask, feature].to_numpy(dtype=float)
            if feature in BINARY_FEATURES:
                center, scale = 0.0, 1.0
                transformed_train, transformed_test = raw_train, raw_test
            else:
                center = float(np.mean(raw_train))
                scale = float(np.std(raw_train, ddof=0))
                _require(scale > 1e-12, f"{feature} degenerate in {env} training fold")
                transformed_train = np.clip((raw_train - center) / scale, -8.0, 8.0)
                transformed_test = np.clip((raw_test - center) / scale, -8.0, 8.0)
            train_values[train_mask] = transformed_train
            test_values[test_mask] = transformed_test
            scaling[feature][env] = {"center": center, "scale": scale}
        train_columns.append(train_values)
        test_columns.append(test_values)
        names.append(feature)
    return (np.column_stack(train_columns), np.column_stack(test_columns), names, scaling)


def crossfit_candidate(
    data: pd.DataFrame,
    features: list[str],
    ridge_lambda: float,
) -> tuple[np.ndarray, list[dict]]:
    prediction = np.full(len(data), np.nan, dtype=float)
    fold_fits: list[dict] = []
    for fold in range(5):
        train = data.loc[data["fold"] != fold]
        test = data.loc[data["fold"] == fold]
        train_design, test_design, names, scaling = make_design(train, test, features)
        target = train["u"].to_numpy(dtype=float)
        weights = np.zeros(len(train), dtype=float)
        initial = np.zeros(train_design.shape[1], dtype=float)
        for index, env in enumerate(ENVIRONMENTS):
            mask = (train["environment"] == env).to_numpy()
            _require(bool(mask.any()), f"empty training environment {env}")
            weights[mask] = 0.25 / float(mask.sum())
            initial[index] = _logit(float(np.mean(target[mask])))
        penalty_mask = np.r_[np.zeros(len(ENVIRONMENTS)), np.ones(len(features))]
        beta, diagnostics = fit_fractional_logit(
            train_design, target, weights, penalty_mask, ridge_lambda, initial)
        heldout = expit(test_design @ beta)
        prediction[test.index.to_numpy()] = heldout
        fold_fits.append({
            "fold": fold,
            "heldout_batches": list(range(fold, 100, 5)),
            "train_rows": int(len(train)),
            "test_rows": int(len(test)),
            "coefficients": dict(zip(names, beta.tolist())),
            "scaling": scaling,
            **diagnostics,
        })
    _require(bool(np.isfinite(prediction).all()), "cross-fit prediction incomplete")
    return prediction, fold_fits


def cluster_metric(values: np.ndarray, data: pd.DataFrame) -> dict:
    """Conditional score dispersion for the four environment means.

    A size-batch is one cluster and contributes a two-orientation influence
    vector.  The resulting 4x4 matrix is block diagonal only conditional on
    the already fitted fold models.  Because slopes are shared across sizes,
    this is not a full nuisance-refit sampling covariance.  No all-size
    meta-estimate is formed.
    """
    values = np.asarray(values, dtype=float)
    estimates = np.zeros(len(ENVIRONMENTS), dtype=float)
    covariance = np.zeros((len(ENVIRONMENTS), len(ENVIRONMENTS)), dtype=float)
    counts = np.zeros(len(ENVIRONMENTS), dtype=int)
    for index, env in enumerate(ENVIRONMENTS):
        mask = (data["environment"] == env).to_numpy()
        estimates[index] = float(np.mean(values[mask]))
        counts[index] = int(mask.sum())

    for size_index, size in enumerate(SIZES):
        positions = [2 * size_index, 2 * size_index + 1]
        influence = np.zeros((100, 2), dtype=float)
        for batch in range(100):
            for local, position in enumerate(positions):
                env = ENVIRONMENTS[position]
                mask = ((data["environment"] == env) &
                        (data["batch"] == batch)).to_numpy()
                influence[batch, local] = (
                    float(np.sum(values[mask])) - estimates[position] * int(mask.sum())
                )
        block = (100.0 / 99.0) * (influence.T @ influence)
        block /= np.outer(counts[positions], counts[positions])
        covariance[np.ix_(positions, positions)] = block

    dispersion_se = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    by_environment = {
        env: {
            "estimate": float(estimates[index]),
            "conditional_score_dispersion_se": float(dispersion_se[index]),
            "conditional_normal_reference_range": [
                float(estimates[index] - 1.96 * dispersion_se[index]),
                float(estimates[index] + 1.96 * dispersion_se[index]),
            ],
            "rows": int(counts[index]),
        }
        for index, env in enumerate(ENVIRONMENTS)
    }
    by_size = {}
    for size_index, size in enumerate(SIZES):
        positions = [2 * size_index, 2 * size_index + 1]
        weight = np.asarray([0.5, 0.5])
        block = covariance[np.ix_(positions, positions)]
        estimate = float(weight @ estimates[positions])
        se = float(math.sqrt(max(float(weight @ block @ weight), 0.0)))
        by_size[size] = {
            "equal_orientation_mean": estimate,
            "conditional_score_dispersion_se_with_direction_covariance": se,
            "conditional_normal_reference_range": [
                estimate - 1.96 * se, estimate + 1.96 * se],
            "conditional_orientation_score_covariance": block.tolist(),
        }
    return {
        "vector_order": list(ENVIRONMENTS),
        "vector": estimates.tolist(),
        "conditional_on_fold_models_batch_score_covariance": covariance.tolist(),
        "conditional_score_dispersion_se": dispersion_se.tolist(),
        "by_environment": by_environment,
        "by_size": by_size,
    }


def feature_support(data: pd.DataFrame, features: Iterable[str]) -> dict:
    output = {}
    for feature in features:
        values = data[feature].to_numpy(dtype=float)
        output[feature] = {
            "min": float(np.min(values)),
            "max": float(np.max(values)),
            "mean": float(np.mean(values)),
            "sd": float(np.std(values, ddof=0)),
            "nonzero_fraction": float(np.mean(values != 0)),
        }
    return output


def coefficient_summary(fold_fits: list[dict], features: list[str]) -> dict:
    output = {}
    for feature in features:
        values = np.asarray([fit["coefficients"][feature] for fit in fold_fits])
        output[feature] = {
            "mean_standardized_log_odds": float(np.mean(values)),
            "min_across_folds": float(np.min(values)),
            "max_across_folds": float(np.max(values)),
            "same_sign_all_folds": bool(np.all(values > 0) or np.all(values < 0)),
        }
    return output


def score(
    contract: dict,
    data: pd.DataFrame,
    audit: dict,
    contract_path: Path,
    analyzer_path: Path,
) -> tuple[dict, dict]:
    candidates = contract["candidate_states_in_fixed_nested_order"]
    cumulative: list[str] = []
    candidate_features: dict[str, list[str]] = {}
    for candidate in candidates:
        cumulative.extend(candidate["add"])
        candidate_features[candidate["id"]] = list(cumulative)
    all_features = candidate_features[candidates[-1]["id"]]
    ridge_lambda = float(contract["model"]["ridge_lambda"])

    y1 = data["branch_clone1_survives"].to_numpy(dtype=float)
    y2 = data["branch_clone2_survives"].to_numpy(dtype=float)
    fitted: dict = {}
    arrays: dict = {}
    for candidate in candidates:
        identifier = candidate["id"]
        features = candidate_features[identifier]
        prediction, fold_fits = crossfit_candidate(data, features, ridge_lambda)
        prediction_for_loss = np.clip(prediction, 1e-8, 1.0 - 1e-8)
        residual = (y1 - prediction) * (y2 - prediction)
        log_loss = -0.5 * (
            y1 * np.log(prediction_for_loss) + (1.0 - y1) * np.log1p(-prediction_for_loss)
            + y2 * np.log(prediction_for_loss) + (1.0 - y2) * np.log1p(-prediction_for_loss)
        )
        brier = 0.5 * ((y1 - prediction) ** 2 + (y2 - prediction) ** 2)
        arrays[identifier] = {
            "prediction": prediction,
            "residual": residual,
            "log_loss": log_loss,
            "brier": brier,
        }
        fitted[identifier] = {
            "features": features,
            "heldout_residual_dependence": cluster_metric(residual, data),
            "heldout_log_loss": cluster_metric(log_loss, data),
            "heldout_brier": cluster_metric(brier, data),
            "heldout_prediction": cluster_metric(prediction, data),
            "coefficient_summary": coefficient_summary(fold_fits, features),
            "fold_fits": fold_fits,
        }

    baseline_id = candidates[0]["id"]
    previous_id = baseline_id
    for candidate in candidates:
        identifier = candidate["id"]
        absorbed = arrays[baseline_id]["residual"] - arrays[identifier]["residual"]
        loss_improvement = arrays[baseline_id]["log_loss"] - arrays[identifier]["log_loss"]
        brier_improvement = arrays[baseline_id]["brier"] - arrays[identifier]["brier"]
        fitted[identifier]["absorbed_vs_intercept"] = cluster_metric(absorbed, data)
        fitted[identifier]["log_loss_improvement_vs_intercept"] = cluster_metric(
            loss_improvement, data)
        fitted[identifier]["brier_improvement_vs_intercept"] = cluster_metric(
            brier_improvement, data)
        block_absorbed = arrays[previous_id]["residual"] - arrays[identifier]["residual"]
        block_loss = arrays[previous_id]["log_loss"] - arrays[identifier]["log_loss"]
        fitted[identifier]["increment_vs_previous"] = {
            "previous": previous_id,
            "absorbed_dependence": cluster_metric(block_absorbed, data),
            "log_loss_improvement": cluster_metric(block_loss, data),
        }
        previous_id = identifier

    baseline = fitted[baseline_id]["heldout_residual_dependence"]
    for candidate in candidates:
        identifier = candidate["id"]
        for size in SIZES:
            base_value = baseline["by_size"][size]["equal_orientation_mean"]
            absorbed_value = fitted[identifier]["absorbed_vs_intercept"]["by_size"][size][
                "equal_orientation_mean"]
            fitted[identifier].setdefault("absorbed_fraction_by_size", {})[size] = (
                absorbed_value / base_value if base_value != 0 else None)

    outcome_audit = {}
    for env in ENVIRONMENTS:
        members = data.loc[data["environment"] == env]
        mean1 = float(members["branch_clone1_survives"].mean())
        mean2 = float(members["branch_clone2_survives"].mean())
        both = float(members["branch_both_survive"].mean())
        outcome_audit[env] = {
            "clone1_mean": mean1,
            "clone2_mean": mean2,
            "pair_mean": 0.5 * (mean1 + mean2),
            "conventional_covariance": both - mean1 * mean2,
            "symmetry_locked_covariance": both - (0.5 * (mean1 + mean2)) ** 2,
            "symmetry_lock_difference": -0.25 * (mean1 - mean2) ** 2,
        }

    parent_labels = audit["parent_score"]["secondary_vector_order"]
    parent_values = audit["parent_score"]["secondary_vector"]
    parent_map = dict(zip(parent_labels, parent_values))
    parent_comparison = {"by_environment": {}, "by_size": {}}
    for env in ENVIRONMENTS:
        size, orientation = env.split("_", 1)
        parent_value = float(parent_map[f"conditional_gap:{size}:{orientation}"])
        observed = outcome_audit[env]["conventional_covariance"]
        _require(abs(parent_value - observed) < 1e-14,
                 f"{env} common-safe covariance no longer reproduces parent score")
        parent_comparison["by_environment"][env] = {
            "parent_conventional_covariance": parent_value,
            "recomputed_conventional_covariance": observed,
            "difference": observed - parent_value,
        }
    for size in SIZES:
        parent = audit["parent_score"]["secondary_size_common_gap"][size]
        crossfit = fitted[baseline_id]["heldout_residual_dependence"]["by_size"][size]
        parent_comparison["by_size"][size] = {
            "parent_GLS_conventional_gap": parent,
            "crossfit_symmetry_locked_equal_orientation_gap": crossfit,
            "comparison_boundary": (
                "The parent uses full-sample clone-specific means and GLS direction weights; "
                "the reanalysis uses fold-held symmetry-locked means and equal direction weights."
            ),
        }

    payload = {
        "schema": "matching-one/p429-common-safe-reanalysis-score/v1",
        "contract": {
            "path": str(contract_path.resolve()),
            "sha256": sha256(contract_path),
        },
        "analyzer": {
            "path": str(analyzer_path.resolve()),
            "sha256": sha256(analyzer_path),
        },
        "contract_status": contract["status"],
        "analysis_unit": contract["analysis_unit"],
        "dependency_contract": contract["dependency_contract"],
        "input_audit": audit,
        "safe_filter_semantics": (
            "orientation-specific branch_common_safe==1; paired orientations share batch/replica "
            "streams but need not both be at risk or common-safe"
        ),
        "feature_support": feature_support(data, all_features),
        "degeneracy_audit": contract["degeneracy_audit"],
        "outcome_audit": outcome_audit,
        "parent_secondary_reproduction": parent_comparison,
        "candidates": fitted,
        "candidate_order": [candidate["id"] for candidate in candidates],
        "dispersion_boundary": (
            "Reported batch quantities are held-out score dispersions conditional on the fitted "
            "fold models. The nuisance models share slopes across sizes and were not refit in a "
            "cluster bootstrap or delete-one loop, so these are not full sampling standard errors "
            "or confidence intervals. Stored cross-size covariance blocks are zero only under "
            "that conditional view. No clone, orientation, or size is counted as an independent "
            "study, and no all-size significance claim or meta-estimate is formed."
        ),
        "metric_identity": (
            "With one shared prediction for both binary clones, pair-average Brier loss minus "
            "the residual cross-product equals mean[(y1+y2)/2-y1*y2], which is model-invariant. "
            "Therefore Brier gain equals absorbed dependence algebraically and is an audit, not "
            "independent corroboration. Log-loss gain is the separate held-out loss diagnostic."
        ),
        "claim_boundary": contract["claim_boundary"],
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
        },
    }
    return payload, arrays


def _fmt(value: float) -> str:
    return f"{value:.7f}"


def render_report(payload: dict) -> str:
    lines = [
        "# P429 common-safe dependence remainder: blocked cross-fit reanalysis",
        "",
        "## Outcome first",
        "",
        "This is a zero-new-sample, post-reveal predictive allocation on the locked P429",
        "production files.  Five-fold cross-fitting holds out whole batches simultaneously",
        "across both orientations and both sizes.  Each row remains one paired-clone unit.",
        "",
        "The tables below report the symmetry-locked held-out residual product, the amount",
        "predictively absorbed relative to the environment-intercept model, and held-out loss.",
        "Direction means use equal weights with the measured within-batch direction covariance;",
        "the two sizes are deliberately not pooled into a single evidence number.",
        "",
        "## Analysis rows",
        "",
        "| size | at-risk rows | common-safe rows | first / second safe | jointly safe replicas |",
        "|---|---:|---:|---:|---:|",
    ]
    for size in SIZES:
        item = payload["input_audit"]["sizes"][size]
        by_o = item["common_safe_rows_by_orientation"]
        lines.append(
            f"| {size} | {item['at_risk_rows']:,} | {item['common_safe_rows']:,} | "
            f"{by_o['first']:,} / {by_o['second']:,} | "
            f"{item['joint_orientation_common_safe_replicas']:,} |"
        )

    baseline = payload["candidates"]["intercept_only"]
    lines.extend([
        "",
        "## Parent-secondary reproduction",
        "",
        "| size | parent conventional gap (GLS SE) | cross-fit baseline (conditional score dispersion) |",
        "|---|---:|---:|",
    ])
    for size in SIZES:
        parent = payload["parent_secondary_reproduction"]["by_size"][size][
            "parent_GLS_conventional_gap"]
        crossfit = baseline["heldout_residual_dependence"]["by_size"][size]
        lines.append(
            f"| {size} | {_fmt(parent['estimate'])} ({_fmt(parent['se'])}) | "
            f"{_fmt(crossfit['equal_orientation_mean'])} "
            f"({_fmt(crossfit['conditional_score_dispersion_se_with_direction_covariance'])}) |"
        )
    lines.extend([
        "",
        "The environment-level conventional covariances reproduce the locked parent score exactly.",
        "The small size-summary differences above are expected: the parent uses full-sample,",
        "clone-specific means and GLS direction weights, whereas this reanalysis uses fold-held,",
        "symmetry-locked means and equal direction weights.",
    ])
    for size in SIZES:
        lines.extend([
            "",
            f"## {size}: fixed nested candidate states",
            "",
            "| candidate state | features | residual dependence (conditional dispersion) | predictively absorbed | absorbed fraction point | log-loss gain | Brier gain audit |",
            "|---|---:|---:|---:|---:|---:|---:|",
        ])
        for identifier in payload["candidate_order"]:
            candidate = payload["candidates"][identifier]
            residual = candidate["heldout_residual_dependence"]["by_size"][size]
            absorbed = candidate["absorbed_vs_intercept"]["by_size"][size]
            log_gain = candidate["log_loss_improvement_vs_intercept"]["by_size"][size]
            brier_gain = candidate["brier_improvement_vs_intercept"]["by_size"][size]
            fraction = candidate["absorbed_fraction_by_size"][size]
            lines.append(
                f"| `{identifier}` | {len(candidate['features'])} | "
                f"{_fmt(residual['equal_orientation_mean'])} "
                f"({_fmt(residual['conditional_score_dispersion_se_with_direction_covariance'])}) | "
                f"{_fmt(absorbed['equal_orientation_mean'])} | {fraction:.1%} | "
                f"{_fmt(log_gain['equal_orientation_mean'])} | "
                f"{_fmt(brier_gain['equal_orientation_mean'])} |"
            )

        lines.extend([
            "",
            f"### {size}: direction audit (point estimates)",
            "",
            "| candidate state | first residual | second residual |",
            "|---|---:|---:|",
        ])
        for identifier in payload["candidate_order"]:
            candidate = payload["candidates"][identifier]
            by_env = candidate["heldout_residual_dependence"]["by_environment"]
            lines.append(
                f"| `{identifier}` | {_fmt(by_env[f'{size}_first']['estimate'])} | "
                f"{_fmt(by_env[f'{size}_second']['estimate'])} |"
            )

        lines.extend([
            "",
            f"### {size}: incremental contribution of each added block",
            "",
            "| added block endpoint | predictive residual reduction beyond previous state | log-loss gain beyond previous state |",
            "|---|---:|---:|",
        ])
        for identifier in payload["candidate_order"][1:]:
            increment = payload["candidates"][identifier]["increment_vs_previous"]
            absorbed = increment["absorbed_dependence"]["by_size"][size]
            log_gain = increment["log_loss_improvement"]["by_size"][size]
            lines.append(
                f"| `{identifier}` | {_fmt(absorbed['equal_orientation_mean'])} "
                f"(conditional dispersion "
                f"{_fmt(absorbed['conditional_score_dispersion_se_with_direction_covariance'])}) | "
                f"{_fmt(log_gain['equal_orientation_mean'])} |"
            )

    lines.extend([
        "",
        "## Interpretation and boundaries",
        "",
        "- The baseline and all candidate scores use one shared prediction for the two suffix",
        "  streams.  The exact audit in `score.json` reports how little this symmetry lock differs",
        "  from the conventional product-of-separate-means covariance.",
        "- `H2_theta` is exactly `H2` in these rows; `H2_figure8` and `H2_separate` are zero.",
        "  `checkpoint_b1_safe_count` is algebraically redundant with `H2`.  None is presented as",
        "  extra evidence.",
        "- `branch_common_safe=1` is filtered within each orientation.  The jointly safe replica",
        "  counts above are overlap diagnostics only; a both-orientations-safe subset is a different",
        "  post-hoc estimand and is not substituted here.",
        "- The cooperative block contains boundary multicontact/contact-pair proxies.  It is not",
        "  the exact microscopic two-step cooperative-pair count.",
        "- H2 is the exact pre-common count of one-step rank-two hazard sites.  Its dominant",
        "  absorption is therefore primarily risk-set accounting, not a newly identified memory",
        "  mechanism; the common update can still change the unobserved successor H2.",
        "- Every covariate is measured before the shared common update.  The analysis therefore",
        "  tests how much pre-update state predicts successor heterogeneity; it does not observe",
        "  the complete successor state.",
        "- Cross-fitting addresses in-sample prediction bias only.  These already revealed rows",
        "  cannot provide independent confirmation, causal mediation, full-state sufficiency,",
        "  Markov closure/nonclosure at scale, a continuum memory field, or a scale exponent.",
        "- Batch quantities are conditional score dispersions, not nuisance-refit sampling SEs or",
        "  confidence intervals.  Shared cross-size slopes induce unconditional coupling that the",
        "  stored zero cross-size blocks do not estimate.  No significance claim is made.",
        "- Absorbed fractions are point estimates without ratio intervals.  Values above 100% and",
        "  negative residuals are model/noise-scale overcorrection, not exact full absorption or",
        "  evidence of negative dependence.  At N425, direction averaging also hides opposite-sign",
        "  first/second residual point estimates, which are shown explicitly above.",
        "- Brier gain is algebraically identical to absorbed dependence under the shared clone",
        "  prediction, so it is only an internal audit.  Log-loss gain is the separate held-out",
        "  predictive-loss diagnostic.",
        "",
        "Full feature transforms, fold fits, four-environment covariance matrices, source hashes,",
        "and support audits are in `score.json`.",
        "",
    ])
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--contract", type=Path,
        default=ROOT / "analysis/p429_common_safe_reanalysis_contract.json")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()

    contract = json.loads(args.contract.read_text())
    _require(contract["schema"] == "matching-one/p429-common-safe-reanalysis-contract/v1",
             "unexpected contract schema")
    data, audit = load_inputs(contract, ROOT)
    payload, _ = score(
        contract, data, audit, args.contract, Path(__file__).resolve())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(render_report(payload))
    concise = {
        "safe_rows": {size: audit["sizes"][size]["common_safe_rows"] for size in SIZES},
        "absorbed_fraction_by_size": {
            identifier: payload["candidates"][identifier]["absorbed_fraction_by_size"]
            for identifier in payload["candidate_order"]
        },
        "final_residual_by_size": payload["candidates"][payload["candidate_order"][-1]][
            "heldout_residual_dependence"]["by_size"],
    }
    print(json.dumps(concise, indent=2))


if __name__ == "__main__":
    main()

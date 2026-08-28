#!/usr/bin/env python3
"""Blind final-tail holdout experiments H2 / H3 / H4.

Test-tail observations are withheld from exponent-set, n_min, model-order,
and ensemble-weight selection. Selection uses training-only rolling-origin
prediction error from finite_size_audit.summarize_model.
"""

from __future__ import annotations

import json
import sys
from dataclasses import asdict
from pathlib import Path

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import finite_size_audit as fsa  # noqa: E402

CSV = ROOT / "data" / "jacobsen_2015_square_site_cylinder.csv"
OUT_DIR = ROOT / "results" / "issue-5"
DPS = 100
MODELS = ((4,), (4, 6), (4, 6, 8), (4, 6, 8, 10), (4, 6, 8, 10, 12))
N_MIN_VALUES = (5, 6, 7, 8, 9, 10)
EXPERIMENTS = (
    {"name": "H2", "holdout": 2, "train_max_n": 19, "test_n": (20, 21)},
    {"name": "H3", "holdout": 3, "train_max_n": 18, "test_n": (19, 20, 21)},
    {"name": "H4", "holdout": 4, "train_max_n": 17, "test_n": (18, 19, 20, 21)},
)


def model_label(powers: tuple[int, ...]) -> str:
    return ",".join(map(str, powers))


def fmt(value: mp.mpf, digits: int = 30) -> str:
    return fsa.fmt(value, digits=digits)


def evaluate_candidate(
    train_obs: list[fsa.Observation],
    powers: tuple[int, ...],
    n_min: int,
    holdout: int,
) -> dict | None:
    try:
        folds, summary = fsa.summarize_model(train_obs, powers, n_min, holdout)
    except ValueError as exc:
        return {
            "model": model_label(powers),
            "n_min": n_min,
            "skipped": True,
            "reason": str(exc),
        }
    subset = [obs for obs in train_obs if obs.n >= n_min]
    coefficients = fsa.fit_linear(subset, powers)
    return {
        "model": model_label(powers),
        "n_min": n_min,
        "skipped": False,
        "folds": summary.folds,
        "median_rmse": summary.median_rmse,
        "worst_rmse": summary.worst_rmse,
        "median_max_abs": summary.median_max_abs,
        "intercept_median": summary.intercept_median,
        "intercept_range": summary.intercept_range,
        "full_fit_intercept": summary.full_fit_intercept,
        "score": summary.score,
        "training_fold_records": [asdict(fold) for fold in folds],
        "coefficients0": fmt(coefficients[0]),
        "_powers": powers,
        "_coefficients": coefficients,
        "_median_rmse": mp.mpf(summary.median_rmse),
        "_score": mp.mpf(summary.score),
    }


def predict_points(
    coefficients: mp.matrix,
    powers: tuple[int, ...],
    test_obs: list[fsa.Observation],
) -> list[dict]:
    rows: list[dict] = []
    errors: list[mp.mpf] = []
    for obs in test_obs:
        predicted = fsa.predict(obs.n, coefficients, powers)
        signed = predicted - obs.value
        errors.append(signed)
        rows.append(
            {
                "n": obs.n,
                "predicted": fmt(predicted),
                "true": fmt(obs.value),
                "signed_error": fmt(signed),
                "absolute_error": fmt(abs(signed)),
            }
        )
    rmse = fsa.rms(errors)
    max_abs = max(abs(error) for error in errors)
    return rows, rmse, max_abs


def ensemble_combine(selected: list[dict], test_obs: list[fsa.Observation]) -> dict:
    """Equal-weight and inverse-RMSE ensembles from training-only scores."""
    intercepts = [mp.mpf(item["full_fit_intercept"]) for item in selected]
    inv_weights = [1 / item["_median_rmse"] for item in selected]
    weight_sum = mp.fsum(inv_weights)
    inv_weights = [weight / weight_sum for weight in inv_weights]
    equal_w = mp.mpf(1) / len(selected)

    equal_mean = mp.fsum(intercepts) / len(intercepts)
    inv_mean = mp.fsum(weight * value for weight, value in zip(inv_weights, intercepts))
    intercepts_sorted = sorted(intercepts)
    mid = len(intercepts_sorted) // 2
    if len(intercepts_sorted) % 2:
        intercept_median = intercepts_sorted[mid]
    else:
        intercept_median = (intercepts_sorted[mid - 1] + intercepts_sorted[mid]) / 2

    equal_rows: list[dict] = []
    inv_rows: list[dict] = []
    equal_errors: list[mp.mpf] = []
    inv_errors: list[mp.mpf] = []
    for obs in test_obs:
        preds = [
            fsa.predict(obs.n, item["_coefficients"], item["_powers"]) for item in selected
        ]
        equal_pred = mp.fsum(preds) / len(preds)
        inv_pred = mp.fsum(weight * value for weight, value in zip(inv_weights, preds))
        equal_err = equal_pred - obs.value
        inv_err = inv_pred - obs.value
        equal_errors.append(equal_err)
        inv_errors.append(inv_err)
        equal_rows.append(
            {
                "n": obs.n,
                "predicted": fmt(equal_pred),
                "true": fmt(obs.value),
                "signed_error": fmt(equal_err),
                "absolute_error": fmt(abs(equal_err)),
            }
        )
        inv_rows.append(
            {
                "n": obs.n,
                "predicted": fmt(inv_pred),
                "true": fmt(obs.value),
                "signed_error": fmt(inv_err),
                "absolute_error": fmt(abs(inv_err)),
            }
        )
    return {
        "n_members": len(selected),
        "members": [
            {
                "model": item["model"],
                "n_min": item["n_min"],
                "training_median_rmse": item["median_rmse"],
                "equal_weight": fmt(equal_w),
                "inverse_rmse_weight": fmt(weight),
                "full_fit_intercept": item["full_fit_intercept"],
            }
            for item, weight in zip(selected, inv_weights)
        ],
        "ensemble_mean_intercept": fmt(equal_mean),
        "ensemble_median_intercept": fmt(intercept_median),
        "ensemble_inverse_rmse_mean_intercept": fmt(inv_mean),
        "minimum_intercept": fmt(min(intercepts)),
        "maximum_intercept": fmt(max(intercepts)),
        "model_spread_exploratory_range": fmt(max(intercepts) - min(intercepts)),
        "note": (
            "min/max intercept is model spread / exploratory range, "
            "not a statistical confidence interval"
        ),
        "equal_weight_predictions": equal_rows,
        "equal_weight_rmse": fmt(fsa.rms(equal_errors)),
        "equal_weight_max_abs": fmt(max(abs(error) for error in equal_errors)),
        "inverse_rmse_predictions": inv_rows,
        "inverse_rmse_rmse": fmt(fsa.rms(inv_errors)),
        "inverse_rmse_max_abs": fmt(max(abs(error) for error in inv_errors)),
    }


def public_candidate(item: dict) -> dict:
    return {key: value for key, value in item.items() if not key.startswith("_")}


def run_experiment(all_obs: list[fsa.Observation], spec: dict) -> dict:
    train_obs = [obs for obs in all_obs if obs.n <= spec["train_max_n"]]
    test_obs = [obs for obs in all_obs if obs.n in spec["test_n"]]
    if len(test_obs) != len(spec["test_n"]):
        raise SystemExit(f"{spec['name']}: missing withheld observations")
    if train_obs and test_obs and train_obs[-1].n >= test_obs[0].n:
        raise SystemExit(f"{spec['name']}: train_max is not strictly below test_min")

    candidates: list[dict] = []
    skipped: list[dict] = []
    for powers in MODELS:
        for n_min in N_MIN_VALUES:
            record = evaluate_candidate(train_obs, powers, n_min, spec["holdout"])
            if record is None or record["skipped"]:
                skipped.append(public_candidate(record))
            else:
                candidates.append(record)

    if not candidates:
        raise SystemExit(f"{spec['name']}: no training-only candidate produced folds")

    candidates.sort(key=lambda item: (item["_median_rmse"], item["_score"], item["n_min"]))
    best = candidates[0]
    point_rows, rmse, max_abs = predict_points(
        best["_coefficients"], best["_powers"], test_obs
    )
    top_k = candidates[: min(5, len(candidates))]
    ensemble = ensemble_combine(top_k, test_obs)

    return {
        "experiment": spec["name"],
        "dps": DPS,
        "holdout": spec["holdout"],
        "train_max_n": spec["train_max_n"],
        "test_n": list(spec["test_n"]),
        "n_training_points": len(train_obs),
        "selection_protocol": (
            "All exponent sets, n_min values, model orders, and ensemble weights "
            "are chosen from training-only rolling-origin prediction error "
            f"(holdout={spec['holdout']}) on n <= {spec['train_max_n']}. "
            "Withheld tail values are read only after the selected configuration "
            "and ensemble are frozen."
        ),
        "n_candidates_evaluated": len(candidates),
        "n_candidates_skipped": len(skipped),
        "skipped": skipped,
        "selected_model": best["model"],
        "selected_n_min": best["n_min"],
        "training_only_score": best["score"],
        "training_only_median_rmse": best["median_rmse"],
        "training_only_worst_rmse": best["worst_rmse"],
        "training_only_folds": best["folds"],
        "full_fit_intercept": best["full_fit_intercept"],
        "predictions": point_rows,
        "rmse": fmt(rmse),
        "maximum_absolute_error": fmt(max_abs),
        "signed_errors": [row["signed_error"] for row in point_rows],
        "absolute_errors": [row["absolute_error"] for row in point_rows],
        "ranking_by_training_median_rmse": [
            {
                "rank": index + 1,
                "model": item["model"],
                "n_min": item["n_min"],
                "training_median_rmse": item["median_rmse"],
                "training_score": item["score"],
                "folds": item["folds"],
                "full_fit_intercept": item["full_fit_intercept"],
            }
            for index, item in enumerate(candidates)
        ],
        "ensemble_top5_training_only": ensemble,
        "all_evaluated_candidates": [public_candidate(item) for item in candidates],
    }


def main() -> int:
    mp.mp.dps = DPS
    observations = fsa.load_observations(CSV)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for spec in EXPERIMENTS:
        payload = run_experiment(observations, spec)
        path = OUT_DIR / f"final_holdout_{spec['name'].lower()}.json"
        path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
        print(
            f"{spec['name']}: selected {payload['selected_model']} "
            f"n_min={payload['selected_n_min']} "
            f"RMSE={mp.nstr(mp.mpf(payload['rmse']), 8)} "
            f"max_abs={mp.nstr(mp.mpf(payload['maximum_absolute_error']), 8)}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

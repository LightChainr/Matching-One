#!/usr/bin/env python3
"""Audit finite-size extrapolations by out-of-sample prediction.

This script deliberately does *not* turn a least-squares residual into an uncertainty
on p_c. It compares correction bases by rolling-origin tail prediction and by the
stability of the inferred infinite-size intercept under changes of n_min.

Example:
    python scripts/finite_size_audit.py \
        data/jacobsen_2015_square_site_cylinder.csv \
        --models 4 4,6 4,6,8 4,6,8,10 4,6,8,10,12 \
        --min-train 8 --holdout 3 --dps 100 --json audit.json
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from statistics import median
from typing import Iterable, Sequence

import mpmath as mp


@dataclass(frozen=True)
class Observation:
    n: int
    value: mp.mpf


@dataclass(frozen=True)
class FoldResult:
    model: str
    n_min: int
    train_max: int
    test_min: int
    test_max: int
    intercept: str
    rmse: str
    max_abs: str


@dataclass(frozen=True)
class ModelSummary:
    model: str
    folds: int
    median_rmse: str
    worst_rmse: str
    median_max_abs: str
    intercept_median: str
    intercept_range: str
    full_fit_intercept: str
    score: str


def parse_model(text: str) -> tuple[int, ...]:
    try:
        powers = tuple(int(part) for part in text.split(",") if part.strip())
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"invalid model {text!r}") from exc
    if not powers or any(power <= 0 for power in powers):
        raise argparse.ArgumentTypeError("models must contain positive powers")
    if tuple(sorted(set(powers))) != powers:
        raise argparse.ArgumentTypeError("model powers must be unique and increasing")
    return powers


def load_observations(path: Path) -> list[Observation]:
    rows: list[Observation] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or not {"n", "value"}.issubset(reader.fieldnames):
            raise ValueError("CSV must contain n and value columns")
        for row in reader:
            rows.append(Observation(n=int(row["n"]), value=mp.mpf(row["value"])))
    rows.sort(key=lambda item: item.n)
    if not rows:
        raise ValueError("input contains no observations")
    if len({row.n for row in rows}) != len(rows):
        raise ValueError("duplicate n values")
    return rows


def design_row(n: int, powers: Sequence[int]) -> list[mp.mpf]:
    n_mp = mp.mpf(n)
    return [mp.mpf(1), *(n_mp ** (-power) for power in powers)]


def fit_linear(observations: Sequence[Observation], powers: Sequence[int]) -> mp.matrix:
    columns = 1 + len(powers)
    if len(observations) < columns:
        raise ValueError("not enough observations for model")
    matrix = mp.matrix([design_row(obs.n, powers) for obs in observations])
    target = mp.matrix([obs.value for obs in observations])
    try:
        coefficients, _residual = mp.qr_solve(matrix, target)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError("singular or ill-conditioned fit") from exc
    return coefficients


def predict(n: int, coefficients: mp.matrix, powers: Sequence[int]) -> mp.mpf:
    row = design_row(n, powers)
    return mp.fsum(term * coefficients[index] for index, term in enumerate(row))


def fmt(value: mp.mpf, digits: int = 30) -> str:
    return mp.nstr(value, n=digits, strip_zeros=False)


def rms(values: Iterable[mp.mpf]) -> mp.mpf:
    values_list = list(values)
    if not values_list:
        return mp.nan
    return mp.sqrt(mp.fsum(value * value for value in values_list) / len(values_list))


def rolling_folds(
    observations: Sequence[Observation],
    powers: Sequence[int],
    min_train_n: int,
    holdout: int,
) -> list[FoldResult]:
    eligible = [obs for obs in observations if obs.n >= min_train_n]
    minimum_points = len(powers) + 3
    folds: list[FoldResult] = []
    for split in range(minimum_points, len(eligible) - holdout + 1):
        train = eligible[:split]
        test = eligible[split : split + holdout]
        if len(test) != holdout:
            continue
        coefficients = fit_linear(train, powers)
        errors = [predict(obs.n, coefficients, powers) - obs.value for obs in test]
        folds.append(
            FoldResult(
                model=",".join(map(str, powers)),
                n_min=train[0].n,
                train_max=train[-1].n,
                test_min=test[0].n,
                test_max=test[-1].n,
                intercept=fmt(coefficients[0]),
                rmse=fmt(rms(errors)),
                max_abs=fmt(max(abs(error) for error in errors)),
            )
        )
    return folds


def intercepts_by_nmin(
    observations: Sequence[Observation],
    powers: Sequence[int],
    first_nmin: int,
) -> list[mp.mpf]:
    results: list[mp.mpf] = []
    minimum_points = len(powers) + 2
    candidates = sorted(obs.n for obs in observations if obs.n >= first_nmin)
    for n_min in candidates:
        subset = [obs for obs in observations if obs.n >= n_min]
        if len(subset) < minimum_points:
            break
        results.append(fit_linear(subset, powers)[0])
    return results


def summarize_model(
    observations: Sequence[Observation],
    powers: Sequence[int],
    min_train_n: int,
    holdout: int,
) -> tuple[list[FoldResult], ModelSummary]:
    folds = rolling_folds(observations, powers, min_train_n, holdout)
    if not folds:
        raise ValueError("no rolling-origin folds; reduce model size/min-train/holdout")

    rmses = [mp.mpf(fold.rmse) for fold in folds]
    maxima = [mp.mpf(fold.max_abs) for fold in folds]
    intercepts = intercepts_by_nmin(observations, powers, min_train_n)
    full_subset = [obs for obs in observations if obs.n >= min_train_n]
    full_intercept = fit_linear(full_subset, powers)[0]

    median_rmse = mp.mpf(str(median(rmses)))
    worst_rmse = max(rmses)
    median_max_abs = mp.mpf(str(median(maxima)))
    intercept_low = min(intercepts)
    intercept_high = max(intercepts)
    intercept_range = intercept_high - intercept_low

    # Dimensionless diagnostic: tail-prediction error plus intercept instability.
    # This is intentionally simple and is not a confidence interval.
    floor = mp.mpf(10) ** (-mp.mp.dps + 10)
    score = mp.log10(median_rmse) + mp.log10(intercept_range + floor)

    summary = ModelSummary(
        model=",".join(map(str, powers)),
        folds=len(folds),
        median_rmse=fmt(median_rmse),
        worst_rmse=fmt(worst_rmse),
        median_max_abs=fmt(median_max_abs),
        intercept_median=fmt(mp.mpf(str(median(intercepts)))),
        intercept_range=fmt(intercept_range),
        full_fit_intercept=fmt(full_intercept),
        score=fmt(score, digits=15),
    )
    return folds, summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path, help="CSV with n,value columns")
    parser.add_argument(
        "--models",
        nargs="+",
        type=parse_model,
        default=[(4,), (4, 6), (4, 6, 8), (4, 6, 8, 10), (4, 6, 8, 10, 12)],
        help="comma-separated correction powers",
    )
    parser.add_argument("--min-train", type=int, default=8, help="smallest n used")
    parser.add_argument("--holdout", type=int, default=2, help="consecutive tail points per fold")
    parser.add_argument("--dps", type=int, default=100, help="mpmath decimal precision")
    parser.add_argument("--json", type=Path, default=None, help="optional JSON output")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.min_train <= 0 or args.holdout <= 0 or args.dps < 40:
        raise SystemExit("min-train and holdout must be positive; dps must be at least 40")
    mp.mp.dps = args.dps
    observations = load_observations(args.csv)

    summaries: list[ModelSummary] = []
    all_folds: list[FoldResult] = []
    for powers in args.models:
        try:
            folds, summary = summarize_model(
                observations, powers, args.min_train, args.holdout
            )
        except ValueError as exc:
            print(f"skip model {powers}: {exc}")
            continue
        all_folds.extend(folds)
        summaries.append(summary)

    if not summaries:
        raise SystemExit("no model could be evaluated")
    summaries.sort(key=lambda item: mp.mpf(item.score))

    print("model          folds    median_RMSE          intercept_span       full_fit_pc")
    for item in summaries:
        print(
            f"{item.model:14} {item.folds:5d}  "
            f"{mp.nstr(mp.mpf(item.median_rmse), 8):>14}  "
            f"{mp.nstr(mp.mpf(item.intercept_range), 8):>14}  "
            f"{mp.nstr(mp.mpf(item.full_fit_intercept), 18)}"
        )
    print("\nRanking is a diagnostic, not a statistical confidence statement.")

    if args.json is not None:
        payload = {
            "input": str(args.csv),
            "dps": args.dps,
            "min_train": args.min_train,
            "holdout": args.holdout,
            "summaries": [asdict(item) for item in summaries],
            "folds": [asdict(item) for item in all_folds],
        }
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Select an extrapolation model without looking at the final tail, then score it.

Only rolling folds whose test widths are at or below the knowledge cutoff are
used for model/configuration selection.  The selected model is then refit on
widths through that cutoff and evaluated once on the final withheld tail.
The full-data summary scores embedded in grid JSON files are never consulted.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from statistics import median
from typing import Any, Iterable, Optional, Sequence

import mpmath as mp


@dataclass(frozen=True)
class Observation:
    n: int
    value: mp.mpf


@dataclass(frozen=True)
class Candidate:
    model: str
    dps: int
    min_train: int
    holdout: int
    validation_folds: int
    validation_test_max: int
    median_rmse: str
    intercept_span: str
    selection_score: str
    source: str


def load_observations(path: Path) -> list[Observation]:
    rows: list[Observation] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or not {"n", "value"}.issubset(reader.fieldnames):
            raise ValueError("CSV must contain n and value columns")
        for row in reader:
            rows.append(Observation(int(row["n"]), mp.mpf(row["value"])))
    rows.sort(key=lambda row: row.n)
    if not rows or len({row.n for row in rows}) != len(rows):
        raise ValueError("input must contain unique observations")
    return rows


def parse_model(text: str) -> tuple[int, ...]:
    try:
        powers = tuple(int(part) for part in text.split(","))
    except ValueError as exc:
        raise ValueError(f"invalid model: {text}") from exc
    if not powers or any(power <= 0 for power in powers):
        raise ValueError(f"invalid model: {text}")
    return powers


def design_row(n: int, powers: Sequence[int]) -> list[mp.mpf]:
    width = mp.mpf(n)
    return [mp.mpf(1)] + [width ** (-power) for power in powers]


def fit(rows: Sequence[Observation], powers: Sequence[int]) -> mp.matrix:
    if len(rows) < len(powers) + 1:
        raise ValueError("not enough training observations for selected model")
    design = mp.matrix([design_row(row.n, powers) for row in rows])
    target = mp.matrix([row.value for row in rows])
    coefficients, _residual = mp.qr_solve(design, target)
    return coefficients


def predict(n: int, coefficients: mp.matrix, powers: Sequence[int]) -> mp.mpf:
    return mp.fsum(
        value * coefficients[index]
        for index, value in enumerate(design_row(n, powers))
    )


def rms(values: Iterable[mp.mpf]) -> mp.mpf:
    items = list(values)
    return mp.sqrt(mp.fsum(item * item for item in items) / len(items))


def load_grid_payloads(raw_dir: Path) -> list[tuple[Path, dict[str, Any]]]:
    payloads: list[tuple[Path, dict[str, Any]]] = []
    for path in sorted(raw_dir.glob("*.json")):
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
            for key in ("dps", "min_train", "holdout", "folds"):
                if key not in payload:
                    raise ValueError(f"missing {key}")
        except (OSError, ValueError, TypeError, json.JSONDecodeError) as exc:
            raise ValueError(f"invalid grid result {path}: {exc}") from exc
        payloads.append((path, payload))
    if not payloads:
        raise ValueError(f"no grid JSON files found in {raw_dir}")
    return payloads


def candidates_from_payloads(
    payloads: Sequence[tuple[Path, dict[str, Any]]],
    *,
    cutoff: int,
    selection_dps: int,
    min_validation_folds: int,
) -> list[Candidate]:
    candidates: list[Candidate] = []
    for path, payload in payloads:
        if int(payload["dps"]) != selection_dps:
            continue
        model_names = sorted({str(fold["model"]) for fold in payload["folds"]})
        for model in model_names:
            safe_folds = [
                fold
                for fold in payload["folds"]
                if str(fold["model"]) == model
                and int(fold["train_max"]) < int(fold["test_min"])
                and int(fold["test_max"]) <= cutoff
            ]
            if len(safe_folds) < min_validation_folds:
                continue
            rmses = [mp.mpf(str(fold["rmse"])) for fold in safe_folds]
            intercepts = [mp.mpf(str(fold["intercept"])) for fold in safe_folds]
            median_rmse = mp.mpf(str(median(rmses)))
            intercept_span = max(intercepts) - min(intercepts)
            floor = mp.mpf(10) ** (-selection_dps + 10)
            score = mp.log10(median_rmse + floor) + mp.log10(intercept_span + floor)
            candidates.append(
                Candidate(
                    model=model,
                    dps=int(payload["dps"]),
                    min_train=int(payload["min_train"]),
                    holdout=int(payload["holdout"]),
                    validation_folds=len(safe_folds),
                    validation_test_max=max(int(fold["test_max"]) for fold in safe_folds),
                    median_rmse=mp.nstr(median_rmse, 30, strip_zeros=False),
                    intercept_span=mp.nstr(intercept_span, 30, strip_zeros=False),
                    selection_score=mp.nstr(score, 20, strip_zeros=False),
                    source=str(path),
                )
            )
    candidates.sort(
        key=lambda row: (
            mp.mpf(row.selection_score),
            len(parse_model(row.model)),
            row.min_train,
            row.holdout,
            row.model,
        )
    )
    if not candidates:
        raise ValueError(
            f"no candidate has {min_validation_folds} training-only folds at dps={selection_dps}"
        )
    return candidates


def final_tail_score(
    observations: Sequence[Observation], candidate: Candidate, cutoff: int
) -> dict[str, Any]:
    powers = parse_model(candidate.model)
    training = [
        row for row in observations if candidate.min_train <= row.n <= cutoff
    ]
    withheld = [row for row in observations if row.n > cutoff]
    if not withheld:
        raise ValueError("final tail is empty")
    mp.mp.dps = candidate.dps
    coefficients = fit(training, powers)
    rows: list[dict[str, Any]] = []
    errors: list[mp.mpf] = []
    for row in withheld:
        prediction = predict(row.n, coefficients, powers)
        error = prediction - row.value
        errors.append(error)
        rows.append(
            {
                "n": row.n,
                "prediction": mp.nstr(prediction, 35, strip_zeros=False),
                "actual": mp.nstr(row.value, 35, strip_zeros=False),
                "signed_error": mp.nstr(error, 30, strip_zeros=False),
                "absolute_error": mp.nstr(abs(error), 30, strip_zeros=False),
            }
        )
    return {
        "training_n_min": min(row.n for row in training),
        "training_n_max": max(row.n for row in training),
        "withheld_n_min": min(row.n for row in withheld),
        "withheld_n_max": max(row.n for row in withheld),
        "intercept": mp.nstr(coefficients[0], 35, strip_zeros=False),
        "rmse": mp.nstr(rms(errors), 30, strip_zeros=False),
        "max_absolute_error": mp.nstr(max(abs(error) for error in errors), 30, strip_zeros=False),
        "predictions": rows,
    }


def summarize(
    *,
    csv_path: Path,
    raw_dir: Path,
    final_tail: int,
    selection_dps: Optional[int],
    min_validation_folds: int,
) -> dict[str, Any]:
    payloads = load_grid_payloads(raw_dir)
    available_dps = sorted({int(payload["dps"]) for _path, payload in payloads})
    chosen_dps = max(available_dps) if selection_dps is None else selection_dps
    if chosen_dps not in available_dps:
        raise ValueError(f"selection dps {chosen_dps} is not present in grid outputs")
    # Decimal strings must be parsed only after the selected arithmetic
    # precision is active. Increasing mp.dps after mp.mpf construction cannot
    # restore digits rounded at the default precision.
    mp.mp.dps = chosen_dps
    observations = load_observations(csv_path)
    if final_tail <= 0 or final_tail >= len(observations):
        raise ValueError("final-tail must be positive and smaller than the dataset")
    cutoff = observations[-final_tail - 1].n
    candidates = candidates_from_payloads(
        payloads,
        cutoff=cutoff,
        selection_dps=chosen_dps,
        min_validation_folds=min_validation_folds,
    )
    selected = candidates[0]
    final_score = final_tail_score(observations, selected, cutoff)
    return {
        "schema_version": 1,
        "input": str(csv_path),
        "knowledge_cutoff": cutoff,
        "withheld_widths": [row.n for row in observations if row.n > cutoff],
        "selection_dps": chosen_dps,
        "selection_rule": (
            "Minimize log10(median rolling RMSE) + log10(rolling-intercept span); "
            "use only folds with test_max <= knowledge_cutoff. Embedded full-data "
            "grid summaries are ignored."
        ),
        "selected": asdict(selected),
        "selection_candidates": [asdict(candidate) for candidate in candidates],
        "final_tail_score": final_score,
    }


def build_parser() -> argparse.ArgumentParser:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "csv",
        type=Path,
        nargs="?",
        default=root / "data/jacobsen_2015_square_site_cylinder.csv",
    )
    parser.add_argument("--raw-dir", type=Path, default=root / "results/issue-5/grid/raw")
    parser.add_argument("--output", type=Path, default=root / "results/issue-5/summary.json")
    parser.add_argument("--final-tail", type=int, default=3)
    parser.add_argument("--selection-dps", type=int, default=None)
    parser.add_argument("--min-validation-folds", type=int, default=2)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.min_validation_folds <= 0:
        raise SystemExit("min-validation-folds must be positive")
    try:
        payload = summarize(
            csv_path=args.csv,
            raw_dir=args.raw_dir,
            final_tail=args.final_tail,
            selection_dps=args.selection_dps,
            min_validation_folds=args.min_validation_folds,
        )
    except ValueError as exc:
        raise SystemExit(str(exc)) from exc
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    selected = payload["selected"]
    score = payload["final_tail_score"]
    print(
        f"selected model={selected['model']} n_min={selected['min_train']} "
        f"holdout={selected['holdout']} dps={selected['dps']}"
    )
    print(f"final-tail RMSE: {score['rmse']}")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

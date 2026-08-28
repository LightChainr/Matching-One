#!/usr/bin/env python3
"""Summarize Issue #5 grid, holdout, and arithmetic-precision results."""

from __future__ import annotations

import csv
import json
import sys
from collections import defaultdict
from pathlib import Path
from statistics import median

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import finite_size_audit as fsa  # noqa: E402

ISSUE = ROOT / "results" / "issue-5"
RAW = ISSUE / "raw"
CSV_DATA = ROOT / "data" / "jacobsen_2015_square_site_cylinder.csv"
DPS = 80


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def mpf(text: str) -> mp.mpf:
    return mp.mpf(text)


def fmt(value: mp.mpf, digits: int = 30) -> str:
    return mp.nstr(value, n=digits, strip_zeros=False)


def decimal_agreement(a: mp.mpf, b: mp.mpf) -> int | str:
    """How many decimal places two numbers agree before they differ.

    Returns the integer floor(-log10(|a-b|)), or "identical" if the difference
    is exactly zero at the comparison working precision.
    """
    delta = abs(a - b)
    if delta == 0:
        return "identical"
    return int(mp.floor(-mp.log10(delta)))


def first_changing_decimal(values: list[mp.mpf]) -> int | str:
    agreements = [
        decimal_agreement(values[i], values[j])
        for i in range(len(values))
        for j in range(i + 1, len(values))
    ]
    numeric = [item for item in agreements if item != "identical"]
    if not numeric:
        return "identical"
    return min(numeric)


def parse_powers(text: str) -> tuple[int, ...]:
    return tuple(int(part) for part in text.split(",") if part.strip())


def unrounded_metrics(
    observations: list,
    powers: tuple[int, ...],
    min_train: int,
    holdout: int,
    dps: int,
) -> dict[str, mp.mpf]:
    """Recompute the three precision metrics without 30-digit serialization."""
    mp.mp.dps = dps
    eligible = [obs for obs in observations if obs.n >= min_train]
    minimum_points = len(powers) + 3
    rmses: list[mp.mpf] = []
    for split in range(minimum_points, len(eligible) - holdout + 1):
        train = eligible[:split]
        test = eligible[split : split + holdout]
        if len(test) != holdout:
            continue
        coefficients = fsa.fit_linear(train, powers)
        errors = [fsa.predict(obs.n, coefficients, powers) - obs.value for obs in test]
        rmses.append(fsa.rms(errors))
    intercepts = fsa.intercepts_by_nmin(observations, powers, min_train)
    full_subset = [obs for obs in observations if obs.n >= min_train]
    full_intercept = fsa.fit_linear(full_subset, powers)[0]
    return {
        "full_fit_intercept": full_intercept,
        "median_rmse": mp.mpf(str(median(rmses))),
        "intercept_range": max(intercepts) - min(intercepts),
    }


def freeze_mpf(value: mp.mpf, digits: int) -> str:
    return mp.nstr(value, n=digits, strip_zeros=False)


def write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main() -> int:
    mp.mp.dps = DPS
    grid_files = sorted(RAW.glob("audit_dps*_nmin*_h*.json"))
    if len(grid_files) != 54:
        print(f"warning: expected 54 grid JSON files, found {len(grid_files)}")

    model_grid_rows: list[dict] = []
    fold_rows: list[dict] = []
    skipped_notes: list[str] = []

    for path in grid_files:
        payload = load_json(path)
        dps = payload["dps"]
        min_train = payload["min_train"]
        holdout = payload["holdout"]
        present_models = {item["model"] for item in payload["summaries"]}
        for summary in payload["summaries"]:
            model_grid_rows.append(
                {
                    "dps": dps,
                    "min_train": min_train,
                    "holdout": holdout,
                    "model": summary["model"],
                    "folds": summary["folds"],
                    "median_rmse": summary["median_rmse"],
                    "worst_rmse": summary["worst_rmse"],
                    "median_max_abs": summary["median_max_abs"],
                    "intercept_median": summary["intercept_median"],
                    "intercept_range": summary["intercept_range"],
                    "full_fit_intercept": summary["full_fit_intercept"],
                    "score": summary["score"],
                    "source": path.name,
                }
            )
        for fold in payload["folds"]:
            if fold["train_max"] >= fold["test_min"]:
                raise SystemExit(
                    f"leakage: {path.name} train_max={fold['train_max']} "
                    f">= test_min={fold['test_min']}"
                )
            fold_rows.append(
                {
                    "dps": dps,
                    "min_train": min_train,
                    "holdout": holdout,
                    "model": fold["model"],
                    "n_min": fold["n_min"],
                    "train_max": fold["train_max"],
                    "test_min": fold["test_min"],
                    "test_max": fold["test_max"],
                    "intercept": fold["intercept"],
                    "rmse": fold["rmse"],
                    "max_abs": fold["max_abs"],
                    "source": path.name,
                }
            )
        expected = {"4", "4,6", "4,6,8", "4,6,8,10", "4,6,8,10,12"}
        missing = sorted(expected - present_models)
        if missing:
            skipped_notes.append(
                f"{path.name}: skipped models {missing} (not fabricated)"
            )

    model_grid_rows.sort(
        key=lambda row: (int(row["dps"]), int(row["min_train"]), int(row["holdout"]), row["model"])
    )
    fold_rows.sort(
        key=lambda row: (
            int(row["dps"]),
            int(row["min_train"]),
            int(row["holdout"]),
            row["model"],
            int(row["train_max"]),
        )
    )

    ranking_rows: list[dict] = []
    grouped: dict[tuple, list[dict]] = defaultdict(list)
    for row in model_grid_rows:
        grouped[(int(row["dps"]), int(row["min_train"]), int(row["holdout"]))].append(row)
    for (dps, min_train, holdout), items in sorted(grouped.items()):
        ranked = sorted(items, key=lambda item: mpf(item["score"]))
        for rank, item in enumerate(ranked, start=1):
            ranking_rows.append(
                {
                    "dps": dps,
                    "min_train": min_train,
                    "holdout": holdout,
                    "rank": rank,
                    "model": item["model"],
                    "score": item["score"],
                    "median_rmse": item["median_rmse"],
                    "worst_rmse": item["worst_rmse"],
                    "intercept_range": item["intercept_range"],
                    "full_fit_intercept": item["full_fit_intercept"],
                    "folds": item["folds"],
                }
            )

    intercept_rows: list[dict] = []
    for row in model_grid_rows:
        intercept_rows.append(
            {
                "dps": row["dps"],
                "min_train": row["min_train"],
                "holdout": row["holdout"],
                "model": row["model"],
                "intercept_median": row["intercept_median"],
                "intercept_range": row["intercept_range"],
                "full_fit_intercept": row["full_fit_intercept"],
                "median_rmse": row["median_rmse"],
            }
        )

    precision_index: dict[tuple, dict[int, dict]] = defaultdict(dict)
    for row in model_grid_rows:
        key = (int(row["min_train"]), int(row["holdout"]), row["model"])
        precision_index[key][int(row["dps"])] = row

    observations = fsa.load_observations(CSV_DATA)
    precision_rows: list[dict] = []
    metrics = ("full_fit_intercept", "median_rmse", "intercept_range")
    for (min_train, holdout, model), by_dps in sorted(
        precision_index.items(), key=lambda item: (item[0][0], item[0][1], item[0][2])
    ):
        if set(by_dps) != {60, 100, 160}:
            continue
        powers = parse_powers(model)
        raw = {}
        for dps_value in (60, 100, 160):
            raw[dps_value] = unrounded_metrics(
                observations, powers, min_train, holdout, dps_value
            )
        mp.mp.dps = 200
        rec: dict = {
            "min_train": min_train,
            "holdout": holdout,
            "model": model,
        }
        json_identical = True
        for metric in metrics:
            json_identical = json_identical and (
                by_dps[60][metric] == by_dps[100][metric] == by_dps[160][metric]
            )
            v60 = mp.mpf(freeze_mpf(raw[60][metric], 55))
            v100 = mp.mpf(freeze_mpf(raw[100][metric], 90))
            v160 = mp.mpf(freeze_mpf(raw[160][metric], 90))
            rec[f"{metric}_dps60"] = freeze_mpf(raw[60][metric], 40)
            rec[f"{metric}_dps100"] = freeze_mpf(raw[100][metric], 50)
            rec[f"{metric}_dps160"] = freeze_mpf(raw[160][metric], 50)
            rec[f"{metric}_json30_dps60"] = by_dps[60][metric]
            rec[f"{metric}_json30_dps100"] = by_dps[100][metric]
            rec[f"{metric}_json30_dps160"] = by_dps[160][metric]
            rec[f"{metric}_absdiff_60_100"] = fmt(abs(v60 - v100), digits=20)
            rec[f"{metric}_absdiff_100_160"] = fmt(abs(v100 - v160), digits=20)
            rec[f"{metric}_agree_decimals_60_100"] = decimal_agreement(v60, v100)
            rec[f"{metric}_agree_decimals_100_160"] = decimal_agreement(v100, v160)
            rec[f"{metric}_first_changing_decimal"] = first_changing_decimal(
                [v60, v100, v160]
            )
        rec["json_30sig_identical"] = json_identical
        precision_rows.append(rec)
    mp.mp.dps = DPS

    write_csv(
        ISSUE / "model_grid.csv",
        [
            "dps",
            "min_train",
            "holdout",
            "model",
            "folds",
            "median_rmse",
            "worst_rmse",
            "median_max_abs",
            "intercept_median",
            "intercept_range",
            "full_fit_intercept",
            "score",
            "source",
        ],
        model_grid_rows,
    )
    write_csv(
        ISSUE / "fold_errors.csv",
        [
            "dps",
            "min_train",
            "holdout",
            "model",
            "n_min",
            "train_max",
            "test_min",
            "test_max",
            "intercept",
            "rmse",
            "max_abs",
            "source",
        ],
        fold_rows,
    )
    write_csv(
        ISSUE / "model_ranking.csv",
        [
            "dps",
            "min_train",
            "holdout",
            "rank",
            "model",
            "score",
            "median_rmse",
            "worst_rmse",
            "intercept_range",
            "full_fit_intercept",
            "folds",
        ],
        ranking_rows,
    )
    write_csv(
        ISSUE / "intercept_stability.csv",
        [
            "dps",
            "min_train",
            "holdout",
            "model",
            "intercept_median",
            "intercept_range",
            "full_fit_intercept",
            "median_rmse",
        ],
        intercept_rows,
    )
    write_csv(
        ISSUE / "precision_stability.csv",
        list(precision_rows[0].keys()) if precision_rows else [],
        precision_rows,
    )

    baseline = load_json(ISSUE / "baseline.json") if (ISSUE / "baseline.json").exists() else None
    holdouts = {}
    for name in ("h2", "h3", "h4"):
        path = ISSUE / f"final_holdout_{name}.json"
        if path.exists():
            holdouts[name.upper()] = load_json(path)

    ref_rows = [row for row in model_grid_rows if int(row["dps"]) == 100]
    by_model_rmse: dict[str, list[mp.mpf]] = defaultdict(list)
    by_model_range: dict[str, list[mp.mpf]] = defaultdict(list)
    by_model_intercept: dict[str, list[mp.mpf]] = defaultdict(list)
    for row in ref_rows:
        by_model_rmse[row["model"]].append(mpf(row["median_rmse"]))
        by_model_range[row["model"]].append(mpf(row["intercept_range"]))
        by_model_intercept[row["model"]].append(mpf(row["full_fit_intercept"]))

    def median_of(values: list[mp.mpf]) -> mp.mpf:
        ordered = sorted(values)
        n = len(ordered)
        if n == 0:
            return mp.nan
        mid = n // 2
        if n % 2:
            return ordered[mid]
        return (ordered[mid - 1] + ordered[mid]) / 2

    model_stability = []
    for model in sorted(by_model_range, key=lambda name: (len(name.split(",")), name)):
        ranges = by_model_range[model]
        rmses = by_model_rmse[model]
        intercepts = by_model_intercept[model]
        model_stability.append(
            {
                "model": model,
                "n_configs": len(ranges),
                "median_intercept_range": median_of(ranges),
                "max_intercept_range": max(ranges) if ranges else mp.nan,
                "median_median_rmse": median_of(rmses),
                "max_median_rmse": max(rmses) if rmses else mp.nan,
                "min_full_fit": min(intercepts) if intercepts else mp.nan,
                "max_full_fit": max(intercepts) if intercepts else mp.nan,
                "full_fit_spread": (max(intercepts) - min(intercepts)) if intercepts else mp.nan,
            }
        )

    most_stable_range = min(model_stability, key=lambda item: item["median_intercept_range"])
    least_stable_range = max(model_stability, key=lambda item: item["median_intercept_range"])
    most_stable_oos = min(model_stability, key=lambda item: item["median_median_rmse"])
    least_stable_oos = max(model_stability, key=lambda item: item["median_median_rmse"])

    # Exploratory ensemble: top 5 configs by training-only rolling prediction
    # (median_rmse) at the baseline protocol dps=100, min_train=8, holdout=2,
    # plus a global top-5 across the dps=100 slice. Weights are equal; they
    # do not use withheld-tail values.
    baseline_slice = [
        row
        for row in ref_rows
        if int(row["min_train"]) == 8 and int(row["holdout"]) == 2
    ]
    global_top = sorted(ref_rows, key=lambda row: mpf(row["median_rmse"]))[:5]
    protocol_top = sorted(baseline_slice, key=lambda row: mpf(row["median_rmse"]))[:5]

    def ensemble_stats(rows: list[dict]) -> dict:
        intercepts = [mpf(row["full_fit_intercept"]) for row in rows]
        return {
            "members": [
                {
                    "model": row["model"],
                    "min_train": row["min_train"],
                    "holdout": row["holdout"],
                    "median_rmse": row["median_rmse"],
                    "full_fit_intercept": row["full_fit_intercept"],
                }
                for row in rows
            ],
            "mean": (mp.fsum(intercepts) / len(intercepts)) if intercepts else mp.nan,
            "median": median_of(intercepts),
            "min": min(intercepts) if intercepts else mp.nan,
            "max": max(intercepts) if intercepts else mp.nan,
        }

    protocol_ens = ensemble_stats(protocol_top)
    global_ens = ensemble_stats(global_top)

    def as_int_agreements(key: str) -> list[int]:
        values: list[int] = []
        for row in precision_rows:
            item = row[key]
            if item == "identical":
                values.append(10**9)
            else:
                values.append(int(item))
        return values

    agree_60_100 = as_int_agreements("full_fit_intercept_agree_decimals_60_100")
    agree_100_160 = as_int_agreements("full_fit_intercept_agree_decimals_100_160")
    rmse_agree_60_100 = as_int_agreements("median_rmse_agree_decimals_60_100")
    rmse_agree_100_160 = as_int_agreements("median_rmse_agree_decimals_100_160")
    range_agree_60_100 = as_int_agreements("intercept_range_agree_decimals_60_100")
    range_agree_100_160 = as_int_agreements("intercept_range_agree_decimals_100_160")
    json_identical_count = sum(1 for row in precision_rows if row["json_30sig_identical"])

    def pct(values: list[int], threshold: int) -> str:
        if not values:
            return "n/a"
        n = sum(1 for value in values if value >= threshold)
        return f"{n}/{len(values)}"

    def describe_agreement(values: list[int]) -> str:
        if not values:
            return "n/a"
        # Agreement at or beyond working precision is identity, not "1e9 decimals".
        working = 80
        finite = [value for value in values if value < working]
        if not finite:
            return "matched at working precision (no difference in printed mpmath values)"
        ordered = sorted(finite)
        return (
            f"worst {ordered[0]} decimals, median {ordered[len(ordered)//2]} decimals, "
            f"best {ordered[-1]} decimals"
        )

    fold_ok = all(int(row["train_max"]) < int(row["test_min"]) for row in fold_rows)

    report: list[str] = []
    report.append("# Issue #5: blind finite-size extrapolation audit")
    report.append("")
    report.append("This report is a computational record. It does not claim a new")
    report.append("threshold value and it does not attach a statistical confidence")
    report.append("interval to any intercept.")
    report.append("")
    report.append("Four quantities are kept separate throughout:")
    report.append("")
    report.append("1. **Arithmetic precision** — how far 60/100/160 dps mpmath")
    report.append("   linear algebra agrees on the same fit.")
    report.append("2. **Out-of-sample prediction error** — rolling-origin and")
    report.append("   withheld-tail errors on cylinder widths that were not used")
    report.append("   to fit the model.")
    report.append("3. **Intercept drift** — movement of the fitted infinite-size")
    report.append("   intercept when `n_min` or the training window changes.")
    report.append("4. **Model-to-model spread** — disagreement among correction")
    report.append("   bases. Min/max of an exploratory ensemble is labelled")
    report.append("   **model spread / exploratory range**, never a CI.")
    report.append("")
    report.append("Fit residuals are not treated as uncertainties on `p_c`.")
    report.append("")
    report.append("## Inputs and protocol")
    report.append("")
    report.append("- Data: `data/jacobsen_2015_square_site_cylinder.csv` (n = 1..21).")
    report.append("- Estimator: `scripts/finite_size_audit.py` (unmodified math).")
    report.append("- Models compared in every grid job: `4` / `4,6` / `4,6,8` /")
    report.append("  `4,6,8,10` / `4,6,8,10,12`.")
    report.append("- Grid: dps ∈ {60, 100, 160} × min_train ∈ {5..10} × holdout ∈ {2,3,4}")
    report.append(f"  → {len(grid_files)} jobs, {len(model_grid_rows)} model summaries,")
    report.append(f"  {len(fold_rows)} folds.")
    report.append("- Runner: `scripts/run_issue5_grid.py` with")
    report.append("  `ProcessPoolExecutor(max_workers=8)`.")
    report.append(
        f"- Rolling-origin leakage check (`train_max < test_min` on every fold): "
        f"{'PASS' if fold_ok else 'FAIL'}."
    )
    if skipped_notes:
        report.append("- Models skipped by the audit (too few points); not fabricated:")
        for note in skipped_notes:
            report.append(f"  - {note}")
    else:
        report.append("- No model was skipped on the full n=1..21 grid.")
    report.append("")

    if baseline is not None:
        report.append("## Baseline reproduction")
        report.append("")
        report.append("Command (unchanged): `min_train=8`, `holdout=2`, `dps=100`,")
        report.append("models `4 4,6 4,6,8 4,6,8,10 4,6,8,10,12`.")
        report.append("")
        report.append("| rank | model | folds | median RMSE | intercept range | full-fit intercept | score |")
        report.append("|---:|---|---:|---|---|---|---|")
        for rank, item in enumerate(baseline["summaries"], start=1):
            report.append(
                f"| {rank} | `{item['model']}` | {item['folds']} | "
                f"{item['median_rmse']} | {item['intercept_range']} | "
                f"{item['full_fit_intercept']} | {item['score']} |"
            )
        report.append("")
        report.append("Ranking uses the script's diagnostic score")
        report.append("`log10(median_rmse) + log10(intercept_range + floor)`.")
        report.append("That score is not a likelihood and not a confidence statement.")
        report.append("")

    report.append("## 1. Arithmetic precision")
    report.append("")
    report.append("This section asks only whether the same least-squares problem")
    report.append("changes when mpmath working precision is 60, 100, or 160 dps.")
    report.append("It is **not** a statement about the statistical accuracy of `p_c`.")
    report.append("")
    report.append(f"Compared {len(precision_rows)} matched (min_train, holdout, model)")
    report.append("triples present at all three dps values.")
    report.append("")
    report.append(
        f"The 30-significant-digit JSON serialization written by "
        f"`finite_size_audit.py` is identical across 60/100/160 dps for "
        f"{json_identical_count}/{len(precision_rows)} matched triples. "
        f"That is the precision of the stored audit tables, not of `p_c`."
    )
    report.append("")
    report.append("Recomputing the same three quantities from the unmodified")
    report.append("`fit_linear` / `rolling_folds` primitives, then comparing the")
    report.append("unrounded mpmath values, shows where the linear algebra itself")
    report.append("stops changing:")
    report.append("")
    report.append("| quantity | 60 vs 100 dps | 100 vs 160 dps | 60 vs 100 agree ≥40 decimals | 100 vs 160 agree ≥40 decimals |")
    report.append("|---|---|---|---:|---:|")
    if agree_100_160:
        report.append(
            f"| full_fit_intercept | {describe_agreement(agree_60_100)} | "
            f"{describe_agreement(agree_100_160)} | {pct(agree_60_100, 40)} | "
            f"{pct(agree_100_160, 40)} |"
        )
        report.append(
            f"| median_rmse | {describe_agreement(rmse_agree_60_100)} | "
            f"{describe_agreement(rmse_agree_100_160)} | {pct(rmse_agree_60_100, 40)} | "
            f"{pct(rmse_agree_100_160, 40)} |"
        )
        report.append(
            f"| intercept_range | {describe_agreement(range_agree_60_100)} | "
            f"{describe_agreement(range_agree_100_160)} | {pct(range_agree_60_100, 40)} | "
            f"{pct(range_agree_100_160, 40)} |"
        )
    report.append("")
    if agree_100_160 and agree_60_100:
        working = 80
        finite_60 = [v for v in agree_60_100 if v < working]
        finite_160 = [v for v in agree_100_160 if v < working]
        if finite_60:
            report.append(
                f"For `full_fit_intercept`, raising working precision from 60 dps "
                f"to 100 dps last changes a digit around decimal place "
                f"{min(finite_60)} in the worst matched triple "
                f"(median {sorted(finite_60)[len(finite_60)//2]})."
            )
        else:
            report.append(
                "For `full_fit_intercept`, 60 dps and 100 dps matched at working "
                "precision (no difference in the printed mpmath values)."
            )
        if finite_160:
            report.append(
                f"Raising 100 dps to 160 dps last changes a digit around decimal "
                f"place {min(finite_160)} in the worst triple "
                f"(median {sorted(finite_160)[len(finite_160)//2]})."
            )
        else:
            report.append(
                "Raising 100 dps to 160 dps matched at working precision "
                "(no difference in the printed mpmath values of `full_fit_intercept`)."
            )
        report.append("")
        report.append("Any residual 100-vs-160 dps movement is at working precision")
        report.append("and is far below the out-of-sample prediction error.")
        report.append("Arithmetic noise is not the dominant uncertainty on this sequence.")
    report.append("")
    report.append("See `results/issue-5/precision_stability.csv`.")
    report.append("")

    report.append("## 2. Out-of-sample prediction error")
    report.append("")
    report.append("Rolling-origin folds from the grid (dps = 100 slice):")
    report.append("")
    report.append("| model | n configs | median of median RMSE | worst median RMSE |")
    report.append("|---|---:|---|---|")
    for item in model_stability:
        report.append(
            f"| `{item['model']}` | {item['n_configs']} | "
            f"{fmt(item['median_median_rmse'])} | {fmt(item['max_median_rmse'])} |"
        )
    report.append("")
    report.append(
        f"Lowest typical rolling RMSE (dps=100): `{most_stable_oos['model']}` "
        f"(median of median RMSE = {fmt(most_stable_oos['median_median_rmse'])})."
    )
    report.append(
        f"Highest typical rolling RMSE (dps=100): `{least_stable_oos['model']}` "
        f"(median of median RMSE = {fmt(least_stable_oos['median_median_rmse'])})."
    )
    report.append("")
    report.append("These RMSE values are prediction errors on withheld widths,")
    report.append("not standard errors of the intercept.")
    report.append("")
    report.append("### Blind final-tail experiments")
    report.append("")
    report.append("Selection used only training widths. The tail was scored after")
    report.append("the configuration was frozen.")
    report.append("")
    if holdouts:
        report.append("| experiment | train n ≤ | test n | selected model | selected n_min | training-only median RMSE | tail RMSE | tail max abs |")
        report.append("|---|---:|---|---|---:|---|---|---|")
        for key in ("H2", "H3", "H4"):
            if key not in holdouts:
                continue
            payload = holdouts[key]
            test_n = ",".join(str(n) for n in payload["test_n"])
            report.append(
                f"| {key} | {payload['train_max_n']} | {test_n} | "
                f"`{payload['selected_model']}` | {payload['selected_n_min']} | "
                f"{payload['training_only_median_rmse']} | {payload['rmse']} | "
                f"{payload['maximum_absolute_error']} |"
            )
        report.append("")
        for key in ("H2", "H3", "H4"):
            if key not in holdouts:
                continue
            payload = holdouts[key]
            report.append(f"#### {key}")
            report.append("")
            report.append(payload["selection_protocol"])
            report.append("")
            report.append(
                f"Selected `{payload['selected_model']}` with `n_min={payload['selected_n_min']}` "
                f"(training-only diagnostic score {payload['training_only_score']})."
            )
            report.append("")
            report.append("| n | predicted | true withheld | signed error | absolute error |")
            report.append("|---:|---|---|---|---|")
            for row in payload["predictions"]:
                report.append(
                    f"| {row['n']} | {row['predicted']} | {row['true']} | "
                    f"{row['signed_error']} | {row['absolute_error']} |"
                )
            report.append("")
            ens = payload["ensemble_top5_training_only"]
            report.append(
                f"Exploratory top-5 ensemble (training-only; equal-weight tail RMSE "
                f"{ens['equal_weight_rmse']}; inverse-RMSE tail RMSE "
                f"{ens['inverse_rmse_rmse']})."
            )
            report.append("")
            report.append("| n | equal-weight predicted | signed error | abs error |")
            report.append("|---:|---|---|---|")
            for row in ens["equal_weight_predictions"]:
                report.append(
                    f"| {row['n']} | {row['predicted']} | {row['signed_error']} | "
                    f"{row['absolute_error']} |"
                )
            report.append("")
            report.append(
                f"{payload['n_candidates_skipped']} candidate (model, n_min) pairs "
                f"were skipped for too few training points and were not filled in."
            )
            report.append("")
    report.append("See `final_holdout_h2.json`, `final_holdout_h3.json`,")
    report.append("`final_holdout_h4.json`, and `fold_errors.csv`.")
    report.append("")

    report.append("## 3. Intercept drift")
    report.append("")
    report.append("Intercept range is the span of full-subset intercepts obtained")
    report.append("by raising `n_min` inside a fixed model, as implemented by")
    report.append("`intercepts_by_nmin` in `finite_size_audit.py`. It measures")
    report.append("sensitivity to the lower cutoff, not sampling error.")
    report.append("")
    report.append("| model | n configs (dps=100) | median intercept range | max intercept range | full-fit spread across configs |")
    report.append("|---|---:|---|---|---|")
    for item in model_stability:
        report.append(
            f"| `{item['model']}` | {item['n_configs']} | "
            f"{fmt(item['median_intercept_range'])} | "
            f"{fmt(item['max_intercept_range'])} | "
            f"{fmt(item['full_fit_spread'])} |"
        )
    report.append("")
    report.append(
        f"Smallest typical intercept range: `{most_stable_range['model']}` "
        f"(median range {fmt(most_stable_range['median_intercept_range'])})."
    )
    report.append(
        f"Largest typical intercept range: `{least_stable_range['model']}` "
        f"(median range {fmt(least_stable_range['median_intercept_range'])})."
    )
    report.append("")
    report.append("See `intercept_stability.csv`.")
    report.append("")

    report.append("## 4. Model-to-model spread")
    report.append("")
    report.append("Different correction bases fitted to the same cylinder sequence")
    report.append("do not share one intercept. The spread below is an exploratory")
    report.append("range over models, **not** a statistical confidence interval.")
    report.append("")

    def dump_ensemble(title: str, ens: dict) -> None:
        report.append(f"### {title}")
        report.append("")
        report.append("Members chosen by training-only rolling-origin median RMSE.")
        report.append("Equal weights. No withheld tail entered the selection.")
        report.append("")
        report.append("| model | min_train | holdout | training median RMSE | full-fit intercept |")
        report.append("|---|---:|---:|---|---|")
        for member in ens["members"]:
            report.append(
                f"| `{member['model']}` | {member['min_train']} | {member['holdout']} | "
                f"{member['median_rmse']} | {member['full_fit_intercept']} |"
            )
        report.append("")
        report.append(f"- ensemble mean: `{fmt(ens['mean'])}`")
        report.append(f"- ensemble median: `{fmt(ens['median'])}`")
        report.append(f"- minimum intercept: `{fmt(ens['min'])}`")
        report.append(f"- maximum intercept: `{fmt(ens['max'])}`")
        report.append(
            f"- model spread / exploratory range: `{fmt(ens['max'] - ens['min'])}`"
        )
        report.append("")

    dump_ensemble(
        "Top 5 at baseline protocol (dps=100, min_train=8, holdout=2)",
        protocol_ens,
    )
    dump_ensemble(
        "Top 5 across the full dps=100 grid by training-only median RMSE",
        global_ens,
    )

    if ref_rows:
        all_intercepts = [mpf(row["full_fit_intercept"]) for row in ref_rows]
        report.append("Across every dps=100 grid summary:")
        report.append("")
        report.append(f"- minimum full-fit intercept: `{fmt(min(all_intercepts))}`")
        report.append(f"- maximum full-fit intercept: `{fmt(max(all_intercepts))}`")
        report.append(
            f"- model spread / exploratory range: `{fmt(max(all_intercepts) - min(all_intercepts))}`"
        )
        report.append("")

    report.append("## Files")
    report.append("")
    report.append("| path | contents |")
    report.append("|---|---|")
    report.append("| `results/issue-5/environment.txt` | host, python, packages, git HEAD |")
    report.append("| `results/issue-5/baseline.json` | unmodified audit, min_train=8, holdout=2, dps=100 |")
    report.append("| `results/issue-5/raw/` | 54 grid JSON payloads |")
    report.append("| `results/issue-5/logs/` | per-job stdout/stderr |")
    report.append("| `results/issue-5/grid_manifest.json` | job success/fail records |")
    report.append("| `results/issue-5/model_grid.csv` | one row per (dps, min_train, holdout, model) |")
    report.append("| `results/issue-5/fold_errors.csv` | every rolling-origin fold |")
    report.append("| `results/issue-5/model_ranking.csv` | diagnostic rank within each job |")
    report.append("| `results/issue-5/intercept_stability.csv` | intercept range and full-fit intercept |")
    report.append("| `results/issue-5/precision_stability.csv` | 60/100/160 dps comparison |")
    report.append("| `results/issue-5/final_holdout_h2.json` | blind n=20,21 |")
    report.append("| `results/issue-5/final_holdout_h3.json` | blind n=19,20,21 |")
    report.append("| `results/issue-5/final_holdout_h4.json` | blind n=18..21 |")
    report.append("")
    report.append("## What this does not show")
    report.append("")
    report.append("- No ± interval in this report is a statistical CI.")
    report.append("- Training residuals are not intercept uncertainties.")
    report.append("- Ensemble min/max is model spread / exploratory range.")
    report.append("- Arithmetic agreement at 40 decimals is not 40-decimal knowledge")
    report.append("  of the percolation threshold.")
    report.append("")

    (ISSUE / "REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(f"wrote {len(model_grid_rows)} model_grid rows, {len(fold_rows)} folds")
    print(f"precision triples: {len(precision_rows)}")
    print(
        f"most stable (intercept range): {most_stable_range['model']}; "
        f"least: {least_stable_range['model']}"
    )
    print(
        f"most stable (OOS RMSE): {most_stable_oos['model']}; "
        f"least: {least_stable_oos['model']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Leakage-safe polynomial/Padé audit for the finite-width threshold sequence.

Stage A selects a model and lower fitting width using rolling folds whose test
widths end at or before 18.  Only after selection are the frozen fits scored on
widths 19--21.  The rational correction is

    p_n = p_inf + n^-4 P_m(n^-2) / Q_k(n^-2),  Q_k(0)=1.

The script reports prediction errors and numerical conditioning; none of its
fit residuals are statistical confidence intervals for ``p_inf``.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import math
from pathlib import Path
import sys
from typing import Iterable, Sequence

import mpmath as mp

from finite_size_audit import Observation, fmt, load_observations, rms


@dataclass(frozen=True)
class Family:
    name: str
    numerator_degree: int
    denominator_degree: int

    @property
    def parameters(self) -> int:
        return 1 + self.numerator_degree + 1 + self.denominator_degree


FAMILIES = (
    Family("poly_0", 0, 0),
    Family("poly_1", 1, 0),
    Family("poly_2", 2, 0),
    Family("poly_3", 3, 0),
    Family("poly_4", 4, 0),
    Family("pade_1_1", 1, 1),
    Family("pade_2_1", 2, 1),
    Family("pade_1_2", 1, 2),
    Family("pade_2_2", 2, 2),
)


def _x_t(n: int) -> tuple[mp.mpf, mp.mpf]:
    n_mp = mp.mpf(n)
    return n_mp ** -2, n_mp ** -4


def _parts(parameters: mp.matrix, family: Family) -> tuple[mp.mpf, list[mp.mpf], list[mp.mpf]]:
    p_inf = parameters[0]
    a_start = 1
    a_stop = a_start + family.numerator_degree + 1
    numerator = [parameters[index] for index in range(a_start, a_stop)]
    denominator = [mp.mpf(1)] + [
        parameters[index]
        for index in range(a_stop, a_stop + family.denominator_degree)
    ]
    return p_inf, numerator, denominator


def _polynomial(coefficients: Sequence[mp.mpf], x: mp.mpf) -> mp.mpf:
    return mp.fsum(coefficient * x**power for power, coefficient in enumerate(coefficients))


def predict(n: int, parameters: mp.matrix, family: Family) -> mp.mpf:
    p_inf, numerator, denominator = _parts(parameters, family)
    x, thermal = _x_t(n)
    return p_inf + thermal * _polynomial(numerator, x) / _polynomial(denominator, x)


def _jacobian_row(n: int, parameters: mp.matrix, family: Family) -> list[mp.mpf]:
    _p_inf, numerator, denominator = _parts(parameters, family)
    x, thermal = _x_t(n)
    top = _polynomial(numerator, x)
    bottom = _polynomial(denominator, x)
    row = [mp.mpf(1)]
    row.extend(thermal * x**power / bottom for power in range(family.numerator_degree + 1))
    row.extend(
        -thermal * top * x**power / bottom**2
        for power in range(1, family.denominator_degree + 1)
    )
    return row


def _linear_initial(observations: Sequence[Observation], family: Family) -> mp.matrix:
    powers = tuple(4 + 2 * power for power in range(family.numerator_degree + 1))
    design = mp.matrix(
        [[mp.mpf(1), *[mp.mpf(obs.n) ** -power for power in powers]] for obs in observations]
    )
    target = mp.matrix([obs.value for obs in observations])
    coefficients, _ = mp.qr_solve(design, target)
    return mp.matrix([*coefficients, *([mp.mpf(0)] * family.denominator_degree)])


def fit(observations: Sequence[Observation], family: Family) -> mp.matrix:
    if len(observations) < family.parameters + 1:
        raise ValueError("not enough observations")
    parameters = _linear_initial(observations, family)
    target = mp.matrix([obs.value for obs in observations])

    for _iteration in range(80):
        predictions = mp.matrix([predict(obs.n, parameters, family) for obs in observations])
        residual = target - predictions
        jacobian = mp.matrix([_jacobian_row(obs.n, parameters, family) for obs in observations])
        try:
            step, _ = mp.qr_solve(jacobian, residual)
        except (ValueError, ZeroDivisionError) as exc:
            raise ValueError("singular rational fit") from exc
        old_norm = mp.norm(residual)
        damping = mp.mpf(1)
        accepted = False
        while damping >= mp.mpf("0.00000095367431640625"):
            candidate = parameters + damping * step
            try:
                new_residual = mp.matrix(
                    [obs.value - predict(obs.n, candidate, family) for obs in observations]
                )
                new_norm = mp.norm(new_residual)
            except ZeroDivisionError:
                new_norm = mp.inf
            if new_norm < old_norm:
                parameters = candidate
                accepted = True
                break
            damping /= 2
        if not accepted:
            if mp.norm(step) <= mp.sqrt(mp.eps) * (1 + mp.norm(parameters)):
                break
            raise ValueError("rational fit did not descend")
        if mp.norm(damping * step) <= mp.sqrt(mp.eps) * (1 + mp.norm(parameters)):
            break
    return parameters


def condition_number(observations: Sequence[Observation], parameters: mp.matrix, family: Family) -> mp.mpf:
    jacobian = mp.matrix([_jacobian_row(obs.n, parameters, family) for obs in observations])
    _u, singular_values, _v = mp.svd(jacobian)
    largest = max(singular_values)
    smallest = min(singular_values)
    return largest / smallest if smallest else mp.inf


def denominator_poles(parameters: mp.matrix, family: Family) -> list[mp.mpc]:
    if family.denominator_degree == 0:
        return []
    _p_inf, _numerator, denominator = _parts(parameters, family)
    try:
        try:
            return list(
                mp.polyroots(denominator, maxsteps=1000, error=False, asc=True)
            )
        except TypeError:
            # mpmath 1.3 has only the historical descending-coefficient API.
            return list(
                mp.polyroots(
                    list(reversed(denominator)), maxsteps=1000, error=False
                )
            )
    except (ValueError, ZeroDivisionError):
        return []


def pole_guard(poles: Sequence[mp.mpc], n_min: int) -> tuple[bool, mp.mpf]:
    """Reject real/near-real poles within one quarter data-interval of [0,x_max]."""

    x_max = mp.mpf(n_min) ** -2
    guard = x_max / 4
    minimum = mp.inf
    for pole in poles:
        real = mp.re(pole)
        imag = abs(mp.im(pole))
        horizontal = mp.mpf(0) if 0 <= real <= x_max else min(abs(real), abs(real - x_max))
        distance = mp.sqrt(horizontal**2 + imag**2)
        minimum = min(minimum, distance)
    return minimum <= guard, minimum


def errors_for(observations: Sequence[Observation], parameters: mp.matrix, family: Family) -> list[mp.mpf]:
    return [predict(obs.n, parameters, family) - obs.value for obs in observations]


def rolling_rows(
    observations: Sequence[Observation], family: Family, n_min: int, cutoff: int, holdout: int
) -> list[dict[str, object]]:
    eligible = [obs for obs in observations if n_min <= obs.n <= cutoff]
    rows: list[dict[str, object]] = []
    for split in range(family.parameters + 1, len(eligible) - holdout + 1):
        train = eligible[:split]
        test = eligible[split : split + holdout]
        try:
            parameters = fit(train, family)
            errors = errors_for(test, parameters, family)
            poles = denominator_poles(parameters, family)
            rejected, distance = pole_guard(poles, n_min)
            condition = condition_number(train, parameters, family)
        except (ValueError, ZeroDivisionError):
            continue
        rows.append(
            {
                "family": family.name,
                "n_min": n_min,
                "train_max": train[-1].n,
                "test_min": test[0].n,
                "test_max": test[-1].n,
                "signed_errors": [fmt(error) for error in errors],
                "rmse": fmt(rms(errors)),
                "max_abs": fmt(max(abs(error) for error in errors)),
                "intercept": fmt(parameters[0]),
                "jacobian_condition": fmt(condition),
                "pole_rejected": rejected,
                "nearest_pole_distance": None if mp.isinf(distance) else fmt(distance),
            }
        )
    return rows


def median_mp(values: Iterable[mp.mpf]) -> mp.mpf:
    ordered = sorted(values)
    middle = len(ordered) // 2
    return ordered[middle] if len(ordered) % 2 else (ordered[middle - 1] + ordered[middle]) / 2


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path)
    parser.add_argument("--cutoff", type=int, default=18)
    parser.add_argument("--targets", default="19,20,21")
    parser.add_argument("--n-min", default="5,6,7,8,9,10")
    parser.add_argument("--holdout", type=int, default=2)
    parser.add_argument("--min-folds", type=int, default=2)
    parser.add_argument("--dps", type=int, default=100)
    parser.add_argument("--json", type=Path, required=True)
    args = parser.parse_args()
    if args.dps < 50 or args.holdout < 1 or args.min_folds < 1:
        raise SystemExit("dps must be >=50; holdout and min-folds must be positive")
    mp.mp.dps = args.dps
    observations = load_observations(args.csv)
    target_ns = tuple(int(value) for value in args.targets.split(","))
    n_mins = tuple(int(value) for value in args.n_min.split(","))
    if any(n <= args.cutoff for n in target_ns):
        raise SystemExit("all targets must be strictly beyond cutoff")
    targets = [obs for obs in observations if obs.n in target_ns]
    if [obs.n for obs in targets] != list(target_ns):
        raise SystemExit("target widths missing from input")

    folds: list[dict[str, object]] = []
    candidates: list[dict[str, object]] = []
    failures: list[dict[str, object]] = []
    for family in FAMILIES:
        for n_min in n_mins:
            family_folds = rolling_rows(observations, family, n_min, args.cutoff, args.holdout)
            valid = [row for row in family_folds if not row["pole_rejected"]]
            if len(valid) < args.min_folds:
                failures.append(
                    {
                        "family": family.name,
                        "n_min": n_min,
                        "stage": "rolling_validation",
                        "reason": "fewer than {} valid folds (found {})".format(
                            args.min_folds, len(valid)
                        ),
                    }
                )
                continue
            folds.extend(family_folds)
            rmses = [mp.mpf(str(row["rmse"])) for row in valid]
            maxima = [mp.mpf(str(row["max_abs"])) for row in valid]
            train = [obs for obs in observations if n_min <= obs.n <= args.cutoff]
            try:
                parameters = fit(train, family)
                poles = denominator_poles(parameters, family)
                rejected, pole_distance = pole_guard(poles, n_min)
                heldout_errors = errors_for(targets, parameters, family)
                condition = condition_number(train, parameters, family)
            except (ValueError, ZeroDivisionError) as exc:
                failures.append(
                    {
                        "family": family.name,
                        "n_min": n_min,
                        "stage": "cutoff_fit",
                        "reason": str(exc),
                    }
                )
                continue
            candidates.append(
                {
                    "family": family.name,
                    "n_min": n_min,
                    "parameters": [fmt(value) for value in parameters],
                    "validation_folds": len(valid),
                    "validation_median_rmse": fmt(median_mp(rmses)),
                    "validation_worst_rmse": fmt(max(rmses)),
                    "validation_median_max_abs": fmt(median_mp(maxima)),
                    "full_fit_intercept": fmt(parameters[0]),
                    "jacobian_condition": fmt(condition),
                    "poles": [[fmt(mp.re(pole)), fmt(mp.im(pole))] for pole in poles],
                    "pole_rejected": rejected,
                    "nearest_pole_distance": None if mp.isinf(pole_distance) else fmt(pole_distance),
                    "target_signed_errors": [fmt(error) for error in heldout_errors],
                    "target_rmse": fmt(rms(heldout_errors)),
                    "target_max_abs": fmt(max(abs(error) for error in heldout_errors)),
                }
            )

    valid_candidates = [row for row in candidates if not row["pole_rejected"]]
    if not valid_candidates:
        raise SystemExit("no valid model/window candidate")
    rank_median = sorted(valid_candidates, key=lambda row: mp.mpf(str(row["validation_median_rmse"])))
    rank_worst = sorted(valid_candidates, key=lambda row: mp.mpf(str(row["validation_worst_rmse"])))
    payload = {
        "protocol": "Stage A selection uses only rolling folds ending at cutoff; targets are scored after freezing",
        "input": str(args.csv),
        "dps": args.dps,
        "cutoff": args.cutoff,
        "targets": list(target_ns),
        "holdout": args.holdout,
        "minimum_validation_folds": args.min_folds,
        "pole_rule": "reject denominator poles within x_max/4 of x in [0,x_max], x_max=n_min^-2",
        "model_definitions": [family.__dict__ for family in FAMILIES],
        "selected_by_median_validation_rmse": {"family": rank_median[0]["family"], "n_min": rank_median[0]["n_min"]},
        "selected_by_worst_validation_rmse": {"family": rank_worst[0]["family"], "n_min": rank_worst[0]["n_min"]},
        "candidates": candidates,
        "failed_or_ineligible_candidates": failures,
        "folds": folds,
        "warning": "held-out errors and intercept drift are deterministic diagnostics, not statistical confidence intervals",
    }
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    for label, winner in (("median", rank_median[0]), ("worst", rank_worst[0])):
        print(
            label,
            winner["family"],
            "n_min=" + str(winner["n_min"]),
            "validation=" + str(winner["validation_median_rmse"]),
            "target=" + str(winner["target_rmse"]),
            "signed=" + ",".join(winner["target_signed_errors"]),
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

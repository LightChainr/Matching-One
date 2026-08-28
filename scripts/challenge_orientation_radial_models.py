#!/usr/bin/env python3
"""Frozen P32 radial/logarithmic challenge for same-N orientation effects.

Input is one row per size and independent seed.  Required logical fields are
``N``, ``seed``, ``delta_M``, ``delta_M_se`` and ``delta_cos4``; documented
aliases are accepted below.  An optional covariance edge list can replace the
default diagonal covariance derived from the reported standard errors.

The training sizes (65, 85, 130), held-out sizes (145, 170), exponent models,
free-alpha bounds, and power-correction candidates are frozen in this file.
Missing frozen sizes produce a machine-readable ``NOT_READY`` result rather
than a partial fit.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from fractions import Fraction
import json
import math
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple


TRAINING_SIZES = (65, 85, 130)
HELDOUT_SIZES = (145, 170)
FIXED_ALPHA = 13.0 / 8.0
OMEGA_CANDIDATES = (1, 2, 3)
FREE_ALPHA_BOUNDS = (0.0, 4.0)

Vector = List[float]
Matrix = List[List[float]]


@dataclass(frozen=True)
class Observation:
    row_id: str
    n: int
    seed: str
    delta_m: float
    se: float
    delta_cos4: float
    delta_cos8: float
    delta_cos4_exact: str
    delta_cos8_exact: str


def _harmonics(a: int, b: int) -> Tuple[Fraction, Fraction]:
    n = a * a + b * b
    cos4 = Fraction(a**4 - 6 * a * a * b * b + b**4, n * n)
    return cos4, 2 * cos4 * cos4 - 1


def _fraction_text(value: Fraction) -> str:
    return "{}/{}".format(value.numerator, value.denominator)


@dataclass
class Fit:
    model: str
    parameters: Dict[str, float]
    train_indices: List[int]
    train_chi_square: float
    train_degrees_of_freedom: int
    condition_number: float
    parameter_covariance: Matrix
    influence: Matrix
    predict: Callable[[Observation], float]
    jacobian: Callable[[Observation], Vector]
    details: Dict[str, object]


def _field(raw: Dict[str, str], names: Sequence[str]) -> str:
    for name in names:
        if name in raw and raw[name] != "":
            return raw[name]
    raise ValueError("missing input field; expected one of " + ", ".join(names))


def read_observations(path: Path) -> List[Observation]:
    output: List[Observation] = []
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for index, raw in enumerate(reader):
            n = int(_field(raw, ("N", "n")))
            seed = _field(raw, ("seed", "rng_seed"))
            delta_m = float(_field(raw, ("delta_M", "deltaM", "delta_m")))
            se = float(_field(raw, ("delta_M_se", "deltaM_se", "delta_m_se", "se")))
            delta_cos4 = float(_field(raw, (
                "delta_cos4", "delta_cos4theta", "delta_cos4_first_minus_second",
            )))
            if raw.get("delta_cos8") not in (None, ""):
                delta_cos8 = float(raw["delta_cos8"])
                delta_cos8_exact = raw.get("delta_cos8_exact") or str(raw["delta_cos8"])
                delta_cos4_exact = raw.get("delta_cos4_exact") or str(delta_cos4)
            elif all(raw.get(name) not in (None, "") for name in ("a1", "b1", "a2", "b2")):
                first4, first8 = _harmonics(int(raw["a1"]), int(raw["b1"]))
                second4, second8 = _harmonics(int(raw["a2"]), int(raw["b2"]))
                exact4 = first4 - second4
                exact8 = first8 - second8
                delta_cos8 = float(exact8)
                delta_cos4_exact = _fraction_text(exact4)
                delta_cos8_exact = _fraction_text(exact8)
                if not math.isclose(delta_cos4, float(exact4), rel_tol=1e-12, abs_tol=1e-15):
                    raise ValueError("delta_cos4 disagrees with the declared Gaussian pair")
            else:
                raise ValueError("delta_cos8 or all four Gaussian representation fields are required")
            row_id = raw.get("row_id") or "{}:{}:{}".format(n, seed, index)
            if n <= 0 or not all(
                math.isfinite(value) for value in (delta_m, se, delta_cos4, delta_cos8)
            ):
                raise ValueError("non-finite or invalid observation in row " + str(index + 2))
            if se <= 0 or delta_cos4 == 0:
                raise ValueError("SE must be positive and delta_cos4 nonzero")
            output.append(Observation(
                row_id, n, seed, delta_m, se, delta_cos4, delta_cos8,
                delta_cos4_exact, delta_cos8_exact,
            ))
    if not output:
        raise ValueError("input CSV contains no observations")
    if len({item.row_id for item in output}) != len(output):
        raise ValueError("row_id values must be unique")
    return output


def read_covariance(path: Optional[Path], observations: Sequence[Observation]) -> Tuple[Matrix, str]:
    covariance = [[0.0] * len(observations) for _ in observations]
    for index, item in enumerate(observations):
        covariance[index][index] = item.se * item.se
    if path is None:
        return covariance, "diagonal_from_delta_M_se"
    indices = {item.row_id: index for index, item in enumerate(observations)}
    seen: Dict[Tuple[int, int], float] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"row_id_i", "row_id_j", "covariance"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("covariance CSV requires row_id_i,row_id_j,covariance")
        for raw in reader:
            if raw["row_id_i"] not in indices or raw["row_id_j"] not in indices:
                raise ValueError("covariance row references an unknown observation")
            i = indices[raw["row_id_i"]]
            j = indices[raw["row_id_j"]]
            value = float(raw["covariance"])
            if not math.isfinite(value):
                raise ValueError("covariance entries must be finite")
            key = (min(i, j), max(i, j))
            if key in seen and not math.isclose(seen[key], value, rel_tol=1e-12, abs_tol=1e-300):
                raise ValueError("inconsistent duplicate covariance entry")
            seen[key] = value
            covariance[i][j] = covariance[j][i] = value
    for index in range(len(observations)):
        if covariance[index][index] <= 0:
            raise ValueError("covariance diagonal must be positive")
    _require_positive_definite(covariance)
    return covariance, "full_edge_list:" + str(path)


def _solve(matrix: Matrix, vector: Vector) -> Vector:
    n = len(vector)
    if len(matrix) != n or any(len(row) != n for row in matrix):
        raise ValueError("linear system must be square")
    augmented = [list(map(float, matrix[i])) + [float(vector[i])] for i in range(n)]
    scale = max((abs(value) for row in matrix for value in row), default=0.0)
    tolerance = max(scale * 1e-13, 1e-300)
    for column in range(n):
        pivot = max(range(column, n), key=lambda row: abs(augmented[row][column]))
        if abs(augmented[pivot][column]) <= tolerance:
            raise ArithmeticError("singular or ill-resolved matrix")
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        divisor = augmented[column][column]
        for entry in range(column, n + 1):
            augmented[column][entry] /= divisor
        for row in range(n):
            if row == column:
                continue
            factor = augmented[row][column]
            for entry in range(column, n + 1):
                augmented[row][entry] -= factor * augmented[column][entry]
    return [augmented[row][-1] for row in range(n)]


def _inverse(matrix: Matrix) -> Matrix:
    n = len(matrix)
    columns = [_solve(matrix, [1.0 if row == column else 0.0 for row in range(n)])
               for column in range(n)]
    return [[columns[column][row] for column in range(n)] for row in range(n)]


def _require_positive_definite(matrix: Matrix) -> None:
    """Validate a covariance matrix by an unpivoted Cholesky factorization."""

    factor = [[0.0] * len(matrix) for _ in matrix]
    for i in range(len(matrix)):
        for j in range(i + 1):
            value = matrix[i][j] - math.fsum(
                factor[i][k] * factor[j][k] for k in range(j)
            )
            if i == j:
                if value <= 0 or not math.isfinite(value):
                    raise ValueError("covariance matrix must be positive definite")
                factor[i][j] = math.sqrt(value)
            else:
                factor[i][j] = value / factor[j][j]


def _transpose(matrix: Matrix) -> Matrix:
    return [list(column) for column in zip(*matrix)]


def _matmul(first: Matrix, second: Matrix) -> Matrix:
    second_t = _transpose(second)
    return [[math.fsum(a * b for a, b in zip(row, column)) for column in second_t]
            for row in first]


def _matvec(matrix: Matrix, vector: Vector) -> Vector:
    return [math.fsum(a * b for a, b in zip(row, vector)) for row in matrix]


def _subtract(first: Matrix, second: Matrix) -> Matrix:
    return [[first[i][j] - second[i][j] for j in range(len(first[i]))]
            for i in range(len(first))]


def _add(first: Matrix, second: Matrix) -> Matrix:
    return [[first[i][j] + second[i][j] for j in range(len(first[i]))]
            for i in range(len(first))]


def _quadratic(vector: Vector, inverse_covariance: Matrix) -> float:
    return math.fsum(a * b for a, b in zip(vector, _matvec(inverse_covariance, vector)))


def _subset(matrix: Matrix, rows: Sequence[int], columns: Optional[Sequence[int]] = None) -> Matrix:
    if columns is None:
        columns = rows
    return [[matrix[i][j] for j in columns] for i in rows]


def _condition_number(matrix: Matrix) -> float:
    inverse = _inverse(matrix)
    norm = max(math.fsum(abs(value) for value in row) for row in matrix)
    inverse_norm = max(math.fsum(abs(value) for value in row) for row in inverse)
    return norm * inverse_norm


def _gls_linear(
    model: str,
    observations: Sequence[Observation],
    covariance: Matrix,
    indices: Sequence[int],
    basis: Callable[[Observation], Vector],
    parameter_names: Sequence[str],
    details: Optional[Dict[str, object]] = None,
) -> Fit:
    selected = list(indices)
    design = [basis(observations[index]) for index in selected]
    y = [observations[index].delta_m for index in selected]
    cov = _subset(covariance, selected)
    cov_inv = _inverse(cov)
    design_t = _transpose(design)
    normal = _matmul(_matmul(design_t, cov_inv), design)
    normal_inv = _inverse(normal)
    rhs = _matvec(_matmul(design_t, cov_inv), y)
    beta = _matvec(normal_inv, rhs)
    residual = [value - math.fsum(x * coefficient for x, coefficient in zip(row, beta))
                for value, row in zip(y, design)]
    influence = _matmul(_matmul(normal_inv, design_t), cov_inv)
    parameters = {name: value for name, value in zip(parameter_names, beta)}
    return Fit(
        model=model,
        parameters=parameters,
        train_indices=selected,
        train_chi_square=_quadratic(residual, cov_inv),
        train_degrees_of_freedom=len(selected) - len(beta),
        condition_number=_condition_number(normal),
        parameter_covariance=normal_inv,
        influence=influence,
        predict=lambda item: math.fsum(x * coefficient for x, coefficient in zip(basis(item), beta)),
        jacobian=basis,
        details=dict(details or {}),
    )


def _base(item: Observation, alpha: float = FIXED_ALPHA) -> float:
    return item.delta_cos4 * item.n ** (-alpha)


def fit_fixed(observations: Sequence[Observation], covariance: Matrix, indices: Sequence[int]) -> Fit:
    return _gls_linear(
        "fixed_13_8", observations, covariance, indices,
        lambda item: [_base(item)], ("A",), {"alpha": FIXED_ALPHA},
    )


def fit_power(
    observations: Sequence[Observation], covariance: Matrix, indices: Sequence[int], omega: int
) -> Fit:
    return _gls_linear(
        "fixed_13_8_power_correction", observations, covariance, indices,
        lambda item: [_base(item), _base(item) * item.n ** (-omega / 2.0)],
        ("A", "B"), {"alpha": FIXED_ALPHA, "omega": omega},
    )


def fit_log(observations: Sequence[Observation], covariance: Matrix, indices: Sequence[int]) -> Fit:
    return _gls_linear(
        "fixed_13_8_log", observations, covariance, indices,
        lambda item: [_base(item), _base(item) * math.log(item.n)],
        ("A", "B"), {"alpha": FIXED_ALPHA},
    )


def fit_h4_h8(
    observations: Sequence[Observation], covariance: Matrix, indices: Sequence[int]
) -> Fit:
    return _gls_linear(
        "fixed_13_8_h4_h8", observations, covariance, indices,
        lambda item: [
            item.delta_cos4 * item.n ** (-FIXED_ALPHA),
            item.delta_cos8 * item.n ** (-FIXED_ALPHA),
        ],
        ("A4", "A8"), {"alpha": FIXED_ALPHA},
    )


def fit_free_alpha(
    observations: Sequence[Observation], covariance: Matrix, indices: Sequence[int]
) -> Fit:
    selected = list(indices)

    def at_alpha(alpha: float) -> Fit:
        return _gls_linear(
            "free_alpha_profile", observations, covariance, selected,
            lambda item: [_base(item, alpha)], ("A",), {"alpha": alpha},
        )

    left, right = FREE_ALPHA_BOUNDS
    ratio = (math.sqrt(5.0) - 1.0) / 2.0
    x1 = right - ratio * (right - left)
    x2 = left + ratio * (right - left)
    f1 = at_alpha(x1).train_chi_square
    f2 = at_alpha(x2).train_chi_square
    for _ in range(100):
        if f1 <= f2:
            right, x2, f2 = x2, x1, f1
            x1 = right - ratio * (right - left)
            f1 = at_alpha(x1).train_chi_square
        else:
            left, x1, f1 = x1, x2, f2
            x2 = left + ratio * (right - left)
            f2 = at_alpha(x2).train_chi_square
    alpha = (left + right) / 2.0
    profiled = at_alpha(alpha)
    amplitude = profiled.parameters["A"]
    train = [observations[index] for index in selected]
    design = [[_base(item, alpha), -amplitude * _base(item, alpha) * math.log(item.n)]
              for item in train]
    cov = _subset(covariance, selected)
    cov_inv = _inverse(cov)
    design_t = _transpose(design)
    normal = _matmul(_matmul(design_t, cov_inv), design)
    parameter_covariance = _inverse(normal)
    influence = _matmul(_matmul(parameter_covariance, design_t), cov_inv)

    def predict(item: Observation) -> float:
        return amplitude * _base(item, alpha)

    def jacobian(item: Observation) -> Vector:
        base = _base(item, alpha)
        return [base, -amplitude * base * math.log(item.n)]

    return Fit(
        model="free_alpha",
        parameters={"A": amplitude, "alpha": alpha},
        train_indices=selected,
        train_chi_square=profiled.train_chi_square,
        train_degrees_of_freedom=len(selected) - 2,
        condition_number=_condition_number(normal),
        parameter_covariance=parameter_covariance,
        influence=influence,
        predict=predict,
        jacobian=jacobian,
        details={"alpha_bounds": list(FREE_ALPHA_BOUNDS)},
    )


def residual_covariance(
    fit: Fit, observations: Sequence[Observation], covariance: Matrix, test_indices: Sequence[int]
) -> Matrix:
    test = list(test_indices)
    train = fit.train_indices
    test_cov = _subset(covariance, test)
    cross = _subset(covariance, test, train)
    jacobian = [fit.jacobian(observations[index]) for index in test]
    prediction_cov = _matmul(_matmul(jacobian, fit.parameter_covariance), _transpose(jacobian))
    cross_term = _matmul(_matmul(cross, _transpose(fit.influence)), _transpose(jacobian))
    return _subtract(_add(test_cov, prediction_cov),
                     _add(cross_term, _transpose(cross_term)))


def score(fit: Fit, observations: Sequence[Observation], covariance: Matrix,
          indices: Sequence[int]) -> Tuple[Dict[str, object], List[Dict[str, object]]]:
    selected = list(indices)
    errors = [observations[index].delta_m - fit.predict(observations[index]) for index in selected]
    residual_cov = residual_covariance(fit, observations, covariance, selected)
    chi_square = _quadratic(errors, _inverse(residual_cov))
    rows: List[Dict[str, object]] = []
    for position, index in enumerate(selected):
        item = observations[index]
        error_se = math.sqrt(max(residual_cov[position][position], 0.0))
        predicted = fit.predict(item)
        rows.append({
            "model": fit.model,
            "row_id": item.row_id,
            "N": item.n,
            "seed": item.seed,
            "observed_delta_M": item.delta_m,
            "predicted_delta_M": predicted,
            "signed_error_observed_minus_predicted": item.delta_m - predicted,
            "prediction_residual_se": error_se,
            "standardized_signed_error": (item.delta_m - predicted) / error_se if error_se else None,
        })
    return {
        "chi_square": chi_square,
        "degrees_of_freedom": len(selected),
        "residual_covariance": residual_cov,
    }, rows


def fit_models(observations: Sequence[Observation], covariance: Matrix,
               training_indices: Sequence[int]) -> Tuple[List[Fit], Dict[str, object]]:
    power_candidates = [fit_power(observations, covariance, training_indices, omega)
                        for omega in OMEGA_CANDIDATES]
    selected_power = min(power_candidates, key=lambda fit: (fit.train_chi_square, fit.details["omega"]))
    selection = {
        "criterion": "minimum covariance-aware training chi-square; ties choose smaller omega",
        "candidates": [
            {"omega": fit.details["omega"], "training_chi_square": fit.train_chi_square,
             "condition_number": fit.condition_number}
            for fit in power_candidates
        ],
        "selected_omega": selected_power.details["omega"],
    }
    return [
        fit_fixed(observations, covariance, training_indices),
        selected_power,
        fit_log(observations, covariance, training_indices),
        fit_h4_h8(observations, covariance, training_indices),
        fit_free_alpha(observations, covariance, training_indices),
    ], selection


def amplitude_rows(observations: Sequence[Observation], covariance: Matrix,
                   fixed_amplitude: float) -> List[Dict[str, object]]:
    transformed = [item.delta_m * item.n ** FIXED_ALPHA / item.delta_cos4 for item in observations]
    factors = [item.n ** FIXED_ALPHA / item.delta_cos4 for item in observations]
    transformed_cov = [[covariance[i][j] * factors[i] * factors[j]
                        for j in range(len(observations))] for i in range(len(observations))]
    output: List[Dict[str, object]] = []
    for n in sorted({item.n for item in observations}):
        indices = [index for index, item in enumerate(observations) if item.n == n]
        cov = _subset(transformed_cov, indices)
        inverse = _inverse(cov)
        ones = [1.0] * len(indices)
        inverse_ones = _matvec(inverse, ones)
        normalizer = math.fsum(inverse_ones)
        weights = [value / normalizer for value in inverse_ones]
        estimate = math.fsum(weights[pos] * transformed[index]
                             for pos, index in enumerate(indices))
        se = math.sqrt(1.0 / normalizer)
        output.append({
            "N": n,
            "seed_count": len(indices),
            "scaled_amplitude": estimate,
            "scaled_amplitude_se": se,
            "drift_from_fixed_training_A": estimate - fixed_amplitude,
            "drift_z": (estimate - fixed_amplitude) / se,
        })
    return output


def write_csv(path: Path, rows: Sequence[Dict[str, object]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fields), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def fit_summary(fit: Fit) -> Dict[str, object]:
    return {
        "model": fit.model,
        "parameters": fit.parameters,
        "details": fit.details,
        "training_chi_square": fit.train_chi_square,
        "training_degrees_of_freedom": fit.train_degrees_of_freedom,
        "condition_number": fit.condition_number,
        "parameter_covariance": fit.parameter_covariance,
    }


def run(args: argparse.Namespace) -> int:
    observations = read_observations(args.input)
    present = sorted({item.n for item in observations})
    required = set(TRAINING_SIZES + HELDOUT_SIZES)
    missing = sorted(required - set(present))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = args.output_dir / "challenge.json"
    if missing:
        payload = {
            "status": "NOT_READY",
            "reason": "frozen P32 sizes are missing; no model was fit",
            "training_sizes": list(TRAINING_SIZES),
            "heldout_sizes": list(HELDOUT_SIZES),
            "present_sizes": present,
            "missing_sizes": missing,
            "input": str(args.input),
        }
        summary_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
        print("NOT_READY: missing frozen sizes " + ",".join(map(str, missing)))
        print("wrote " + str(summary_path))
        return 0

    covariance, covariance_source = read_covariance(args.covariance, observations)
    training = [index for index, item in enumerate(observations) if item.n in TRAINING_SIZES]
    heldout = [index for index, item in enumerate(observations) if item.n in HELDOUT_SIZES]
    fits, power_selection = fit_models(observations, covariance, training)
    heldout_covariance = _subset(covariance, heldout)
    heldout_zero_errors = [observations[index].delta_m for index in heldout]
    zero_effect_chi_square = _quadratic(
        heldout_zero_errors, _inverse(heldout_covariance)
    )
    zero_effect = {
        "model": "zero_effect",
        "heldout_chi_square": zero_effect_chi_square,
        "heldout_degrees_of_freedom": len(heldout),
    }
    heldout_rows: List[Dict[str, object]] = []
    model_results: List[Dict[str, object]] = []
    for fit in fits:
        heldout_score, rows = score(fit, observations, covariance, heldout)
        heldout_rows.extend(rows)
        result = fit_summary(fit)
        result["heldout"] = heldout_score
        result["heldout"]["chi_square_improvement_over_zero"] = (
            zero_effect_chi_square / heldout_score["chi_square"]
            if heldout_score["chi_square"] > 0 else None
        )
        model_results.append(result)

    loso_rows: List[Dict[str, object]] = []
    for omitted_n in TRAINING_SIZES:
        fold_training = [index for index in training if observations[index].n != omitted_n]
        fold_test = [index for index in training if observations[index].n == omitted_n]
        fold_power_candidates = [
            fit_power(observations, covariance, fold_training, omega)
            for omega in OMEGA_CANDIDATES
        ]
        fold_power = min(
            fold_power_candidates,
            key=lambda fit: (fit.train_chi_square, fit.details["omega"]),
        )
        fold_fits = [
            fit_fixed(observations, covariance, fold_training),
            fold_power,
            fit_log(observations, covariance, fold_training),
            fit_h4_h8(observations, covariance, fold_training),
            fit_free_alpha(observations, covariance, fold_training),
        ]
        for fit in fold_fits:
            _, rows = score(fit, observations, covariance, fold_test)
            for row in rows:
                row["omitted_training_size"] = omitted_n
                row["power_omega_selected_within_fold"] = (
                    fit.details["omega"] if fit.model == "fixed_13_8_power_correction" else ""
                )
                loso_rows.append(row)

    fixed_amplitude = fits[0].parameters["A"]
    amplitudes = amplitude_rows(observations, covariance, fixed_amplitude)
    angular_design = []
    seen_design_sizes = set()
    for item in observations:
        if item.n in seen_design_sizes:
            continue
        seen_design_sizes.add(item.n)
        angular_design.append({
            "N": item.n,
            "delta_cos4": item.delta_cos4,
            "delta_cos8": item.delta_cos8,
            "delta_cos4_exact": item.delta_cos4_exact,
            "delta_cos8_exact": item.delta_cos8_exact,
        })
    payload = {
        "status": "READY",
        "protocol": {
            "training_sizes": list(TRAINING_SIZES),
            "heldout_sizes": list(HELDOUT_SIZES),
            "fixed_alpha": FIXED_ALPHA,
            "power_omega_candidates": list(OMEGA_CANDIDATES),
            "free_alpha_bounds": list(FREE_ALPHA_BOUNDS),
            "orientation_order": "first_minus_second",
            "selection_uses_heldout": False,
        },
        "input": str(args.input),
        "covariance_source": covariance_source,
        "observation_count": len(observations),
        "angular_design": angular_design,
        "power_correction_selection": power_selection,
        "zero_effect_benchmark": zero_effect,
        "models": model_results,
        "caveats": [
            "held-out chi-square includes the supplied observation covariance and linearized training-fit uncertainty",
            "free-alpha uncertainty and held-out covariance use a local Jacobian approximation",
            "LOSO is a training-set diagnostic and does not reveal or select on held-out sizes",
            "model uncertainty is not folded into Monte Carlo standard errors",
        ],
    }
    summary_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    write_csv(
        args.output_dir / "heldout_signed_errors.csv", heldout_rows,
        ("model", "row_id", "N", "seed", "observed_delta_M", "predicted_delta_M",
         "signed_error_observed_minus_predicted", "prediction_residual_se",
         "standardized_signed_error"),
    )
    write_csv(
        args.output_dir / "amplitude_drift.csv", amplitudes,
        ("N", "seed_count", "scaled_amplitude", "scaled_amplitude_se",
         "drift_from_fixed_training_A", "drift_z"),
    )
    write_csv(
        args.output_dir / "loso_predictions.csv", loso_rows,
        ("model", "omitted_training_size", "row_id", "N", "seed", "observed_delta_M",
         "predicted_delta_M", "signed_error_observed_minus_predicted",
         "prediction_residual_se", "standardized_signed_error",
         "power_omega_selected_within_fold"),
    )
    print("wrote " + str(summary_path))
    print("wrote heldout_signed_errors.csv, amplitude_drift.csv, loso_predictions.csv")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="deltaM_by_size_seed CSV")
    parser.add_argument("--covariance", type=Path, help="optional row-id covariance edge-list CSV")
    parser.add_argument("--output-dir", type=Path, required=True)
    return run(parser.parse_args())


if __name__ == "__main__":
    raise SystemExit(main())

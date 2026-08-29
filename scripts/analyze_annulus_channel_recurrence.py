#!/usr/bin/env python3
"""Covariance-aware radial recurrence scorer for Issue #253.

The scorer treats the two norm-5 tori as two readouts of one scalar radial
channel.  Radii 2,4,8 define equal log steps.  Radius 7 is kept out of the
recurrence construction and used only as a correlated, post-reveal design
check.  No row in this file is a new evidence block beyond PR #247.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Sequence, Tuple


CHANNELS = ("A_plus", "A_minus")
GEOMETRIES = (325, 425)
SOURCE_RADII = (2, 4, 8)
TARGET_RADIUS = 7
X_RANGE = (-8.0, 8.0)
GAP_RANGE = (0.05, 16.0)
THETA_RANGE = (0.05, math.pi - 0.05)


def solve(matrix: Sequence[Sequence[float]], rhs: Sequence[float]) -> List[float]:
    size = len(rhs)
    a = [list(map(float, matrix[row])) + [float(rhs[row])] for row in range(size)]
    for column in range(size):
        pivot = max(range(column, size), key=lambda row: abs(a[row][column]))
        if abs(a[pivot][column]) < 1e-22:
            raise ArithmeticError("singular matrix")
        a[column], a[pivot] = a[pivot], a[column]
        scale = a[column][column]
        a[column] = [value / scale for value in a[column]]
        for row in range(size):
            if row == column:
                continue
            factor = a[row][column]
            a[row] = [a[row][j] - factor * a[column][j] for j in range(size + 1)]
    return [a[row][-1] for row in range(size)]


def inverse(matrix: Sequence[Sequence[float]]) -> List[List[float]]:
    size = len(matrix)
    columns = [solve(matrix, [1.0 if row == col else 0.0 for row in range(size)])
               for col in range(size)]
    return [[columns[col][row] for col in range(size)] for row in range(size)]


def matvec(matrix: Sequence[Sequence[float]], vector: Sequence[float]) -> List[float]:
    return [sum(a * b for a, b in zip(row, vector)) for row in matrix]


def quadratic(vector: Sequence[float], precision: Sequence[Sequence[float]]) -> float:
    return sum(a * b for a, b in zip(vector, matvec(precision, vector)))


def chi_square_survival(value: float, degrees: int) -> float:
    """Regularized upper incomplete gamma Q(df/2, value/2)."""
    if degrees <= 0:
        raise ValueError("positive degrees of freedom required")
    a, x = degrees / 2.0, value / 2.0
    if x <= 0.0:
        return 1.0
    epsilon, tiny = 3e-14, 1e-300
    if x < a + 1.0:
        term = 1.0 / a
        total = term
        ap = a
        for _ in range(10000):
            ap += 1.0
            term *= x / ap
            total += term
            if abs(term) < abs(total) * epsilon:
                break
        lower = total * math.exp(-x + a * math.log(x) - math.lgamma(a))
        return max(0.0, min(1.0, 1.0 - lower))
    b = x + 1.0 - a
    c = 1.0 / tiny
    d = 1.0 / b
    h = d
    for i in range(1, 10000):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < tiny:
            d = tiny
        c = b + an / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < epsilon:
            break
    return math.exp(-x + a * math.log(x) - math.lgamma(a)) * h


def subvector(order: Sequence[str], point: Sequence[float], labels: Sequence[str]) -> List[float]:
    return [float(point[order.index(label)]) for label in labels]


def submatrix(order: Sequence[str], matrix: Sequence[Sequence[float]],
              labels: Sequence[str]) -> List[List[float]]:
    indices = [order.index(label) for label in labels]
    return [[float(matrix[i][j]) for j in indices] for i in indices]


def labels_for(channels: Sequence[str], radii: Sequence[int]) -> List[str]:
    return [f"N{n}_R{radius}_Delta_{channel}"
            for channel in channels for n in GEOMETRIES for radius in radii]


def basis(model: str, radius: int, parameters: Sequence[float]) -> List[float]:
    step = math.log(radius / 2.0, 2.0)
    if model == "R1":
        eigenvalue = 2.0 ** (-parameters[0])
        return [eigenvalue ** step]
    if model == "J2":
        eigenvalue = 2.0 ** (-parameters[0])
        return [eigenvalue ** step, step * eigenvalue ** step]
    if model == "R2":
        center, gap = parameters
        first = 2.0 ** (-(center - gap / 2.0))
        second = 2.0 ** (-(center + gap / 2.0))
        return [first ** step, second ** step]
    if model == "R2_gap1":
        center = parameters[0]
        first = 2.0 ** (-(center - 0.5))
        second = 2.0 ** (-(center + 0.5))
        return [first ** step, second ** step]
    if model == "C2":
        decay, theta = parameters
        radius_eigenvalue = 2.0 ** (-decay)
        return [radius_eigenvalue ** step * math.cos(theta * step),
                radius_eigenvalue ** step * math.sin(theta * step)]
    raise ValueError(model)


def gls_profile(point: Sequence[float], covariance: Sequence[Sequence[float]],
                channels: Sequence[str], model: str,
                parameters: Sequence[float]) -> Dict[str, object]:
    modes = len(basis(model, 2, parameters))
    outputs = len(channels) * len(GEOMETRIES)
    design = []
    for output in range(outputs):
        for radius in SOURCE_RADII:
            row = [0.0] * (outputs * modes)
            row[output * modes:(output + 1) * modes] = basis(model, radius, parameters)
            design.append(row)
    precision = inverse(covariance)
    normal = [[sum(design[i][a] * precision[i][j] * design[j][b]
                   for i in range(len(point)) for j in range(len(point)))
               for b in range(outputs * modes)] for a in range(outputs * modes)]
    rhs = [sum(design[i][a] * precision[i][j] * point[j]
               for i in range(len(point)) for j in range(len(point)))
           for a in range(outputs * modes)]
    coefficients = solve(normal, rhs)
    fitted = [sum(value * coefficient for value, coefficient in zip(row, coefficients))
              for row in design]
    residual = [observed - expected for observed, expected in zip(point, fitted)]
    return {
        "chi_square": quadratic(residual, precision),
        "coefficients": coefficients,
        "fitted": fitted,
        "residual": residual,
    }


def golden_minimize(function: Callable[[float], float], low: float, high: float,
                    iterations: int = 80) -> Tuple[float, float]:
    ratio = (math.sqrt(5.0) - 1.0) / 2.0
    left = high - ratio * (high - low)
    right = low + ratio * (high - low)
    f_left, f_right = function(left), function(right)
    for _ in range(iterations):
        if f_left <= f_right:
            high, right, f_right = right, left, f_left
            left = high - ratio * (high - low)
            f_left = function(left)
        else:
            low, left, f_left = left, right, f_right
            right = low + ratio * (high - low)
            f_right = function(right)
    value = (low + high) / 2.0
    return value, function(value)


def grid_minimize_2d(function: Callable[[float, float], float],
                     first_range: Tuple[float, float], second_range: Tuple[float, float],
                     rounds: int = 4, width: int = 25) -> Tuple[List[float], float]:
    first_low, first_high = first_range
    second_low, second_high = second_range
    best = [0.0, 0.0]
    best_value = float("inf")
    for _ in range(rounds):
        first_step = (first_high - first_low) / (width - 1)
        second_step = (second_high - second_low) / (width - 1)
        for i in range(width):
            first = first_low + i * first_step
            for j in range(width):
                second = second_low + j * second_step
                value = function(first, second)
                if value < best_value:
                    best, best_value = [first, second], value
        first_low = max(first_range[0], best[0] - first_step)
        first_high = min(first_range[1], best[0] + first_step)
        second_low = max(second_range[0], best[1] - second_step)
        second_high = min(second_range[1], best[1] + second_step)
    return best, best_value


def profile_model(point: Sequence[float], covariance: Sequence[Sequence[float]],
                  channels: Sequence[str], model: str) -> Dict[str, object]:
    objective = lambda parameters: float(
        gls_profile(point, covariance, channels, model, parameters)["chi_square"])
    if model in ("R1", "J2", "R2_gap1"):
        value, _ = golden_minimize(lambda x: objective([x]), *X_RANGE)
        parameters = [value]
        ranges = [X_RANGE]
    elif model == "R2":
        parameters, _ = grid_minimize_2d(
            lambda center, gap: objective([center, gap]), X_RANGE, GAP_RANGE)
        ranges = [X_RANGE, GAP_RANGE]
    elif model == "C2":
        parameters, _ = grid_minimize_2d(
            lambda decay, theta: objective([decay, theta]), X_RANGE, THETA_RANGE)
        ranges = [X_RANGE, THETA_RANGE]
    else:
        raise ValueError(model)
    result = gls_profile(point, covariance, channels, model, parameters)
    nonlinear_count = len(parameters)
    linear_count = len(result["coefficients"])
    degrees = len(point) - nonlinear_count - linear_count
    boundary = [min(abs(value - low), abs(high - value)) < 1e-3 * (high - low)
                for value, (low, high) in zip(parameters, ranges)]
    return {
        "model": model,
        "parameters": parameters,
        "parameter_names": {
            "R1": ["x"], "J2": ["x"], "R2": ["center_x", "gap_x"],
            "R2_gap1": ["center_x_with_gap_fixed_to_1"],
            "C2": ["decay_x", "theta_per_log2_step"],
        }[model],
        "parameter_ranges": ranges,
        "boundary_hit": boundary,
        "chi_square": result["chi_square"],
        "degrees_of_freedom": degrees,
        "chi_square_survival": chi_square_survival(result["chi_square"], degrees)
                               if degrees > 0 else None,
        "linear_parameter_count": linear_count,
        "nonlinear_parameter_count": nonlinear_count,
    }


def recurrence(source: Sequence[float]) -> Dict[str, float]:
    if len(source) != 6:
        raise ValueError("source must be two geometry rows at radii 2,4,8")
    first = source[:3]
    second = source[3:]
    matrix = [[first[1], -first[0]], [second[1], -second[0]]]
    trace, determinant = solve(matrix, [first[2], second[2]])
    discriminant = trace * trace - 4.0 * determinant
    return {"T": trace, "D": determinant, "Delta": discriminant}


def continuous_prediction(source: Sequence[float], radial_class: str,
                          target_radius: float = TARGET_RADIUS) -> List[float]:
    rec = recurrence(source)
    trace, determinant, discriminant = rec["T"], rec["D"], rec["Delta"]
    step = math.log(target_radius / 2.0, 2.0)
    rows = (source[:3], source[3:])
    predictions = []
    if radial_class == "R2":
        if discriminant <= 0.0:
            raise ValueError("R2 requires positive discriminant")
        root = math.sqrt(discriminant)
        first, second = (trace + root) / 2.0, (trace - root) / 2.0
        if first <= 0.0 or second <= 0.0:
            raise ValueError("continuous real radial roots must be positive")
        for g0, g1, _ in rows:
            a = (g1 - second * g0) / (first - second)
            b = (first * g0 - g1) / (first - second)
            predictions.append(a * first ** step + b * second ** step)
    elif radial_class == "C2":
        if discriminant >= 0.0 or determinant <= 0.0:
            raise ValueError("C2 requires D>0 and negative discriminant")
        modulus = math.sqrt(determinant)
        theta = math.acos(trace / (2.0 * modulus))
        for g0, g1, _ in rows:
            sine_amplitude = (g1 / modulus - g0 * math.cos(theta)) / math.sin(theta)
            predictions.append(modulus ** step * (
                g0 * math.cos(theta * step) + sine_amplitude * math.sin(theta * step)))
    elif radial_class == "J2":
        if trace <= 0.0:
            raise ValueError("J2 requires a positive repeated eigenvalue")
        eigenvalue = trace / 2.0
        for g0, g1, _ in rows:
            slope = g1 / eigenvalue - g0
            predictions.append(eigenvalue ** step * (g0 + slope * step))
    else:
        raise ValueError(radial_class)
    return predictions


def finite_jacobian(function: Callable[[Sequence[float]], Sequence[float]],
                    point: Sequence[float], covariance: Sequence[Sequence[float]]) -> List[List[float]]:
    base = list(function(point))
    jacobian = [[0.0] * len(point) for _ in base]
    for column in range(len(point)):
        scale = max(math.sqrt(max(covariance[column][column], 0.0)) * 1e-3,
                    abs(point[column]) * 1e-6, 1e-9)
        plus, minus = list(point), list(point)
        plus[column] += scale
        minus[column] -= scale
        right, left = list(function(plus)), list(function(minus))
        for row in range(len(base)):
            jacobian[row][column] = (right[row] - left[row]) / (2.0 * scale)
    return jacobian


def propagated_covariance(jacobian: Sequence[Sequence[float]],
                          covariance: Sequence[Sequence[float]]) -> List[List[float]]:
    return [[sum(jacobian[i][a] * covariance[a][b] * jacobian[j][b]
                 for a in range(len(covariance)) for b in range(len(covariance)))
             for j in range(len(jacobian))] for i in range(len(jacobian))]


def channel_recurrence(order: Sequence[str], point: Sequence[float],
                       covariance: Sequence[Sequence[float]], channel: str) -> Dict[str, object]:
    all_radii = SOURCE_RADII + (TARGET_RADIUS,)
    labels = labels_for([channel], all_radii)
    values = subvector(order, point, labels)
    matrix = submatrix(order, covariance, labels)
    # Reorder from geometry-major [2,4,8,7] to source [N325...,N425...] + targets.
    source = values[0:3] + values[4:7]
    target = [values[3], values[7]]
    source_indices = [0, 1, 2, 4, 5, 6]
    source_covariance = [[matrix[i][j] for j in source_indices] for i in source_indices]
    rec = recurrence(source)
    jac = finite_jacobian(lambda row: list(recurrence(row).values()), source, source_covariance)
    rec_cov = propagated_covariance(jac, source_covariance)
    delta_se = math.sqrt(max(rec_cov[2][2], 0.0))
    if rec["D"] > 0.0 and rec["Delta"] < 0.0:
        point_class = "C2_principal_branch"
        prediction_class = "C2"
    elif rec["Delta"] > 0.0:
        roots = [(rec["T"] + math.sqrt(rec["Delta"])) / 2.0,
                 (rec["T"] - math.sqrt(rec["Delta"])) / 2.0]
        if min(roots) > 0.0:
            point_class = "R2_positive_roots"
            prediction_class = "R2"
        else:
            point_class = "real_roots_not_both_positive"
            prediction_class = None
    else:
        point_class = "unclassified"
        prediction_class = None
    prediction_record = None
    if prediction_class is not None:
        prediction = continuous_prediction(source, prediction_class)
        residual = [observed - expected for observed, expected in zip(target, prediction)]

        def residual_function(row: Sequence[float]) -> List[float]:
            candidate_source = list(row[0:3]) + list(row[4:7])
            candidate_target = [row[3], row[7]]
            candidate_prediction = continuous_prediction(candidate_source, prediction_class)
            return [observed - expected for observed, expected in
                    zip(candidate_target, candidate_prediction)]

        residual_jac = finite_jacobian(residual_function, values, matrix)
        residual_cov = propagated_covariance(residual_jac, matrix)
        prediction_record = {
            "radius": TARGET_RADIUS,
            "log2_step_from_R2": math.log(TARGET_RADIUS / 2.0, 2.0),
            "branch": "positive-real roots" if prediction_class == "R2" else
                      "principal complex phase theta in (0,pi)",
            "predicted": prediction,
            "observed": target,
            "residual": residual,
            "residual_covariance": residual_cov,
            "chi_square": quadratic(residual, inverse(residual_cov)),
            "degrees_of_freedom": 2,
            "chi_square_survival": chi_square_survival(
                quadratic(residual, inverse(residual_cov)), 2),
            "status": "correlated post-reveal design check; not independent evidence",
        }
    return {
        "channel": channel,
        "source_radii": list(SOURCE_RADII),
        "source_geometry_order": list(GEOMETRIES),
        "recurrence": rec,
        "recurrence_covariance_delta_method": rec_cov,
        "Delta_SE": delta_se,
        "Delta_over_SE": rec["Delta"] / delta_se if delta_se else None,
        "point_class": point_class,
        "identifiability": "T,D are exactly saturated by two outputs x three dyadic radii; sign is descriptive",
        "R7_check": prediction_record,
    }


def synthetic_oracles() -> Dict[str, object]:
    step = math.log(TARGET_RADIUS / 2.0, 2.0)
    records = {}
    definitions = {
        "J2": ([0.7, 0.7], [[1.2, -0.4], [-0.3, 0.8]]),
        "R2": ([0.8, 0.35], [[1.1, -0.5], [0.4, 0.9]]),
        "C2": ([0.82, 0.73], [[1.0, -0.6], [-0.2, 1.1]]),
    }
    for model, (spectrum, amplitudes) in definitions.items():
        rows = []
        target = []
        for amplitude in amplitudes:
            values = []
            for n in (0.0, 1.0, 2.0, step):
                if model == "J2":
                    value = spectrum[0] ** n * (amplitude[0] + amplitude[1] * n)
                elif model == "R2":
                    value = amplitude[0] * spectrum[0] ** n + amplitude[1] * spectrum[1] ** n
                else:
                    modulus, theta = spectrum
                    value = modulus ** n * (amplitude[0] * math.cos(theta * n) +
                                              amplitude[1] * math.sin(theta * n))
                values.append(value)
            rows.extend(values[:3])
            target.append(values[3])
        rec = recurrence(rows)
        predicted = continuous_prediction(rows, model)
        records[model] = {
            "recurrence": rec,
            "target": target,
            "predicted": predicted,
            "max_abs_error": max(abs(a - b) for a, b in zip(target, predicted)),
        }
    return records


def analyze(path: Path) -> Dict[str, object]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    block = payload["contrast_vector"]
    order, point, covariance = block["order"], block["point"], block["covariance"]
    profiles = {}
    for name, channels in (("plus", ["A_plus"]), ("minus", ["A_minus"]),
                           ("joint_common_spectrum", list(CHANNELS))):
        source_labels = labels_for(channels, SOURCE_RADII)
        source_point = subvector(order, point, source_labels)
        source_cov = submatrix(order, covariance, source_labels)
        profiles[name] = {model: profile_model(source_point, source_cov, channels, model)
                          for model in ("R1", "J2", "R2_gap1", "R2", "C2")}
    return {
        "schema": "matching-one/annulus-channel-recurrence/v1",
        "issue": 253,
        "source": str(path),
        "evidence_boundary": "one post-reveal reanalysis of PR247; never count as a new evidence row",
        "radial_coordinate": "n=log2(R/2)",
        "profile_bounds_are_design_not_theory": {
            "x": list(X_RANGE), "real_gap": list(GAP_RANGE),
            "complex_phase": list(THETA_RANGE),
        },
        "model_profiles": profiles,
        "channel_recurrences": {
            channel: channel_recurrence(order, point, covariance, channel)
            for channel in CHANNELS
        },
        "identifiability": {
            "R1": "testable: 6 source values, 3 parameters, df=3 per channel",
            "J2": "weakly testable: 6 source values, 5 parameters, df=1 per channel",
            "R2_gap1": "weakly testable ordinary-real adversary from PR260: x=17/4 versus x=21/4 fixes the relative gap to one before observable normalization",
            "R2": "saturated per channel: 6 source values, 6 parameters; Delta>0 is not a goodness-of-fit result",
            "C2": "saturated per channel: 6 source values, 6 parameters; Delta<0 has phase aliases and is not a goodness-of-fit result",
            "joint_common_spectrum": "plus and minus jointly give df 7/3/2/2 for R1/J2/R2/C2, but this adds a matching-parity sharing assumption",
        },
        "next_acquisition": {
            "geometry": {"N": 365, "first": [14, 13], "second": [19, 2]},
            "radii": [2, 4, 7, 8],
            "reason": "a third primitive orientation-contrast output overidentifies T,D on 2,4,8 and retains R7 as an off-grid propagation target",
            "minimality": "one additional geometry block; no denser radius ladder",
        },
        "mechanism_link": {
            "statement": "the minus point estimate is a complex/rotating two-state realization, consistent in form with the rotating rank-two topology plane found exactly in Issue 256",
            "non_duplication": "Issue 256 and this score are coordinate-design links, not independent votes for one mechanism",
        },
        "synthetic_oracles": synthetic_oracles(),
    }


def main(argv: Iterable[str] = None) -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    result = analyze(args.input)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()

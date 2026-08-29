#!/usr/bin/env python3
"""Continuous radial-spectrum and matching-sector tests for Issue 253.

This is deliberately downstream of the frozen recurrence score.  It uses the
off-dyadic R=7 point to profile named continuous transfer classes, and keeps
the old N325/N425 and new N365 Monte Carlo streams as two dependency groups.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Callable, Sequence

try:
    from scripts.analyze_annulus_channel_recurrence import (
        X_RANGE,
        GAP_RANGE,
        THETA_RANGE,
        basis,
        chi_square_survival,
        golden_minimize,
        grid_minimize_2d,
        inverse,
        quadratic,
        solve,
        submatrix,
        subvector,
    )
    from scripts.analyze_norm5_multiradius_pivotal import (
        contrast_vectors,
        read_rows,
    )
    from scripts.score_p253_n365_heldout import read_n365
except ModuleNotFoundError:  # Direct ``python scripts/...`` execution.
    from analyze_annulus_channel_recurrence import (  # type: ignore
        X_RANGE,
        GAP_RANGE,
        THETA_RANGE,
        basis,
        chi_square_survival,
        golden_minimize,
        grid_minimize_2d,
        inverse,
        quadratic,
        solve,
        submatrix,
        subvector,
    )
    from analyze_norm5_multiradius_pivotal import (  # type: ignore
        contrast_vectors,
        read_rows,
    )
    from score_p253_n365_heldout import read_n365  # type: ignore


CHANNELS = ("A_plus", "A_minus")
RADII = (2, 4, 7, 8)
OLD_GEOMETRIES = (325, 425)
ALL_GEOMETRIES = (325, 425, 365)
MODELS = ("R1", "J2", "R2_gap1", "R2", "C2")
RANK2_HELDOUT_MODELS = ("J2", "R2_gap1")
LABEL = re.compile(r"N(?P<N>\d+)_R(?P<R>\d+)_Delta_(?P<channel>A_plus|A_minus)$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def add_matrices(first: Sequence[Sequence[float]],
                 second: Sequence[Sequence[float]]) -> list[list[float]]:
    return [[float(a) + float(b) for a, b in zip(row_a, row_b)]
            for row_a, row_b in zip(first, second)]


def jackknife_covariance(replicates: Sequence[Sequence[float]]) -> list[list[float]]:
    count = len(replicates)
    width = len(replicates[0])
    mean = [sum(row[j] for row in replicates) / count for j in range(width)]
    return [[
        (count - 1) / count * sum(
            (row[i] - mean[i]) * (row[j] - mean[j]) for row in replicates)
        for j in range(width)] for i in range(width)]


def leave_one_means(rows: Sequence[Sequence[float]]) -> list[list[float]]:
    count = len(rows)
    totals = [sum(row[j] for row in rows) for j in range(len(rows[0]))]
    return [[(totals[j] - row[j]) / (count - 1) for j in range(len(totals))]
            for row in rows]


def canonical_labels(geometries: Sequence[int] = ALL_GEOMETRIES,
                     channels: Sequence[str] = CHANNELS) -> list[str]:
    return [f"N{geometry}_R{radius}_Delta_{channel}"
            for channel in channels for geometry in geometries for radius in RADII]


def parse_label(label: str) -> tuple[str, int, int]:
    match = LABEL.fullmatch(label)
    if match is None:
        raise ValueError(f"bad annulus label {label}")
    return match.group("channel"), int(match.group("N")), int(match.group("R"))


def ranges_for(model: str) -> list[tuple[float, float]]:
    if model in ("R1", "J2", "R2_gap1"):
        return [X_RANGE]
    if model == "R2":
        return [X_RANGE, GAP_RANGE]
    if model == "C2":
        return [X_RANGE, THETA_RANGE]
    raise ValueError(model)


def parameter_names(model: str) -> list[str]:
    return {
        "R1": ["x"],
        "J2": ["x"],
        "R2_gap1": ["center_x_with_gap_fixed_to_1"],
        "R2": ["center_x", "gap_x"],
        "C2": ["decay_x", "theta_per_log2_step"],
    }[model]


def transfer_signature(model: str, parameters: Sequence[float]) -> dict:
    if model in ("R1", "J2"):
        eigenvalue = 2.0 ** (-parameters[0])
        return {
            "eigenvalues_per_radius_doubling": [eigenvalue] * (2 if model == "J2" else 1),
            "jordan_candidate": model == "J2",
        }
    if model in ("R2", "R2_gap1"):
        center = parameters[0]
        gap = parameters[1] if model == "R2" else 1.0
        return {
            "eigenvalues_per_radius_doubling": [
                2.0 ** (-(center - gap / 2.0)),
                2.0 ** (-(center + gap / 2.0)),
            ],
            "dimension_gap": gap,
        }
    if model == "C2":
        modulus = 2.0 ** (-parameters[0])
        return {
            "complex_pair_modulus_per_radius_doubling": modulus,
            "phase_per_log2_step": parameters[1],
            "eigenvalues_real_imag": [
                [modulus * math.cos(parameters[1]), modulus * math.sin(parameters[1])],
                [modulus * math.cos(parameters[1]), -modulus * math.sin(parameters[1])],
            ],
        }
    raise ValueError(model)


class GLSContext:
    """Profile geometry-specific amplitudes at fixed sector spectra."""

    def __init__(self, labels: Sequence[str], covariance: Sequence[Sequence[float]],
                 model: str, shared_spectrum: bool):
        self.labels = list(labels)
        self.entries = [parse_label(label) for label in labels]
        self.precision = inverse(covariance)
        self.model = model
        self.shared = shared_spectrum
        self.channels = tuple(channel for channel in CHANNELS
                              if any(row[0] == channel for row in self.entries))
        self.groups = []
        for channel, geometry, _ in self.entries:
            group = (channel, geometry)
            if group not in self.groups:
                self.groups.append(group)
        self.modes = len(basis(model, 2, [0.0] * len(ranges_for(model))))
        self.nonlinear_count = len(ranges_for(model)) * (1 if shared_spectrum else len(self.channels))

    def spectra(self, parameters: Sequence[float]) -> dict[str, list[float]]:
        width = len(ranges_for(self.model))
        if self.shared:
            return {channel: list(parameters[:width]) for channel in self.channels}
        return {channel: list(parameters[i * width:(i + 1) * width])
                for i, channel in enumerate(self.channels)}

    def profile(self, point: Sequence[float], parameters: Sequence[float]) -> dict:
        spectra = self.spectra(parameters)
        columns = len(self.groups) * self.modes
        design = []
        for channel, geometry, radius in self.entries:
            row = [0.0] * columns
            offset = self.groups.index((channel, geometry)) * self.modes
            row[offset:offset + self.modes] = basis(self.model, radius, spectra[channel])
            design.append(row)
        precision_design = [[sum(self.precision[i][j] * design[j][a]
                                 for j in range(len(point)))
                             for a in range(columns)] for i in range(len(point))]
        normal = [[sum(design[i][a] * precision_design[i][b]
                       for i in range(len(point)))
                   for b in range(columns)] for a in range(columns)]
        precision_point = [sum(self.precision[i][j] * point[j]
                               for j in range(len(point))) for i in range(len(point))]
        rhs = [sum(design[i][a] * precision_point[i] for i in range(len(point)))
               for a in range(columns)]
        coefficients = solve(normal, rhs)
        fitted = [sum(a * b for a, b in zip(row, coefficients)) for row in design]
        residual = [observed - expected for observed, expected in zip(point, fitted)]
        return {
            "chi_square": quadratic(residual, self.precision),
            "coefficients": coefficients,
            "fitted": fitted,
            "residual": residual,
        }

    def objective(self, point: Sequence[float]) -> Callable[[Sequence[float]], float]:
        return lambda parameters: float(self.profile(point, parameters)["chi_square"])


def coordinate_minimize(function: Callable[[Sequence[float]], float],
                        initial: Sequence[float], ranges: Sequence[tuple[float, float]],
                        rounds: int = 5) -> tuple[list[float], float]:
    """Deterministic profiled coordinate descent, initialized by a global fit."""
    point = list(initial)
    best = function(point)
    for _ in range(rounds):
        previous = best
        for index, (low, high) in enumerate(ranges):
            def scalar(value: float) -> float:
                candidate = list(point)
                candidate[index] = value
                return function(candidate)
            value, score = global_minimize_1d(scalar, low, high, width=31)
            if score < best:
                point[index], best = value, score
        if previous - best < 1e-10:
            break
    return point, best


def global_minimize_1d(function: Callable[[float], float], low: float, high: float,
                       width: int = 65) -> tuple[float, float]:
    """Grid-bracketed scalar minimization, including both profile boundaries."""
    step = (high - low) / (width - 1)
    grid = [(low + i * step, function(low + i * step)) for i in range(width)]
    index = min(range(width), key=lambda i: grid[i][1])
    best = grid[index]
    if 0 < index < width - 1:
        candidate = golden_minimize(function, grid[index - 1][0], grid[index + 1][0],
                                    iterations=55)
        if candidate[1] < best[1]:
            best = candidate
    return best


def fit_context(context: GLSContext, point: Sequence[float],
                warm_start: Sequence[float] | None = None) -> dict:
    objective = context.objective(point)
    single_ranges = ranges_for(context.model)
    all_ranges = single_ranges * (1 if context.shared else len(context.channels))
    if warm_start is not None:
        parameters, _ = coordinate_minimize(objective, warm_start, all_ranges, rounds=3)
    elif len(all_ranges) == 1:
        value, _ = global_minimize_1d(lambda x: objective([x]), *all_ranges[0])
        parameters = [value]
    elif len(all_ranges) == 2:
        parameters, _ = grid_minimize_2d(
            lambda a, b: objective([a, b]), all_ranges[0], all_ranges[1])
    else:
        # Four-dimensional sector-separated R2/C2: seed with independent
        # channel fits, then restore the full cross-sector precision.
        parameters = []
        for channel in context.channels:
            indices = [i for i, row in enumerate(context.entries) if row[0] == channel]
            labels = [context.labels[i] for i in indices]
            covariance = [[inverse(context.precision)[i][j] for j in indices]
                          for i in indices]
            subpoint = [point[i] for i in indices]
            subcontext = GLSContext(labels, covariance, context.model, True)
            parameters.extend(fit_context(subcontext, subpoint)["parameters"])
        starts = [parameters]
        shared_seed = []
        for i in range(len(single_ranges)):
            average = sum(parameters[j * len(single_ranges) + i]
                          for j in range(len(context.channels))) / len(context.channels)
            shared_seed.extend([average] * len(context.channels))
        # Reorder [p0,p0,p1,p1] to channel-major [p0,p1,p0,p1].
        if len(single_ranges) == 2:
            shared_seed = [shared_seed[0], shared_seed[2], shared_seed[1], shared_seed[3]]
        starts.append(shared_seed)
        candidates = [coordinate_minimize(objective, start, all_ranges, rounds=6)
                      for start in starts]
        parameters, _ = min(candidates, key=lambda row: row[1])
    result = context.profile(point, parameters)
    degrees = len(point) - len(result["coefficients"]) - context.nonlinear_count
    boundary = [min(abs(value - low), abs(high - value)) < 1e-3 * (high - low)
                for value, (low, high) in zip(parameters, all_ranges)]
    names = parameter_names(context.model)
    if not context.shared:
        names = [f"{channel}:{name}" for channel in context.channels for name in names]
    spectra = context.spectra(parameters)
    return {
        "model": context.model,
        "spectrum_scope": "shared" if context.shared else "sector_separated",
        "parameters": parameters,
        "parameter_names": names,
        "parameter_ranges": all_ranges,
        "transfer_spectrum": {channel: transfer_signature(context.model, spectrum)
                              for channel, spectrum in spectra.items()},
        "boundary_hit": boundary,
        "chi_square": result["chi_square"],
        "degrees_of_freedom": degrees,
        "chi_square_survival": chi_square_survival(result["chi_square"], degrees),
        "linear_parameter_count": len(result["coefficients"]),
        "nonlinear_parameter_count": context.nonlinear_count,
    }


def combine_blocks(old: dict, new: dict) -> tuple[list[str], list[float], list[list[float]]]:
    labels = canonical_labels()
    point = []
    covariance = [[0.0] * len(labels) for _ in labels]
    for label in labels:
        _, geometry, _ = parse_label(label)
        source = old if geometry in OLD_GEOMETRIES else new
        point.append(float(source["point"][source["order"].index(label)]))
    for i, left in enumerate(labels):
        _, first_geometry, _ = parse_label(left)
        first_source = old if first_geometry in OLD_GEOMETRIES else new
        for j, right in enumerate(labels):
            _, second_geometry, _ = parse_label(right)
            if (first_geometry in OLD_GEOMETRIES) != (second_geometry in OLD_GEOMETRIES):
                continue
            source = first_source
            covariance[i][j] = float(source["covariance"][source["order"].index(left)]
                                                       [source["order"].index(right)])
    return labels, point, covariance


def spectrum_map(labels: Sequence[str], point: Sequence[float],
                 covariance: Sequence[Sequence[float]]) -> dict:
    per_sector = {}
    for channel in CHANNELS:
        selected = [label for label in labels if parse_label(label)[0] == channel]
        values = subvector(labels, point, selected)
        matrix = submatrix(labels, covariance, selected)
        per_sector[channel] = {
            model: fit_context(GLSContext(selected, matrix, model, True), values)
            for model in MODELS
        }
    shared = {}
    separated = {}
    sharing_tests = {}
    for model in MODELS:
        shared[model] = fit_context(GLSContext(labels, covariance, model, True), point)
        separated[model] = fit_context(GLSContext(labels, covariance, model, False), point)
        delta = shared[model]["chi_square"] - separated[model]["chi_square"]
        degrees = (separated[model]["nonlinear_parameter_count"] -
                   shared[model]["nonlinear_parameter_count"])
        sharing_tests[model] = {
            "null": "matching-even and matching-odd sectors share the continuous radial spectrum",
            "delta_chi_square": delta,
            "degrees_of_freedom": degrees,
            "chi_square_survival": chi_square_survival(delta, degrees),
            "nested_models": ["shared_spectrum", "sector_separated_spectra"],
            "calibration_note": "asymptotic Wilks reference; profile-boundary hits make this descriptive",
        }
    return {
        "per_sector": per_sector,
        "joint_shared_spectrum": shared,
        "joint_sector_separated_spectra": separated,
        "sector_sharing_likelihood_ratio": sharing_tests,
    }


def read_old_jackknife(batch_path: Path, analysis: dict,
                       labels: Sequence[str]) -> tuple[list[float], list[list[float]]]:
    rows, batches = read_rows(batch_path)
    full = contrast_vectors(rows, batches)[3]
    raw_order = contrast_vectors(rows, batches)[2]
    if raw_order != analysis["order"]:
        raise ValueError("old raw/analysis label order mismatch")
    delete = [contrast_vectors(rows, batches, omitted)[3] for omitted in batches]
    return subvector(raw_order, full, labels), [subvector(raw_order, row, labels) for row in delete]


def read_new_batch_vectors(metadata_path: Path, labels: Sequence[str]) -> tuple[list[float], list[list[float]]]:
    # Reuse the frozen reader for validation and canonical point/covariance.
    frozen = read_n365(metadata_path)
    import csv
    with Path(frozen["batch_csv"]).open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    indexed = {(row["label"], int(row["radius"]), int(row["batch"])): row for row in rows}
    batches = sorted({int(row["batch"]) for row in rows})
    batch_vectors = []
    for batch in batches:
        vector = []
        for label in labels:
            channel, _, radius = parse_label(label)
            column = {"A_plus": "h4_plus", "A_minus": "h4_minus"}[channel]
            first = indexed["n365_first", radius, batch]
            second = indexed["n365_second", radius, batch]
            vector.append((int(first[column]) - int(second[column])) / int(first["samples"]))
        batch_vectors.append(vector)
    return subvector(frozen["order"], frozen["point"], labels), leave_one_means(batch_vectors)


def rank2_residual(model: str, spectra: dict[str, list[float]],
                   new_point: Sequence[float]) -> list[float]:
    result = []
    for index, channel in enumerate(CHANNELS):
        values = list(new_point[index * len(RADII):(index + 1) * len(RADII)])
        source_design = [basis(model, radius, spectra[channel]) for radius in RADII[:2]]
        amplitudes = solve(source_design, values[:2])
        predicted = [sum(a * b for a, b in zip(basis(model, radius, spectra[channel]), amplitudes))
                     for radius in RADII[2:]]
        result.extend([values[2] - predicted[0], values[3] - predicted[1]])
    return result


def heldout_spectrum_score(model: str, shared: bool, old_labels: Sequence[str],
                           old_point: Sequence[float], old_delete: Sequence[Sequence[float]],
                           old_covariance: Sequence[Sequence[float]],
                           new_point: Sequence[float], new_delete: Sequence[Sequence[float]]) -> dict:
    context = GLSContext(old_labels, old_covariance, model, shared)
    fit = fit_context(context, old_point)
    spectra = context.spectra(fit["parameters"])
    residual = rank2_residual(model, spectra, new_point)

    old_replicates = []
    for row in old_delete:
        candidate = fit_context(context, row, warm_start=fit["parameters"])
        old_replicates.append(rank2_residual(model, context.spectra(candidate["parameters"]), new_point))
    new_replicates = [rank2_residual(model, spectra, row) for row in new_delete]
    old_component = jackknife_covariance(old_replicates)
    new_component = jackknife_covariance(new_replicates)
    covariance = add_matrices(old_component, new_component)
    chi_square = quadratic(residual, inverse(covariance))
    standard_errors = [math.sqrt(max(covariance[i][i], 0.0))
                       for i in range(len(residual))]
    return {
        "model": model,
        "old_training_spectrum_scope": "shared" if shared else "sector_separated",
        "old_training_fit": fit,
        "N365_calibration_radii": list(RADII[:2]),
        "N365_heldout_radii": list(RADII[2:]),
        "residual_order": [f"{channel}_R{radius}" for channel in CHANNELS for radius in RADII[2:]],
        "residual": residual,
        "standard_error": standard_errors,
        "marginal_z": [value / error for value, error in zip(residual, standard_errors)],
        "covariance": covariance,
        "covariance_components": {
            "G_old_spectrum_training_jackknife": old_component,
            "G_n365_calibration_and_target_jackknife": new_component,
        },
        "chi_square": chi_square,
        "degrees_of_freedom": len(residual),
        "chi_square_survival": chi_square_survival(chi_square, len(residual)),
        "uncertainty_method": "fixed-weight delete-one jackknife within each dependency group; independent group covariances added",
    }


def archive_inventory(root: Path) -> list[dict]:
    rows = [
        ("P43 heldout full curve", "results/server-20260828/P43-heldout-fullcurve-500m/raw/n185.metadata.json"),
        ("P49 full curve", "results/server-20260828/P49-fullcurve-doubling-100m/raw/n130.metadata.json"),
        ("P50 third lineage full curve", "results/server-20260829/P50-n145-n290-fullcurve/raw/n145_100m.metadata.json"),
        ("P57 norm5 full curve", "results/server-20260829/P57-norm5-500m/raw/n325_500m.metadata.json"),
        ("P225 C4 marked-pivotal pilot", "results/local-20260829/P225-multiradius-pivotal/raw/n130_n170_20k.metadata.json"),
    ]
    result = []
    for name, relative in rows:
        path = root / relative
        result.append({
            "name": name,
            "path": relative,
            "sha256": sha256(path),
            "eligible_for_numeric_merge": False,
            "reason": ("global threshold/rank/full-curve archive has no annulus radius and marked-insertion label"
                       if "marked-pivotal" not in name else
                       "different C4 geometry, size and radius grid; no aligned norm5 transfer readout"),
            "allowed_use": "operator-context metadata only",
        })
    return result


def render(root: Path, old_analysis_path: Path, old_batch_path: Path,
           n365_metadata_path: Path) -> dict:
    old_payload = json.loads(old_analysis_path.read_text(encoding="utf-8"))["contrast_vector"]
    new_payload = read_n365(n365_metadata_path)
    labels, point, covariance = combine_blocks(old_payload, new_payload)
    models = spectrum_map(labels, point, covariance)

    old_labels = canonical_labels(OLD_GEOMETRIES)
    new_labels = canonical_labels((365,))
    old_point, old_delete = read_old_jackknife(old_batch_path, old_payload, old_labels)
    new_point, new_delete = read_new_batch_vectors(n365_metadata_path, new_labels)
    old_covariance = submatrix(old_payload["order"], old_payload["covariance"], old_labels)
    heldout = {
        model: {
            scope: heldout_spectrum_score(
                model, scope == "shared", old_labels, old_point, old_delete,
                old_covariance, new_point, new_delete)
            for scope in ("shared", "sector_separated")
        }
        for model in RANK2_HELDOUT_MODELS
    }
    return {
        "schema": "matching-one/p253-annulus-continuous-spectrum/v1",
        "issue": 253,
        "status": "post_reveal_continuous_spectrum_and_independent_geometry_holdout",
        "radial_coordinate": "n=log2(R/2); R7 is an off-grid continuous-transfer constraint",
        "sources": {
            "old_analysis": {"path": str(old_analysis_path.relative_to(root)), "sha256": sha256(old_analysis_path)},
            "old_batches": {"path": str(old_batch_path.relative_to(root)), "sha256": sha256(old_batch_path)},
            "n365_metadata": {"path": str(n365_metadata_path.relative_to(root)), "sha256": sha256(n365_metadata_path)},
            "n365_batches": {"path": str(Path(new_payload["batch_csv"]).relative_to(root)),
                             "sha256": sha256(Path(new_payload["batch_csv"]))},
        },
        "dependency_groups": {
            "G_old": {
                "members": "N325/N425 x all four radii x both sectors",
                "batches": len(old_delete),
                "covariance": "full 16x16 delete-one covariance retained",
            },
            "G_n365": {
                "members": "N365 x all four radii x both sectors",
                "batches": len(new_delete),
                "covariance": "full 8x8 batch covariance retained",
            },
            "between_groups": "independent counter domains; block-diagonal only across G_old/G_n365",
        },
        "combined_order": labels,
        "continuous_spectrum_profiles": models,
        "heldout_N365_spectrum_transfer": heldout,
        "archive_eligibility": archive_inventory(root),
        "interpretation_boundary": {
            "exact": [
                "All-data profiles use 24 observations and the complete within-group covariance.",
                "Heldout scores train the named spectrum only on G_old, calibrate N365 amplitudes at R2/R4, and score N365 R7/R8.",
                "J2 and R2_gap1 are nonnested continuous-transfer classes; their chi-square difference is descriptive, not Wilks-calibrated.",
            ],
            "mechanism_inference": "sector-sharing likelihood ratios test one radial generator versus matching-parity-specific generators",
            "exploratory_conjecture": "a finite complex phase in A_minus would be a rotating two-state sector; a theta-boundary fit is not such evidence",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--old-analysis", type=Path,
                        default=Path("results/server-20260829/P225-norm5-multiradius/analysis.json"))
    parser.add_argument("--old-batches", type=Path,
                        default=Path("results/server-20260829/P225-norm5-multiradius/raw/norm5_200k.batches.csv"))
    parser.add_argument("--n365-metadata", type=Path,
                        default=Path("results/server-20260829/P253-n365-annulus/raw/n365_200k.metadata.json"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    resolve = lambda path: path if path.is_absolute() else root / path
    payload = render(root, resolve(args.old_analysis), resolve(args.old_batches),
                     resolve(args.n365_metadata))
    output = resolve(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {output}")


if __name__ == "__main__":
    main()

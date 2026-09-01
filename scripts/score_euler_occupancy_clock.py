#!/usr/bin/env python3
"""Score the exact occupancy-clock and fixed-K source residual in P267.

The input is the committed per-batch, per-occupation path aggregate.  No
Monte Carlo samples are generated.  ``n_occ`` below is the pre-insertion
occupation count and is deliberately not called K, to avoid confusion with
the activation thresholds K1/K2.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import subprocess
from typing import Any, Iterable

import mpmath as mp
import numpy as np
import yaml


SOURCE_COLUMNS = (
    "sum_q",
    "sum_O_ext",
    "sum_O_near",
    "sum_J_S_re",
    "sum_J_S_im",
    "sum_J_D_re",
    "sum_J_D_im",
    "sum_O_ext_J_S_re",
    "sum_O_ext_J_S_im",
    "sum_O_ext_J_D_re",
    "sum_O_ext_J_D_im",
    "sum_O_near_J_S_re",
    "sum_O_near_J_S_im",
    "sum_O_near_J_D_re",
    "sum_O_near_J_D_im",
    "sum_J_D_conj_J_S_re",
    "sum_J_D_conj_J_S_im",
    "sum_abs_J_S2",
)

COMPLEX_METRICS = (
    "P4_ext_raw_D",
    "P4_ext_clock_D",
    "P4_ext_within_n_occ_D",
    "P4_far_raw_D",
    "P4_far_clock_D",
    "P4_far_within_n_occ_D",
    "P4_far_within_n_occ_S",
    "P4_far_n_occ_JS_residual",
)

SCALAR_METRICS = (
    "beta_first_re",
    "beta_first_im",
    "beta_second_re",
    "beta_second_im",
    "ext_closure_abs",
    "far_closure_abs",
)

METRIC_ORDER = tuple(
    item
    for name in COMPLEX_METRICS
    for item in (name + "_re", name + "_im")
) + SCALAR_METRICS


def _git_bytes(commit: str, path: str) -> bytes:
    proc = subprocess.run(
        ["git", "show", f"{commit}:{path}"],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    return proc.stdout


def _input_bytes(source: dict[str, Any]) -> bytes:
    local = source.get("local_path")
    if local and Path(local).is_file():
        payload = Path(local).read_bytes()
    else:
        payload = _git_bytes(source["commit"], source["path"])
    digest = hashlib.sha256(payload).hexdigest()
    if digest != source["sha256"]:
        raise ValueError(
            f"input SHA256 mismatch for {source['path']}: {digest}"
        )
    return payload


def _falling(x: int, degree: int) -> int:
    if degree < 0:
        raise ValueError("negative falling-factorial degree")
    value = 1
    for offset in range(degree):
        value *= x - offset
    return value


def _ratio_falling(x: int, population: int, degree: int) -> float:
    if degree == 0:
        return 1.0
    if x < degree or population < degree:
        return 0.0
    return _falling(x, degree) / _falling(population, degree)


def mu_ext(n: int, n_occ: int) -> float:
    """Exact E[O_ext | n_occ] for O_ext=V-E_NN+F0."""
    return (
        n_occ
        - 2.0 * n * _ratio_falling(n_occ, n, 2)
        + n * _ratio_falling(n_occ, n, 4)
    )


def mu_near(n: int, n_occ: int) -> float:
    """Exact root-absent mean of the frozen Chebyshev-R2 nuisance.

    The 5x5 southwest-anchor stencil contains 24 non-root vertex terms, 23
    non-root horizontal and 23 non-root vertical occupied edges, four empty
    faces containing the known-absent next site, and 21 other empty faces.
    The production geometries are separately checked to be injective on the
    required local coordinate patch.
    """
    population = n - 1
    empty = population - n_occ
    return (
        24.0 * _ratio_falling(n_occ, population, 1)
        - 46.0 * _ratio_falling(n_occ, population, 2)
        + 4.0 * _ratio_falling(empty, population, 3)
        + 21.0 * _ratio_falling(empty, population, 4)
    )


def _is_lattice_vector(dx: int, dy: int, matrix: list[list[int]]) -> bool:
    a, b = matrix[0]
    c, d = matrix[1]
    det = a * d - b * c
    if det == 0:
        raise ValueError("singular period matrix")
    first = d * dx - b * dy
    second = -c * dx + a * dy
    return first % det == 0 and second % det == 0


def validate_local_injectivity(matrix: list[list[int]]) -> None:
    # O_near anchors range from -2..2 and their faces reach +1.  Checking the
    # complete -2..3 patch is stronger than the individual term requirement.
    points = [(x, y) for x in range(-2, 4) for y in range(-2, 4)]
    for left_index, left in enumerate(points):
        for right in points[left_index + 1 :]:
            if _is_lattice_vector(
                left[0] - right[0], left[1] - right[1], matrix
            ):
                raise ValueError(
                    f"local Euler patch aliases under period matrix {matrix}"
                )


def binomial_weights(size: int, p: float) -> np.ndarray:
    if not 0.0 < p < 1.0:
        raise ValueError("intrinsic root must be in (0,1)")
    weights = np.empty(size + 1, dtype=float)
    weights[0] = (1.0 - p) ** size
    ratio = p / (1.0 - p)
    for k in range(size):
        weights[k + 1] = weights[k] * (size - k) / (k + 1) * ratio
    weights /= weights.sum()
    return weights


def cos4(a: int, b: int) -> float:
    norm = a * a + b * b
    return (a**4 - 6 * a * a * b * b + b**4) / (norm * norm)


def read_archive(source: dict[str, Any], n: int) -> dict[str, Any]:
    payload = gzip.decompress(_input_bytes(source))
    text = io.TextIOWrapper(io.BytesIO(payload), encoding="utf-8", newline="")
    reader = csv.DictReader(text)
    required = {
        "n", "a", "b", "orientation", "batch", "samples", "k", *SOURCE_COLUMNS
    }
    if reader.fieldnames is None or not required.issubset(reader.fieldnames):
        raise ValueError(f"missing path columns: {sorted(required - set(reader.fieldnames or []))}")
    rows: dict[tuple[str, int], list[dict[str, str]]] = {}
    for raw in reader:
        if int(raw["n"]) != n:
            raise ValueError("manifest N differs from path N")
        rows.setdefault((raw["orientation"], int(raw["batch"])), []).append(raw)
    orientations = sorted({key[0] for key in rows})
    batches = sorted({key[1] for key in rows})
    if orientations != ["first", "second"]:
        raise ValueError(f"unexpected orientations: {orientations}")
    if len(batches) < 2:
        raise ValueError("aligned delete-one requires at least two batches")
    values = np.zeros((2, len(batches), n, len(SOURCE_COLUMNS)), dtype=float)
    samples = np.zeros((2, len(batches), n), dtype=float)
    reps: dict[str, tuple[int, int]] = {}
    for side_index, side in enumerate(orientations):
        for batch_index, batch in enumerate(batches):
            block = sorted(rows[(side, batch)], key=lambda row: int(row["k"]))
            if [int(row["k"]) for row in block] != list(range(n)):
                raise ValueError(f"incomplete n_occ path for {side} batch {batch}")
            reps[side] = (int(block[0]["a"]), int(block[0]["b"]))
            for k, raw in enumerate(block):
                samples[side_index, batch_index, k] = float(raw["samples"])
                values[side_index, batch_index, k, :] = [
                    float(raw[column]) for column in SOURCE_COLUMNS
                ]
    if np.any(samples <= 0):
        raise ValueError("non-positive microcanonical sample count")
    return {
        "orientations": orientations,
        "batches": batches,
        "reps": reps,
        "values": values,
        "samples": samples,
        "source_sha256": source["sha256"],
    }


def _combined(archive: dict[str, Any], omitted: int | None) -> tuple[np.ndarray, np.ndarray]:
    values = archive["values"]
    samples = archive["samples"]
    if omitted is None:
        return values.sum(axis=1), samples.sum(axis=1)
    return values.sum(axis=1) - values[:, omitted, :, :], samples.sum(axis=1) - samples[:, omitted, :]


def _means(values: np.ndarray, samples: np.ndarray) -> dict[str, np.ndarray]:
    return {
        name: values[:, :, index] / samples
        for index, name in enumerate(SOURCE_COLUMNS)
    }


def matching_curve(q: np.ndarray, p: float) -> float:
    coefficients = np.concatenate([q, np.array([1.0])])
    return float(np.dot(binomial_weights(len(q), p), coefficients))


def intrinsic_center(means: dict[str, np.ndarray]) -> float:
    lower, upper = 0.0, 1.0
    for _ in range(80):
        middle = (lower + upper) / 2.0
        value = (
            matching_curve(means["sum_q"][0], middle)
            + matching_curve(means["sum_q"][1], middle)
        ) / 2.0
        if value < 0:
            lower = middle
        else:
            upper = middle
    return (lower + upper) / 2.0


def _complex(means: dict[str, np.ndarray], stem: str, side: int) -> np.ndarray:
    return means[stem + "_re"][side] + 1j * means[stem + "_im"][side]


def orientation_statistics(
    means: dict[str, np.ndarray], side: int, n: int, p: float
) -> dict[str, complex | float]:
    weights = binomial_weights(n - 1, p)
    scale = n / np.arange(n, 0, -1, dtype=float)
    o_ext = means["sum_O_ext"][side]
    o_near = means["sum_O_near"][side]
    o_far = o_ext - o_near
    mu_e = np.array([mu_ext(n, k) for k in range(n)])
    mu_n = np.array([mu_near(n, k) for k in range(n)])
    mu_f = mu_e - mu_n

    j_s = _complex(means, "sum_J_S", side)
    j_d = _complex(means, "sum_J_D", side)
    oe_js = _complex(means, "sum_O_ext_J_S", side)
    oe_jd = _complex(means, "sum_O_ext_J_D", side)
    on_js = _complex(means, "sum_O_near_J_S", side)
    on_jd = _complex(means, "sum_O_near_J_D", side)
    of_js, of_jd = oe_js - on_js, oe_jd - on_jd

    def connected_raw(o: np.ndarray, oj: np.ndarray, j: np.ndarray) -> complex:
        o_mean = float(np.dot(weights, o))
        j_mean = complex(np.dot(weights, scale * j))
        return complex(np.dot(weights, scale * oj)) - o_mean * j_mean

    def clock(mu: np.ndarray, j: np.ndarray) -> complex:
        mu_mean = float(np.dot(weights, mu))
        j_mean = complex(np.dot(weights, scale * j))
        return complex(np.dot(weights, mu * scale * j)) - mu_mean * j_mean

    def within(o: np.ndarray, oj: np.ndarray, j: np.ndarray) -> complex:
        # Sum of fixed-n_occ conditional covariances.  This uses the archived
        # conditional means, not a parametric fit.
        return complex(np.dot(weights, scale * (oj - o * j)))

    raw_ext_d = connected_raw(o_ext, oe_jd, j_d)
    raw_far_d = connected_raw(o_far, of_jd, j_d)
    clock_ext_d = clock(mu_e, j_d)
    clock_far_d = clock(mu_f, j_d)
    within_ext_d = within(o_ext, oe_jd, j_d)
    within_far_d = within(o_far, of_jd, j_d)
    within_far_s = within(o_far, of_js, j_s)

    gram_ds = (
        means["sum_J_D_conj_J_S_re"][side]
        + 1j * means["sum_J_D_conj_J_S_im"][side]
    )
    gram_ss = means["sum_abs_J_S2"][side]
    centered_ds = np.dot(weights, gram_ds - j_d * np.conjugate(j_s))
    centered_ss = float(np.dot(weights, gram_ss - np.abs(j_s) ** 2))
    if not centered_ss > 0:
        raise ValueError("non-positive fixed-n_occ JS Gram denominator")
    beta = complex(centered_ds / centered_ss)
    residual = within_far_d - beta * within_far_s

    # Raw = exact-clock + within-n_occ in expectation.  The finite archive
    # closure uses the sampled O mean; reporting it catches semantic drift.
    ext_closure = raw_ext_d - clock_ext_d - within_ext_d
    far_closure = raw_far_d - clock_far_d - within_far_d
    return {
        "ext_raw_D": raw_ext_d,
        "ext_clock_D": clock_ext_d,
        "ext_within_n_occ_D": within_ext_d,
        "far_raw_D": raw_far_d,
        "far_clock_D": clock_far_d,
        "far_within_n_occ_D": within_far_d,
        "far_within_n_occ_S": within_far_s,
        "far_n_occ_JS_residual": residual,
        "beta": beta,
        "ext_closure": ext_closure,
        "far_closure": far_closure,
        "sample_minus_exact_mu_ext": float(np.dot(weights, o_ext - mu_e)),
        "sample_minus_exact_mu_near": float(np.dot(weights, o_near - mu_n)),
    }


def projected(
    values: np.ndarray,
    samples: np.ndarray,
    archive: dict[str, Any],
    n: int,
) -> tuple[float, dict[str, float], dict[str, Any]]:
    means = _means(values, samples)
    p = intrinsic_center(means)
    stats = [orientation_statistics(means, side, n, p) for side in range(2)]
    first_rep = archive["reps"]["first"]
    second_rep = archive["reps"]["second"]
    leverage = cos4(*first_rep) - cos4(*second_rep)
    if leverage == 0:
        raise ValueError("zero H4 leverage")
    points: dict[str, complex] = {}
    source_names = {
        "P4_ext_raw_D": "ext_raw_D",
        "P4_ext_clock_D": "ext_clock_D",
        "P4_ext_within_n_occ_D": "ext_within_n_occ_D",
        "P4_far_raw_D": "far_raw_D",
        "P4_far_clock_D": "far_clock_D",
        "P4_far_within_n_occ_D": "far_within_n_occ_D",
        "P4_far_within_n_occ_S": "far_within_n_occ_S",
        "P4_far_n_occ_JS_residual": "far_n_occ_JS_residual",
    }
    for output_name, source_name in source_names.items():
        points[output_name] = complex(
            (stats[0][source_name] - stats[1][source_name]) / leverage
        )
    flattened: dict[str, float] = {}
    for name, value in points.items():
        flattened[name + "_re"] = value.real
        flattened[name + "_im"] = value.imag
    for side_name, side in zip(("first", "second"), stats):
        beta = complex(side["beta"])
        flattened[f"beta_{side_name}_re"] = beta.real
        flattened[f"beta_{side_name}_im"] = beta.imag
    flattened["ext_closure_abs"] = abs(
        points["P4_ext_raw_D"]
        - points["P4_ext_clock_D"]
        - points["P4_ext_within_n_occ_D"]
    )
    flattened["far_closure_abs"] = abs(
        points["P4_far_raw_D"]
        - points["P4_far_clock_D"]
        - points["P4_far_within_n_occ_D"]
    )
    details = {
        "leverage": leverage,
        "representatives": {"first": list(first_rep), "second": list(second_rep)},
        "orientation": {
            name: {
                key: ([complex(value).real, complex(value).imag] if isinstance(value, complex) else value)
                for key, value in row.items()
            }
            for name, row in zip(("first", "second"), stats)
        },
    }
    return p, flattened, details


def jackknife_covariance(rows: Iterable[dict[str, float]]) -> np.ndarray:
    matrix = np.array([[row[name] for name in METRIC_ORDER] for row in rows])
    mean = matrix.mean(axis=0)
    count = matrix.shape[0]
    return (count - 1.0) / count * (matrix - mean).T @ (matrix - mean)


def chi2_tail(value: float, degrees: int) -> dict[str, float | str]:
    probability = mp.gammainc(
        degrees / 2.0, value / 2.0, mp.inf
    ) / mp.gamma(degrees / 2.0)
    return {
        "p_value": mp.nstr(probability, 12),
        "log10_p_value": float(mp.log10(probability)),
    }


def mahalanobis_2d(point: complex, covariance: np.ndarray) -> float:
    vector = np.array([point.real, point.imag])
    return float(vector @ np.linalg.pinv(covariance, rcond=1e-12) @ vector)


def _complex_point(point: dict[str, float], name: str) -> complex:
    return complex(point[name + "_re"], point[name + "_im"])


def score_run(run: dict[str, Any], source_commit: str) -> dict[str, Any]:
    n = int(run["N"])
    for matrix in run["period_matrices"]:
        validate_local_injectivity(matrix)
    source = {
        "commit": source_commit,
        "path": run["path"],
        "sha256": run["sha256"],
        "local_path": run.get("local_path"),
    }
    archive = read_archive(source, n)
    values, samples = _combined(archive, None)
    center, point, details = projected(values, samples, archive, n)
    delete_rows = []
    delete_centers = []
    for omitted in range(len(archive["batches"])):
        leave_values, leave_samples = _combined(archive, omitted)
        leave_center, leave_point, _ = projected(
            leave_values, leave_samples, archive, n
        )
        delete_centers.append(leave_center)
        delete_rows.append(leave_point)
    covariance = jackknife_covariance(delete_rows)
    standard_errors = {
        name: math.sqrt(max(0.0, covariance[index, index]))
        for index, name in enumerate(METRIC_ORDER)
    }
    residual_name = "P4_far_n_occ_JS_residual"
    residual = _complex_point(point, residual_name)
    indices = [METRIC_ORDER.index(residual_name + suffix) for suffix in ("_re", "_im")]
    residual_covariance = covariance[np.ix_(indices, indices)]
    chi2 = mahalanobis_2d(residual, residual_covariance)
    raw_far = _complex_point(point, "P4_far_raw_D")
    clock_far = _complex_point(point, "P4_far_clock_D")
    within_far = _complex_point(point, "P4_far_within_n_occ_D")
    tail = chi2_tail(chi2, 2)
    return {
        "N": n,
        "source": {
            "commit": source_commit,
            "path": run["path"],
            "sha256": run["sha256"],
        },
        "batches": len(archive["batches"]),
        "intrinsic_center": center,
        "intrinsic_center_delete_one_se": math.sqrt(
            (len(delete_centers) - 1.0) / len(delete_centers)
            * sum((value - sum(delete_centers) / len(delete_centers)) ** 2 for value in delete_centers)
        ),
        "metric_order": list(METRIC_ORDER),
        "point": point,
        "standard_error": standard_errors,
        "delete_one_covariance": covariance.tolist(),
        "primary_residual": {
            "complex": [residual.real, residual.imag],
            "covariance_re_im": residual_covariance.tolist(),
            "chi2_2d_zero": chi2,
            **tail,
            "clock_magnitude_fraction_of_raw_far_D": abs(clock_far) / abs(raw_far),
            "fixed_n_occ_magnitude_fraction_of_raw_far_D": abs(within_far) / abs(raw_far),
            "magnitude_retained_from_raw_far_D": abs(residual) / abs(raw_far),
            "magnitude_retained_after_fixed_n_occ_before_JS": abs(residual) / abs(within_far),
        },
        "details": details,
    }


def build(manifest: dict[str, Any]) -> dict[str, Any]:
    runs = [score_run(run, manifest["source_commit"]) for run in manifest["runs"]]
    joint_chi2 = sum(run["primary_residual"]["chi2_2d_zero"] for run in runs)
    joint_tail = chi2_tail(joint_chi2, 2 * len(runs))
    return {
        "schema": "matching-one.euler-occupancy-clock.v1",
        "status": "existing_data_reanalysis_no_new_samples",
        "dependency_group": manifest["dependency_group"],
        "source_commit": manifest["source_commit"],
        "definitions": {
            "occupation_coordinate": "n_occ is pre-insertion occupied-site count; it is not activation K1 or K2",
            "O_ext": "V-E_NN+F0",
            "mu_ext": "n_occ-2N(n_occ)_2/(N)_2+N(n_occ)_4/(N)_4",
            "O_near": "frozen Chebyshev-R2 next-site nuisance with exact root-absent conditional mean",
            "source_projection": "within each n_occ center JD and JS in the same-next-site Horvitz Gram, then JD_res=JD-beta*JS",
            "primary": "P4 Cov(O_far-mu_far(n_occ), JD_res) with every root, conditional mean and beta recomputed delete-one",
        },
        "runs": runs,
        "joint_primary_zero": {
            "chi2": joint_chi2,
            "degrees_of_freedom": 2 * len(runs),
            **joint_tail,
            "covariance_rule": "sizes use disjoint committed counter ranges; orientation and all model views within a size share one aligned delete-one covariance",
        },
        "claim_boundary": [
            "a nonzero residual escapes the declared sigma(n_occ)+span(JS) nuisance in this observer/source metric only",
            "it is not by itself a new microscopic field, Q4 identification, continuum exponent, or independent evidence block",
            "clock, within-n_occ, JS projection, O_ext and O_far are correlated coordinates of the same N325/N425 production archives",
        ],
    }


def render_markdown(result: dict[str, Any]) -> str:
    lines = [
        "# Euler occupancy-clock and source residual",
        "",
        "This is a zero-new-sample reanalysis of the committed N325/N425 path",
        "aggregates. `n_occ` is the pre-insertion occupation count, not K1/K2.",
        "",
        "## Decision",
        "",
    ]
    joint = result["joint_primary_zero"]
    lines.append(
        f"The residual after the exact occupancy clock and the fixed-`n_occ` JS "
        f"source projection has joint `chi2={joint['chi2']:.6g}/{joint['degrees_of_freedom']}`, "
        f"`p={joint['p_value']}` (`log10 p={joint['log10_p_value']:.3f}`)."
    )
    lines.extend([
        "",
        "| N | residual complex | zero p | clock / raw | residual / raw | residual / fixed-`n_occ` |",
        "|---:|---:|---:|---:|---:|---:|",
    ])
    for run in result["runs"]:
        primary = run["primary_residual"]
        z = primary["complex"]
        lines.append(
            f"| {run['N']} | `{z[0]:.8g}{z[1]:+.8g}i` | "
            f"`{primary['p_value']}` | `{primary['clock_magnitude_fraction_of_raw_far_D']:.6g}` | "
            f"`{primary['magnitude_retained_from_raw_far_D']:.6g}` | "
            f"`{primary['magnitude_retained_after_fixed_n_occ_before_JS']:.6g}` |"
        )
    lines.extend(
        [
            "",
            "## Exact population decomposition",
            "",
            "For each orientation and delete-one replicate:",
            "",
            "```text",
            "O_ext = mu_ext(n_occ) + fixed-n_occ residual,",
            "mu_ext(k) = k - 2N (k)_2/(N)_2 + N (k)_4/(N)_4.",
            "```",
            "",
            "The frozen radius-2 local nuisance also has an exact root-absent",
            "conditional mean on these locally injective period quotients. The report",
            "then centers JD and JS within each `n_occ`, recomputes their same-next-site",
            "Gram coefficient, and scores the far-Euler coupling to the remaining source.",
            "",
            "The finite archive records the small sample closure residual between the",
            "analytic clock and the empirical fixed-occupation covariance in the JSON.",
            "The clock, within-occupation and source-projected rows are coordinates of",
            "the same production block, not independent votes. A surviving residual",
            "escapes only the declared `sigma(n_occ)+span(JS)` nuisance. It is not a",
            "Q4 field identification, exponent, or proof of a second microscopic source.",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-md", type=Path, required=True)
    args = parser.parse_args()
    manifest = yaml.safe_load(args.manifest.read_text(encoding="utf-8"))
    result = build(manifest)
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    args.output_md.parent.mkdir(parents=True, exist_ok=True)
    args.output_md.write_text(render_markdown(result), encoding="utf-8")


if __name__ == "__main__":
    main()

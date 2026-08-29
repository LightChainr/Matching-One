#!/usr/bin/env python3
"""Score the frozen seven-quantile shape contract on threshold-rank archives."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

import mpmath as mp


GRID = (0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)
TAIL_INDEX = (0, 1, 5, 6)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


@dataclass
class Archive:
    path: Path
    n: int
    batches: tuple[int, ...]
    per_batch: dict[int, list[int]]
    samples_per_component_batch: int

    @classmethod
    def read(cls, path: Path) -> "Archive":
        required = {"n", "orientation", "batch", "samples", "kind", "k", "count"}
        components: dict[tuple[int, str, str], list[int]] = {}
        samples: dict[tuple[int, str, str], int] = {}
        n_value: int | None = None
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            missing = required - set(reader.fieldnames or ())
            if missing:
                raise ValueError(f"{path} lacks columns {sorted(missing)}")
            for raw in reader:
                n = int(raw["n"])
                if n_value is None:
                    n_value = n
                if n != n_value:
                    raise ValueError(f"mixed N in {path}")
                orientation = raw["orientation"]
                kind = raw["kind"]
                if orientation not in ("first", "second") or kind not in ("minus", "plus"):
                    raise ValueError(f"invalid component in {path}")
                batch = int(raw["batch"])
                k = int(raw["k"])
                count = int(raw["count"])
                if not 1 <= k <= n or count <= 0:
                    raise ValueError(f"invalid threshold row in {path}")
                key = (batch, orientation, kind)
                histogram = components.setdefault(key, [0] * (n + 1))
                histogram[k] += count
                sample_count = int(raw["samples"])
                if key in samples and samples[key] != sample_count:
                    raise ValueError(f"inconsistent sample count in {path}")
                samples[key] = sample_count
        if n_value is None:
            raise ValueError(f"empty archive {path}")
        batches = tuple(sorted({key[0] for key in components}))
        expected = {
            (batch, orientation, kind)
            for batch in batches
            for orientation in ("first", "second")
            for kind in ("minus", "plus")
        }
        if set(components) != expected:
            raise ValueError(f"incomplete batch/orientation/kind grid in {path}")
        distinct_samples = set(samples.values())
        if len(distinct_samples) != 1:
            raise ValueError(f"unequal component batch sizes in {path}")
        per_batch: dict[int, list[int]] = {}
        for batch in batches:
            mixed = [0] * (n_value + 1)
            for orientation in ("first", "second"):
                for kind in ("minus", "plus"):
                    component = components[(batch, orientation, kind)]
                    if sum(component) != samples[(batch, orientation, kind)]:
                        raise ValueError(f"histogram total mismatch in {path}, batch {batch}")
                    for k, count in enumerate(component):
                        mixed[k] += count
            per_batch[batch] = mixed
        return cls(path, n_value, batches, per_batch, distinct_samples.pop())

    def mixture(self, omitted: int | None = None) -> tuple[list[int], int]:
        selected = [batch for batch in self.batches if batch != omitted]
        if not selected:
            raise ValueError("cannot omit the only batch")
        histogram = [0] * (self.n + 1)
        for batch in selected:
            for k, count in enumerate(self.per_batch[batch]):
                histogram[k] += count
        total = 4 * self.samples_per_component_batch * len(selected)
        if sum(histogram) != total:
            raise ValueError("mixture total mismatch")
        return histogram, total


def cdf(histogram: Sequence[int], total: int, p: float) -> float:
    n = len(histogram) - 1
    if p <= 0.0:
        return 0.0
    if p >= 1.0:
        return 1.0
    q = 1.0 - p
    binomial = q**n
    cumulative = 0
    value = 0.0
    for occupied in range(n + 1):
        if occupied:
            cumulative += histogram[occupied]
        value += cumulative * binomial
        if occupied < n:
            binomial *= (n - occupied) * p / ((occupied + 1) * q)
    return value / total


def quantile(histogram: Sequence[int], total: int, probability: float) -> float:
    lower, upper = 0.0, 1.0
    for _ in range(64):
        midpoint = (lower + upper) / 2.0
        if cdf(histogram, total, midpoint) < probability:
            lower = midpoint
        else:
            upper = midpoint
    return (lower + upper) / 2.0


def profile(histogram: Sequence[int], total: int) -> dict[str, Any]:
    quantiles = [quantile(histogram, total, probability) for probability in GRID]
    center = quantiles[3]
    scale = quantiles[4] - quantiles[2]
    if not scale > 0.0:
        raise ValueError("nonpositive IQR")
    standardized = [(value - center) / scale for value in quantiles]
    return {"center": center, "scale": scale, "quantiles": quantiles, "standardized": standardized}


def jackknife_covariance(rows: Sequence[Sequence[float]]) -> list[list[float]]:
    count = len(rows)
    means = [math.fsum(row[index] for row in rows) / count for index in range(len(rows[0]))]
    factor = (count - 1) / count
    return [
        [
            factor * math.fsum((row[i] - means[i]) * (row[j] - means[j]) for row in rows)
            for j in range(len(means))
        ]
        for i in range(len(means))
    ]


def add_covariances(left: Sequence[Sequence[float]], right: Sequence[Sequence[float]]) -> list[list[float]]:
    return [[left[i][j] + right[i][j] for j in range(len(left))] for i in range(len(left))]


def quadratic(residual: Sequence[float], covariance: Sequence[Sequence[float]]) -> tuple[float, int, list[float]]:
    matrix = mp.matrix(covariance)
    eigenvalues, eigenvectors = mp.eigsy(matrix)
    values = [float(eigenvalues[index]) for index in range(len(residual))]
    largest = max(values)
    keep = [index for index, value in enumerate(values) if value > largest * 1e-10]
    if not keep:
        raise ValueError("shape covariance has zero numerical rank")
    vector = mp.matrix(residual)
    transformed = eigenvectors.T * vector
    score = mp.fsum(transformed[index] ** 2 / eigenvalues[index] for index in keep)
    return float(score), len(keep), values


def compare_archives(reference: Archive, target: Archive, paired_batches: bool) -> dict[str, Any]:
    ref_hist, ref_total = reference.mixture()
    target_hist, target_total = target.mixture()
    ref = profile(ref_hist, ref_total)
    child = profile(target_hist, target_total)
    residual_all = [child["standardized"][i] - ref["standardized"][i] for i in range(len(GRID))]
    residual = [residual_all[i] for i in TAIL_INDEX]

    ref_deletes = []
    for batch in reference.batches:
        histogram, total = reference.mixture(batch)
        row = profile(histogram, total)["standardized"]
        ref_deletes.append([row[i] for i in TAIL_INDEX])
    target_deletes = []
    for batch in target.batches:
        histogram, total = target.mixture(batch)
        row = profile(histogram, total)["standardized"]
        target_deletes.append([row[i] for i in TAIL_INDEX])
    if paired_batches:
        if reference.batches != target.batches:
            raise ValueError("paired comparison requires identical batch labels")
        delete_residuals = [
            [target_deletes[row][column] - ref_deletes[row][column] for column in range(len(TAIL_INDEX))]
            for row in range(len(ref_deletes))
        ]
        covariance = jackknife_covariance(delete_residuals)
        dependency = "common seed/counter batches; covariance formed from aligned delete-one residuals"
    else:
        covariance = add_covariances(jackknife_covariance(ref_deletes), jackknife_covariance(target_deletes))
        dependency = "independent archives; reference and target jackknife covariances added"
    chi2, rank, eigenvalues = quadratic(residual, covariance)
    p_mp = mp.gammainc(mp.mpf(rank) / 2, chi2 / 2, mp.inf, regularized=True)
    p_value = float(p_mp)
    log10_p = float(mp.log(p_mp, 10))

    return {
        "reference": {"N": reference.n, **ref},
        "target": {"N": target.n, **child},
        "location_shift": child["center"] - ref["center"],
        "scale_ratio": child["scale"] / ref["scale"],
        "shape_residuals_all_grid": residual_all,
        "shape_sse_all_grid": math.fsum(value * value for value in residual_all),
        "shape_max_abs_all_grid": max(abs(value) for value in residual_all),
        "tail_test_indices": list(TAIL_INDEX),
        "tail_covariance": covariance,
        "tail_covariance_eigenvalues": eigenvalues,
        "tail_chi_square": chi2,
        "tail_covariance_rank": rank,
        "tail_p_value": p_value,
        "tail_log10_p_value": log10_p,
        "dependency": dependency,
    }


def cosine(left: Sequence[float], right: Sequence[float]) -> float:
    numerator = math.fsum(a * b for a, b in zip(left, right))
    denominator = math.sqrt(math.fsum(a * a for a in left) * math.fsum(b * b for b in right))
    if denominator == 0.0:
        raise ValueError("zero shape residual cannot define a direction")
    return numerator / denominator


def projection_scale(source: Sequence[float], target: Sequence[float]) -> float:
    denominator = math.fsum(value * value for value in source)
    if denominator == 0.0:
        raise ValueError("zero source residual cannot define a projection scale")
    return math.fsum(left * right for left, right in zip(source, target)) / denominator


def effective_power_from_cover_ratio(ratio: float) -> float:
    """Solve (1-5^-q)/(1-2^-q)=ratio for the leading correction power q."""
    lower, upper = 1e-8, 8.0
    limiting = math.log(5.0) / math.log(2.0)
    if not 1.0 < ratio < limiting:
        raise ValueError("cover ratio lies outside the positive-power range")
    for _ in range(100):
        midpoint = (lower + upper) / 2.0
        value = (1.0 - 5.0 ** (-midpoint)) / (1.0 - 2.0 ** (-midpoint))
        if value > ratio:
            lower = midpoint
        else:
            upper = midpoint
    return (lower + upper) / 2.0


def render_markdown(artifact: dict[str, Any]) -> str:
    lines = [
        "# Production threshold-shape collapse score",
        "",
        "The frozen q05/q10/q25/q50/q75/q90/q95 contract is applied to the equal mixture of",
        "Kminus/Kplus and both Gaussian orientations. Central standardized coordinates are",
        "identically fixed, so the covariance score uses the four tail coordinates.",
        "",
        "| lineage | location shift | scale ratio | shape SSE | max | tail chi-square | log10 p |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for name, comparison in artifact["comparisons"].items():
        lines.append(
            "| %s | %.10g | %.10g | %.10g | %.10g | %.6g/%d | %.6g |"
            % (
                name,
                comparison["location_shift"],
                comparison["scale_ratio"],
                comparison["shape_sse_all_grid"],
                comparison["shape_max_abs_all_grid"],
                comparison["tail_chi_square"],
                comparison["tail_covariance_rank"],
                comparison["tail_log10_p_value"],
            )
        )
    diagnostic = artifact["direction_diagnostic"]
    lines += [
        "",
        "## Post-reveal shape-flow direction",
        "",
        "- same-parent norm-2 versus norm-5 cosine: N65 `%.9f`, N85 `%.9f`;"
        % (
            diagnostic["same_parent_N65_norm2_vs_norm5_cosine"],
            diagnostic["same_parent_N85_norm2_vs_norm5_cosine"],
        ),
        "- norm-5/norm-2 projection scale: N65 `%.6f`, N85 `%.6f`;"
        % (
            diagnostic["norm5_over_norm2_projection_scale"]["N65"],
            diagnostic["norm5_over_norm2_projection_scale"]["N85"],
        ),
        "- effective positive correction power: N65 `%.6f`, N85 `%.6f` (5/8 predicts ratio `%.6f`)."
        % (
            diagnostic["effective_positive_correction_power"]["N65"],
            diagnostic["effective_positive_correction_power"]["N85"],
            diagnostic["power_5_over_8_predicted_cover_ratio"],
        ),
        "- independent N185->N265 direction cosine range across the four cover vectors: `%.6f` to `%.6f`."
        % (
            min(diagnostic["independent_N185_to_N265_cosines"].values()),
            max(diagnostic["independent_N185_to_N265_cosines"].values()),
        ),
        "",
        "## Interpretation",
        "",
        artifact["interpretation"],
        "",
        "This is retrospective production evidence from existing archives. It is not a",
        "prospective universality test and does not compare distinct microscopic models.",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n65", type=Path, required=True)
    parser.add_argument("--n85", type=Path, required=True)
    parser.add_argument("--n130", type=Path, required=True)
    parser.add_argument("--n170", type=Path, required=True)
    parser.add_argument("--n185", type=Path, required=True)
    parser.add_argument("--n265", type=Path, required=True)
    parser.add_argument("--n325", type=Path, required=True)
    parser.add_argument("--n425", type=Path, required=True)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args()
    paths = {
        name: getattr(args, name)
        for name in ("n65", "n85", "n130", "n170", "n185", "n265", "n325", "n425")
    }
    archives = {name: Archive.read(path) for name, path in paths.items()}
    comparisons = {
        "N65_to_N130_norm2": compare_archives(archives["n65"], archives["n130"], True),
        "N85_to_N170_norm2": compare_archives(archives["n85"], archives["n170"], True),
        "N65_to_N325_norm5": compare_archives(archives["n65"], archives["n325"], False),
        "N85_to_N425_norm5": compare_archives(archives["n85"], archives["n425"], False),
        "N185_to_N265_independent": compare_archives(archives["n185"], archives["n265"], False),
    }
    n65_scale = projection_scale(
        comparisons["N65_to_N130_norm2"]["shape_residuals_all_grid"],
        comparisons["N65_to_N325_norm5"]["shape_residuals_all_grid"],
    )
    n85_scale = projection_scale(
        comparisons["N85_to_N170_norm2"]["shape_residuals_all_grid"],
        comparisons["N85_to_N425_norm5"]["shape_residuals_all_grid"],
    )
    direction_diagnostic = {
        "same_parent_N65_norm2_vs_norm5_cosine": cosine(
            comparisons["N65_to_N130_norm2"]["shape_residuals_all_grid"],
            comparisons["N65_to_N325_norm5"]["shape_residuals_all_grid"],
        ),
        "same_parent_N85_norm2_vs_norm5_cosine": cosine(
            comparisons["N85_to_N170_norm2"]["shape_residuals_all_grid"],
            comparisons["N85_to_N425_norm5"]["shape_residuals_all_grid"],
        ),
        "norm5_over_norm2_projection_scale": {"N65": n65_scale, "N85": n85_scale},
        "effective_positive_correction_power": {
            "N65": effective_power_from_cover_ratio(n65_scale),
            "N85": effective_power_from_cover_ratio(n85_scale),
        },
        "power_5_over_8_predicted_cover_ratio": (1.0 - 5.0 ** (-0.625)) / (1.0 - 2.0 ** (-0.625)),
        "independent_N185_to_N265_cosines": {
            name: cosine(
                comparisons["N185_to_N265_independent"]["shape_residuals_all_grid"],
                comparisons[name]["shape_residuals_all_grid"],
            )
            for name in (
                "N65_to_N130_norm2",
                "N85_to_N170_norm2",
                "N65_to_N325_norm5",
                "N85_to_N425_norm5",
            )
        },
        "boundary": "post-reveal Euclidean direction diagnostic; no covariance-weighted mechanism score",
    }
    artifact = {
        "schema": "matching-one/threshold-shape-collapse-production/v1",
        "issue": 122,
        "status": "retrospective_existing_production_archives",
        "quantile_convention": "continuous generalized inverse by monotone bisection",
        "threshold_distribution": "equal mixture of Kminus/Kplus beta order statistics and both orientations",
        "grid": list(GRID),
        "input_sha256": {name: sha256(path) for name, path in paths.items()},
        "comparisons": comparisons,
        "direction_diagnostic": direction_diagnostic,
        "interpretation": (
            "A small raw standardized residual is not automatically statistical collapse: the high-statistic "
            "archives resolve finite-size tail deformation. Compare the two independent norm-5 lineages as "
            "one mechanism pattern, not as additive confirmations."
        ),
    }
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.markdown.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    args.markdown.write_text(render_markdown(artifact), encoding="utf-8")


if __name__ == "__main__":
    main()

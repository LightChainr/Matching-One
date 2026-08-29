#!/usr/bin/env python3
"""Fail-closed location/scale/shape contract for weighted threshold distributions."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable, Sequence


QUANTILE_GRID = tuple(Fraction(value, 100) for value in (5, 10, 25, 50, 75, 90, 95))


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


@dataclass(frozen=True)
class WeightedDistribution:
    values: tuple[Fraction, ...]
    counts: tuple[int, ...]

    def __post_init__(self) -> None:
        if not self.values or len(self.values) != len(self.counts):
            raise ValueError("values and counts must be nonempty and aligned")
        if any(right <= left for left, right in zip(self.values, self.values[1:])):
            raise ValueError("values must be strictly increasing")
        if any(not isinstance(count, int) or count <= 0 for count in self.counts):
            raise ValueError("counts must be positive integers")

    @property
    def total(self) -> int:
        return sum(self.counts)

    def quantile(self, probability: Fraction) -> Fraction:
        if not Fraction(0) <= probability <= Fraction(1):
            raise ValueError("probability must lie in [0,1]")
        if probability == 0:
            return self.values[0]
        cumulative = 0
        for value, count in zip(self.values, self.counts):
            cumulative += count
            if cumulative * probability.denominator >= probability.numerator * self.total:
                return value
        return self.values[-1]

    def affine(self, offset: Fraction, scale: Fraction) -> "WeightedDistribution":
        if scale <= 0:
            raise ValueError("affine scale must be positive")
        return WeightedDistribution(tuple(offset + scale * value for value in self.values), self.counts)


def standardized_profile(
    distribution: WeightedDistribution,
    grid: Sequence[Fraction] = QUANTILE_GRID,
) -> dict[str, Any]:
    center = distribution.quantile(Fraction(1, 2))
    lower = distribution.quantile(Fraction(1, 4))
    upper = distribution.quantile(Fraction(3, 4))
    scale = upper - lower
    if scale <= 0:
        raise ValueError("interquartile scale must be positive")
    quantiles = [distribution.quantile(probability) for probability in grid]
    standardized = [(value - center) / scale for value in quantiles]
    return {
        "center": center,
        "scale": scale,
        "quantiles": quantiles,
        "standardized_quantiles": standardized,
    }


def compare_shapes(
    reference: WeightedDistribution,
    target: WeightedDistribution,
    grid: Sequence[Fraction] = QUANTILE_GRID,
) -> dict[str, Any]:
    first = standardized_profile(reference, grid)
    second = standardized_profile(target, grid)
    residuals = [right - left for left, right in zip(
        first["standardized_quantiles"], second["standardized_quantiles"])]
    return {
        "location_shift": second["center"] - first["center"],
        "scale_ratio": second["scale"] / first["scale"],
        "standardized_shape_residuals": residuals,
        "shape_sse": sum((value * value for value in residuals), Fraction(0)),
        "shape_max_abs": max((abs(value) for value in residuals), default=Fraction(0)),
    }


def serialize(value: Any) -> Any:
    if isinstance(value, Fraction):
        return fraction_text(value)
    if isinstance(value, dict):
        return {key: serialize(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [serialize(item) for item in value]
    return value


def fixture(values: Iterable[int], counts: Iterable[int]) -> WeightedDistribution:
    return WeightedDistribution(tuple(Fraction(value) for value in values), tuple(counts))


def build_artifact() -> dict[str, Any]:
    reference = fixture(range(7), (1, 2, 4, 6, 4, 2, 1))
    affine = reference.affine(Fraction(5), Fraction(3))
    deformed = fixture(range(7), (3, 1, 2, 6, 4, 2, 2))
    affine_comparison = compare_shapes(reference, affine)
    deformation_comparison = compare_shapes(reference, deformed)
    assert affine_comparison["shape_sse"] == 0
    assert deformation_comparison["shape_sse"] > 0
    return serialize({
        "schema": "matching-one/threshold-shape-collapse-contract/v1",
        "issue": 122,
        "data_class": "synthetic_contract_fixtures_only",
        "quantile_convention": "left-continuous generalized inverse: min{x:F(x)>=q}",
        "center": "q50",
        "scale": "q75-q25",
        "quantile_grid": QUANTILE_GRID,
        "fixtures": {
            "reference": {"values": reference.values, "counts": reference.counts},
            "positive_affine": {"values": affine.values, "counts": affine.counts},
            "tail_deformation": {"values": deformed.values, "counts": deformed.counts},
        },
        "comparisons": {
            "positive_affine": affine_comparison,
            "tail_deformation": deformation_comparison,
        },
        "decision_contract": {
            "location": "target q50 - reference q50",
            "scale": "target IQR / reference IQR",
            "shape": "standardized quantile residual vector on the frozen grid",
            "scalar_shape_summaries": ["sum of squared residuals", "maximum absolute residual"],
            "fail_closed": ["nonpositive counts", "unsorted support", "zero IQR", "nonpositive affine scale"],
        },
        "boundary": (
            "The contract proves positive-affine invariance and deformation sensitivity on synthetic "
            "fixtures only. It provides no covariance, p-value, conformal map, or universality result."
        ),
    })


def render_markdown(artifact: dict[str, Any]) -> str:
    affine = artifact["comparisons"]["positive_affine"]
    deformed = artifact["comparisons"]["tail_deformation"]
    return "\n".join([
        "# Threshold-distribution shape-collapse contract", "",
        "All values below are synthetic contract fixtures; they are not simulation results.", "",
        "- quantile: left-continuous generalized inverse;",
        "- center: `q50`;",
        "- scale: `q75-q25`;",
        "- frozen grid: `%s`." % artifact["quantile_grid"], "",
        "| comparison | location shift | scale ratio | shape SSE | max absolute shape residual |",
        "|---|---:|---:|---:|---:|",
        "| positive affine | `%s` | `%s` | `%s` | `%s` |" % (
            affine["location_shift"], affine["scale_ratio"], affine["shape_sse"], affine["shape_max_abs"]),
        "| tail deformation | `%s` | `%s` | `%s` | `%s` |" % (
            deformed["location_shift"], deformed["scale_ratio"], deformed["shape_sse"],
            deformed["shape_max_abs"]), "", "## Interpretation boundary", "",
        artifact["boundary"], "",
    ])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--format", choices=("json", "markdown"), default="json")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    artifact = build_artifact()
    rendered = json.dumps(artifact, indent=2, sort_keys=True) + "\n" if args.format == "json" else render_markdown(artifact)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()

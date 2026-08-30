#!/usr/bin/env python3
"""Dependency-free exact algebra for finite derivative jets."""

from __future__ import annotations

from fractions import Fraction
import math
import re
from typing import Any, Sequence


def canonical_fraction(value: Any, label: str = "value") -> Fraction:
    if isinstance(value, bool):
        raise ValueError(f"{label} must not be boolean")
    if isinstance(value, int):
        return Fraction(value)
    if not isinstance(value, str) or re.fullmatch(r"-?\d+(?:/[1-9]\d*)?", value) is None:
        raise ValueError(f"{label} must be a canonical exact fraction")
    result = Fraction(value)
    if str(result) != value:
        raise ValueError(f"{label} must be reduced")
    return result


def parse_derivatives(values: Sequence[Any], order: int, label: str) -> list[Fraction]:
    if len(values) != order + 1:
        raise ValueError(f"{label} must contain derivatives 0..{order}")
    return [canonical_fraction(value, f"{label}[{index}]") for index, value in enumerate(values)]


def derivatives_to_series(derivatives: Sequence[Fraction]) -> list[Fraction]:
    return [value / math.factorial(order) for order, value in enumerate(derivatives)]


def series_to_derivatives(series: Sequence[Fraction]) -> list[Fraction]:
    return [value * math.factorial(order) for order, value in enumerate(series)]


def multiply_series(
    left: Sequence[Fraction], right: Sequence[Fraction], order: int
) -> list[Fraction]:
    output = [Fraction(0) for _ in range(order + 1)]
    for i, left_value in enumerate(left[: order + 1]):
        for j, right_value in enumerate(right[: order + 1 - i]):
            output[i + j] += left_value * right_value
    return output


def compose_series(
    outer: Sequence[Fraction], inner: Sequence[Fraction], order: int
) -> list[Fraction]:
    """Return ``outer(inner(t))`` through ``t**order`` for ``inner(0)=0``."""

    if not inner or inner[0] != 0:
        raise ValueError("inner series must have zero constant term")
    result = [Fraction(0) for _ in range(order + 1)]
    power = [Fraction(1)] + [Fraction(0) for _ in range(order)]
    for coefficient in outer[: order + 1]:
        for index in range(order + 1):
            result[index] += coefficient * power[index]
        power = multiply_series(power, inner, order)
    return result


def inverse_series(series: Sequence[Fraction], order: int) -> list[Fraction]:
    """Invert ``series(t)`` through ``order`` with exact arithmetic."""

    if len(series) < order + 1 or series[0] != 0:
        raise ValueError("coordinate series must contain orders 0..maximum with zero constant")
    if series[1] == 0:
        raise ValueError("coordinate first derivative must be nonzero")
    inverse = [Fraction(0) for _ in range(order + 1)]
    inverse[1] = 1 / series[1]
    for degree in range(2, order + 1):
        known = compose_series(series, inverse, degree)[degree]
        inverse[degree] = -known / series[1]
    identity = compose_series(series, inverse, order)
    expected = [Fraction(0), Fraction(1)] + [Fraction(0) for _ in range(order - 1)]
    if identity != expected:
        raise ArithmeticError("exact series inversion failed")
    return inverse


def compose_derivatives(
    outer_derivatives: Sequence[Fraction],
    inner_delta_derivatives: Sequence[Fraction],
    order: int,
) -> list[Fraction]:
    """Compose a jet with a centered inner-coordinate delta jet."""

    if len(outer_derivatives) < order + 1 or len(inner_delta_derivatives) < order + 1:
        raise ValueError("both derivative jets must cover the requested order")
    outer = derivatives_to_series(outer_derivatives[: order + 1])
    inner = derivatives_to_series(inner_delta_derivatives[: order + 1])
    return series_to_derivatives(compose_series(outer, inner, order))


def parametric_derivatives(
    matching_derivatives: Sequence[Fraction],
    reference_derivatives: Sequence[Fraction],
    order: int,
) -> list[Fraction]:
    """Eliminate the bare coordinate and return derivatives of ``M(U)``."""

    if len(matching_derivatives) < order + 1 or len(reference_derivatives) < order + 1:
        raise ValueError("both derivative jets must cover the requested order")
    matching = derivatives_to_series(matching_derivatives[: order + 1])
    reference = derivatives_to_series(reference_derivatives[: order + 1])
    reference[0] = Fraction(0)
    bare_delta = inverse_series(reference, order)
    return series_to_derivatives(compose_series(matching, bare_delta, order))


def normalized_odd_invariants(
    derivatives: Sequence[Fraction], orders: Sequence[int] = (3, 5)
) -> dict[str, Fraction]:
    if len(derivatives) < 2 or derivatives[1] == 0:
        raise ValueError("parametric first derivative must be nonzero")
    output: dict[str, Fraction] = {}
    for order in orders:
        if order < 3 or order % 2 == 0 or order >= len(derivatives):
            raise ValueError("normalized orders must be available odd orders >=3")
        output[str(order)] = derivatives[order] / derivatives[1] ** order
    return output

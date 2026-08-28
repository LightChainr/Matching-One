#!/usr/bin/env python3
"""Generate exact Gaussian-integer orientation genealogies and scale predictions."""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from fractions import Fraction
from typing import Dict, List, Tuple


@dataclass(frozen=True)
class Gaussian:
    a: int
    b: int

    @property
    def norm(self) -> int:
        return self.a * self.a + self.b * self.b

    @property
    def content(self) -> int:
        return math.gcd(abs(self.a), abs(self.b))

    def smith_invariants(self) -> Tuple[int, int]:
        if self.norm == 0:
            raise ValueError("zero Gaussian integer does not define a torus")
        d1 = self.content
        return d1, self.norm // d1

    @property
    def cyclic_translation_group(self) -> bool:
        return self.smith_invariants()[0] == 1

    def multiply(self, other: "Gaussian") -> "Gaussian":
        return Gaussian(
            self.a * other.a - self.b * other.b,
            self.a * other.b + self.b * other.a,
        )

    def canonical_d4(self) -> "Gaussian":
        first, second = abs(self.a), abs(self.b)
        if second > first:
            first, second = second, first
        return Gaussian(first, second)

    def cos4(self) -> Fraction:
        n = self.norm
        if n == 0:
            raise ValueError("zero Gaussian integer has no orientation")
        return Fraction(
            self.a**4 - 6 * self.a * self.a * self.b * self.b + self.b**4,
            n * n,
        )

    def sin4(self) -> Fraction:
        n = self.norm
        if n == 0:
            raise ValueError("zero Gaussian integer has no orientation")
        return Fraction(
            4 * self.a * self.b * (self.a * self.a - self.b * self.b),
            n * n,
        )

    def cos4m(self, m: int) -> Fraction:
        if m < 0:
            raise ValueError("harmonic multiplier m must be nonnegative")
        if m == 0:
            return Fraction(1)
        value = self.cos4()
        if m == 1:
            return value
        previous, current = Fraction(1), value
        for _degree in range(2, m + 1):
            previous, current = current, 2 * value * current - previous
        return current

    def as_pair(self) -> List[int]:
        return [self.a, self.b]


def fraction_payload(value: Fraction) -> Dict[str, object]:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": float(value),
    }


def orientation_payload(value: Gaussian) -> Dict[str, object]:
    smith = value.smith_invariants()
    return {
        "pair": value.as_pair(),
        "N": value.norm,
        "translation_group": {
            "smith_invariants": list(smith),
            "cyclic": smith[0] == 1,
        },
        "cos4": fraction_payload(value.cos4m(1)),
        "cos8": fraction_payload(value.cos4m(2)),
        "cos12": fraction_payload(value.cos4m(3)),
        "sin4_raw": fraction_payload(value.sin4()),
    }


def _same_translation_group(first: Gaussian, second: Gaussian, label: str) -> None:
    if first.smith_invariants() != second.smith_invariants():
        raise ValueError(
            "{} orientations have different finite translation groups: {} versus {}".format(
                label, first.smith_invariants(), second.smith_invariants()
            )
        )


def _target_expression(
    angular_ratio: Fraction,
    norm_ratio: int,
    alpha_num: int,
    alpha_den: int,
) -> str:
    return "({}/{})*{}^(-{}/{})".format(
        angular_ratio.numerator,
        angular_ratio.denominator,
        norm_ratio,
        alpha_num,
        alpha_den,
    )


def lineage_payload(
    first: Gaussian,
    second: Gaussian,
    multiplier: Gaussian,
    *,
    alpha_num: int = 13,
    alpha_den: int = 8,
) -> Dict[str, object]:
    if first.norm != second.norm:
        raise ValueError("the two parent orientations must have the same norm")
    _same_translation_group(first, second, "parent")
    if multiplier.norm == 0:
        raise ValueError("multiplier must be nonzero")
    if alpha_den <= 0:
        raise ValueError("alpha denominator must be positive")

    child_first_raw = first.multiply(multiplier)
    child_second_raw = second.multiply(multiplier)
    if child_first_raw.norm != child_second_raw.norm:
        raise AssertionError("Gaussian multiplication did not preserve paired norm")
    _same_translation_group(child_first_raw, child_second_raw, "child")

    radial_factor = multiplier.norm ** (-alpha_num / alpha_den)
    radial_expression = "{}^(-{}/{})".format(
        multiplier.norm, alpha_num, alpha_den
    )
    harmonic_predictions: Dict[str, object] = {}
    for m in (1, 2, 3):
        parent_delta = first.cos4m(m) - second.cos4m(m)
        child_delta = child_first_raw.cos4m(m) - child_second_raw.cos4m(m)
        if parent_delta == 0:
            harmonic_predictions["H{}".format(4 * m)] = {
                "parent_delta": fraction_payload(parent_delta),
                "child_delta": fraction_payload(child_delta),
                "angular_ratio": None,
                "target_expression": None,
                "target_delta_M_ratio": None,
            }
            continue
        angular_ratio = child_delta / parent_delta
        harmonic_predictions["H{}".format(4 * m)] = {
            "parent_delta": fraction_payload(parent_delta),
            "child_delta": fraction_payload(child_delta),
            "angular_ratio": fraction_payload(angular_ratio),
            "target_expression": _target_expression(
                angular_ratio,
                multiplier.norm,
                alpha_num,
                alpha_den,
            ),
            "target_delta_M_ratio": radial_factor * float(angular_ratio),
        }

    return {
        "parent": {
            "first": orientation_payload(first),
            "second": orientation_payload(second),
        },
        "multiplier": orientation_payload(multiplier),
        "child": {
            "first_raw": orientation_payload(child_first_raw),
            "second_raw": orientation_payload(child_second_raw),
            "first_canonical": orientation_payload(child_first_raw.canonical_d4()),
            "second_canonical": orientation_payload(child_second_raw.canonical_d4()),
        },
        "pair_translation_group_contract": {
            "parent_smith_invariants": list(first.smith_invariants()),
            "child_smith_invariants": list(child_first_raw.smith_invariants()),
            "parent_pair_matches": True,
            "child_pair_matches": True,
        },
        "radial_exponent_in_N": {
            "numerator": alpha_num,
            "denominator": alpha_den,
        },
        "norm_ratio": multiplier.norm,
        "radial_factor_expression": radial_expression,
        "harmonic_predictions": harmonic_predictions,
        "formula": "q^(-alpha) * DeltaCos(4m)_child / DeltaCos(4m)_parent",
    }


def default_catalog() -> Dict[str, object]:
    doubling = Gaussian(1, 1)
    doubling_lineages = {
        "65_to_130": lineage_payload(Gaussian(8, 1), Gaussian(7, 4), doubling),
        "85_to_170": lineage_payload(Gaussian(9, 2), Gaussian(7, 6), doubling),
        "145_to_290": lineage_payload(Gaussian(12, 1), Gaussian(9, 8), doubling),
    }
    norm5_lineages = {
        "65_to_325": lineage_payload(
            Gaussian(8, 1), Gaussian(7, 4), Gaussian(2, -1)
        ),
        "85_to_425": lineage_payload(
            Gaussian(9, 2), Gaussian(7, 6), Gaussian(2, 1)
        ),
    }
    commuting_square = [
        ("65_x_17_plus", Gaussian(8, 1), Gaussian(7, 4), Gaussian(4, 1)),
        ("65_x_17_minus", Gaussian(8, 1), Gaussian(7, 4), Gaussian(4, -1)),
        ("85_x_13_plus", Gaussian(9, 2), Gaussian(7, 6), Gaussian(3, 2)),
        ("85_x_13_minus", Gaussian(9, 2), Gaussian(7, 6), Gaussian(3, -2)),
        ("221_x_5_plus", Gaussian(14, 5), Gaussian(11, 10), Gaussian(2, 1)),
        ("221_x_5_minus", Gaussian(14, 5), Gaussian(11, 10), Gaussian(2, -1)),
    ]
    return {
        "schema_version": 3,
        "purpose": "outcome-free exact Gaussian semigroup design",
        "doubling_lineages": doubling_lineages,
        "norm5_harmonic_discrimination": norm5_lineages,
        "N1105_edges": {
            name: lineage_payload(first, second, multiplier)
            for name, first, second, multiplier in commuting_square
        },
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    catalog = subparsers.add_parser("catalog")
    catalog.add_argument("--output")
    lineage = subparsers.add_parser("lineage")
    lineage.add_argument("--first", nargs=2, type=int, required=True)
    lineage.add_argument("--second", nargs=2, type=int, required=True)
    lineage.add_argument("--multiplier", nargs=2, type=int, required=True)
    lineage.add_argument("--alpha", nargs=2, type=int, default=(13, 8))
    lineage.add_argument("--output")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.command == "catalog":
        payload = default_catalog()
    else:
        payload = lineage_payload(
            Gaussian(*args.first),
            Gaussian(*args.second),
            Gaussian(*args.multiplier),
            alpha_num=args.alpha[0],
            alpha_den=args.alpha[1],
        )
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output:
        from pathlib import Path
        path = Path(args.output)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
    else:
        print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Cluster-sign-integrated sufficient statistics for Issue #234.

The Camia--Feng spin is zero on a white site and is an independent symmetric
sign shared by every site of one black cluster.  Signs are integrated out
analytically.  A tiny 6x2 periodic triangular quotient supplies an exhaustive
4096-configuration oracle; it is an algebra/convention test, not a scaling
experiment.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import json
from pathlib import Path
from typing import Iterable, Mapping, Sequence


Pair = tuple[int, int]
Coordinate = tuple[int, int]
MOMENT_ORDER = ("LL", "LD", "DD")
PAIR_ORDER = ("L1", "L2", "D1", "D2")
PRODUCT_ORDER = ("LL", "L1D2", "D1L2", "DD")


class DisjointSet:
    def __init__(self, size: int) -> None:
        self.parent = list(range(size))
        self.rank = [0] * size

    def find(self, item: int) -> int:
        root = item
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[item] != item:
            parent = self.parent[item]
            self.parent[item] = root
            item = parent
        return root

    def union(self, first: int, second: int) -> None:
        first = self.find(first)
        second = self.find(second)
        if first == second:
            return
        if self.rank[first] < self.rank[second]:
            first, second = second, first
        self.parent[second] = first
        if self.rank[first] == self.rank[second]:
            self.rank[first] += 1


def vertex_id(width: int, height: int, coordinate: Coordinate) -> int:
    x, y = coordinate
    return (y % height) * width + (x % width)


def triangular_edges(width: int, height: int) -> tuple[Pair, ...]:
    if width < 2 or height < 2:
        raise ValueError("periodic triangular quotient must be at least 2x2")
    edges = set()
    for y in range(height):
        for x in range(width):
            first = vertex_id(width, height, (x, y))
            for dx, dy in ((1, 0), (0, 1), (1, -1)):
                second = vertex_id(width, height, (x + dx, y + dy))
                if first != second:
                    edges.add(tuple(sorted((first, second))))
    return tuple(sorted(edges))


def black_cluster_roots(
    width: int, height: int, mask: int, edges: Sequence[Pair] | None = None
) -> list[int | None]:
    count = width * height
    if mask < 0 or mask >= 1 << count:
        raise ValueError("occupation mask is outside the quotient")
    union_find = DisjointSet(count)
    for first, second in edges or triangular_edges(width, height):
        if (mask >> first) & 1 and (mask >> second) & 1:
            union_find.union(first, second)
    return [union_find.find(site) if (mask >> site) & 1 else None for site in range(count)]


def sign_integrated_pair(roots: Sequence[int | None], pair: Pair) -> int:
    """E_sigma[S_x S_y | omega]: black-connectivity indicator."""

    first, second = pair
    root = roots[first]
    return int(root is not None and root == roots[second])


def sign_integrated_four_spin(
    roots: Sequence[int | None], first_pair: Pair, second_pair: Pair
) -> int:
    """E_sigma[prod of four endpoint spins | omega].

    It is one exactly when all endpoints are black and every touched black
    cluster occurs with even multiplicity among the four endpoint slots.
    Endpoint slots, rather than distinct sites, are counted so the identity
    also covers coincident endpoints.
    """

    parity: dict[int, int] = {}
    for site in (*first_pair, *second_pair):
        root = roots[site]
        if root is None:
            return 0
        parity[root] = parity.get(root, 0) ^ 1
    return int(not any(parity.values()))


def brute_force_sign_average(
    roots: Sequence[int | None], first_pair: Pair, second_pair: Pair
) -> Fraction:
    """Tiny independent oracle that explicitly enumerates touched signs."""

    endpoint_roots = [roots[site] for site in (*first_pair, *second_pair)]
    if any(root is None for root in endpoint_roots):
        return Fraction(0)
    distinct = sorted(set(endpoint_roots))  # type: ignore[arg-type]
    total = 0
    for sign_mask in range(1 << len(distinct)):
        signs = {
            root: (1 if (sign_mask >> index) & 1 else -1)
            for index, root in enumerate(distinct)
        }
        product = 1
        for root in endpoint_roots:
            product *= signs[root]  # type: ignore[index]
        total += product
    return Fraction(total, 1 << len(distinct))


@dataclass
class SufficientSums:
    samples: int
    pair_sums: dict[str, int]
    within_products: dict[str, int]
    four_spin_sums: dict[str, int]

    @classmethod
    def empty(cls) -> "SufficientSums":
        return cls(
            samples=0,
            pair_sums={name: 0 for name in PAIR_ORDER},
            within_products={name: 0 for name in PRODUCT_ORDER},
            four_spin_sums={name: 0 for name in PRODUCT_ORDER},
        )

    def add(self, pair_values: Mapping[str, int], four_values: Mapping[str, int]) -> None:
        self.samples += 1
        for name in PAIR_ORDER:
            self.pair_sums[name] += int(pair_values[name])
        products = {
            "LL": pair_values["L1"] * pair_values["L2"],
            "L1D2": pair_values["L1"] * pair_values["D2"],
            "D1L2": pair_values["D1"] * pair_values["L2"],
            "DD": pair_values["D1"] * pair_values["D2"],
        }
        for name in PRODUCT_ORDER:
            self.within_products[name] += int(products[name])
            self.four_spin_sums[name] += int(four_values[name])


def configuration_values(
    roots: Sequence[int | None], pairs: Mapping[str, Pair]
) -> tuple[dict[str, int], dict[str, int]]:
    pair_values = {
        name: sign_integrated_pair(roots, pairs[name]) for name in PAIR_ORDER
    }
    four_values = {
        "LL": sign_integrated_four_spin(roots, pairs["L1"], pairs["L2"]),
        "L1D2": sign_integrated_four_spin(roots, pairs["L1"], pairs["D2"]),
        "D1L2": sign_integrated_four_spin(roots, pairs["D1"], pairs["L2"]),
        "DD": sign_integrated_four_spin(roots, pairs["D1"], pairs["D2"]),
    }
    return pair_values, four_values


def translation_averaged_configuration_sums(
    roots: Sequence[int | None],
    width: int,
    height: int,
    delta_radius: int,
    center_displacement: Coordinate,
) -> tuple[dict[str, int], dict[str, int]]:
    """Integer spatial sums for one production configuration.

    The returned values are sums over all ``width*height`` translations.
    Keeping them integral lets a production engine archive exact sufficient
    statistics; pass that translation count as ``placements`` when forming
    the U-statistic.
    """

    if delta_radius <= 1 or 2 * delta_radius >= width:
        raise ValueError("bilocal radius must exceed one and fit the x period")
    pair_sums = {name: 0 for name in PAIR_ORDER}
    four_sums = {name: 0 for name in PRODUCT_ORDER}

    def axis_pair(center: Coordinate, radius: int) -> Pair:
        x, y = center
        return (
            vertex_id(width, height, (x - radius, y)),
            vertex_id(width, height, (x + radius, y)),
        )

    dx, dy = center_displacement
    for y in range(height):
        for x in range(width):
            first_center = (x, y)
            second_center = (x + dx, y + dy)
            pairs = {
                "L1": axis_pair(first_center, 1),
                "L2": axis_pair(second_center, 1),
                "D1": axis_pair(first_center, delta_radius),
                "D2": axis_pair(second_center, delta_radius),
            }
            pair_values, four_values = configuration_values(roots, pairs)
            for name in PAIR_ORDER:
                pair_sums[name] += pair_values[name]
            for name in PRODUCT_ORDER:
                four_sums[name] += four_values[name]
    return pair_sums, four_sums


def _unbiased_centered(
    sums: SufficientSums,
    first: str,
    second: str,
    product: str,
    placements: int,
) -> Fraction:
    """Unbiased Eh_AB-Eu_A Eu_B using a cross-configuration U-statistic."""

    count = sums.samples
    if count < 2:
        raise ValueError("at least two configurations are required")
    if placements <= 0:
        raise ValueError("placements must be positive")
    raw = Fraction(sums.four_spin_sums[product], count * placements)
    cross_product = Fraction(
        sums.pair_sums[first] * sums.pair_sums[second]
        - sums.within_products[product],
        count * (count - 1) * placements * placements,
    )
    return raw - cross_product


def unbiased_moment_vector(
    sums: SufficientSums, *, placements: int = 1
) -> tuple[Fraction, ...]:
    ll = _unbiased_centered(sums, "L1", "L2", "LL", placements)
    ld = (
        _unbiased_centered(sums, "L1", "D2", "L1D2", placements)
        + _unbiased_centered(sums, "D1", "L2", "D1L2", placements)
    ) / 2
    dd = _unbiased_centered(sums, "D1", "D2", "DD", placements)
    return ll, ld, dd


def covariance_of_mean(rows: Sequence[Sequence[float]]) -> list[list[float]]:
    """Covariance of a mean across independent block-level U-statistics."""

    count = len(rows)
    if count < 2 or any(len(row) != 3 for row in rows):
        raise ValueError("need at least two three-coordinate block estimates")
    means = [sum(row[index] for row in rows) / count for index in range(3)]
    return [
        [
            sum(
                (row[first] - means[first]) * (row[second] - means[second])
                for row in rows
            )
            / (count * (count - 1))
            for second in range(3)
        ]
        for first in range(3)
    ]


def tiny_oracle_pairs(width: int = 6, height: int = 3) -> dict[str, Pair]:
    centers = ((0, 0), (3, 1))

    def axis_pair(center: Coordinate, radius: int) -> Pair:
        x, y = center
        return (
            vertex_id(width, height, (x - radius, y)),
            vertex_id(width, height, (x + radius, y)),
        )

    return {
        "L1": axis_pair(centers[0], 1),
        "L2": axis_pair(centers[1], 1),
        "D1": axis_pair(centers[0], 2),
        "D2": axis_pair(centers[1], 2),
    }


def _fraction_payload(value: Fraction) -> dict[str, object]:
    return {
        "fraction": f"{value.numerator}/{value.denominator}",
        "float": float(value),
    }


def exact_oracle(width: int = 6, height: int = 3) -> dict[str, object]:
    count = width * height
    if count > 20:
        raise ValueError("tiny exact oracle is capped at 20 sites")
    edges = triangular_edges(width, height)
    pairs = tiny_oracle_pairs(width, height)
    sums = SufficientSums.empty()
    explicit_sign_checks = 0
    for mask in range(1 << count):
        roots = black_cluster_roots(width, height, mask, edges)
        pair_values, four_values = configuration_values(roots, pairs)
        sums.add(pair_values, four_values)
        for product, first, second in (
            ("LL", "L1", "L2"),
            ("L1D2", "L1", "D2"),
            ("D1L2", "D1", "L2"),
            ("DD", "D1", "D2"),
        ):
            brute = brute_force_sign_average(roots, pairs[first], pairs[second])
            if brute != four_values[product]:
                raise AssertionError("analytic cluster-sign integration failed")
            explicit_sign_checks += 1

    configurations = 1 << count
    pair_means = {
        name: Fraction(sums.pair_sums[name], configurations) for name in PAIR_ORDER
    }
    raw_means = {
        name: Fraction(sums.four_spin_sums[name], configurations)
        for name in PRODUCT_ORDER
    }
    exact_centered = (
        raw_means["LL"] - pair_means["L1"] * pair_means["L2"],
        (
            raw_means["L1D2"] - pair_means["L1"] * pair_means["D2"]
            + raw_means["D1L2"] - pair_means["D1"] * pair_means["L2"]
        )
        / 2,
        raw_means["DD"] - pair_means["D1"] * pair_means["D2"],
    )
    return {
        "schema": "matching-one/p234-triangular-energy-logpair-tiny-oracle/v1",
        "issue": 234,
        "role": "algebra_and_periodic_geometry_oracle_not_scaling_evidence",
        "geometry": {
            "width": width,
            "height": height,
            "sites": count,
            "edges_after_periodic_deduplication": len(edges),
            "neighbor_steps": [[1, 0], [0, 1], [1, -1]],
            "configurations": configurations,
            "p": "1/2",
        },
        "pairs": {name: list(pair) for name, pair in pairs.items()},
        "pair_connection_means": {
            name: _fraction_payload(value) for name, value in pair_means.items()
        },
        "four_spin_means": {
            name: _fraction_payload(value) for name, value in raw_means.items()
        },
        "centered_moment_order": list(MOMENT_ORDER),
        "exact_centered_moments": [
            _fraction_payload(value) for value in exact_centered
        ],
        "explicit_cluster_sign_checks": explicit_sign_checks,
        "sign_integration_passed": True,
        "u_statistic_identity": (
            "mean(h_AB)-[sum(u_A)sum(u_B)-sum(u_Au_B)]/[n(n-1)] "
            "is unbiased for E[h_AB]-E[u_A]E[u_B] under iid configurations"
        ),
        "boundary": (
            "This small quotient has deliberately non-scaling local and bilocal "
            "radii; its rational values must not be interpreted as continuum amplitudes."
        ),
    }


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args(list(argv) if argv is not None else None)
    payload = exact_oracle()
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(text, encoding="utf-8")
    else:
        print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

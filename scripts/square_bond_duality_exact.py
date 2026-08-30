#!/usr/bin/env python3
"""Exact tiny-torus conventions for the square-bond self-dual control.

Issue #42 asks whether the two-sector matching-even / matching-odd picture
survives an exactly self-dual square model. The primal/dual wrapping channels
and even/odd combinations are fixed here by exhaustive identities, so future
production or exploratory runs can share one convention.

On an L x L square torus, each of the 2 L^2 primal bonds is paired with the
unique dual bond that crosses it. Under the vertex identification already
used by ``square_bond_pairs``, every dual edge is equal as a displaced edge
to some primal edge. Geometric dual transport is the bijection

    T(mask)_j = 1  iff  the dual edge of a vacant primal bond i equals
                        the primal edge of bond j.

Naive bit-complement of the same mask is **not** this map: it occupies the
vacant primal edges rather than the identified dual edges, so it does not
swap primal and dual wrapping channels.

Define, for each wrapping channel,

    S = (R_primal + R_dual) / 2
    D =  R_primal - R_dual

T is a bijection of configurations, sends D -> -D, and leaves S invariant.
At the exact self-dual point p = 1/2 every configuration is equiprobable,
so

    E[D(p=1/2)] = 0

identically, for every finite L and every channel. Equivalently, primal and
dual wrapping are equal in law because the dual graph is the same square
grid and dual occupation is i.i.d. Bernoulli(1/2).

This is a C5 finite identity, not a large-N amplitude test. It does **not**
measure the orientation projector P4[S] ~ L^{-2} or the derivative sector
P4[D']. Those remain empirical questions for #42.

This script enumerates L=2 (8 bonds) and L=3 (18 bonds) and records exact
rational channel expectations. It is a convention/oracle freeze.
"""
from __future__ import annotations

import argparse
import json
from fractions import Fraction
from typing import Dict, Iterable, List, Sequence, Tuple

from square_bond_kappa3 import BondPair, square_bond_pairs
from torus_homology import (
    ComponentHomology,
    HomologyUnionFind,
    WrappingChannels,
    wrapping_channels,
)


CHANNELS: Tuple[str, ...] = (
    "cross",
    "both",
    "either",
    "direction_0",
    "direction_1",
)


def wrapping_from_union_find(union_find: HomologyUnionFind) -> WrappingChannels:
    components: List[ComponentHomology] = []
    for vertex in range(len(union_find.parent)):
        root, _, _ = union_find.find(vertex)
        if root != vertex:
            continue
        components.append(
            ComponentHomology(root, union_find.size[root], tuple(union_find.basis[root]))
        )
    return wrapping_channels(components)


def primal_dual_wrapping(
    length: int, mask: int, pairs: Sequence[BondPair]
) -> Tuple[WrappingChannels, WrappingChannels]:
    vertex_count = length * length
    primal = HomologyUnionFind(vertex_count, (length, length))
    dual = HomologyUnionFind(vertex_count, (length, length))
    for bond_index, pair in enumerate(pairs):
        if (mask >> bond_index) & 1:
            primal.add_edge(*pair.primal)
        else:
            dual.add_edge(*pair.dual)
    return wrapping_from_union_find(primal), wrapping_from_union_find(dual)


def geometric_dual_mask(mask: int, pairs: Sequence[BondPair]) -> int:
    primal_index = {pair.primal: index for index, pair in enumerate(pairs)}
    transported = 0
    for bond_index, pair in enumerate(pairs):
        if (mask >> bond_index) & 1:
            continue
        transported |= 1 << primal_index[pair.dual]
    return transported


def channel_indicators(channels: WrappingChannels) -> Dict[str, int]:
    payload = channels.as_dict()
    return {name: int(bool(payload[name])) for name in CHANNELS}


def enumerate_exact(length: int) -> dict[str, object]:
    if length < 2:
        raise ValueError("L must be at least 2")
    pairs = square_bond_pairs(length)
    bond_count = len(pairs)
    if bond_count > 24:
        raise ValueError(
            f"exact enumeration of L={length} requires 2^{bond_count} states"
        )

    configuration_count = 1 << bond_count
    full_mask = configuration_count - 1
    primal_sum = {name: 0 for name in CHANNELS}
    dual_sum = {name: 0 for name in CHANNELS}
    even_sum = {name: 0 for name in CHANNELS}
    odd_sum = {name: 0 for name in CHANNELS}
    naive_complement_failures = 0
    dual_transport_failures = 0
    even_odd_failures = 0
    transported_images = set()

    for mask in range(configuration_count):
        primal, dual = primal_dual_wrapping(length, mask, pairs)
        complement_primal, complement_dual = primal_dual_wrapping(
            length, full_mask ^ mask, pairs
        )
        transported = geometric_dual_mask(mask, pairs)
        transported_images.add(transported)
        transported_primal, transported_dual = primal_dual_wrapping(
            length, transported, pairs
        )
        primal_bits = channel_indicators(primal)
        dual_bits = channel_indicators(dual)
        complement_primal_bits = channel_indicators(complement_primal)
        complement_dual_bits = channel_indicators(complement_dual)
        transported_primal_bits = channel_indicators(transported_primal)
        transported_dual_bits = channel_indicators(transported_dual)
        if primal_bits != complement_dual_bits or dual_bits != complement_primal_bits:
            naive_complement_failures += 1
        if dual_bits != transported_primal_bits or primal_bits != transported_dual_bits:
            dual_transport_failures += 1
        for name in CHANNELS:
            r_p = primal_bits[name]
            r_d = dual_bits[name]
            even = r_p + r_d
            odd = r_p - r_d
            transported_even = transported_primal_bits[name] + transported_dual_bits[name]
            transported_odd = transported_primal_bits[name] - transported_dual_bits[name]
            if even != transported_even or odd != -transported_odd:
                even_odd_failures += 1
            primal_sum[name] += r_p
            dual_sum[name] += r_d
            even_sum[name] += even
            odd_sum[name] += odd

    channels: dict[str, object] = {}
    denominator = configuration_count
    for name in CHANNELS:
        odd_total = odd_sum[name]
        even_total = even_sum[name]
        channels[name] = {
            "E_R_primal": _fraction_payload(Fraction(primal_sum[name], denominator)),
            "E_R_dual": _fraction_payload(Fraction(dual_sum[name], denominator)),
            "E_S": _fraction_payload(Fraction(even_total, 2 * denominator)),
            "E_D": _fraction_payload(Fraction(odd_total, denominator)),
            "D_vanishes_at_half": odd_total == 0,
            "S_identically_zero": even_total == 0,
            "primal_dual_equidistributed": primal_sum[name] == dual_sum[name],
        }

    passed = (
        dual_transport_failures == 0
        and even_odd_failures == 0
        and len(transported_images) == configuration_count
        and all(row["D_vanishes_at_half"] for row in channels.values())  # type: ignore[index]
        and all(row["primal_dual_equidistributed"] for row in channels.values())  # type: ignore[index]
    )
    return {
        "model": "square_bond_square_torus",
        "L": length,
        "N_vertices": length * length,
        "N_bonds": bond_count,
        "configurations": configuration_count,
        "p": "1/2",
        "naive_complement_swap_failures": naive_complement_failures,
        "geometric_dual_transport_failures": dual_transport_failures,
        "even_odd_involution_failures": even_odd_failures,
        "geometric_dual_map_bijective": len(transported_images) == configuration_count,
        "channels": channels,
        "passed": passed,
        "scientific_boundary": (
            "Exact geometric dual-transport identity and E[D](p=1/2)=0. "
            "Naive bit-complement is not the duality map. "
            "Does not measure orientation amplitudes or L^{-13/4} scaling."
        ),
    }


def _fraction_payload(value: Fraction) -> dict[str, int | str]:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "fraction": str(value),
    }


def json_ready(value: object) -> object:
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    return value


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--L", type=int, default=2)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args(list(argv) if argv is not None else None)
    result = enumerate_exact(args.L)
    if args.json:
        print(json.dumps(json_ready(result), indent=2, sort_keys=True))
    else:
        print(
            f"L={result['L']} bonds={result['N_bonds']} "
            f"configs={result['configurations']} passed={result['passed']}"
        )
        print(
            "  naive_complement_failures="
            f"{result['naive_complement_swap_failures']} "
            "geometric_dual_failures="
            f"{result['geometric_dual_transport_failures']} "
            f"T_bijective={result['geometric_dual_map_bijective']}"
        )
        channels = result["channels"]
        assert isinstance(channels, dict)
        for name in CHANNELS:
            row = channels[name]
            print(
                f"  {name:12} E[S]={row['E_S']['fraction']:>12} "
                f"E[D]={row['E_D']['fraction']}"
            )
    if not result["passed"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

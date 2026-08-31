#!/usr/bin/env python3
"""Exact left-right reversal anti-involution for typed terminal partitions."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

try:
    from scripts.terminal_partition_canonical import RGS, apply_permutation, enumerate_rgs, validate_rgs
    from scripts.terminal_partition_serial_category import serial_compose
except ModuleNotFoundError:
    from terminal_partition_canonical import RGS, apply_permutation, enumerate_rgs, validate_rgs
    from terminal_partition_serial_category import serial_compose


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "analysis" / "terminal_partition_serial_reversal_certificate.json"
SCHEMA = "matching-one/terminal-partition-serial-reversal/v1"
TERMINAL_ORDER = ("L0", "L1", "R0", "R1")
PORT_REVERSAL = (2, 3, 0, 1)


def reverse_ports(state: Sequence[int]) -> RGS:
    """Swap the typed left and right port pairs without swapping lanes."""

    return apply_permutation(validate_rgs(state, 4), PORT_REVERSAL)


def reverse_index_vector(values: Sequence[Any], partitions: Sequence[RGS]) -> tuple[Any, ...]:
    """Transport a state-indexed vector through port reversal."""

    states = tuple(validate_rgs(state, 4) for state in partitions)
    if len(states) != 15 or len(set(states)) != 15 or set(states) != set(enumerate_rgs(4)):
        raise ValueError("partition catalog must contain every four-terminal state exactly once")
    if len(values) != len(states):
        raise ValueError("state vector length does not match partition catalog")
    index = {state: position for position, state in enumerate(states)}
    return tuple(values[index[reverse_ports(state)]] for state in states)


def build_artifact() -> dict[str, Any]:
    partitions = enumerate_rgs(4)
    index = {state: position for position, state in enumerate(partitions)}
    reversal = [index[reverse_ports(state)] for state in partitions]
    involution_failures = sum(reversal[reversal[position]] != position for position in range(15))
    anti_failures = 0
    table_failures = 0
    for left in partitions:
        for right in partitions:
            expected = reverse_ports(serial_compose(left, right))
            observed = serial_compose(reverse_ports(right), reverse_ports(left))
            anti_failures += expected != observed
            table_failures += index[expected] != index[observed]

    fixed = [position for position, target in enumerate(reversal) if position == target]
    two_cycles = [
        [position, target]
        for position, target in enumerate(reversal)
        if position < target
    ]
    identity = (0, 1, 0, 1)
    return {
        "schema": SCHEMA,
        "issue": 13,
        "status": "exact_typed_serial_reversal_anti_involution",
        "terminal_order": list(TERMINAL_ORDER),
        "permutation_convention": "old_terminal_to_new_terminal",
        "port_reversal_permutation": list(PORT_REVERSAL),
        "partition_catalog": [list(state) for state in partitions],
        "reversal_index_map": reversal,
        "orbit_decomposition": {
            "fixed_indices": fixed,
            "fixed_partitions": [list(partitions[position]) for position in fixed],
            "two_cycles": two_cycles,
            "fixed_count": len(fixed),
            "two_cycle_count": len(two_cycles),
        },
        "exhaustive_pairs": {
            "cases": len(partitions) ** 2,
            "anti_homomorphism_failures": anti_failures,
            "cayley_table_anti_automorphism_failures": table_failures,
        },
        "exact_checks": {
            "reversal_is_involutive_on_all_states": involution_failures == 0,
            "all_225_products_reverse_order": anti_failures == 0,
            "full_cayley_table_has_anti_automorphism": table_failures == 0,
            "wire_identity_is_fixed": reverse_ports(identity) == identity,
            "orbits_partition_all_states": len(fixed) + 2 * len(two_cycles) == len(partitions),
        },
        "claim_boundary": {
            "included": "typed left-right port reversal, all-state involution, full serial-table anti-automorphism, and orbit classification",
            "excluded": "planar duality, complement duality, periodic gluing, reliability polynomials, thresholds, or bounds",
            "parent_issue": "remain open",
        },
    }


def validate_artifact(artifact: Mapping[str, Any]) -> dict[str, Any]:
    expected = build_artifact()
    if artifact != expected:
        raise ValueError("serial-reversal artifact does not exactly reproduce")
    if set(expected["exact_checks"].values()) != {True}:
        raise ValueError("all exact checks must pass")
    return {
        "schema": SCHEMA,
        "status": "valid_exact_typed_serial_reversal",
        "states": len(expected["partition_catalog"]),
        "pairs": expected["exhaustive_pairs"]["cases"],
        "fixed_states": expected["orbit_decomposition"]["fixed_count"],
        "two_cycles": expected["orbit_decomposition"]["two_cycle_count"],
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--validate", type=Path)
    args = parser.parse_args(argv)
    if args.validate:
        print(json.dumps(validate_artifact(json.loads(args.validate.read_text())), indent=2, sort_keys=True))
        return 0
    rendered = json.dumps(build_artifact(), indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

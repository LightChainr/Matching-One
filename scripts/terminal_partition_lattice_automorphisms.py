#!/usr/bin/env python3
"""Enumerate every automorphism of the labelled four-terminal partition lattice."""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import permutations
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

try:
    from scripts.terminal_partition_canonical import (
        RGS,
        apply_permutation,
        enumerate_rgs,
        full_symmetric_group,
        rgs_to_blocks,
    )
    from scripts.terminal_partition_gluing_algebra import partition_join
except ModuleNotFoundError:
    from terminal_partition_canonical import RGS, apply_permutation, enumerate_rgs, full_symmetric_group, rgs_to_blocks
    from terminal_partition_gluing_algebra import partition_join


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "analysis" / "terminal_partition_lattice_automorphisms_certificate.json"
SCHEMA = "matching-one/terminal-partition-lattice-automorphisms/v1"


def block_count(partition: Sequence[int]) -> int:
    return len(rgs_to_blocks(partition))


def atom_catalog(partitions: Sequence[RGS]) -> tuple[RGS, ...]:
    atoms = tuple(value for value in partitions if block_count(value) == 3)
    if len(atoms) != 6:
        raise ValueError("four-terminal lattice must have six atoms")
    return atoms


def atom_support(partition: Sequence[int], atoms: Sequence[RGS]) -> tuple[int, ...]:
    blocks = tuple(set(block) for block in rgs_to_blocks(partition))
    result = []
    for index, atom in enumerate(atoms):
        pair = next(set(block) for block in rgs_to_blocks(atom) if len(block) == 2)
        if any(pair.issubset(block) for block in blocks):
            result.append(index)
    return tuple(result)


def join_atoms(indices: Sequence[int], atoms: Sequence[RGS]) -> RGS:
    result = (0, 1, 2, 3)
    for index in indices:
        if isinstance(index, bool) or not isinstance(index, int) or not 0 <= index < len(atoms):
            raise ValueError("atom index out of range")
        result = partition_join(result, atoms[index])
    return result


def induced_map_from_atom_permutation(
    atom_permutation: Sequence[int], partitions: Sequence[RGS], atoms: Sequence[RGS]
) -> tuple[int, ...] | None:
    atom_permutation = tuple(atom_permutation)
    if len(atom_permutation) != len(atoms) or set(atom_permutation) != set(range(len(atoms))):
        raise ValueError("atom permutation must be a bijection")
    index = {partition: position for position, partition in enumerate(partitions)}
    mapping = tuple(
        index[join_atoms(tuple(atom_permutation[item] for item in atom_support(partition, atoms)), atoms)]
        for partition in partitions
    )
    if len(set(mapping)) != len(partitions):
        return None
    for left_index, left in enumerate(partitions):
        for right_index, right in enumerate(partitions):
            joined = index[partition_join(left, right)]
            image_joined = index[partition_join(partitions[mapping[left_index]], partitions[mapping[right_index]])]
            if mapping[joined] != image_joined:
                return None
    return mapping


def enumerate_lattice_automorphisms(partitions: Sequence[RGS]) -> tuple[tuple[int, ...], ...]:
    atoms = atom_catalog(partitions)
    accepted = []
    for atom_permutation in permutations(range(len(atoms))):
        mapping = induced_map_from_atom_permutation(atom_permutation, partitions, atoms)
        if mapping is not None:
            accepted.append(mapping)
    return tuple(sorted(set(accepted)))


def terminal_relabeling_maps(partitions: Sequence[RGS]) -> tuple[tuple[int, ...], ...]:
    index = {partition: position for position, partition in enumerate(partitions)}
    return tuple(
        sorted(
            {
                tuple(index[apply_permutation(partition, permutation)] for partition in partitions)
                for permutation in full_symmetric_group(4)
            }
        )
    )


def orbit_catalog(maps: Sequence[Sequence[int]], size: int) -> list[list[int]]:
    unseen = set(range(size))
    orbits = []
    while unseen:
        start = min(unseen)
        orbit = sorted({mapping[start] for mapping in maps})
        unseen.difference_update(orbit)
        orbits.append(orbit)
    return orbits


def build_artifact() -> dict[str, Any]:
    partitions = enumerate_rgs(4)
    atoms = atom_catalog(partitions)
    lattice_maps = enumerate_lattice_automorphisms(partitions)
    relabeling_maps = terminal_relabeling_maps(partitions)
    orbits = orbit_catalog(lattice_maps, len(partitions))
    orbit_rows = []
    for orbit in orbits:
        representative = partitions[orbit[0]]
        orbit_rows.append(
            {
                "representative": list(representative),
                "block_sizes": sorted((len(block) for block in rgs_to_blocks(representative)), reverse=True),
                "orbit_size": len(orbit),
                "stabilizer_size": len(lattice_maps) // len(orbit),
                "members": [list(partitions[index]) for index in orbit],
            }
        )
    atom_candidate_count = 1
    for value in range(2, len(atoms) + 1):
        atom_candidate_count *= value
    return {
        "schema": SCHEMA,
        "issue": 13,
        "status": "exact_four_terminal_partition_lattice_automorphisms",
        "partition_catalog": [list(value) for value in partitions],
        "atom_catalog": [list(value) for value in atoms],
        "enumeration": {
            "atom_permutations_checked": atom_candidate_count,
            "lattice_automorphisms": len(lattice_maps),
            "terminal_relabeling_maps": len(relabeling_maps),
            "rejected_atom_permutations": atom_candidate_count - len(lattice_maps),
        },
        "automorphism_maps": [list(mapping) for mapping in lattice_maps],
        "partition_type_orbits": orbit_rows,
        "orbit_size_histogram": {
            str(size): count for size, count in sorted(Counter(len(orbit) for orbit in orbits).items())
        },
        "exact_checks": {
            "all_lattice_automorphisms_are_terminal_relabelings": lattice_maps == relabeling_maps,
            "automorphism_group_has_order_24": len(lattice_maps) == 24,
            "all_720_atom_permutations_classified": atom_candidate_count == 720,
            "five_integer_partition_type_orbits": len(orbits) == 5,
            "orbit_stabilizer_products_are_24": all(
                row["orbit_size"] * row["stabilizer_size"] == 24 for row in orbit_rows
            ),
        },
        "claim_boundary": {
            "included": "automorphism group of the labelled four-terminal partition lattice",
            "excluded": "D4 embedding choice, noncrossing states, planar duality, reliability, thresholds, or bounds",
            "parent_issue": "remain open",
        },
    }


def validate_artifact(artifact: Mapping[str, Any]) -> dict[str, Any]:
    expected = build_artifact()
    if artifact != expected:
        raise ValueError("partition-automorphism artifact does not exactly reproduce")
    if set(expected["exact_checks"].values()) != {True}:
        raise ValueError("all exact checks must pass")
    return {
        "schema": SCHEMA,
        "status": "valid_exact_four_terminal_partition_lattice_automorphisms",
        "atom_permutations_checked": 720,
        "automorphisms": 24,
        "partition_type_orbits": 5,
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

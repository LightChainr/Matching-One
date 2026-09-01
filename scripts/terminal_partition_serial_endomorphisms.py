#!/usr/bin/env python3
"""Complete unit-preserving endomorphism census of the typed serial monoid."""

from __future__ import annotations

import argparse
from collections import Counter, deque
from itertools import combinations, product
import json
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

try:
    from scripts.terminal_partition_canonical import enumerate_rgs
    from scripts.terminal_partition_serial_category import serial_compose
except ModuleNotFoundError:
    from terminal_partition_canonical import enumerate_rgs
    from terminal_partition_serial_category import serial_compose

SCHEMA = "matching-one/terminal-partition-serial-endomorphisms/v1"
IDENTITY = 6


def multiplication_table() -> tuple[tuple[int, ...], ...]:
    states = enumerate_rgs(4)
    index = {state: i for i, state in enumerate(states)}
    return tuple(tuple(index[serial_compose(a, b)] for b in states) for a in states)


def closure(seed: Sequence[int], table: Sequence[Sequence[int]]) -> frozenset[int]:
    result = {IDENTITY, *seed}
    changed = True
    while changed:
        changed = False
        for a in tuple(result):
            for b in tuple(result):
                if table[a][b] not in result:
                    result.add(table[a][b])
                    changed = True
    return frozenset(result)


def word_representatives(seed: Sequence[int], table: Sequence[Sequence[int]]) -> dict[int, tuple[int, ...]]:
    words = {IDENTITY: ()}
    queue = deque([IDENTITY])
    while queue:
        source = queue.popleft()
        for position, generator in enumerate(seed):
            target = table[source][generator]
            if target not in words:
                words[target] = words[source] + (position,)
                queue.append(target)
    if len(words) != len(table):
        raise ValueError("seed does not generate every state")
    return words


def evaluate_word(word: Sequence[int], images: Sequence[int], table: Sequence[Sequence[int]]) -> int:
    value = IDENTITY
    for position in word:
        value = table[value][images[position]]
    return value


def compose_maps(left: Sequence[int], right: Sequence[int]) -> tuple[int, ...]:
    if len(left) != len(right):
        raise ValueError("map widths differ")
    return tuple(left[right[i]] for i in range(len(left)))


def kernel_labels(mapping: Sequence[int]) -> list[int]:
    labels: dict[int, int] = {}
    return [labels.setdefault(value, len(labels)) for value in mapping]


def build_artifact() -> dict[str, Any]:
    states = enumerate_rgs(4)
    table = multiplication_table()
    generating_sets = [seed for seed in combinations(range(15), 3) if closure(seed, table) == frozenset(range(15))]
    base = generating_sets[0]
    words = word_representatives(base, table)
    endomorphisms = set()
    for images in product(range(15), repeat=3):
        mapping = tuple(evaluate_word(words[i], images, table) for i in range(15))
        if all(mapping[table[a][b]] == table[mapping[a]][mapping[b]] for a in range(15) for b in range(15)):
            endomorphisms.add(mapping)
    endomorphisms = sorted(endomorphisms)
    automorphisms = [value for value in endomorphisms if len(set(value)) == 15]
    idempotents = [value for value in endomorphisms if compose_maps(value, value) == value]
    unseen = set(endomorphisms)
    orbits = []
    while unseen:
        value = min(unseen)
        orbit = {
            compose_maps(compose_maps(auto, value), auto)
            for auto in automorphisms
        }
        unseen -= orbit
        orbits.append(sorted(endomorphisms.index(item) for item in orbit))
    records = []
    for value in endomorphisms:
        image = sorted(set(value))
        records.append({
            "mapping": list(value),
            "image": image,
            "image_size": len(image),
            "kernel_labels": kernel_labels(value),
            "idempotent_retraction": value in idempotents and all(value[x] == x for x in image),
        })
    return {
        "schema": SCHEMA,
        "issue": 13,
        "status": "complete_unit_preserving_endomorphism_census",
        "partition_catalog": [list(state) for state in states],
        "base_generators": list(base),
        "generator_image_candidates": 15 ** 3,
        "endomorphisms": records,
        "endomorphism_count": len(endomorphisms),
        "image_size_histogram": {str(k): v for k, v in sorted(Counter(record["image_size"] for record in records).items())},
        "automorphism_indices": [endomorphisms.index(value) for value in automorphisms],
        "idempotent_retraction_indices": [endomorphisms.index(value) for value in idempotents],
        "automorphism_conjugacy_orbits": orbits,
        "exact_checks": {
            "all_generator_images_exhausted": 15 ** len(base) == 3375,
            "all_maps_preserve_identity": all(value[IDENTITY] == IDENTITY for value in endomorphisms),
            "all_maps_preserve_every_product": all(value[table[a][b]] == table[value[a]][value[b]] for value in endomorphisms for a in range(15) for b in range(15)),
            "automorphisms_are_exactly_bijective_endomorphisms": len(automorphisms) == 2,
            "idempotents_are_retractions_onto_images": all(all(value[x] == x for x in set(value)) for value in idempotents),
            "conjugacy_orbits_partition_census": sorted(i for orbit in orbits for i in orbit) == list(range(len(endomorphisms))),
        },
        "claim_boundary": {
            "included": "complete unit-preserving endomorphism census with images, kernels, retractions, and automorphism-conjugacy orbits",
            "excluded": "non-unital maps, graph-gadget maps, planar duality, periodic gluing, reliability, or thresholds",
            "parent_issue": "remain open",
        },
    }


def validate_artifact(artifact: Mapping[str, Any]) -> dict[str, Any]:
    expected = build_artifact()
    if artifact != expected:
        raise ValueError("endomorphism artifact does not exactly reproduce")
    if set(expected["exact_checks"].values()) != {True}:
        raise ValueError("all exact checks must pass")
    return {"schema": SCHEMA, "status": "valid", "endomorphisms": expected["endomorphism_count"]}


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

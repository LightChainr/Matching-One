#!/usr/bin/env python3
"""Canonical set-partition encoding under an explicit terminal group."""

from __future__ import annotations

import argparse
from itertools import permutations
import json
from pathlib import Path
from typing import Any, Iterable, Mapping, Optional, Sequence, Tuple


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "analysis" / "terminal_partition_canonical_manifest.json"
FORBIDDEN_KEYS = frozenset(
    {
        "edges",
        "gadget",
        "graphs",
        "planarity",
        "polynomial",
        "probabilities",
        "replacement_claim",
    }
)
RGS = Tuple[int, ...]
Permutation = Tuple[int, ...]


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _walk_forbidden(value: Any, path: str = "$") -> None:
    if isinstance(value, Mapping):
        bad = sorted(FORBIDDEN_KEYS.intersection(value))
        _require(not bad, "%s contains out-of-scope fields: %s" % (path, ",".join(bad)))
        for key, child in value.items():
            _walk_forbidden(child, "%s.%s" % (path, key))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            _walk_forbidden(child, "%s[%d]" % (path, index))


def validate_rgs(value: Sequence[int], expected_n: Optional[int] = None) -> RGS:
    """Validate and normalize a restricted-growth string."""

    rgs = tuple(value)
    _require(bool(rgs), "RGS must not be empty")
    _require(all(type(label) is int for label in rgs), "RGS labels must be integers")
    _require(rgs[0] == 0, "RGS must start with 0")
    maximum = 0
    for index, label in enumerate(rgs[1:], start=1):
        _require(0 <= label <= maximum + 1, "RGS growth rule fails at index %d" % index)
        maximum = max(maximum, label)
    if expected_n is not None:
        _require(len(rgs) == expected_n, "RGS length does not match terminal count")
    return rgs


def rgs_to_blocks(value: Sequence[int]) -> Tuple[Tuple[int, ...], ...]:
    rgs = validate_rgs(value)
    blocks = [[] for _ in range(max(rgs) + 1)]
    for terminal, label in enumerate(rgs):
        blocks[label].append(terminal)
    return tuple(tuple(block) for block in blocks)


def blocks_to_rgs(blocks: Iterable[Iterable[int]], n: int) -> RGS:
    """Encode a partition, independent of its input block ordering."""

    _require(type(n) is int and n > 0, "terminal count must be positive")
    owner = {}
    block_count = 0
    for block_index, source_block in enumerate(blocks):
        block = tuple(source_block)
        _require(bool(block), "partition blocks must not be empty")
        block_count += 1
        for terminal in block:
            _require(type(terminal) is int, "terminal labels must be integers")
            _require(0 <= terminal < n, "terminal label out of range")
            _require(terminal not in owner, "terminal occurs in multiple blocks")
            owner[terminal] = block_index
    _require(block_count > 0 and set(owner) == set(range(n)), "blocks must partition every terminal exactly once")

    canonical_labels = {}
    encoded = []
    for terminal in range(n):
        source_label = owner[terminal]
        if source_label not in canonical_labels:
            canonical_labels[source_label] = len(canonical_labels)
        encoded.append(canonical_labels[source_label])
    return validate_rgs(encoded, n)


def enumerate_rgs(n: int) -> Tuple[RGS, ...]:
    _require(type(n) is int and n > 0, "terminal count must be positive")
    values = [(0,)]
    for _ in range(1, n):
        values = [prefix + (label,) for prefix in values for label in range(max(prefix) + 2)]
    return tuple(values)


def validate_permutation(value: Sequence[int], n: int) -> Permutation:
    permutation = tuple(value)
    _require(len(permutation) == n, "permutation length does not match terminal count")
    _require(all(type(label) is int for label in permutation), "permutation labels must be integers")
    _require(set(permutation) == set(range(n)), "permutation must be a bijection")
    return permutation


def compose(p: Sequence[int], q: Sequence[int]) -> Permutation:
    """Compose old-to-new maps, applying q first and p second."""

    _require(len(p) == len(q), "cannot compose permutations of different sizes")
    p_valid = validate_permutation(p, len(p))
    q_valid = validate_permutation(q, len(q))
    return tuple(p_valid[q_valid[index]] for index in range(len(p_valid)))


def inverse(permutation: Sequence[int]) -> Permutation:
    valid = validate_permutation(permutation, len(permutation))
    result = [0] * len(valid)
    for old, new in enumerate(valid):
        result[new] = old
    return tuple(result)


def validate_group(group: Iterable[Sequence[int]], n: int) -> Tuple[Permutation, ...]:
    normalized = tuple(validate_permutation(permutation, n) for permutation in group)
    _require(bool(normalized), "permutation group must not be empty")
    _require(len(normalized) == len(set(normalized)), "permutation group contains duplicates")
    members = set(normalized)
    identity = tuple(range(n))
    _require(identity in members, "permutation group is missing the identity")
    for member in normalized:
        _require(inverse(member) in members, "permutation group is not closed under inverses")
        for other in normalized:
            _require(compose(member, other) in members, "permutation group is not closed under composition")
    return tuple(sorted(normalized))


def full_symmetric_group(n: int) -> Tuple[Permutation, ...]:
    _require(type(n) is int and n > 0, "terminal count must be positive")
    return tuple(permutations(range(n)))


def apply_permutation(value: Sequence[int], permutation: Sequence[int]) -> RGS:
    rgs = validate_rgs(value)
    valid = validate_permutation(permutation, len(rgs))
    moved_blocks = [[valid[terminal] for terminal in block] for block in rgs_to_blocks(rgs)]
    return blocks_to_rgs(moved_blocks, len(rgs))


def canonical_orbit(value: Sequence[int], group: Iterable[Sequence[int]]) -> RGS:
    rgs = validate_rgs(value)
    valid_group = validate_group(group, len(rgs))
    return min(apply_permutation(rgs, permutation) for permutation in valid_group)


def orbit_catalog(n: int, group: Iterable[Sequence[int]]) -> Mapping[RGS, Tuple[RGS, ...]]:
    valid_group = validate_group(group, n)
    buckets = {}
    for rgs in enumerate_rgs(n):
        representative = min(apply_permutation(rgs, permutation) for permutation in valid_group)
        buckets.setdefault(representative, []).append(rgs)
    return {representative: tuple(members) for representative, members in sorted(buckets.items())}


def validate_manifest(manifest: Mapping[str, Any]) -> Mapping[str, Any]:
    _walk_forbidden(manifest)
    _require(manifest.get("schema") == "matching-one/terminal-partition-canonical/v1", "unknown schema")
    _require(manifest.get("issue") == 13, "wrong issue")
    _require(manifest.get("status") == "encoding_primitive_only", "scope status drift")
    scope = manifest.get("scope", {})
    _require(scope.get("terminal_counts") == [3, 4], "terminal-count scope drift")
    _require(scope.get("encoding") == "restricted_growth_string", "encoding drift")
    group_contract = manifest.get("group_contract", {})
    _require(group_contract.get("permutation_convention") == "old_terminal_to_new_terminal", "group action convention drift")
    _require(group_contract.get("composition") == "p_after_q[i] = p[q[i]]", "composition convention drift")
    _require(
        group_contract.get("validation") == ["bijection", "identity", "inverse", "composition_closure"],
        "group validation contract drift",
    )
    _require(
        group_contract.get("canonical_representative") == "lexicographic_minimum_rgs_over_explicit_group",
        "canonicalization rule drift",
    )

    audited = {}
    reference_counts = manifest.get("reference_counts", {})
    for n in (3, 4):
        reference = reference_counts.get(str(n), {})
        partitions = enumerate_rgs(n)
        catalog = orbit_catalog(n, full_symmetric_group(n))
        representatives = [list(value) for value in catalog]
        _require(reference.get("bell_number") == len(partitions), "Bell count drift for n=%d" % n)
        _require(reference.get("full_symmetric_orbits") == len(catalog), "orbit count drift for n=%d" % n)
        _require(reference.get("orbit_representatives") == representatives, "orbit representatives drift for n=%d" % n)
        audited[str(n)] = {
            "partitions": len(partitions),
            "full_symmetric_orbits": len(catalog),
            "representatives": representatives,
        }
    _require(manifest.get("claim_boundary", {}).get("parent_issue") == "remain open", "parent boundary changed")
    return {
        "schema": manifest["schema"],
        "status": "valid_encoding_primitive_only",
        "terminal_counts": audited,
        "parent_issue": "remain open",
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    args = parser.parse_args(argv)
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    print(json.dumps(validate_manifest(manifest), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

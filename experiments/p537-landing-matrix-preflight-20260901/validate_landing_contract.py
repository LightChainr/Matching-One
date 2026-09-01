#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
import json
import os
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SOURCE = Path(
    os.environ.get(
        "MATCHING_ONE_THERMAL_GATE",
        ROOT / "experiments/p337-thermal-gate-audit-20260901/thermal_gate.py",
    )
)
spec = importlib.util.spec_from_file_location("thermal_gate", SOURCE)
tg = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(tg)


def clean(pattern):
    shared = set(pattern[:4]) & set(pattern[4:])
    incidence = sorted(
        (pattern[:4].count(value), pattern[4:].count(value)) for value in shared
    )
    return incidence == [(1, 1), (1, 1)] and max(pattern) + 1 == 6


def rotate(pattern, first, second):
    return tg.canon(
        tuple(pattern[(i + first) % 4] for i in range(4))
        + tuple(pattern[4 + (i + second) % 4] for i in range(4))
    )


def swap(pattern):
    return tg.canon(pattern[4:] + pattern[:4])


states = {pattern for pattern in tg.partitions(8) if clean(pattern)}
orbits = []
while states:
    representative = min(states)
    orbit = {
        rotate(representative, first, second)
        for first in range(4)
        for second in range(4)
    }
    orbit |= {swap(pattern) for pattern in tuple(orbit)}
    assert len({tg.g16(pattern) for pattern in orbit}) == 1
    orbits.append(
        {
            "representative": "".join(map(str, representative)),
            "size": len(orbit),
            "g16": tg.g16(representative),
        }
    )
    states -= orbit

assert len(orbits) == 4
assert Counter((row["size"], row["g16"]) for row in orbits) == Counter(
    {(16, 4): 2, (32, 0): 1, (8, 0): 1}
)
payload = {
    "Bell8_clean_six_block_partition_count": sum(row["size"] for row in orbits),
    "C4xC4_and_site_exchange_orbits": orbits,
    "nonzero_partition_count": sum(row["size"] for row in orbits if row["g16"]),
    "nonzero_orbits": 2,
    "all_nonzero_g16": 4,
}
(Path(__file__).resolve().parent / "landing_contract_validation.json").write_text(
    json.dumps(payload, indent=2) + "\n"
)
print(json.dumps(payload, indent=2))

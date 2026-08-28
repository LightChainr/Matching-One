#!/usr/bin/env python3
"""Stable entry point for the adjacent-axis annihilator scorer.

This module keeps the research scorer in ``score_axis_pair_annihilator`` as the
single implementation of the numerical fits, while fixing the batch-reader
ordering contract at the entry point used by tests and production.  The
underlying ``calculate`` function resolves ``read_pair_histograms`` through its
module globals, so installing the corrected reader here also fixes CLI use via
``main``.
"""
from __future__ import annotations

import csv
from pathlib import Path
from typing import Dict, Sequence, Tuple

import score_axis_pair_annihilator as _base


PairHistogram = _base.PairHistogram


def read_pair_histograms(paths: Sequence[Path]) -> Dict[Tuple[int, str, int], PairHistogram]:
    required = {"pair_L", "n", "L", "role", "batch", "samples", "kind", "k", "count"}
    records: Dict[Tuple[int, str, int], PairHistogram] = {}
    for path in paths:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            missing = required - set(reader.fieldnames or ())
            if missing:
                raise ValueError(f"{path}: missing fields {sorted(missing)}")
            for raw in reader:
                pair_L = int(raw["pair_L"])
                L = int(raw["L"])
                n = int(raw["n"])
                role = raw["role"]
                batch = int(raw["batch"])
                samples = int(raw["samples"])
                kind = raw["kind"]
                rank = int(raw["k"])
                count = int(raw["count"])
                if role not in ("upper", "lower") or kind not in ("minus", "plus"):
                    raise ValueError("unknown role/kind")
                if pair_L < 3 or L not in (pair_L, pair_L - 1) or n != L * L:
                    raise ValueError("inconsistent adjacent-axis geometry")
                if (role == "upper") != (L == pair_L):
                    raise ValueError("role does not match L")
                key = (pair_L, role, batch)
                if key not in records:
                    records[key] = PairHistogram(
                        pair_L,
                        n,
                        L,
                        role,
                        batch,
                        samples,
                        [0] * (n + 1),
                        [0] * (n + 1),
                    )
                row = records[key]
                if (row.n, row.L, row.samples) != (n, L, samples):
                    raise ValueError("metadata changed within batch")
                getattr(row, kind)[rank] += count
    if not records:
        raise ValueError("no histogram rows")

    for row in records.values():
        if sum(row.minus) != row.samples or sum(row.plus) != row.samples:
            raise ValueError("histogram total differs from samples")

    pair_sizes = sorted({key[0] for key in records})
    for pair_L in pair_sizes:
        signatures = []
        for role in ("upper", "lower"):
            selected = sorted(
                (row for key, row in records.items() if key[:2] == (pair_L, role)),
                key=lambda row: row.batch,
            )
            if not selected:
                raise ValueError(f"missing role {role} at pair L={pair_L}")
            batches = [row.batch for row in selected]
            if batches != list(range(len(batches))) or len(batches) < 2:
                raise ValueError("batches must be contiguous and at least two")
            signatures.append((batches, [row.samples for row in selected]))
        if signatures[0] != signatures[1]:
            raise ValueError("upper/lower batch alignment is absent")
    return records


# Install the corrected reader in the implementation module.  All functions
# below then use the same numerical code and the corrected batch contract.
_base.read_pair_histograms = read_pair_histograms

fit_f_shape = _base.fit_f_shape
fit_root_power = _base.fit_root_power
pair_jackknife = _base.pair_jackknife
calculate = _base.calculate
report = _base.report


def main() -> int:
    return _base.main()


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Run the frozen Issue #43 scorer with the exact production allocation locked.

This is a drop-in preflight wrapper around score_issue43_full_curve.py.  The
underlying reconstruction and statistical score are unchanged; this layer only
strengthens the provenance checks to match the append-only production freeze in
Issue #43 before N=185/265 target values are revealed.

If a genuinely necessary production-code fix ever requires a new clean source
commit, update the Issue #43 freeze *before* target execution/reveal and change
FROZEN_SOURCE_COMMIT here in the same pre-target change.  Do not relax this
validator after seeing a target.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Mapping

import score_issue43_full_curve as frozen


FROZEN_SOURCE_COMMIT = "302464c3a08bdf74a8cea079a50cfebd7fc8843f"
FROZEN_SAMPLES = 500_000_000
FROZEN_BATCHES = 100
FROZEN_SEED = 2_026_104_301
FROZEN_THREADS = 8
FROZEN_COUNTERS = {
    185: (7_000_000_000, 7_500_000_000),
    265: (7_500_000_000, 8_000_000_000),
}


_base_validate_metadata = frozen.validate_metadata
_base_validate_moments = frozen.validate_moments


def validate_metadata(path: Path) -> dict:
    """Run the frozen validator and then enforce the exact production freeze."""
    run = _base_validate_metadata(path)
    raw = json.loads(path.read_text(encoding="utf-8"))
    n = int(run["N"])

    expected_first, expected_last = FROZEN_COUNTERS[n]
    checks = {
        "source commit": str(run["commit"]).lower() == FROZEN_SOURCE_COMMIT,
        "samples": int(run["samples"]) == FROZEN_SAMPLES,
        "batches": int(run["batches"]) == FROZEN_BATCHES,
        "seed": int(run["seed"]) == FROZEN_SEED,
        "counter first": int(run["counter_first"]) == expected_first,
        "counter last": int(run["counter_last"]) == expected_last,
        "OpenMP": raw.get("openmp") is True,
        "threads_requested": int(raw.get("threads_requested", -1)) == FROZEN_THREADS,
        "compiler recorded": bool(str(raw.get("compiler", "")).strip()),
    }
    failed = [name for name, ok in checks.items() if not ok]
    if failed:
        raise ValueError("Issue #43 frozen production mismatch: " + ", ".join(failed))

    command = str(raw.get("command", ""))
    if int(frozen.command_option(command, "--threads")) != FROZEN_THREADS:
        raise ValueError("Issue #43 frozen production mismatch: command --threads")

    return run


def validate_moments(
    path: Path,
    run: Mapping[str, object],
    records: Mapping[tuple[str, int], Mapping[str, object]],
) -> None:
    """Retain the original checks and add the exact joint gap-square identity."""
    _base_validate_moments(path, run, records)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        for raw in reader:
            kminus2 = int(raw["sum_kminus2"])
            kplus2 = int(raw["sum_kplus2"])
            product = int(raw["sum_product"])
            gap2 = int(raw["sum_gap2"])
            expected_gap2 = kplus2 + kminus2 - 2 * product
            if gap2 != expected_gap2:
                raise ValueError("moment squared-gap identity failed")
            if gap2 < 0 or product < 0:
                raise ValueError("joint moment fields must be nonnegative")


def main() -> int:
    # Preserve one implementation of parsing, reconstruction, and scoring.
    frozen.validate_metadata = validate_metadata
    frozen.validate_moments = validate_moments
    return frozen.main()


if __name__ == "__main__":
    raise SystemExit(main())

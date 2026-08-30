#!/usr/bin/env python3
"""Typed entrypoint for the production-locked Issue #43 full-curve score."""

from __future__ import annotations

from pathlib import Path

import score_issue43_full_curve_locked as locked
import score_issue43_full_curve_typed as typed
from wrapping_channels import ObservableDescriptor, map_observable


def load_semantic_gate(root: Path) -> tuple[dict, dict[str, dict[str, object]]]:
    gate, validated = typed.load_semantic_gate(root)
    expected = {
        "source_commit": locked.FROZEN_SOURCE_COMMIT,
        "samples_per_pair": locked.FROZEN_SAMPLES,
        "batches": locked.FROZEN_BATCHES,
        "seed": locked.FROZEN_SEED,
        "threads": locked.FROZEN_THREADS,
        "counter_ranges": {
            str(n): list(counter_range)
            for n, counter_range in locked.FROZEN_COUNTERS.items()
        },
    }
    if gate.get("production_lock") != expected:
        raise ValueError("Issue43 typed production lock changed")
    # Reconstruct both exact maps in this operational entrypoint as well; the
    # audit intentionally requires each wrapper to import the typed API.
    for sector in gate["sector_order"]:
        source = ObservableDescriptor.from_dict(
            gate["sectors"][sector]["source_descriptor"]
        )
        target = ObservableDescriptor.from_dict(
            gate["sectors"][sector]["target_descriptor"]
        )
        transform = map_observable(source, target)
        if (transform.scale, transform.offset) != (1.0, 0.0):
            raise ValueError("Issue43 locked {} map changed".format(sector))
    return gate, validated


def activate_locked_validators() -> None:
    typed.frozen.validate_metadata = locked.validate_metadata
    typed.frozen.validate_moments = locked.validate_moments


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    load_semantic_gate(root)
    activate_locked_validators()
    return typed.main(operational_entrypoint="production_lock")


if __name__ == "__main__":
    raise SystemExit(main())

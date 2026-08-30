#!/usr/bin/env python3
"""Type-safe entrypoint for the frozen rank-gap boundary-target score."""

from __future__ import annotations

import argparse
from contextlib import redirect_stdout
import io
import json
from pathlib import Path
import sys
import tempfile
from typing import Callable, Sequence

from wrapping_channels import ObservableDescriptor, map_observable


SEMANTIC_GATE = "predictions/rank_gap_boundary_targets_semantic_gate_20260830.json"
SOURCE_SIZES = [65, 85, 130, 170, 185, 265]
TARGET_SIZES = [325, 425]


def load_semantic_gate(root: Path) -> tuple[dict, ObservableDescriptor, ObservableDescriptor, object]:
    gate = json.loads((root / SEMANTIC_GATE).read_text(encoding="utf-8"))
    if gate.get("status") != "semantic_gate_added_after_frozen_rank_gap_boundary_score":
        raise ValueError("rank-gap boundary semantic gate status changed")
    if gate.get("frozen_kernel_git_blob") != "8f9d41b5ea77828b7683c3a48d7f0517654b980e":
        raise ValueError("rank-gap boundary frozen kernel identity changed")
    if gate.get("source_size_order") != SOURCE_SIZES:
        raise ValueError("rank-gap boundary source size order changed")
    if gate.get("target_size_order") != TARGET_SIZES:
        raise ValueError("rank-gap boundary target size order changed")
    if (gate.get("paired_quantity"), gate.get("units")) != (
        "K_plus_minus_K_minus", "rank"
    ):
        raise ValueError("rank-gap paired quantity or units changed")
    if gate.get("orientation_pooling") != "equal_mean_of_first_and_second":
        raise ValueError("rank-gap orientation pooling changed")
    if gate.get("model") != {
        "expression": "E[G]=A*N^(5/8)+B",
        "exponent_in_N": {"numerator": 5, "denominator": 8, "fitted": False},
        "fit_parameters_in_order": ["A", "B"],
    }:
        raise ValueError("rank-gap frozen model changed")
    source = ObservableDescriptor.from_dict(gate["source_descriptor"])
    target = ObservableDescriptor.from_dict(gate["target_descriptor"])
    transform = map_observable(source, target)
    expected = gate["exact_registered_map"]
    if (transform.scale, transform.offset) != (
        float(expected["scale"]), float(expected["offset"])
    ) or (transform.scale, transform.offset) != (1.0, 0.0):
        raise ValueError("rank-gap registered map changed")
    return gate, source, target, transform


def _run_frozen_kernel(
    manifest: Path,
    source_score: Path,
    target_runs: Sequence[str],
    dps: int,
) -> dict:
    import score_rank_gap_boundary_targets as frozen_kernel

    with tempfile.TemporaryDirectory() as directory:
        output = Path(directory) / "score.json"
        argv = [
            "score_rank_gap_boundary_targets.py",
            "--manifest", str(manifest),
            "--source-score", str(source_score),
            "--dps", str(dps),
            "--output", str(output),
        ]
        for specification in target_runs:
            argv.extend(["--target-run", specification])
        old_argv = sys.argv
        try:
            sys.argv = argv
            with redirect_stdout(io.StringIO()):
                status = frozen_kernel.main()
        finally:
            sys.argv = old_argv
        if status != 0:
            raise RuntimeError(f"frozen rank-gap boundary kernel returned {status}")
        return json.loads(output.read_text(encoding="utf-8"))


def score_typed(
    root: Path,
    manifest: Path,
    source_score: Path,
    target_runs: Sequence[str] = (),
    dps: int = 50,
    *,
    runner: Callable[[Path, Path, Sequence[str], int], dict] = _run_frozen_kernel,
) -> dict:
    gate, source, target, transform = load_semantic_gate(root)
    result = runner(manifest, source_score, target_runs, dps)
    if result.get("model") != "E[G]=A*N^(5/8)+B; exponent fixed, not fitted":
        raise ValueError("frozen rank-gap model differs from semantic gate")
    if result.get("target_order") != gate["target_size_order"]:
        raise ValueError("frozen rank-gap target order differs from semantic gate")
    if result.get("status") not in {
        "frozen_source_fit_targets_unseen",
        "targets_revealed_and_scored_against_frozen_source_fit",
    }:
        raise ValueError("frozen rank-gap score status differs from semantic gate")
    result["observable_semantics"] = {
        "semantic_gate": SEMANTIC_GATE,
        "semantic_gate_status": gate["status"],
        "source_descriptor": source.to_dict(),
        "target_descriptor": target.to_dict(),
        "applied_transform": transform.to_dict(),
        "paired_quantity": gate["paired_quantity"],
        "units": gate["units"],
        "orientation_pooling": gate["orientation_pooling"],
        "source_target_relation": gate["source_target_relation"],
        "validation_order": "semantic_map_before_frozen_source_or_target_score",
        "evidence_boundary": gate["evidence_boundary"],
    }
    return result


def main(argv: Sequence[str] | None = None) -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--manifest", type=Path,
        default=root / "predictions/rank_gap_boundary_correction_targets_20260829.yaml",
    )
    parser.add_argument("--source-score", type=Path, required=True)
    parser.add_argument("--target-run", action="append", default=[])
    parser.add_argument("--dps", type=int, default=50)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    payload = score_typed(
        root, args.manifest, args.source_score, args.target_run, args.dps
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

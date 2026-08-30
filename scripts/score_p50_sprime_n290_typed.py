#!/usr/bin/env python3
"""Type-safe entrypoint for the frozen P50 N=290 S-prime scorer."""

from __future__ import annotations

import argparse
from contextlib import redirect_stdout
import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
from typing import Callable, Sequence

from wrapping_channels import ObservableDescriptor, map_observable


SEMANTIC_GATE = "predictions/p50_sprime_n290_semantic_gate_20260830.json"


def git_blob_sha(data: bytes) -> str:
    header = f"blob {len(data)}\0".encode("ascii")
    return hashlib.sha1(header + data).hexdigest()


def load_semantic_gate(root: Path) -> tuple[dict, ObservableDescriptor, ObservableDescriptor, object]:
    gate = json.loads((root / SEMANTIC_GATE).read_text(encoding="utf-8"))
    if gate.get("status") != "semantic_gate_added_after_frozen_p50_sprime_score":
        raise ValueError("P50 S-prime semantic gate status changed")
    if gate.get("frozen_kernel_git_blob") != "1d6d044f35758088e67b959c0a33a432a8c43b5b":
        raise ValueError("P50 S-prime frozen kernel identity changed")
    if gate.get("target_size") != 290 or gate.get("observable") != "P4_S_prime":
        raise ValueError("P50 S-prime target contract changed")
    if float(gate.get("leading_power_in_N")) != 1.25:
        raise ValueError("P50 S-prime leading power changed")
    if gate.get("models_in_scoring_order") != [
        "q2_even_scalar_correction",
        "rank2_jordan_log",
    ]:
        raise ValueError("P50 S-prime model order changed")
    source = ObservableDescriptor.from_dict(gate["source_descriptor"])
    target = ObservableDescriptor.from_dict(gate["target_descriptor"])
    transform = map_observable(source, target)
    expected = gate["exact_registered_map"]
    if (transform.scale, transform.offset) != (
        float(expected["scale"]),
        float(expected["offset"]),
    ):
        raise ValueError("registered P50 S-prime map differs from semantic gate")
    if (transform.scale, transform.offset) != (1.0, 0.0):
        raise ValueError("P50 S-prime descriptor map must be exact identity")
    return gate, source, target, transform


def validate_prediction_files(gate: dict, q2: Path, jordan: Path) -> None:
    supplied = (q2, jordan)
    for path, expected in zip(supplied, gate["prediction_files"]):
        if git_blob_sha(path.read_bytes()) != expected["git_blob_sha"]:
            raise ValueError(f"prediction identity changed: {path}")


def _run_frozen_kernel(child_hist: Path, q2: Path, jordan: Path) -> dict:
    import score_p50_sprime_n290 as frozen_kernel

    with tempfile.TemporaryDirectory() as directory:
        output = Path(directory) / "score.json"
        argv = [
            "score_p50_sprime_n290.py",
            "--child-hist", str(child_hist),
            "--q2", str(q2),
            "--jordan", str(jordan),
            "--output", str(output),
        ]
        old_argv = sys.argv
        try:
            sys.argv = argv
            with redirect_stdout(io.StringIO()):
                status = frozen_kernel.main()
        finally:
            sys.argv = old_argv
        if status != 0:
            raise RuntimeError(f"frozen P50 S-prime kernel returned {status}")
        return json.loads(output.read_text(encoding="utf-8"))


def score_typed(
    root: Path,
    child_hist: Path,
    q2: Path,
    jordan: Path,
    *,
    runner: Callable[[Path, Path, Path], dict] = _run_frozen_kernel,
) -> dict:
    gate, source, target, transform = load_semantic_gate(root)
    validate_prediction_files(gate, q2, jordan)
    result = runner(child_hist, q2, jordan)
    if result.get("observable") != gate["observable"] or result.get("N") != gate["target_size"]:
        raise ValueError("frozen P50 S-prime result differs from semantic gate")
    if list(result.get("models", {})) != gate["models_in_scoring_order"]:
        raise ValueError("frozen P50 S-prime result model order differs from semantic gate")
    result["observable_semantics"] = {
        "semantic_gate": SEMANTIC_GATE,
        "semantic_gate_status": gate["status"],
        "source_descriptor": source.to_dict(),
        "target_descriptor": target.to_dict(),
        "applied_transform": transform.to_dict(),
        "validation_order": "semantic_map_and_prediction_identity_before_frozen_kernel_score",
        "evidence_boundary": gate["evidence_boundary"],
    }
    return result


def main(argv: Sequence[str] | None = None) -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--child-hist", type=Path, required=True)
    parser.add_argument(
        "--q2",
        type=Path,
        default=root / "predictions/p48_sprime_q2_correction_20260828.yaml",
    )
    parser.add_argument(
        "--jordan",
        type=Path,
        default=root / "predictions/p48_sprime_jordan_log_20260828.yaml",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    payload = score_typed(root, args.child_hist, args.q2, args.jordan)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Type-safe entrypoint for the frozen intrinsic quantile-center score."""

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


SEMANTIC_GATE = "predictions/intrinsic_quantile_center_n145_n290_semantic_gate_20260830.json"
FEATURE_ORDER = ["Q", "w_0.025_scaled", "w_0.05_scaled", "c_0.025", "c_0.05"]
RESIDUAL_ORDER = [
    "Q290_minus_frozen_ratio_Q145",
    "scaled_width_drift_0.025",
    "scaled_width_drift_0.05",
]


def load_semantic_gate(root: Path) -> tuple[dict, ObservableDescriptor, ObservableDescriptor, object]:
    gate = json.loads((root / SEMANTIC_GATE).read_text(encoding="utf-8"))
    if gate.get("status") != "semantic_gate_added_after_frozen_intrinsic_quantile_center_score":
        raise ValueError("intrinsic quantile-center semantic gate status changed")
    if gate.get("frozen_kernel_git_blob") != "837652a580ef915eb5649ea63e5b1d1ba73a3e7e":
        raise ValueError("intrinsic quantile-center frozen kernel identity changed")
    if (gate.get("source_size"), gate.get("target_size")) != (145, 290):
        raise ValueError("intrinsic quantile-center size order changed")
    if gate.get("frozen_u") != [0.025, 0.05]:
        raise ValueError("intrinsic quantile-center levels changed")
    if gate.get("size_local_feature_order") != FEATURE_ORDER:
        raise ValueError("intrinsic quantile-center feature order changed")
    if gate.get("joint_residual_order") != RESIDUAL_ORDER:
        raise ValueError("intrinsic quantile-center residual order changed")
    if gate.get("normalization_powers_in_N") != {
        "Q": {"numerator": 3, "denominator": 4},
        "width": {"numerator": 3, "denominator": 8},
    }:
        raise ValueError("intrinsic quantile-center normalization changed")
    if gate.get("cross_size_covariance") != "zero_by_independent_rng_domains":
        raise ValueError("intrinsic quantile-center covariance boundary changed")
    source = ObservableDescriptor.from_dict(gate["source_descriptor"])
    target = ObservableDescriptor.from_dict(gate["target_descriptor"])
    transform = map_observable(source, target)
    expected = gate["exact_registered_map"]
    if (transform.scale, transform.offset) != (
        float(expected["scale"]), float(expected["offset"])
    ) or (transform.scale, transform.offset) != (1.0, 0.0):
        raise ValueError("intrinsic quantile-center registered map changed")
    return gate, source, target, transform


def _run_frozen_kernel(
    parent_hist: Path,
    parent_metadata: Path,
    child_hist: Path,
    child_metadata: Path,
    freeze: Path,
) -> dict:
    import score_intrinsic_quantile_center_n145_n290 as frozen_kernel

    with tempfile.TemporaryDirectory() as directory:
        output = Path(directory) / "score.json"
        argv = [
            "score_intrinsic_quantile_center_n145_n290.py",
            "--parent-hist", str(parent_hist),
            "--parent-metadata", str(parent_metadata),
            "--child-hist", str(child_hist),
            "--child-metadata", str(child_metadata),
            "--freeze", str(freeze),
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
            raise RuntimeError(f"frozen intrinsic quantile-center kernel returned {status}")
        return json.loads(output.read_text(encoding="utf-8"))


def score_typed(
    root: Path,
    parent_hist: Path,
    parent_metadata: Path,
    child_hist: Path,
    child_metadata: Path,
    freeze: Path,
    *,
    runner: Callable[[Path, Path, Path, Path, Path], dict] = _run_frozen_kernel,
) -> dict:
    gate, source, target, transform = load_semantic_gate(root)
    result = runner(parent_hist, parent_metadata, child_hist, child_metadata, freeze)
    if result.get("size_local_feature_order") != gate["size_local_feature_order"]:
        raise ValueError("frozen intrinsic quantile-center feature order differs from gate")
    if result.get("joint_residual_order") != gate["joint_residual_order"]:
        raise ValueError("frozen intrinsic quantile-center residual order differs from gate")
    if list(result.get("observations", {})) != ["N145", "N290"]:
        raise ValueError("frozen intrinsic quantile-center size order differs from gate")
    if result.get("frozen", {}).get("u") != gate["frozen_u"]:
        raise ValueError("frozen intrinsic quantile-center levels differ from gate")
    result["observable_semantics"] = {
        "semantic_gate": SEMANTIC_GATE,
        "semantic_gate_status": gate["status"],
        "source_descriptor": source.to_dict(),
        "target_descriptor": target.to_dict(),
        "applied_transform": transform.to_dict(),
        "normalization_powers_in_N": gate["normalization_powers_in_N"],
        "cross_size_covariance": gate["cross_size_covariance"],
        "validation_order": "semantic_map_before_frozen_coordinate_score",
        "evidence_boundary": gate["evidence_boundary"],
    }
    return result


def main(argv: Sequence[str] | None = None) -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parent-hist", type=Path, required=True)
    parser.add_argument("--parent-metadata", type=Path, required=True)
    parser.add_argument("--child-hist", type=Path, required=True)
    parser.add_argument("--child-metadata", type=Path, required=True)
    parser.add_argument(
        "--freeze", type=Path,
        default=root / "predictions/intrinsic_quantile_center_20260829.yaml",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    payload = score_typed(
        root, args.parent_hist, args.parent_metadata,
        args.child_hist, args.child_metadata, args.freeze,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

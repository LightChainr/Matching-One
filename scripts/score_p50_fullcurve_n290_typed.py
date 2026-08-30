#!/usr/bin/env python3
"""Type-safe entrypoint for the frozen P50 N=145 to N=290 full-curve score."""

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
from validate_p50_fullcurve_quantity_contract import (
    CONTRACT,
    git_blob,
    load_contract,
    validate_repository_files,
)


FROZEN_SCHEMA = "matching-one/P50-N145-N290-fullcurve-score/v1"
FROZEN_STATUS = "frozen primary full-curve score; independent parent/child streams"
FROZEN_COVARIANCE_RULE = (
    "N145 and N290 are jackknifed internally and treated as independent; "
    "residual covariance is Cov_child + ratio^2 Cov_parent"
)
P4_ORDER = ["P4_S", "P4_D", "P4_S_prime", "P4_D_prime"]


def _run_frozen_kernel(
    parent_hist: Path,
    parent_metadata: Path,
    child_hist: Path,
    child_metadata: Path,
    prediction: Path,
) -> dict:
    import score_p50_fullcurve_n290 as frozen_kernel

    with tempfile.TemporaryDirectory() as directory:
        output = Path(directory) / "score.json"
        argv = [
            "score_p50_fullcurve_n290.py",
            "--parent-hist", str(parent_hist),
            "--parent-metadata", str(parent_metadata),
            "--child-hist", str(child_hist),
            "--child-metadata", str(child_metadata),
            "--prediction", str(prediction),
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
            raise RuntimeError(f"frozen P50 full-curve kernel returned {status}")
        return json.loads(output.read_text(encoding="utf-8"))


def validate_frozen_result(result: dict, contract: dict) -> None:
    if result.get("schema") != FROZEN_SCHEMA or result.get("status") != FROZEN_STATUS:
        raise ValueError("frozen P50 full-curve result identity changed")
    if result.get("scoring_order") != contract["scoring_order"]:
        raise ValueError("frozen P50 full-curve scoring order changed")
    observations = result.get("observations", {})
    if list(observations) != ["N145", "N290"]:
        raise ValueError("frozen P50 full-curve size order changed")
    expected_features = [*contract["feature_order"], "p0"]
    if any(list(observations[name]) != expected_features for name in observations):
        raise ValueError("frozen P50 full-curve feature order changed")
    if result.get("primary_deltaM_transfer", {}).get("levels") != [0.0, 0.025, 0.05]:
        raise ValueError("frozen P50 full-curve DeltaM levels changed")
    diagnostics = result.get("p48_diagnostics", [])
    if [row.get("metric") for row in diagnostics] != P4_ORDER:
        raise ValueError("frozen P50 full-curve P4 order changed")
    if result.get("covariance_rule") != FROZEN_COVARIANCE_RULE:
        raise ValueError("frozen P50 full-curve covariance rule changed")


def score_typed(
    root: Path,
    parent_hist: Path,
    parent_metadata: Path,
    child_hist: Path,
    child_metadata: Path,
    prediction: Path,
    *,
    runner: Callable[[Path, Path, Path, Path, Path], dict] = _run_frozen_kernel,
) -> dict:
    contract, validated = load_contract(root)
    for row in validated.values():
        source = row["source_descriptor"]
        target = row["target_descriptor"]
        if (
            not isinstance(source, ObservableDescriptor)
            or map_observable(source, target) != row["transform"]
        ):
            raise ValueError("P50 wrapper topology-map validation changed")
    validate_repository_files(root, contract)
    expected_prediction = contract["frozen_prediction"]
    if git_blob(prediction) != expected_prediction["git_blob"]:
        raise ValueError("runtime P50 prediction identity changed")
    result = runner(
        parent_hist, parent_metadata, child_hist, child_metadata, prediction
    )
    validate_frozen_result(result, contract)
    result["observable_semantics"] = {
        "quantity_contract": CONTRACT,
        "quantity_contract_status": contract["status"],
        "topology_anchors": {
            name: {
                "source_descriptor": row["source_descriptor"].to_dict(),
                "target_descriptor": row["target_descriptor"].to_dict(),
                "applied_transform": row["transform"].to_dict(),
            }
            for name, row in validated.items()
        },
        "response_coordinates": contract["response_coordinates"],
        "sizes_in_order": contract["sizes_in_order"],
        "representations": contract["representations"],
        "lineage_sign": contract["lineage_sign"],
        "rng_relation": contract["rng_relation"],
        "feature_order": contract["feature_order"],
        "scoring_order": contract["scoring_order"],
        "covariance_contract": contract["covariance_contract"],
        "validation_order": (
            "contract_maps_and_repository_blobs_before_frozen_kernel; "
            "frozen_output_contract_before_annotation"
        ),
        "evidence_boundary": contract["evidence_boundary"],
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
        "--prediction",
        type=Path,
        default=root / "predictions/p49_slope_two_sector_145_290_20260828.yaml",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    payload = score_typed(
        root,
        args.parent_hist,
        args.parent_metadata,
        args.child_hist,
        args.child_metadata,
        args.prediction,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

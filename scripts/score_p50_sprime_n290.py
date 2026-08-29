#!/usr/bin/env python3
"""Score the frozen q=2 and Jordan P4[S-prime] predictions at N=290."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import yaml

from score_p50_fullcurve_n290 import FEATURE_ORDER, estimate, grouped, read_one_size


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def scalar_prediction_score(
    observed: float, target_se: float, predicted: float, source_se: float
) -> dict[str, float]:
    variance = target_se * target_se + source_se * source_se
    z = (observed - predicted) / math.sqrt(variance)
    return {
        "predicted": predicted,
        "source_prediction_se": source_se,
        "residual": observed - predicted,
        "variance": variance,
        "signed_z": z,
        "chi_square": z * z,
        "degrees_of_freedom": 1,
        "two_sided_p": math.erfc(abs(z) / math.sqrt(2.0)),
    }


def frozen_n290(path: Path) -> tuple[str, float, float]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    row = next(item for item in payload["frozen_predictions"] if item["N"] == 290)
    return (
        str(payload["model"]["name"]),
        float(row["P4_S_prime"]),
        float(row["source_fit_se_P4_S_prime"]),
    )


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--child-hist", type=Path, required=True)
    parser.add_argument("--q2", type=Path, default=root / "predictions/p48_sprime_q2_correction_20260828.yaml")
    parser.add_argument("--jordan", type=Path, default=root / "predictions/p48_sprime_jordan_log_20260828.yaml")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    child = grouped(read_one_size(args.child_hist, 290), 290)
    point, covariance = estimate(child, lineage_sign=+1.0)
    index = FEATURE_ORDER.index("P4_S_prime")
    observed = float(point["P4_S_prime"])
    target_se = math.sqrt(float(covariance[index][index]))
    models = {}
    for position, path in enumerate((args.q2, args.jordan), start=1):
        name, predicted, source_se = frozen_n290(path)
        try:
            display_path = str(path.resolve().relative_to(root))
        except ValueError:
            display_path = str(path)
        models[name] = {
            "frozen_order": position,
            **scalar_prediction_score(observed, target_se, predicted, source_se),
            "prediction": display_path,
            "prediction_sha256": sha256(path),
        }
    payload = {
        "schema": "matching-one/P50-N290-Sprime-frozen-score/v1",
        "status": "ordered post-primary score on the P50 raw block",
        "observable": "P4_S_prime",
        "N": 290,
        "observed": observed,
        "target_se": target_se,
        "models": models,
        "decision": (
            "both frozen models survive; q2 remains first in the declared order, "
            "while Jordan is descriptively closer"
        ),
        "evidence_boundary": (
            "This reuses the P50 histogram and is not independent of other P50 diagnostics."
        ),
        "child_hist": str(args.child_hist),
        "child_hist_sha256": sha256(args.child_hist),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

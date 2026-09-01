#!/usr/bin/env python3
"""Three-model held-out N145 signed-real phase scorer."""

from __future__ import annotations

import argparse
import importlib.util
import json
from math import sqrt
from pathlib import Path
import sys
from typing import Sequence


ROOT = Path(__file__).resolve().parents[2]
BASE_SCORER = ROOT / "experiments" / "p275-gaussian-c3-phase-20260901" / "score.py"


def _load_base():
    spec = importlib.util.spec_from_file_location("p275_n65_phase_score", BASE_SCORER)
    if spec is None or spec.loader is None:
        raise ImportError(BASE_SCORER)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


BASE = _load_base()
ALPHA = BASE.ALPHA
DELTA = __import__("math").atan2(3.0, 4.0)
COORDINATES = BASE.COORDINATES
INVARIANT_PREDICTIONS = {
    "H0": 1.0,
    "H4": 0.42197248,
    "H8": -0.6438784522452992,
}


def projective_invariant(vector: Sequence[float]) -> float:
    z1 = complex(vector[0], vector[1])
    z2 = complex(vector[2], vector[3])
    product = z2 * z1.conjugate()
    denominator = abs(product) ** 2
    if denominator <= 0.0:
        raise ValueError("projective invariant requires both complex coordinates nonzero")
    return (product * product).real / denominator


def jackknife_invariant(rows: Sequence[Sequence[float]]) -> dict[str, object]:
    count = len(rows)
    means = [sum(row[j] for row in rows) / count for j in range(4)]
    replicates = []
    for omitted in range(count):
        replicate = [
            (count * means[j] - rows[omitted][j]) / (count - 1) for j in range(4)
        ]
        replicates.append(projective_invariant(replicate))
    center = projective_invariant(means)
    replicate_mean = sum(replicates) / count
    variance = (count - 1) / count * sum(
        (value - replicate_mean) ** 2 for value in replicates
    )
    standard_error = sqrt(max(0.0, variance))
    return {
        "value": center,
        "jackknife_standard_error": standard_error,
        "normal_95_interval": [center - 1.96 * standard_error, center + 1.96 * standard_error],
        "predictions": INVARIANT_PREDICTIONS,
        "decision_role": "diagnostic_only",
    }


def score(rows: Sequence[Sequence[float]]) -> dict[str, object]:
    means = [sum(row[j] for row in rows) / len(rows) for j in range(4)]
    covariance = BASE.covariance_of_mean(rows)
    models = {
        "H0": BASE.profile_model(means, covariance, 0.0),
        "H4": BASE.profile_model(means, covariance, +4.0 * DELTA),
        "H8": BASE.profile_model(means, covariance, -8.0 * DELTA),
    }
    passes = {name: model["p_value"] >= ALPHA for name, model in models.items()}
    survivors = [name for name, passed in passes.items() if passed]
    if survivors == ["H8"]:
        decision = "H8_OBSERVER_HARMONIC_SELECTED"
    elif survivors == ["H0"]:
        decision = "H0_EVEN_CHARACTER_SELECTED"
    elif survivors == ["H4"]:
        decision = "H4_SELECTED"
    else:
        decision = "STOP_MIXED_OR_UNRESOLVED"
    return {
        "schema": "matching-one/p275-gaussian-c3-phase-n145-score/v1",
        "batch_count": len(rows),
        "coordinate_order": list(COORDINATES),
        "mean": means,
        "covariance_of_mean": covariance,
        "delta_radians": DELTA,
        "alpha": ALPHA,
        "models": models,
        "projective_invariant": jackknife_invariant(rows),
        "decision": decision,
        "model_passes": passes,
        "surviving_models": survivors,
        "top_up_authorized": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("batches", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError(args.out)
    payload = score(BASE.read_rows(args.batches))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"decision": payload["decision"], "output": str(args.out)}))


if __name__ == "__main__":
    main()

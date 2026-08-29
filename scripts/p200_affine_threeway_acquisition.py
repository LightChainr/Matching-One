#!/usr/bin/env python3
"""Freeze the minimal N580/N650 three-way acquisition for Issue #200."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Sequence


ORDERS = (2, 3, 4, 5, 6)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def width_corrected_jet(state: dict, orders: Sequence[int] = ORDERS) -> list[float]:
    width = float(state["canonical_width"])
    jet = [float(value) for value in state["finite_thermal_jet"]]
    return [width**order * jet[order] for order in orders]


def forecast_from_norm5_pair(
    parent: Sequence[float], norm5: Sequence[float], model: str
) -> list[float]:
    if model == "q2_analytic":
        return [(9.0 * child - base) / 8.0 for base, child in zip(parent, norm5)]
    if model == "rank2_Jordan":
        ratio = math.log(2.0) / math.log(5.0)
        return [child + ratio * (child - base) for base, child in zip(parent, norm5)]
    raise ValueError(model)


def forecast_from_norm2_pair(
    parent: Sequence[float], norm2: Sequence[float], model: str
) -> list[float]:
    if model == "q2_analytic":
        return [(9.0 * child - 4.0 * base) / 5.0 for base, child in zip(parent, norm2)]
    if model == "rank2_Jordan":
        ratio = math.log(5.0) / math.log(2.0)
        return [child + ratio * (child - base) for base, child in zip(parent, norm2)]
    raise ValueError(model)


def cocycle_multiplier(model: str) -> float:
    if model == "q2_analytic":
        return 8.0 / 5.0
    if model == "rank2_Jordan":
        return math.log(5.0) / math.log(2.0)
    raise ValueError(model)


def path_difference_factor(model: str) -> float:
    if model == "q2_analytic":
        return 9.0 / 8.0
    c = cocycle_multiplier(model)
    return 1.0 + 1.0 / c


def scaled_matrix(matrix: Sequence[Sequence[float]], factor: float) -> list[list[float]]:
    return [[factor * factor * float(value) for value in row] for row in matrix]


def render(p180_path: Path, jet_path: Path, parity_path: Path) -> dict:
    p180 = json.loads(p180_path.read_text(encoding="utf-8"))
    jet = json.loads(jet_path.read_text(encoding="utf-8"))
    parity = json.loads(parity_path.read_text(encoding="utf-8"))

    states = jet["states"]
    x1 = width_corrected_jet(states["65"])
    x2 = width_corrected_jet(states["130"])
    x5 = width_corrected_jet(states["325"])
    models = {}
    for model in ("q2_analytic", "rank2_Jordan"):
        source = jet["secondary_multiplier_cocycles"][model]
        residual = [float(value) for value in source["residual"][:5]]
        covariance = [[float(value) for value in row[:5]] for row in source["covariance"][:5]]
        via_norm5 = forecast_from_norm5_pair(x1, x5, model)
        via_norm2 = forecast_from_norm2_pair(x1, x2, model)
        factor = path_difference_factor(model)
        difference = [left - right for left, right in zip(via_norm5, via_norm2)]
        expected_difference = [factor * value for value in residual]
        if max(abs(a - b) for a, b in zip(difference, expected_difference)) > 1e-10:
            raise AssertionError("path difference is not the transformed frozen cocycle residual")
        difference_covariance = scaled_matrix(covariance, factor)
        models[model] = {
            "orders": list(ORDERS),
            "N650_prediction_via_N65_N325": via_norm5,
            "N650_prediction_via_N65_N130": via_norm2,
            "path_difference": difference,
            "path_difference_factor_times_existing_residual": factor,
            "path_difference_covariance": difference_covariance,
            "path_difference_standard_error": [
                math.sqrt(max(0.0, difference_covariance[i][i])) for i in range(5)
            ],
            "existing_source_score": source["score"],
        }

    odd = parity["odd"]
    odd_template = [float(value) for value in odd["residual"]]
    odd_covariance = [[float(value) for value in row] for row in odd["covariance"]]

    return {
        "schema": "matching-one/p200-affine-threeway-acquisition/v1",
        "issue": 200,
        "status": "frozen_post_reveal_before_N580_N650_acquisition",
        "N580_radial_clock": {
            "designs": [[24, 2], [18, 16]],
            "paths": ["145 --(1+i)^2--> 580", "145 --(1+i)--> 290 --(1+i)--> 580"],
            "state_order": p180["state_order"],
            "models": p180["predictions"],
            "decision": "joint four-state affine score selects q2, Jordan, both, or neither",
        },
        "N650_mixed_clock": {
            "parent_designs_N65": [[8, 1], [7, 4]],
            "endpoint_designs_N650": [[23, 11], [17, 19]],
            "path_A": "65 --(2-i)--> 325 --(1+i)--> 650",
            "path_B": "65 --(1+i)--> 130 --(2-i)--> 650",
            "same_endpoint_multiplier": "(2-i)(1+i)=(1+i)(2-i)=3+i",
            "spin4_character": "chi4(3+i)=(7+24i)/25",
            "unmarked_rule": "one endpoint histogram only; the two factorizations are not independent replicas",
            "width_corrected_jet_models": models,
            "marked_commutator": {
                "definition": "C_mark=H_[(2-i) then (1+i)]-H_[(1+i) then (2-i)]",
                "q2_target": "zero",
                "Jordan_target": "zero",
                "morphism_memory_target": "nonzero rank-one vector aligned with the frozen P57 conjugation-odd template",
                "frozen_P57_odd_template_orders": list(ORDERS),
                "frozen_P57_odd_template": odd_template,
                "frozen_P57_odd_template_covariance": odd_covariance,
                "direction_score": "fit one amplitude to the frozen template, then score the four-dimensional covariance-orthogonal residual",
            },
        },
        "three_way_decision": {
            "q2": "N580 obeys (-1/2,3/2); N650 agrees with both q2 path forecasts; marked commutator is zero",
            "Jordan": "N580 obeys (-1,2); N650 agrees with both Jordan path forecasts; marked commutator is zero",
            "morphism_memory": "one unmarked radial law survives but the marked N650 commutator is nonzero and rank-one aligned with the P57 odd template",
            "higher_rank": "neither N580 affine law survives, or marked defect is not rank-one/template-aligned",
        },
        "claim_boundary": {
            "exact": "Gaussian multiplier products, Cayley-Hamilton coefficients, and path-difference/cocycle algebra",
            "source_data": "N145/N290 and norm5/P57 are already revealed and provide predictions/covariances only",
            "high_risk": "a morphism-sensitive N650 defect lies in the one-dimensional P57 odd template direction",
        },
        "provenance": {
            "p180": {"path": str(p180_path), "sha256": sha256(p180_path)},
            "thermal_jet": {"path": str(jet_path), "sha256": sha256(jet_path)},
            "parity": {"path": str(parity_path), "sha256": sha256(parity_path)},
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--p180",
        type=Path,
        default=Path("results/full-curve-transfer/p180_n145_n290_affine_clock.json"),
    )
    parser.add_argument(
        "--thermal-jet",
        type=Path,
        default=Path("results/server-20260829/P57-norm5-500m/thermal_jet_score.json"),
    )
    parser.add_argument(
        "--parity",
        type=Path,
        default=Path("results/server-20260829/P57-norm5-500m/conjugation_parity_diagnostic.json"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("predictions/p200_affine_threeway_acquisition_20260829.json"),
    )
    args = parser.parse_args()
    payload = render(args.p180, args.thermal_jet, args.parity)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

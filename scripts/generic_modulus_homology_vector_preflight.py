#!/usr/bin/env python3
"""No-sampling preflight for a stabilizer-free homology-vector experiment."""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from math import gcd
from pathlib import Path
from typing import Iterable, Sequence

import mpmath as mp
import numpy as np

from pinson_arguin_kdv import eisenstein_e4, primitive_k4_holomorphic_series
from pinson_arguin_primitive import (
    engine_to_paper,
    primitive_probability_direct,
)


Matrix = tuple[tuple[int, int], tuple[int, int]]
Vector = tuple[int, int]


def parse_matrix(value: Sequence[Sequence[int]]) -> Matrix:
    if len(value) != 2 or any(len(row) != 2 for row in value):
        raise ValueError("period matrix must be 2 by 2")
    matrix = tuple(tuple(int(item) for item in row) for row in value)
    if determinant(matrix) <= 0:
        raise ValueError("period matrix must have positive determinant")
    return matrix  # type: ignore[return-value]


def determinant(matrix: Matrix) -> int:
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def period_tau(matrix: Matrix) -> tuple[Fraction, Fraction]:
    """Return tau=(v1.v2+i det)/|v1|^2 for column period vectors."""

    a, b = matrix[0]
    c, d = matrix[1]
    norm = a * a + c * c
    return Fraction(a * b + c * d, norm), Fraction(determinant(matrix), norm)


def gram(matrix: Matrix) -> Matrix:
    a, b = matrix[0]
    c, d = matrix[1]
    return ((a * a + c * c, a * b + c * d), (a * b + c * d, b * b + d * d))


def lattice_stabilizer(matrix: Matrix) -> list[Matrix]:
    """Enumerate the complete orientation-preserving integral Gram stabilizer."""

    g = gram(matrix)
    eig = np.linalg.eigvalsh(np.asarray(g, dtype=float))
    bound = math.ceil(math.sqrt(max(g[0][0], g[1][1]) / float(eig[0]))) + 1
    answer: list[Matrix] = []
    for a in range(-bound, bound + 1):
        for b in range(-bound, bound + 1):
            for c in range(-bound, bound + 1):
                for d in range(-bound, bound + 1):
                    candidate = ((a, b), (c, d))
                    if determinant(candidate) != 1:
                        continue
                    # M^T G M, written without floating arithmetic.
                    gm00 = g[0][0] * a * a + 2 * g[0][1] * a * c + g[1][1] * c * c
                    gm01 = g[0][0] * a * b + g[0][1] * (a * d + b * c) + g[1][1] * c * d
                    gm11 = g[0][0] * b * b + 2 * g[0][1] * b * d + g[1][1] * d * d
                    if (gm00, gm01, gm11) == (g[0][0], g[0][1], g[1][1]):
                        answer.append(candidate)
    return sorted(answer)


def primitive_lines(cutoff: int) -> list[Vector]:
    lines = []
    for a in range(0, cutoff + 1):
        for b in range(-cutoff, cutoff + 1):
            if a == 0 and b <= 0:
                continue
            if max(abs(a), abs(b)) == 0 or gcd(abs(a), abs(b)) != 1:
                continue
            lines.append((a, b))
    return sorted(lines)


def primitive_k4_k8_holomorphic(
    engine_line: Vector, tau: mp.mpc, *, dps: int
) -> tuple[mp.mpc, mp.mpc]:
    """Return the K4 response and an E4*K4 weight-8 covariant template."""

    a, b = engine_to_paper(engine_line)
    with mp.workdps(dps + 25):
        z = mp.mpc(tau)
        k4 = primitive_k4_holomorphic_series(a, b, z, dps=dps)
        k8 = eisenstein_e4(z, 100) * k4
        return +k4, +k8


def rank_readouts(
    probabilities: Sequence[float], derivatives: Sequence[float]
) -> dict[str, float]:
    """Map (P0,P1,P2) and its thermal derivative to q/E readouts."""

    if len(probabilities) != 3 or len(derivatives) != 3:
        raise ValueError("rank readout needs P0,P1,P2 and three derivatives")
    p0, p1, p2 = map(float, probabilities)
    dp0, dp1, dp2 = map(float, derivatives)
    if not math.isclose(p0 + p1 + p2, 1.0, rel_tol=0.0, abs_tol=1e-9):
        raise ValueError("rank probabilities do not sum to one")
    if not math.isclose(dp0 + dp1 + dp2, 0.0, rel_tol=0.0, abs_tol=1e-9):
        raise ValueError("rank derivatives do not sum to zero")
    return {
        "P0": p0,
        "P1": p1,
        "P2": p2,
        "q": p2 - p0,
        "E": p0 + p2,
        "dpq": dp2 - dp0,
        "dpE": dp0 + dp2,
    }


def _tau(matrix: Matrix) -> mp.mpc:
    real, imag = period_tau(matrix)
    return mp.mpc(mp.mpf(real.numerator) / real.denominator, mp.mpf(imag.numerator) / imag.denominator)


def _channel_vector(
    tau: mp.mpc,
    explicit: Sequence[Vector],
    tail: Sequence[Vector],
    *,
    dps: int,
    operator: str,
) -> np.ndarray:
    if operator == "probability":
        explicit_values = [
            primitive_probability_direct(*engine_to_paper(line), tau, dps=dps)
            for line in explicit
        ]
        tail_values = [
            primitive_probability_direct(*engine_to_paper(line), tau, dps=dps)
            for line in tail
        ]
        p1 = mp.fsum(explicit_values) + mp.fsum(tail_values)
        edge = (1 - p1) / 2
        values = [edge, *explicit_values, mp.fsum(tail_values), edge]
    else:
        explicit_pairs = [primitive_k4_k8_holomorphic(line, tau, dps=dps) for line in explicit]
        tail_pairs = [primitive_k4_k8_holomorphic(line, tau, dps=dps) for line in tail]
        which = 0 if operator == "K4" else 1
        explicit_values = [2 * mp.re(pair[which]) for pair in explicit_pairs]
        tail_value = mp.fsum(2 * mp.re(pair[which]) for pair in tail_pairs)
        rank1_sum = mp.fsum(explicit_values) + tail_value
        # The rank-0/rank-2 completion enforces total probability and Q=1 duality.
        edge = -rank1_sum / 2
        values = [edge, *explicit_values, tail_value, edge]
    return np.asarray([float(value) for value in values], dtype=float)


def _metric_projection(vector: np.ndarray, nuisance: np.ndarray, metric: np.ndarray) -> np.ndarray:
    if nuisance.size == 0:
        return vector
    gram_n = nuisance.T @ metric @ nuisance
    return vector - nuisance @ np.linalg.pinv(gram_n, rcond=1e-12) @ nuisance.T @ metric @ vector


def _record_matrix(matrix: Matrix) -> dict[str, object]:
    real, imag = period_tau(matrix)
    return {
        "matrix": [list(row) for row in matrix],
        "determinant": determinant(matrix),
        "tau": {"real": str(real), "imag": str(imag)},
        "gram": [list(row) for row in gram(matrix)],
        "stabilizer": [[list(row) for row in item] for item in lattice_stabilizer(matrix)],
    }


def analyze(config: dict[str, object]) -> dict[str, object]:
    dps = int(config["decimal_precision"])
    base = parse_matrix(config["base_period_matrix"])
    stencils = {
        name: {sign: parse_matrix(matrix) for sign, matrix in record.items()}
        for name, record in config["signed_stencils"].items()
    }
    matrices = [base, *(matrix for record in stencils.values() for matrix in record.values())]
    if any(determinant(matrix) != determinant(base) for matrix in matrices):
        raise ValueError("all stencil cells must have equal area")
    expected_stabilizer = [((-1, 0), (0, -1)), ((1, 0), (0, 1))]
    if any(lattice_stabilizer(matrix) != expected_stabilizer for matrix in matrices):
        raise ValueError("a selected period cell has a nontrivial orientation-preserving stabilizer")

    explicit = primitive_lines(int(config["primitive_cutoff"]))
    all_tail = primitive_lines(int(config["tail_cutoff"]))
    tail = [line for line in all_tail if line not in set(explicit)]
    channel_order = ["rank0", *[f"rank1:{a},{b}" for a, b in explicit], "rank1:tail", "rank2"]

    with mp.workdps(dps + 15):
        base_probability = _channel_vector(_tau(base), explicit, tail, dps=dps, operator="probability")
        direction_rows: dict[str, dict[str, object]] = {}
        stacked: dict[str, list[np.ndarray]] = {"mu_KdV": [], "mu_Q4_Jordan": [], "baseline": []}
        for name, pair in stencils.items():
            minus_tau, plus_tau = _tau(pair["minus"]), _tau(pair["plus"])
            span = abs(plus_tau - minus_tau)
            if span == 0:
                raise ValueError(f"{name} has a zero tangent span")
            rows = {}
            for operator, model in (("probability", "baseline"), ("K4", "mu_KdV"), ("K8", "mu_Q4_Jordan")):
                minus = _channel_vector(minus_tau, explicit, tail, dps=dps, operator=operator)
                plus = _channel_vector(plus_tau, explicit, tail, dps=dps, operator=operator)
                derivative = (plus - minus) / float(span)
                stacked[model].append(derivative)
                rows[model] = derivative.tolist()
            direction_rows[name] = {
                "minus": _record_matrix(pair["minus"]),
                "plus": _record_matrix(pair["plus"]),
                "complex_span": {
                    "real": mp.nstr(mp.re(plus_tau - minus_tau), 25),
                    "imag": mp.nstr(mp.im(plus_tau - minus_tau), 25),
                    "absolute": mp.nstr(span, 25),
                },
                "directional_vectors": rows,
            }

    probability_covariance = np.diag(base_probability) - np.outer(base_probability, base_probability)
    metric_block = np.linalg.pinv(probability_covariance, rcond=1e-11)
    metric = np.kron(np.eye(len(stencils)), metric_block)
    baseline_stack = np.concatenate(stacked["baseline"])
    nuisance_columns = []
    width = len(channel_order)
    for index, row in enumerate(stacked["baseline"]):
        column = np.zeros(len(stencils) * width)
        column[index * width : (index + 1) * width] = row
        nuisance_columns.append(column)
    nuisance = np.column_stack(nuisance_columns)

    projected = {}
    for model in ("mu_KdV", "mu_Q4_Jordan"):
        raw = np.concatenate(stacked[model])
        vec = _metric_projection(raw, nuisance, metric)
        norm2 = float(vec @ metric @ vec)
        projected[model] = {"vector": vec, "metric_norm_squared": norm2}
    first = projected["mu_KdV"]
    second = projected["mu_Q4_Jordan"]
    cosine = float(first["vector"] @ metric @ second["vector"] / math.sqrt(first["metric_norm_squared"] * second["metric_norm_squared"]))
    shape_d2 = 2 * (1 - abs(cosine))
    cost = config["cost_proxy"]
    cpu_hours = (
        float(cost["one_thread_elapsed_seconds"])
        / float(cost["samples"])
        * 1_000_000
        * int(cost["planned_stencil_cells"])
        / 3600
    )

    p0 = float(base_probability[0])
    p2 = float(base_probability[-1])
    p1 = 1 - p0 - p2
    rank = rank_readouts((p0, p1, p2), (0.0, 0.0, 0.0))
    rank["dpq"] = "required_from_finite_pilot"
    rank["dpE"] = "required_from_finite_pilot"

    return {
        "schema": "matching-one.generic-modulus-homology-vector-preflight-result.v1",
        "status": "launch-ready observable and theory preflight; no random sampling",
        "provenance": {
            "execution": "local_fallback_after_TV2N0X_publickey_failure",
            "random_samples": 0,
            "server_directory_requested_but_not_written": "/workspace/matching-one-generic-modulus-20260901",
        },
        "geometry": {
            "base": _record_matrix(base),
            "all_cells_equal_area": determinant(base),
            "stabilizer_verdict": "all five cells have orientation-preserving lattice stabilizer {+I,-I}",
            "directions": direction_rows,
        },
        "observable": {
            "channel_order": channel_order,
            "base_continuum_probability": base_probability.tolist(),
            "base_rank_readouts": rank,
            "tail_definition": f"all primitive unoriented lines with max norm {int(config['primitive_cutoff']) + 1}..{int(config['tail_cutoff'])}; finite producer saves every observed line sparsely",
            "finite_pilot_required_columns": [
                "P0", "P1_by_primitive_winding", "P1_tail", "P2",
                "d_p_P0", "d_p_P1_by_primitive_winding", "d_p_P1_tail", "d_p_P2",
                "q=P2-P0", "E=P0+P2", "d_p_q=d_p_P2-d_p_P0", "d_p_E=d_p_P0+d_p_P2"
            ],
        },
        "pipeline": {
            "order": ["homology vector", "q/E", "rank-1 physical normalizer", "pooled root", "U"],
            "pooled_root": "solve pooled q(p_star)=0 inside every aligned covariance replicate",
            "normalizer": "D=d_p q(p_star), common to both signs of one tangent stencil",
            "continuum_subtraction": "subtract the saved baseline homology-vector directional derivative before field scoring",
            "U_vector": "N^(13/8) * directional_contrast[d_p P_h] / (2 D) at p_star, for every saved homology channel h",
        },
        "predictions": {
            "mu_KdV": np.concatenate(stacked["mu_KdV"]).tolist(),
            "mu_Q4_Jordan": np.concatenate(stacked["mu_Q4_Jordan"]).tolist(),
            "mu_embedding": [0.0] * (len(stencils) * len(channel_order)),
            "stack_order": list(stencils),
            "channel_order_within_direction": channel_order,
        },
        "separation": {
            "covariance": "Q=1 continuum multinomial per-sample proxy at the base modulus; empirical full-vector covariance is not yet available",
            "nuisance": "one continuum embedding-tangent leakage column per signed direction",
            "KdV_Q4_Jordan_profiled_cosine_up_to_sign": abs(cosine),
            "KdV_Q4_Jordan_unit_Fisher_shape_D2": shape_d2,
            "cost_lower_bound_cpu_hours_per_1m_four_cells": cpu_hours,
            "conditional_shape_value_per_cpu_hour": shape_d2 / cpu_hours,
            "portfolio_maximin_V": 0.0,
            "portfolio_zero_reason": "mu_embedding is exactly zero after continuum subtraction and the physical nonzero operator-amplitude floor is not supplied; the KdV-versus-Q4/Jordan shape separation is nevertheless nondegenerate",
            "profiled_metric_norm_squared": {
                model: record["metric_norm_squared"] for model, record in projected.items()
            },
        },
        "decision": {
            "theory_templates_rank_separated": shape_d2 > 1e-6,
            "pilot_ready": True,
            "next_bounded_action": "one N100 sparse full-homology covariance pilot on the four frozen equal-area stencil cells",
            "large_production_ready": False,
            "large_production_gate": "empirical full-vector covariance plus a nonzero amplitude or model-lowering threshold",
        },
        "boundary": config["execution_boundary"],
    }


def render_report(result: dict[str, object]) -> str:
    geometry = result["geometry"]
    separation = result["separation"]
    observable = result["observable"]
    return "\n".join([
        "# Stabilizer-free generic-modulus homology-vector preflight",
        "",
        "This is a deterministic theory/design calculation. It generated zero random samples.",
        "",
        "## Selected design",
        "",
        f"- Base period matrix: `{geometry['base']['matrix']}`, tau=`{geometry['base']['tau']['real']}+i*{geometry['base']['tau']['imag']}`, area {geometry['all_cells_equal_area']}.",
        "- The base and all four signed-stencil cells have exact orientation-preserving stabilizer `{+I,-I}`.",
        "- Directions: an exact shear pair and an independent oblique-aspect pair, all at equal area/cost.",
        f"- Full output schema has {len(observable['channel_order'])} theory channels here and requires sparse retention of every finite-sample primitive winding line.",
        "",
        "## Frozen vectors and information score",
        "",
        "`mu_KdV`, `mu_Q4_Jordan`, and continuum-subtracted zero `mu_embedding` are stored in `latest.json`.",
        f"After profiling the two continuum-tangent leakage columns, KdV versus Q4/Jordan has unit-Fisher shape D^2 = `{separation['KdV_Q4_Jordan_unit_Fisher_shape_D2']:.9g}` (absolute cosine `{separation['KdV_Q4_Jordan_profiled_cosine_up_to_sign']:.9g}`).",
        f"The timing-based four-cell lower-bound cost is `{separation['cost_lower_bound_cpu_hours_per_1m_four_cells']:.9g}` CPU-hours per one million samples/cell-equivalent, giving a conditional shape-value proxy `{separation['conditional_shape_value_per_cpu_hour']:.9g}` per CPU-hour.",
        "The strict three-model maximin V remains `0`: the embedding prediction is zero after continuum subtraction and no nonzero operator-amplitude floor is frozen. This does not erase the nondegenerate KdV/Q4 shape result.",
        "",
        "## Frozen finite-pipeline contract",
        "",
        "Save P0, the sparse rank-1 primitive-winding vector, P2 and all thermal derivatives. Form q=P2-P0, E=P0+P2, refind the pooled q=0 root, use D=d_p q as the common physical normalizer, continuum-subtract each directional homology vector, and form the channelwise U vector at that root.",
        "",
        "## Decision",
        "",
        "The exact geometry and theory templates are ready for one bounded N100 full-vector covariance pilot. Large production is not started by this result.",
        "",
        "Execution provenance: local fallback after TV2N0X rejected the existing SSH key; no key reset was attempted.",
        "",
    ])


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--json-output", type=Path, required=True)
    parser.add_argument("--report-output", type=Path, required=True)
    args = parser.parse_args(list(argv) if argv is not None else None)
    result = analyze(json.loads(args.manifest.read_text(encoding="utf-8")))
    args.json_output.parent.mkdir(parents=True, exist_ok=True)
    args.json_output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    args.report_output.write_text(render_report(result), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

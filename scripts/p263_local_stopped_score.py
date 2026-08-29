#!/usr/bin/env python3
"""Exact identity and tiny variance oracle for the localized #263 Q-score."""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path


CORE_EDGES = ((0, 1), (1, 2), (2, 3), (3, 0))
HIGH_REQUIRED_BITS = ((1, 1), (3, 1), (0, 0), (2, 0))


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def graph_edges(spectator_edges: int) -> tuple[tuple[int, int], ...]:
    return CORE_EDGES + tuple(
        (4 + 2 * index, 5 + 2 * index) for index in range(spectator_edges)
    )


def cluster_score_j(mask: int, spectator_edges: int) -> int:
    edges = graph_edges(spectator_edges)
    vertices = 4 + 2 * spectator_edges
    parent = list(range(vertices))

    def find(vertex: int) -> int:
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    bonds = 0
    for edge_index, (first, second) in enumerate(edges):
        if not ((mask >> edge_index) & 1):
            continue
        bonds += 1
        first_root = find(first)
        second_root = find(second)
        if first_root != second_root:
            parent[second_root] = first_root
    clusters = sum(find(vertex) == vertex for vertex in range(vertices))
    return 2 * clusters + bonds


def high_event(mask: int) -> bool:
    return (mask & 0b1111) == 0b1010


def stopped_transcript(mask: int) -> tuple[tuple[int, int], ...]:
    transcript = []
    for edge, required in HIGH_REQUIRED_BITS:
        value = (mask >> edge) & 1
        transcript.append((edge, value))
        if value != required:
            break
    return tuple(transcript)


def overwrite(mask: int, transcript: tuple[tuple[int, int], ...]) -> int:
    result = mask
    for edge, value in transcript:
        if value:
            result |= 1 << edge
        else:
            result &= ~(1 << edge)
    return result


def _mean(values: list[Fraction]) -> Fraction:
    return sum(values, Fraction()) / len(values)


def _variance(values: list[Fraction]) -> Fraction:
    center = _mean(values)
    return _mean([(value - center) ** 2 for value in values])


def variance_row(spectator_edges: int, inner_replicates: tuple[int, ...]) -> dict:
    edge_count = 4 + spectator_edges
    states = list(range(1 << edge_count))
    score_mean = _mean(
        [Fraction(cluster_score_j(mask, spectator_edges), 2) for mask in states]
    )
    global_values = [
        Fraction(int(high_event(mask)))
        * (Fraction(cluster_score_j(mask, spectator_edges), 2) - score_mean)
        for mask in states
    ]
    target = _mean(global_values)
    global_variance = _variance(global_values)

    stopped_means = []
    stopped_conditional_variances = []
    transcript_lengths = []
    for outer in states:
        transcript = stopped_transcript(outer)
        transcript_lengths.append(len(transcript))
        if not high_event(outer):
            stopped_means.append(Fraction())
            stopped_conditional_variances.append(Fraction())
            continue
        coupled_values = [
            Fraction(
                cluster_score_j(overwrite(base, transcript), spectator_edges)
                - cluster_score_j(base, spectator_edges),
                2,
            )
            for base in states
        ]
        stopped_means.append(_mean(coupled_values))
        stopped_conditional_variances.append(_variance(coupled_values))
    if _mean(stopped_means) != target:
        raise AssertionError("stopped coupled identity failed")
    ideal_variance = _variance(stopped_means)
    completion_noise = _mean(stopped_conditional_variances)
    coupled = {}
    for replicates in inner_replicates:
        variance = ideal_variance + completion_noise / replicates
        coupled[str(replicates)] = {
            "variance": fraction_text(variance),
            "variance_decimal": float(variance),
            "ratio_to_global": float(variance / global_variance),
        }
    return {
        "spectator_edges": spectator_edges,
        "configurations": len(states),
        "target_covariance": fraction_text(target),
        "global_centered_variance": fraction_text(global_variance),
        "global_centered_variance_decimal": float(global_variance),
        "ideal_stopped_variance": fraction_text(ideal_variance),
        "completion_noise_at_one_inner_draw": fraction_text(completion_noise),
        "coupled_inner_replicates": coupled,
        "mean_revealed_core_edges": fraction_text(_mean(
            [Fraction(value) for value in transcript_lengths]
        )),
        "revealed_core_edges_on_high_event": 4,
    }


def render() -> dict:
    rows = [variance_row(spectators, (1, 4, 8, 16)) for spectators in range(5)]
    return {
        "schema": "matching-one.p263-local-stopped-qscore-exact.v1",
        "issue": 263,
        "status": "exact_identity_and_tiny_variance_oracle",
        "score_identity": {
            "global_score": "S=J/2-E[J/2]",
            "stopped_sigma_field": (
                "T is an adaptive p=1/2 edge transcript that determines the event I"
            ),
            "rao_blackwell": "S_T=E[S|T] and Cov(I,S)=E[I*S_T]",
            "coupled_completion": (
                "C is independent p=1/2; C<-T overwrites only transcript edges"
            ),
            "unbiased_estimator": (
                "I(T)*(J(C<-T)-J(C))/2 has expectation Cov(I,J/2)"
            ),
            "local_bound": (
                "if r transcript edges are overwritten, |J(C<-T)-J(C)|/2 <= r/2"
            ),
        },
        "tiny_graph": {
            "core": "open four-cycle with terminals at all four vertices",
            "event": "14|23, core mask 1010 in edge order (12,23,34,41)",
            "adaptive_reveal_order": "23=open,41=open,12=closed,34=closed; stop at first mismatch",
            "spectators": (
                "disjoint Bernoulli edges invisible to the event; they model far-field extensive noise"
            ),
            "rows": rows,
        },
        "variance_statement": {
            "ideal": (
                "I*S_T is the Rao-Blackwell projection of I*S and cannot have larger variance"
            ),
            "finite_inner_replicates": (
                "Var(local_K)=Var(I*S_T)+E[I*Var(delta S|T)]/K"
            ),
            "scaling_conjecture": (
                "far-field area noise cancels exactly; remaining variance follows the marked-arm transcript size and rare-event count"
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = render()
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(encoded, encoding="utf-8")
    else:
        print(encoded, end="")


if __name__ == "__main__":
    main()

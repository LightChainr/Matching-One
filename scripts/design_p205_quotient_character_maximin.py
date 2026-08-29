#!/usr/bin/env python3
"""Design a quotient-character mechanism selector after the P205 10M score.

The search enumerates square Gaussian period matrices at equal area, retains
pairs with different Smith type, and treats H4/H8/H12 as three one-amplitude
character fingerprints.  The maximin target is the smallest bidirectional
noise-whitened distance between any two fingerprint lines.  Noise and runtime
are calibrated from the completed P205 common-field block; no new target data
are read.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Iterable, Mapping, Sequence


HARMONICS = (4, 8, 12)
REFERENCE_ALPHA = Fraction(13, 8)
REFERENCE_AMPLITUDE = 0.7885
BASELINE_SAMPLES = 10_000_000
CHI2_DF2_95 = 5.99146454710798
MODEL_ORDER = ("H4", "H8", "H12")


@dataclass(frozen=True)
class QuotientPair:
    n: int
    first: tuple[int, int]
    second: tuple[int, int]

    @property
    def first_smith(self) -> tuple[int, int]:
        return smith_invariants(self.n, self.first)

    @property
    def second_smith(self) -> tuple[int, int]:
        return smith_invariants(self.n, self.second)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def gaussian_matrix(value: tuple[int, int]) -> tuple[tuple[int, int], tuple[int, int]]:
    a, b = value
    return ((a, -b), (b, a))


def determinant(matrix: Sequence[Sequence[int]]) -> int:
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def smith_invariants(n: int, value: tuple[int, int]) -> tuple[int, int]:
    divisor = math.gcd(abs(value[0]), abs(value[1]))
    return divisor, n // divisor


def cos4(value: tuple[int, int]) -> Fraction:
    a, b = value
    n = a * a + b * b
    return Fraction(a**4 - 6 * a * a * b * b + b**4, n * n)


def harmonic_character(value: tuple[int, int], harmonic: int) -> Fraction:
    x = cos4(value)
    if harmonic == 4:
        return x
    if harmonic == 8:
        return 2 * x * x - 1
    if harmonic == 12:
        return 4 * x * x * x - 3 * x
    raise ValueError("harmonic must be 4, 8, or 12")


def character_difference(pair: QuotientPair, harmonic: int) -> Fraction:
    return harmonic_character(pair.first, harmonic) - harmonic_character(
        pair.second, harmonic
    )


def gaussian_representations(n: int) -> list[tuple[int, int]]:
    rows = []
    for b in range(math.isqrt(n // 2) + 1):
        a = math.isqrt(n - b * b)
        if a >= b and a * a + b * b == n:
            rows.append((a, b))
    return rows


def enumerate_pairs(n_min: int, n_max: int) -> list[QuotientPair]:
    rows = []
    for n in range(n_min, n_max + 1):
        for first, second in itertools.combinations(gaussian_representations(n), 2):
            pair = QuotientPair(n, first, second)
            if pair.first_smith == pair.second_smith:
                continue
            if any(character_difference(pair, harmonic) == 0 for harmonic in HARMONICS):
                continue
            rows.append(pair)
    return rows


def fraction_payload(value: Fraction) -> dict[str, object]:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "exact": str(value),
        "decimal": float(value),
    }


def pair_payload(pair: QuotientPair) -> dict[str, object]:
    return {
        "N": pair.n,
        "first": list(pair.first),
        "second": list(pair.second),
        "first_matrix": [list(row) for row in gaussian_matrix(pair.first)],
        "second_matrix": [list(row) for row in gaussian_matrix(pair.second)],
        "first_smith_invariants": list(pair.first_smith),
        "second_smith_invariants": list(pair.second_smith),
        "character_difference_first_minus_second": {
            f"H{harmonic}": fraction_payload(character_difference(pair, harmonic))
            for harmonic in HARMONICS
        },
        "sign_code": {
            f"H{harmonic}": "+" if character_difference(pair, harmonic) > 0 else "-"
            for harmonic in HARMONICS
        },
    }


def contrast_variance(covariance: Sequence[Sequence[float]], left: int, right: int) -> float:
    return (
        covariance[left][left] + covariance[right][right]
        - covariance[left][right] - covariance[right][left]
    )


def calibrate(score_path: Path, repo_root: Path) -> dict[str, object]:
    score = json.loads(score_path.read_text(encoding="utf-8"))
    if score.get("schema") != "matching-one/p205-norm5-conjugate-coalescence-score/v1":
        raise ValueError("calibration is not the P205 production score")
    variances = []
    for n in ("325", "425"):
        covariance = [
            [float(value) for value in row]
            for row in score["by_size"][n]["raw_M_covariance"]
        ]
        # Geometry order is C,A,B.  The two production jobs are C-A and C-B.
        variances.extend(
            [contrast_variance(covariance, 0, 1), contrast_variance(covariance, 0, 2)]
        )
    pair_variance = math.fsum(variances) / len(variances)

    cpu_per_site = []
    metadata_rows = []
    for row in score["provenance"]["inputs"]:
        path = Path(row["metadata"])
        if not path.is_absolute():
            path = repo_root / path
        metadata = json.loads(path.read_text(encoding="utf-8"))
        n = int(metadata["designs"][0]["N"])
        samples = int(metadata["samples_per_pair"])
        cpu_seconds = float(metadata["elapsed_seconds"]) * int(metadata["threads_requested"])
        cpu_per_site.append(cpu_seconds / (n * samples))
        metadata_rows.append({
            "path": str(path.relative_to(repo_root)),
            "sha256": sha256(path),
            "N": n,
            "samples": samples,
            "elapsed_seconds": float(metadata["elapsed_seconds"]),
            "threads": int(metadata["threads_requested"]),
        })
    model_scores = {
        row["name"]: {
            "joint_chi_square": float(row["joint_chi_square"]),
            "survival": float(row["chi_square_survival_df2"]),
            "delta_chi_square_from_best": float(row["delta_chi_square_from_best"]),
        }
        for row in score["models"]
    }
    return {
        "score_path": str(score_path.relative_to(repo_root)),
        "score_sha256": sha256(score_path),
        "samples": BASELINE_SAMPLES,
        "pair_contrast_variances": variances,
        "mean_pair_contrast_variance": pair_variance,
        "mean_pair_contrast_se": math.sqrt(pair_variance),
        "cpu_seconds_per_site_update": math.fsum(cpu_per_site) / len(cpu_per_site),
        "metadata": metadata_rows,
        "existing_model_scores": model_scores,
    }


def fingerprint(pair: QuotientPair, harmonic: int) -> float:
    return (
        pair.n ** (-float(REFERENCE_ALPHA))
        * float(character_difference(pair, harmonic))
    )


def projection_distance_squared(vector: Sequence[float], line: Sequence[float]) -> float:
    norm_vector = math.fsum(value * value for value in vector)
    norm_line = math.fsum(value * value for value in line)
    if norm_vector == 0.0 or norm_line == 0.0:
        return 0.0
    inner = math.fsum(left * right for left, right in zip(vector, line))
    return max(0.0, norm_vector - inner * inner / norm_line)


def campaign_score(
    pairs: Sequence[QuotientPair], calibration: Mapping[str, object],
    samples: int = BASELINE_SAMPLES,
) -> dict[str, object]:
    variance = float(calibration["mean_pair_contrast_variance"]) * BASELINE_SAMPLES / samples
    vectors = {
        harmonic: [
            REFERENCE_AMPLITUDE * fingerprint(pair, harmonic) / math.sqrt(variance)
            for pair in pairs
        ]
        for harmonic in HARMONICS
    }
    pairwise = {}
    for left, right in itertools.combinations(HARMONICS, 2):
        left_to_right = projection_distance_squared(vectors[left], vectors[right])
        right_to_left = projection_distance_squared(vectors[right], vectors[left])
        pairwise[f"H{left}__H{right}"] = {
            "H_left_true_noncentrality": left_to_right,
            "H_right_true_noncentrality": right_to_left,
            "bidirectional_min_noncentrality": min(left_to_right, right_to_left),
        }
    maximin = min(
        row["bidirectional_min_noncentrality"] for row in pairwise.values()
    )
    cpu_seconds = (
        float(calibration["cpu_seconds_per_site_update"])
        * samples * sum(pair.n for pair in pairs)
    )
    lane_loads = [0, 0]
    for pair in sorted(pairs, key=lambda row: row.n, reverse=True):
        lane = 0 if lane_loads[0] <= lane_loads[1] else 1
        lane_loads[lane] += pair.n
    return {
        "pairs": [pair_payload(pair) for pair in pairs],
        "samples_per_pair": samples,
        "model_fingerprints_N_minus_13_8_times_character": {
            f"H{harmonic}": [fingerprint(pair, harmonic) for pair in pairs]
            for harmonic in HARMONICS
        },
        "noise_whitened_reference_fingerprints": {
            f"H{harmonic}": vectors[harmonic] for harmonic in HARMONICS
        },
        "pairwise_bidirectional_projection": pairwise,
        "maximin_noncentrality": maximin,
        "sqrt_maximin_noncentrality": math.sqrt(maximin),
        "estimated_cpu_seconds": cpu_seconds,
        "estimated_wall_seconds_serial_8_threads": cpu_seconds / 8.0,
        "estimated_wall_seconds_two_8_thread_lanes": (
            float(calibration["cpu_seconds_per_site_update"])
            * samples * max(lane_loads) / 8.0
        ),
        "maximin_noncentrality_per_cpu_second": maximin / cpu_seconds,
        "samples_per_pair_for_df2_95_if_reference_amplitude": (
            samples * CHI2_DF2_95 / maximin if maximin > 0 else math.inf
        ),
    }


def rank_campaigns(
    candidates: Sequence[QuotientPair], count: int,
    calibration: Mapping[str, object], top_n: int = 12,
) -> list[dict[str, object]]:
    scored = []
    for pairs in itertools.combinations(candidates, count):
        row = campaign_score(pairs, calibration)
        if row["maximin_noncentrality"] <= 1e-15:
            continue
        scored.append(row)
    scored.sort(
        key=lambda row: (
            float(row["maximin_noncentrality_per_cpu_second"]),
            float(row["maximin_noncentrality"]),
        ),
        reverse=True,
    )
    return scored[:top_n]


def selected_pair_key(row: Mapping[str, object]) -> tuple[int, tuple[int, int], tuple[int, int]]:
    return int(row["N"]), tuple(row["first"]), tuple(row["second"])


def design_payload(score_path: Path, repo_root: Path, n_max: int = 2000) -> dict[str, object]:
    calibration = calibrate(score_path, repo_root)
    candidates = enumerate_pairs(25, n_max)
    measured = {
        (325, (18, 1), (15, 10)),
        (325, (17, 6), (15, 10)),
        (425, (20, 5), (19, 8)),
        (425, (20, 5), (16, 13)),
    }
    medium_candidates = [
        pair for pair in candidates
        if pair.n >= 325 and (pair.n, pair.first, pair.second) not in measured
    ]
    top_two = rank_campaigns(candidates, 2, calibration)
    top_three = rank_campaigns(candidates, 3, calibration)
    medium_three = rank_campaigns(medium_candidates, 3, calibration)
    selected = top_three[0]

    scale_matched = (
        QuotientPair(400, (20, 0), (16, 12)),
        QuotientPair(450, (21, 3), (15, 15)),
        QuotientPair(500, (22, 4), (20, 10)),
    )
    scale_matched_row = campaign_score(scale_matched, calibration)

    return {
        "schema": "matching-one/p205-quotient-character-maximin-design/v1",
        "status": "design_only_after_P205_10M__before_new_quotient_targets",
        "purpose": (
            "Select equal-area cross-Smith quotient pairs whose H4/H8/H12 "
            "character fingerprints are maximally separated per measured CPU cost."
        ),
        "calibration": calibration,
        "model": {
            "harmonic_order": list(MODEL_ORDER),
            "radial_exponent_in_N": str(REFERENCE_ALPHA),
            "reference_amplitude_for_planning_only": REFERENCE_AMPLITUDE,
            "nuisance_parameter_in_score": "one shared amplitude per harmonic model",
            "observable": "DeltaM_N = M_first(p_ref)-M_second(p_ref)",
            "model_vector": "DeltaM_N = A_s * N^(-13/8) * Delta cos(s theta)",
        },
        "enumeration": {
            "N_min": 25,
            "N_max": n_max,
            "canonical_square_Gaussian_representations_only": True,
            "require_different_Smith_invariants": True,
            "candidate_pair_count": len(candidates),
        },
        "top_two_pair_campaigns": top_two,
        "top_three_pair_campaigns": top_three,
        "selected_three_pair_character_prism": selected,
        "selected_sign_code": {
            f"H{harmonic}": [
                row["sign_code"][f"H{harmonic}"] for row in selected["pairs"]
            ]
            for harmonic in HARMONICS
        },
        "medium_N_325_plus_best_new_three_pair_campaigns": medium_three,
        "scale_matched_N400_450_500_sign_code_bridge": scale_matched_row,
        "existing_P205_limitation": {
            "H4_H8_delta_chi_square": (
                calibration["existing_model_scores"]["H4"]["joint_chi_square"]
                - calibration["existing_model_scores"]["H8"]["joint_chi_square"]
            ),
            "H12_H8_delta_chi_square": (
                calibration["existing_model_scores"]["H12"]["joint_chi_square"]
                - calibration["existing_model_scores"]["H8"]["joint_chi_square"]
            ),
            "reason": (
                "N325/N425 each supplied only one affine-null coordinate, the raw "
                "M noise scale was comparable to the fixed-model residuals, and all "
                "three harmonics survived. The character prism adds orthogonal sign "
                "flips instead of only shrinking the same covariance."
            ),
        },
        "minimal_sufficient_statistics": {
            "per_pair": (
                "For each synchronized batch and each orientation, the marginal "
                "K_plus and K_minus threshold-rank histograms, or equivalently the "
                "two fixed-p M values after a frozen binomial transform."
            ),
            "not_required_for_primary": [
                "root", "derivative", "joint K_plus*K_minus moment",
                "full-curve refit", "fitted exponent", "quotient offset",
            ],
            "covariance": (
                "Common priority field within each equal-N pair; independent seeds "
                "across N25/N50/N125; one 3-vector covariance is therefore diagonal "
                "across sizes after the within-pair DeltaM construction."
            ),
        },
    }


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--p205-score", type=Path,
        default=root / "results/server-20260829/P205-norm5-conjugate-coalescence/analysis/score.json",
    )
    parser.add_argument("--n-max", type=int, default=2000)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = design_payload(args.p205_score, root, args.n_max)
    text = json.dumps(payload, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
        print(args.output)
    else:
        print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

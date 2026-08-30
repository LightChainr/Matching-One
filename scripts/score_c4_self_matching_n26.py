#!/usr/bin/env python3
"""Score the two frozen N=26 Beta laws and report exact structure.

The scoring order is read from the pre-target artifact.  A failed fixed target
is only described; this program contains no generalized-Beta fitting path.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from typing import Sequence


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def vector_sha256(values: Sequence[int]) -> str:
    canonical = json.dumps(list(values), separators=(",", ":"))
    return hashlib.sha256(canonical.encode("ascii")).hexdigest()


def bernstein_to_power(counts: Sequence[int]) -> list[int]:
    n = len(counts) - 1
    power = [0] * (n + 1)
    for k, value in enumerate(counts):
        for degree in range(k, n + 1):
            power[degree] += (
                value
                * (-1) ** (degree - k)
                * math.comb(n - k, degree - k)
            )
    while len(power) > 1 and power[-1] == 0:
        power.pop()
    return power


def _deterministic_enumeration_payload(payload: dict[str, object]) -> dict[str, object]:
    return {
        key: value
        for key, value in payload.items()
        if key not in {"threads", "elapsed_seconds"}
    }


def _factorization(power: Sequence[int]) -> str | None:
    try:
        import sympy as sp
    except ImportError:
        return None
    p = sp.symbols("p")
    expression = sum(sp.Integer(value) * p**degree for degree, value in enumerate(power))
    return str(sp.factor(expression))


def score(
    prediction_path: Path,
    primary_path: Path,
    reproduction_path: Path,
) -> dict[str, object]:
    prediction = json.loads(prediction_path.read_text(encoding="utf-8"))
    primary = json.loads(primary_path.read_text(encoding="utf-8"))
    reproduction = json.loads(reproduction_path.read_text(encoding="utf-8"))
    if prediction["status"] != "FROZEN_BEFORE_N26_ENUMERATION":
        raise ValueError("prediction artifact is not marked pre-target frozen")
    if prediction["post_failure_rule"] != "STOP_WITHOUT_GENERALIZED_BETA_FIT":
        raise ValueError("post-failure no-fit rule changed")
    if primary["geometry"]["N"] != 26 or reproduction["geometry"]["N"] != 26:
        raise ValueError("both target files must be N=26")
    reproduction_identical = (
        _deterministic_enumeration_payload(primary)
        == _deterministic_enumeration_payload(reproduction)
    )
    if not reproduction_identical:
        raise ValueError("independent enumeration coefficients disagree")

    channel = prediction["geometry"]["wrapping_channel"]
    actual = primary["channels"][channel]["M_bernstein_integer_coefficients"]
    actual_power = bernstein_to_power(actual)
    if actual_power != primary["channels"][channel]["M_power_coefficients_ascending"]:
        raise ValueError("independent Bernstein-to-power conversion disagrees")

    hypotheses_by_name = {
        item["name"]: item for item in prediction["hypotheses"]
    }
    hypothesis_scores = []
    for name in prediction["scoring_order"]:
        hypothesis = hypotheses_by_name[name]
        target = hypothesis["bernstein_integer_coefficients"]
        if vector_sha256(target) != hypothesis["bernstein_vector_sha256"]:
            raise ValueError(f"frozen vector hash mismatch for {name}")
        differences = [observed - expected for observed, expected in zip(actual, target)]
        nonzero = [index for index, value in enumerate(differences) if value]
        first = nonzero[0] if nonzero else None
        target_power = bernstein_to_power(target)
        difference_power = bernstein_to_power(differences)
        hypothesis_scores.append(
            {
                "name": name,
                "beta_parameters": hypothesis["beta_parameters"],
                "exact_pass": not nonzero,
                "differing_coefficient_count": len(nonzero),
                "first_difference": None
                if first is None
                else {
                    "occupation_k": first,
                    "observed": actual[first],
                    "target": target[first],
                    "observed_minus_target": differences[first],
                },
                "differences_observed_minus_target": differences,
                "target_power_degree": len(target_power) - 1,
                "difference_power_factorization_over_Q": _factorization(
                    difference_power
                ),
            }
        )

    all_channel_vectors = {
        name: row["M_bernstein_integer_coefficients"]
        for name, row in primary["channels"].items()
    }
    channel_independent = len({tuple(row) for row in all_channel_vectors.values()}) == 1
    n = len(actual) - 1
    f_numerators = [math.comb(n, k) + actual[k] for k in range(n + 1)]
    if any(value % 2 for value in f_numerators):
        raise AssertionError("F=(1+M)/2 has nonintegral Bernstein numerator")
    f_counts = [value // 2 for value in f_numerators]
    first_f_support = next(k for k, value in enumerate(f_counts) if value)
    beta_f_at_first = {}
    for hypothesis in prediction["hypotheses"]:
        vector = hypothesis["bernstein_integer_coefficients"]
        beta_f_at_first[hypothesis["name"]] = (
            math.comb(n, first_f_support) + vector[first_f_support]
        ) // 2

    exact_structure = {
        "M_bernstein_integer_coefficients": actual,
        "M_bernstein_vector_sha256": vector_sha256(actual),
        "M_coefficients_anti_palindromic": all(
            actual[k] == -actual[n - k] for k in range(n + 1)
        ),
        "M_identical_across_all_five_wrapping_channels": channel_independent,
        "M_power_coefficients_ascending": actual_power,
        "M_power_degree": len(actual_power) - 1,
        "M_power_factorization_over_Q": _factorization(actual_power),
        "F_first_nonzero_occupation": first_f_support,
        "F_bernstein_numerator_at_first_support": f_counts[first_f_support],
        "frozen_target_F_numerators_at_same_occupation": beta_f_at_first,
        "raw_either_wrap_count_at_first_support": primary["channels"]["either"][
            "R_bernstein_integer_coefficients"
        ][first_f_support],
        "raw_direction_wrap_counts_at_first_support": {
            "direction_0": primary["channels"]["direction_0"][
                "R_bernstein_integer_coefficients"
            ][first_f_support],
            "direction_1": primary["channels"]["direction_1"][
                "R_bernstein_integer_coefficients"
            ][first_f_support],
        },
    }
    all_pass = all(item["exact_pass"] for item in hypothesis_scores)
    return {
        "schema": "matching-one/c4-self-matching-n26-score/v1",
        "inputs": {
            "prediction": str(prediction_path),
            "prediction_sha256": file_sha256(prediction_path),
            "primary": str(primary_path),
            "primary_sha256": file_sha256(primary_path),
            "reproduction": str(reproduction_path),
            "reproduction_sha256": file_sha256(reproduction_path),
        },
        "independent_reproduction_identical": reproduction_identical,
        "scoring_order": prediction["scoring_order"],
        "hypothesis_scores": hypothesis_scores,
        "all_frozen_hypotheses_pass": all_pass,
        "protocol_conclusion": (
            "FROZEN_EXACT_LAW_PASSED"
            if all_pass
            else "BOTH_FROZEN_LAWS_FAILED_STOP_NO_GENERALIZED_BETA_FIT"
        ),
        "generalized_beta_fit_performed": False,
        "exact_structure_after_fixed_scoring": exact_structure,
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prediction", type=Path, required=True)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--reproduction", type=Path, required=True)
    parser.add_argument("--json", type=Path)
    args = parser.parse_args(argv)
    payload = score(args.prediction, args.primary, args.reproduction)
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

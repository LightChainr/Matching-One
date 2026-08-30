#!/usr/bin/env python3
"""Reconstruct essential-H1 birth observables from a committed rank archive."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path
from typing import Any, Mapping

import rigorous_pc_confidence_gate as exact_probability


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "results/server-20260828/C05/axis_L8_pilot"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_checksums(path: Path) -> dict[str, str]:
    result = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        digest, name = line.split("  ", 1)
        result[name] = digest
    return result


def read_marginal(path: Path, field: str) -> dict[int, int]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = csv.DictReader(handle)
        return {int(row[field]): int(row["count"]) for row in rows}


def read_joint(path: Path) -> dict[tuple[int, int], int]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = csv.DictReader(handle)
        return {
            (int(row["K_minus"]), int(row["K_plus"])): int(row["count"])
            for row in rows
        }


def marginalize_joint(joint: Mapping[tuple[int, int], int]) -> tuple[dict[int, int], dict[int, int]]:
    minus = {}
    plus = {}
    for (k_minus, k_plus), count in joint.items():
        minus[k_minus] = minus.get(k_minus, 0) + count
        plus[k_plus] = plus.get(k_plus, 0) + count
    return minus, plus


def joint_moments(joint: Mapping[tuple[int, int], int]) -> dict[str, int]:
    result = {
        "sum_K_minus": 0,
        "sum_K_plus": 0,
        "sum_K_minus_squared": 0,
        "sum_K_plus_squared": 0,
        "sum_K_minus_times_K_plus": 0,
        "sum_gap": 0,
        "sum_gap_squared": 0,
    }
    for (k_minus, k_plus), count in joint.items():
        gap = k_plus - k_minus
        result["sum_K_minus"] += count * k_minus
        result["sum_K_plus"] += count * k_plus
        result["sum_K_minus_squared"] += count * k_minus * k_minus
        result["sum_K_plus_squared"] += count * k_plus * k_plus
        result["sum_K_minus_times_K_plus"] += count * k_minus * k_plus
        result["sum_gap"] += count * gap
        result["sum_gap_squared"] += count * gap * gap
    return result


def _decimal_text(value: Decimal, digits: int = 50) -> str:
    return format(value, ".%dE" % (digits - 1))


def _binomial_tail_decimal(n: int, cutoff: int, p: Decimal) -> Decimal:
    return sum(
        (
            Decimal(math.comb(n, occupied))
            * p**occupied
            * (1 - p) ** (n - occupied)
            for occupied in range(cutoff, n + 1)
        ),
        Decimal(0),
    )


def _order_density_decimal(n: int, rank: int, p: Decimal) -> Decimal:
    return (
        Decimal(n)
        * math.comb(n - 1, rank - 1)
        * p ** (rank - 1)
        * (1 - p) ** (n - rank)
    )


def _decimal_ln(value: Decimal) -> Decimal:
    """Natural logarithm by the atanh series, sufficient near balance."""

    if value <= 0:
        raise ValueError("logarithm requires a positive value")
    z = (value - 1) / (value + 1)
    z_squared = z * z
    term = z
    total = term
    denominator = 1
    threshold = Decimal(1).scaleb(-70)
    while True:
        term *= z_squared
        denominator += 2
        addition = term / denominator
        total += addition
        if abs(addition) < threshold:
            return 2 * total


def evaluate_decimal(
    n: int,
    sample_count: int,
    minus: Mapping[int, int],
    plus: Mapping[int, int],
    p_text: str,
) -> dict[str, str]:
    with localcontext() as context:
        context.prec = 80
        p = Decimal(p_text)
        f1 = sum(
            (count * _binomial_tail_decimal(n, rank, p) for rank, count in minus.items()),
            Decimal(0),
        ) / sample_count
        f2 = sum(
            (count * _binomial_tail_decimal(n, rank, p) for rank, count in plus.items()),
            Decimal(0),
        ) / sample_count
        density1 = sum(
            (count * _order_density_decimal(n, rank, p) for rank, count in minus.items()),
            Decimal(0),
        ) / sample_count
        density2 = sum(
            (count * _order_density_decimal(n, rank, p) for rank, count in plus.items()),
            Decimal(0),
        ) / sample_count
        probability0 = 1 - f1
        probability1 = f1 - f2
        probability2 = f2
        matching = probability2 - probability0
        mixture_cdf = (f1 + f2) / 2
        rho = (density1 + density2) / 2
        phi = _decimal_ln(probability2 / probability0)
        return {
            "p": p_text,
            "F_first": _decimal_text(f1),
            "F_second": _decimal_text(f2),
            "P_rank_0": _decimal_text(probability0),
            "P_rank_1": _decimal_text(probability1),
            "P_rank_2": _decimal_text(probability2),
            "M": _decimal_text(matching),
            "essential_birth_mixture_cdf": _decimal_text(mixture_cdf),
            "f_first": _decimal_text(density1),
            "f_second": _decimal_text(density2),
            "rho_equals_M_prime_over_2": _decimal_text(rho),
            "M_prime": _decimal_text(density1 + density2),
            "rank_one_derivative": _decimal_text(density1 - density2),
            "Phi_log_P2_over_P0": _decimal_text(phi),
        }


def evaluate_half_exact(
    n: int,
    sample_count: int,
    minus: Mapping[int, int],
    plus: Mapping[int, int],
) -> dict[str, str]:
    p = Fraction(1, 2)
    f1 = sum(
        (count * exact_probability.binomial_tail(n, rank, p) for rank, count in minus.items()),
        Fraction(0),
    ) / sample_count
    f2 = sum(
        (count * exact_probability.binomial_tail(n, rank, p) for rank, count in plus.items()),
        Fraction(0),
    ) / sample_count
    values = {
        "F_first": f1,
        "F_second": f2,
        "P_rank_0": 1 - f1,
        "P_rank_1": f1 - f2,
        "P_rank_2": f2,
        "M": f1 + f2 - 1,
        "essential_birth_mixture_cdf": (f1 + f2) / 2,
    }
    return {key: exact_probability.fraction_text(value) for key, value in values.items()}


def priority_moments(n: int, sample_count: int, moments: Mapping[str, int]) -> dict[str, str]:
    scale1 = sample_count * (n + 1)
    scale2 = sample_count * (n + 1) * (n + 2)
    mean1 = Fraction(moments["sum_K_minus"], scale1)
    mean2 = Fraction(moments["sum_K_plus"], scale1)
    second1 = Fraction(
        moments["sum_K_minus_squared"] + moments["sum_K_minus"], scale2
    )
    second2 = Fraction(
        moments["sum_K_plus_squared"] + moments["sum_K_plus"], scale2
    )
    cross = Fraction(
        moments["sum_K_minus_times_K_plus"] + moments["sum_K_minus"], scale2
    )
    mean_center = (mean1 + mean2) / 2
    mean_lifetime = mean2 - mean1
    second_center = (second1 + 2 * cross + second2) / 4
    second_lifetime = second1 - 2 * cross + second2
    values = {
        "E_tau_first": mean1,
        "E_tau_second": mean2,
        "E_center": mean_center,
        "E_lifetime": mean_lifetime,
        "Var_center": second_center - mean_center * mean_center,
        "Var_lifetime": second_lifetime - mean_lifetime * mean_lifetime,
        "neutral_area_integral": Fraction(moments["sum_gap"], scale1),
    }
    return {
        key: exact_probability.fraction_text(value)
        for key, value in values.items()
    } | {
        key + "_decimal": exact_probability.decimal_text(value)
        for key, value in values.items()
    }


def build_artifact() -> dict[str, Any]:
    metadata = json.loads((ARCHIVE / "metadata.json").read_text(encoding="utf-8"))
    derived = json.loads((ARCHIVE / "derived_summary.json").read_text(encoding="utf-8"))
    expected_checksums = read_checksums(ARCHIVE / "checksums.sha256")
    checksum_rows = {
        name: {"expected": digest, "observed": _sha256(ARCHIVE / name)}
        for name, digest in expected_checksums.items()
    }
    minus = read_marginal(ARCHIVE / "kminus_hist.csv", "K_minus")
    plus = read_marginal(ARCHIVE / "kplus_hist.csv", "K_plus")
    joint = read_joint(ARCHIVE / "joint_hist.csv")
    joint_minus, joint_plus = marginalize_joint(joint)
    moments = joint_moments(joint)
    n = int(metadata["N"])
    sample_count = int(metadata["sample_count"])

    support_violations = [
        {"K_minus": pair[0], "K_plus": pair[1], "count": count}
        for pair, count in joint.items()
        if pair[0] > pair[1]
    ]
    root_evaluation = evaluate_decimal(n, sample_count, minus, plus, derived["root"])
    declared_evaluation = evaluate_decimal(
        n, sample_count, minus, plus, derived["evaluations"][0]["p"]
    )
    with localcontext() as context:
        context.prec = 80
        declared_m_residual = Decimal(declared_evaluation["M"]) - Decimal(
            derived["evaluations"][0]["M"]
        )
        declared_derivative_residual = Decimal(declared_evaluation["M_prime"]) - Decimal(
            derived["evaluations"][0]["M_prime"]
        )

    assert all(row["expected"] == row["observed"] for row in checksum_rows.values())
    assert sum(minus.values()) == sample_count
    assert sum(plus.values()) == sample_count
    assert sum(joint.values()) == sample_count
    assert minus == joint_minus
    assert plus == joint_plus
    assert moments == metadata["first_second_joint_integer_moments"]
    assert not support_violations
    assert abs(Decimal(root_evaluation["M"])) < Decimal("1e-45")
    assert abs(declared_m_residual) < Decimal("1e-45")
    assert abs(declared_derivative_residual) < Decimal("1e-45")

    priority = priority_moments(n, sample_count, moments)
    assert priority["E_lifetime"] == priority["neutral_area_integral"]
    return {
        "schema": "matching-one/essential-birth-histogram-analysis/v1",
        "issue": 269,
        "status": "no_new_production_reinterpretation",
        "source_archive": str(ARCHIVE.relative_to(ROOT)),
        "source_sample_count": sample_count,
        "source_N": n,
        "archive_validation": {
            "checksums": checksum_rows,
            "marginal_totals": {"K_minus": sum(minus.values()), "K_plus": sum(plus.values())},
            "joint_total": sum(joint.values()),
            "joint_reconstructs_both_marginals": minus == joint_minus and plus == joint_plus,
            "moments_match_metadata": moments == metadata["first_second_joint_integer_moments"],
            "K_minus_above_K_plus_support_violations": support_violations,
        },
        "exact_at_p_half": evaluate_half_exact(n, sample_count, minus, plus),
        "archived_root_evaluation": root_evaluation,
        "archived_declared_point_reproduction": {
            "computed": declared_evaluation,
            "declared": derived["evaluations"][0],
            "M_residual": _decimal_text(declared_m_residual),
            "M_prime_residual": _decimal_text(declared_derivative_residual),
        },
        "priority_center_lifetime_moments": priority,
        "identities": {
            "rank_probabilities": "P0=1-F_first, P1=F_first-F_second, P2=F_second",
            "matching_curve": "M=P2-P0=F_first+F_second-1",
            "canonical_cdf": "(1+M)/2=(F_first+F_second)/2",
            "density": "rho=M_prime/2=(f_first+f_second)/2",
            "neutral_area": "integral_0^1 P1(p)dp=E[tau_second-tau_first]=E[(K_plus-K_minus)/(N+1)]",
        },
        "recoverability_boundary": {
            "available_without_new_production": [
                "both marginal birth CDFs and densities",
                "rank-0/1/2 probabilities",
                "joint discrete center/gap distribution and priority center/lifetime moments",
                "topological balance ratio P2/P0",
            ],
            "missing_from_archive": [
                "projective winding line ell",
                "integral saturation index",
                "first/second birth-site local marks",
            ],
        },
        "scientific_boundary": (
            "the 100000 committed pilot samples are reinterpreted, not new independent evidence"
        ),
    }


def render_markdown(artifact: Mapping[str, Any]) -> str:
    validation = artifact["archive_validation"]
    root = artifact["archived_root_evaluation"]
    priority = artifact["priority_center_lifetime_moments"]
    lines = [
        "# Essential-birth reconstruction from the committed axis L=8 pilot",
        "",
        "This is a no-new-production reinterpretation of 100,000 existing samples.",
        "",
        "## Archive audit",
        "",
        "- every declared SHA-256 checksum matches;",
        "- both marginal totals and the joint total equal `%d`;" % artifact["source_sample_count"],
        "- the joint table reconstructs both marginals: `%s`;" % validation["joint_reconstructs_both_marginals"],
        "- all joint integer moments reproduce metadata: `%s`;" % validation["moments_match_metadata"],
        "- support rows with `K_minus>K_plus`: `%d`." % len(validation["K_minus_above_K_plus_support_violations"]),
        "",
        "## Homology-birth interpretation",
        "",
        "The two threshold histograms are the first- and second-essential-birth distributions.",
        "They reconstruct `P0=1-F1`, `P1=F1-F2`, `P2=F2`, and `M=P2-P0` exactly.",
        "At the archived finite root `%s`, the recomputed `M` is `%s` and the equal-weight birth"
        % (root["p"], root["M"]),
        "mixture CDF is `%s`." % root["essential_birth_mixture_cdf"],
        "",
        "The neutral-area identity becomes an exact priority lifetime:",
        "",
        "```text",
        "integral P(R=1) dp = E[tau_second-tau_first]",
        "                    = %s" % priority["E_lifetime"],
        "                    = %s (decimal)" % priority["E_lifetime_decimal"],
        "```",
        "",
        "The joint archive also determines `E[C]`, `Var(C)`, and `Var(W)` for",
        "`C=(tau_first+tau_second)/2` and `W=tau_second-tau_first`.",
        "",
        "## Missing marks",
        "",
        "The historical files do not contain the projective winding line `ell`, integral saturation",
        "index, or first/second birth-site local marks. Those quantities require a future stream; they",
        "cannot be reconstructed from marginal or joint endpoint counts.",
        "",
        "This analysis reuses the same pilot data and is not new independent evidence.",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--format", choices=("json", "markdown"), default="json")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    artifact = build_artifact()
    rendered = (
        json.dumps(artifact, indent=2, sort_keys=True) + "\n"
        if args.format == "json"
        else render_markdown(artifact)
    )
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()

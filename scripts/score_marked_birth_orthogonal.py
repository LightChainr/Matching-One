#!/usr/bin/env python3
"""Batch-Gram orthogonalization of the production marked-birth source.

The same lifted line carries two gate characters,

    J_D = chi4(ell) (I12-I01),   J_S = chi4(ell) (I12+I01).

We use the even line-gate source J_S as the explicit thermal control.  The
coefficient alpha is trained only on N=65 batch fluctuations and N=130 is a
held-out q2 score.  This is an estimator-level Gram projection; it is not
mislabelled as the unavailable configuration-level source Gram matrix.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any, Sequence

import mpmath as mp

from score_marked_birth_path import (
    PathRow,
    combine,
    cos4,
    projected,
    read_path,
)


@dataclass
class ComplexScore:
    p: mp.mpf
    J_D: mp.mpc
    J_S: mp.mpc
    C_D: mp.mpc
    C_S: mp.mpc
    contact_residual_D: mp.mpf
    contact_residual_S: mp.mpf


def _complex(values: dict[str, mp.mpf], stem: str) -> mp.mpc:
    return mp.mpc(values[stem + "_re"], values[stem + "_im"])


def _orientation(values: dict[str, mp.mpf], p: mp.mpf) -> dict[str, mp.mpc]:
    q = values["A_top"]
    jd = _complex(values, "J_D")
    js = _complex(values, "J_S")
    environment_qjd = _complex(values, "environment_q_J_D")
    cd = environment_qjd + p * jd - q * jd
    # Gate algebra: q_before*J_S = -q_before*J_D and S*J_S=J_S.
    cs = -environment_qjd + p * js - q * js
    exact_cd = mp.mpf("0.5") * js + (p - mp.mpf("0.5") - q) * jd
    exact_cs = mp.mpf("0.5") * jd + (p - mp.mpf("0.5") - q) * js
    return {
        "J_D": jd,
        "J_S": js,
        "C_D": cd,
        "C_S": cs,
        "contact_D": cd - exact_cd,
        "contact_S": cs - exact_cs,
    }


def score_pair(first: Sequence[PathRow], second: Sequence[PathRow]) -> ComplexScore:
    p, _, orientations = projected(first, second)
    delta = cos4(first[0].a, first[0].b) - cos4(second[0].a, second[0].b)
    left = _orientation(orientations["first"], p)
    right = _orientation(orientations["second"], p)
    return ComplexScore(
        p=p,
        J_D=(left["J_D"] - right["J_D"]) / delta,
        J_S=(left["J_S"] - right["J_S"]) / delta,
        C_D=(left["C_D"] - right["C_D"]) / delta,
        C_S=(left["C_S"] - right["C_S"]) / delta,
        contact_residual_D=max(abs(left["contact_D"]), abs(right["contact_D"])),
        contact_residual_S=max(abs(left["contact_S"]), abs(right["contact_S"])),
    )


def load_scores(prefix: Path) -> tuple[ComplexScore, list[ComplexScore]]:
    groups = read_path(Path(str(prefix) + ".path.csv"))
    n = next(iter(groups))[0]
    batches = sorted(
        set(key[2] for key in groups if key[1] == "first")
        & set(key[2] for key in groups if key[1] == "second")
    )
    first = combine([groups[(n, "first", batch)] for batch in batches])
    second = combine([groups[(n, "second", batch)] for batch in batches])
    point = score_pair(first, second)
    per_batch = [
        score_pair(groups[(n, "first", batch)], groups[(n, "second", batch)])
        for batch in batches
    ]
    return point, per_batch


def gram_alpha(rows: Sequence[ComplexScore]) -> mp.mpc:
    mean_j = mp.fsum(row.J_D for row in rows) / len(rows)
    mean_t = mp.fsum(row.J_S for row in rows) / len(rows)
    numerator = mp.fsum(
        mp.conj(row.J_S - mean_t) * (row.J_D - mean_j) for row in rows
    )
    denominator = mp.fsum(abs(row.J_S - mean_t) ** 2 for row in rows)
    if denominator == 0:
        raise ValueError("thermal-control batch variance vanished")
    return numerator / denominator


def _mean(values: Sequence[mp.mpc]) -> mp.mpc:
    return mp.fsum(values) / len(values)


def _jackknife_se(values: Sequence[mp.mpc]) -> tuple[mp.mpf, mp.mpf]:
    mean = _mean(values)
    factor = mp.mpf(len(values) - 1) / len(values)
    return (
        mp.sqrt(factor * mp.fsum((mp.re(value - mean)) ** 2 for value in values)),
        mp.sqrt(factor * mp.fsum((mp.im(value - mean)) ** 2 for value in values)),
    )


def _text(value: mp.mpf | mp.mpc) -> str:
    return mp.nstr(value, 24)


def _complex_json(value: mp.mpc) -> dict[str, str]:
    return {"re": _text(mp.re(value)), "im": _text(mp.im(value)), "abs": _text(abs(value))}


def _score_json(point: ComplexScore, alpha: mp.mpc) -> dict[str, Any]:
    return {
        "p": _text(point.p),
        "raw_source_J_D": _complex_json(point.J_D),
        "thermal_control_J_S": _complex_json(point.J_S),
        "raw_connected_C_D": _complex_json(point.C_D),
        "thermal_connected_C_S": _complex_json(point.C_S),
        "orthogonal_source_J_perp": _complex_json(point.J_D - alpha * point.J_S),
        "orthogonal_connected_C_perp": _complex_json(point.C_D - alpha * point.C_S),
        "contact_identity_max_residual": {
            "D": _text(point.contact_residual_D),
            "S": _text(point.contact_residual_S),
        },
    }


def build_report(n65_prefix: Path, n130_prefix: Path, p50_prefix: Path) -> dict[str, Any]:
    n65, batches65 = load_scores(n65_prefix)
    n130, batches130 = load_scores(n130_prefix)
    p50, batches50 = load_scores(p50_prefix)
    alpha = gram_alpha(batches65)
    alpha_delete = [
        gram_alpha([row for index, row in enumerate(batches65) if index != omitted])
        for omitted in range(len(batches65))
    ]
    alpha_se = _jackknife_se(alpha_delete)

    # Each training batch is scored using an alpha that did not see that batch.
    train_crossfit_j = [
        row.J_D - alpha_delete[index] * row.J_S
        for index, row in enumerate(batches65)
    ]
    train_crossfit_c = [
        row.C_D - alpha_delete[index] * row.C_S
        for index, row in enumerate(batches65)
    ]

    centered_t = [row.J_S - _mean([item.J_S for item in batches65]) for row in batches65]
    centered_residual = [
        row.J_D - _mean([item.J_D for item in batches65]) - alpha * t
        for row, t in zip(batches65, centered_t)
    ]
    gram_residual = mp.fsum(
        mp.conj(t) * residual for t, residual in zip(centered_t, centered_residual)
    )

    heldout_j = n130.J_D - alpha * n130.J_S
    heldout_c = n130.C_D - alpha * n130.C_S
    p50_j = p50.J_D - alpha * p50.J_S
    p50_c = p50.C_D - alpha * p50.C_S
    train_j = n65.J_D - alpha * n65.J_S
    train_c = n65.C_D - alpha * n65.C_S

    raw_source_transfer = n130.J_D / n65.J_D
    raw_connected_transfer = n130.C_D / mp.conj(n65.C_D)
    residual_source_transfer = heldout_j / train_j
    residual_connected_transfer = heldout_c / mp.conj(train_c)

    return {
        "schema": "matching-one/marked-birth-batch-gram-orthogonal/v1",
        "issues": [215, 275, 276],
        "training": "N65 only",
        "held_out": "N130 q2 child",
        "external_direction": "P50 N145",
        "thermal_control": "J_S4=chi4(ell)(I01+I12)",
        "gram_definition": "alpha=<delta J_S,delta J_D>/<delta J_S,delta J_S> over 20 aligned N65 batch estimators",
        "claim_boundary": {
            "batch_gram": "estimator-level projection, not a configuration-level source Gram matrix",
            "exact_contact": "C_D=J_S/2+(p-1/2-<q>)J_D and C_S=J_D/2+(p-1/2-<q>)J_S orientation by orientation",
            "missing_for_field_gram": "per-configuration or two-root Horvitz sums of J_D*conj(J_S) and |J_S|^2",
        },
        "alpha": _complex_json(alpha),
        "alpha_delete_one_se": {"re": _text(alpha_se[0]), "im": _text(alpha_se[1])},
        "gram_normal_equation_residual": _complex_json(gram_residual),
        "scores": {
            "N65_training": _score_json(n65, alpha),
            "N130_heldout": _score_json(n130, alpha),
            "P50_N145_external": _score_json(p50, alpha),
        },
        "training_crossfit": {
            "mean_J_perp": _complex_json(_mean(train_crossfit_j)),
            "mean_C_perp": _complex_json(_mean(train_crossfit_c)),
            "J_perp_batch_se": dict(zip(("re", "im"), map(_text, _jackknife_se(train_crossfit_j)))),
            "C_perp_batch_se": dict(zip(("re", "im"), map(_text, _jackknife_se(train_crossfit_c)))),
        },
        "transfers": {
            "raw_source_N130_over_N65": _complex_json(raw_source_transfer),
            "raw_connected_N130_over_conj_N65": _complex_json(raw_connected_transfer),
            "orthogonal_source_N130_over_N65": _complex_json(residual_source_transfer),
            "orthogonal_connected_N130_over_conj_N65": _complex_json(residual_connected_transfer),
            "targets": {
                "H4_mean_source": "-2^(-13/8)=-0.3242098886627524",
                "thermal_connected": "2^(3/8)=1.2968395546510096 after q2 conjugation",
            },
        },
        "conclusion": (
            "The N65 batch-Gram projection does not remove the held-out thermal growth. "
            "More strongly, the connected response is exactly closed by gate-source means, "
            "so it is not an independent Q4 matrix element."
        ),
        "minimal_additional_statistics": [
            "two-root Horvitz sum of J_D4*conj(J_S4) per k/batch",
            "two-root Horvitz sum of |J_S4|^2 per k/batch",
            "an independent global observer O_ext and O_ext*J_D4/O_ext*J_S4 products",
        ],
        "scientific_card": [
            "MECHANISM SPACE: the apparent thermal growth is an exact rank-gate contact channel, not free evidence for a Q4 matrix element.",
            "NOT PROVED: the noisy mean-J_D q2 exponent remains a candidate; batch Gram is not the missing field-level Gram.",
            "OBSERVER-SECTOR-SOURCE-GEOMETRY: A_top | Alexander odd | J_D/J_S gate doublet | q2 Gaussian parent and held-out child.",
            "DEPENDENCY GROUP: orthogonal and raw scores are deterministic transforms of the same 634040d pilot.",
            "UPWEIGHT OBSERVATION: mean-J_D radial transfer or an independent-observer coupling after two-root field orthogonalization.",
        ],
    }


def render_markdown(report: dict[str, Any]) -> str:
    alpha = report["alpha"]
    transfers = report["transfers"]
    return "\n".join(
        [
            "# Held-out batch-Gram orthogonalization",
            "",
            f"N65 trains `alpha={alpha['re']}+({alpha['im']})i` using the centered complex batch Gram normal equation.",
            "N130 is held out; P50 N145 is external direction only.",
            "",
            "The projection does not remove thermal growth:",
            "",
            f"- raw connected q2 transfer: `{transfers['raw_connected_N130_over_conj_N65']}`;",
            f"- orthogonal connected q2 transfer: `{transfers['orthogonal_connected_N130_over_conj_N65']}`.",
            "",
            "The reason is stronger than a weak regression. Gate algebra gives exactly, orientation by orientation,",
            "",
            "```text",
            "Cov(q,J_D)=J_S/2+(p-1/2-<q>)J_D,",
            "Cov(q,J_S)=J_D/2+(p-1/2-<q>)J_S.",
            "```",
            "",
            "Thus the connected rank response contains no independent matrix element beyond source means.",
            "A true field Gram additionally needs two-root `J_D*conj(J_S)` and `|J_S|^2`, and a non-tautological coupling needs an independent global observer.",
            "",
            "## Scientific card",
            "",
            *[f"{index}. {line}" for index, line in enumerate(report["scientific_card"], 1)],
            "",
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n65-prefix", type=Path, required=True)
    parser.add_argument("--n130-prefix", type=Path, required=True)
    parser.add_argument("--p50-prefix", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path)
    parser.add_argument("--dps", type=int, default=50)
    args = parser.parse_args()
    mp.mp.dps = args.dps
    report = build_report(args.n65_prefix, args.n130_prefix, args.p50_prefix)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.markdown:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(render_markdown(report), encoding="utf-8")


if __name__ == "__main__":
    main()

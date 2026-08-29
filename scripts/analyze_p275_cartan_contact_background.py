#!/usr/bin/env python3
"""Post-reveal Cartan-contact analysis of the P275 nine-geometry stream.

This script does not rescore the frozen P275 selector.  It verifies the exact
finite-path Ward identity behind its revealed scale-zero response and uses the
already-revealed full covariance for descriptive (discovery-only) fits.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Sequence

import mpmath as mp
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_p275_atop_field_identity import (  # noqa: E402
    MODULUS_ORDER,
    SIZE_ORDER,
    _expectation,
    _geometry_map,
    _load_prediction,
    _matching_root,
    _number,
    _pool_levels,
    _read_run,
    validate_runs,
)


def exact_transition_oracle() -> list[dict]:
    """Check q_- D=(S-D)/2 and q_+ D=(S+D)/2 transition by transition."""
    transitions = [
        (-1, 0, 1, -1, "01"),
        (0, 1, 1, 1, "12"),
        # A direct 0->2 rank jump has no canonical birth line in the runner;
        # its line-character S and D are therefore both zero.
        (-1, 1, 0, 0, "02_line_null"),
    ]
    result = []
    for q_minus, q_plus, source_s, source_d, kind in transitions:
        lhs_minus = q_minus * source_d
        rhs_minus = (source_s - source_d) / 2
        lhs_plus = q_plus * source_d
        rhs_plus = (source_s + source_d) / 2
        result.append({
            "transition": kind,
            "q_minus": q_minus,
            "q_plus": q_plus,
            "S": source_s,
            "D": source_d,
            "q_minus_D": lhs_minus,
            "half_S_minus_D": rhs_minus,
            "q_plus_D": lhs_plus,
            "half_S_plus_D": rhs_plus,
            "pass": lhs_minus == rhs_minus and lhs_plus == rhs_plus,
        })
    return result


def discover_runs(source_dir: Path) -> list[dict]:
    specs = []
    for metadata_path in sorted((source_dir / "raw").glob("N*/*.metadata.json")):
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        csv_path = metadata_path.with_name(
            metadata_path.name.removesuffix(".metadata.json") + ".batches.csv"
        )
        specs.append(f"{metadata['N']}:{metadata['modulus']}:{csv_path}:{metadata_path}")
    return [_read_run(spec) for spec in specs]


def geometry_contact_estimate(run: dict, target: dict, omitted_batch: int | None = None) -> dict:
    levels = _pool_levels(run["rows"], omitted_batch)
    p = _matching_root(levels)
    mean_q = _expectation(levels, "sum_q", p)
    birth = _expectation(levels, "sum_birth_mass", p)
    transport = target["transport"]
    phase = complex(_number(transport["real"]), _number(transport["imag"]))
    gamma_lab = []
    contact_lab = []
    thermal_lab = []
    relative_lab = []
    ingredients = {}
    for component in ("Re", "Im"):
        mean_d = _expectation(levels, f"sum_{component}_J_D4", p)
        mean_s = _expectation(levels, f"sum_{component}_J_S4", p)
        mean_qd = _expectation(levels, f"sum_q_{component}_J_D4", p)
        gamma_component = (mean_qd - mean_q * mean_d) / birth
        thermal_component = 0.5 * mean_s / birth
        relative_component = (p - 0.5 - mean_q) * mean_d / birth
        contact_component = thermal_component + relative_component
        gamma_lab.append(gamma_component)
        contact_lab.append(contact_component)
        thermal_lab.append(thermal_component)
        relative_lab.append(relative_component)
        ingredients[component] = {
            "mean_J_D4": mean_d,
            "mean_J_S4": mean_s,
            "mean_qJ_D4": mean_qd,
        }
    gamma = phase * complex(*gamma_lab)
    contact = phase * complex(*contact_lab)
    thermal = phase * complex(*thermal_lab)
    relative = phase * complex(*relative_lab)
    residual = gamma - contact
    return {
        "p_matching": p,
        "mean_q_residual": mean_q,
        "birth_mass": birth,
        "gamma_canonical": [gamma.real, gamma.imag],
        "contact_canonical": [contact.real, contact.imag],
        "thermal_contact_canonical": [thermal.real, thermal.imag],
        "relative_contact_canonical": [relative.real, relative.imag],
        "ward_residual_canonical": [residual.real, residual.imag],
        "ingredients_lab": ingredients,
    }


def vectors(
    runs: Sequence[dict], prediction: dict, omitted: tuple[int, int] | None = None
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict]:
    geometry = _geometry_map(prediction)
    run_map = {(run["N"], run["modulus"]): run for run in runs}
    gamma, contact, residual, thermal, relative, details = [], [], [], [], [], {}
    for n in SIZE_ORDER:
        for modulus in MODULUS_ORDER:
            batch = omitted[1] if omitted is not None and omitted[0] == n else None
            item = geometry_contact_estimate(
                run_map[(n, modulus)], geometry[(n, modulus)], batch
            )
            gamma.extend(item["gamma_canonical"])
            contact.extend(item["contact_canonical"])
            residual.extend(item["ward_residual_canonical"])
            thermal.extend(item["thermal_contact_canonical"])
            relative.extend(item["relative_contact_canonical"])
            details[f"N{n}:{modulus}"] = item
    return (
        np.asarray(gamma),
        np.asarray(contact),
        np.asarray(residual),
        np.asarray(thermal),
        np.asarray(relative),
        details,
    )


def jackknife(runs: Sequence[dict], prediction: dict) -> dict:
    full_gamma, full_contact, full_residual, full_thermal, full_relative, details = vectors(
        runs, prediction
    )
    covariance = np.zeros((18, 18), dtype=float)
    part_covariance = np.zeros((36, 36), dtype=float)
    residual_deleted = []
    batches = int(prediction["phase1_microcanonical_matching_root"]["batches"])
    for size_index, n in enumerate(SIZE_ORDER):
        block = slice(6 * size_index, 6 * (size_index + 1))
        deleted_gamma = []
        deleted_parts = []
        for batch in range(batches):
            gamma, _, residual, thermal, relative, _ = vectors(
                runs, prediction, (n, batch)
            )
            deleted_gamma.append(gamma[block])
            deleted_parts.append(np.concatenate((thermal[block], relative[block])))
            residual_deleted.append(residual[block])
        deleted_gamma = np.asarray(deleted_gamma)
        centered = deleted_gamma - deleted_gamma.mean(axis=0)
        covariance[block, block] = (batches - 1) / batches * (centered.T @ centered)
        deleted_parts = np.asarray(deleted_parts)
        centered_parts = deleted_parts - deleted_parts.mean(axis=0)
        local_part_covariance = (batches - 1) / batches * (
            centered_parts.T @ centered_parts
        )
        thermal_block = slice(6 * size_index, 6 * (size_index + 1))
        relative_block = slice(18 + 6 * size_index, 18 + 6 * (size_index + 1))
        indices = np.r_[
            np.arange(thermal_block.start, thermal_block.stop),
            np.arange(relative_block.start, relative_block.stop),
        ]
        part_covariance[np.ix_(indices, indices)] = local_part_covariance
    return {
        "gamma": full_gamma,
        "contact": full_contact,
        "residual": full_residual,
        "thermal": full_thermal,
        "relative": full_relative,
        "details": details,
        "covariance": covariance,
        "part_covariance": part_covariance,
        "max_delete_one_ward_residual": float(
            np.max(np.abs(np.concatenate(residual_deleted)))
        ),
    }


def _survival(chi_square: float, dof: int) -> float | None:
    if dof <= 0:
        return None
    return float(mp.gammainc(dof / 2, chi_square / 2, mp.inf, regularized=True))


def discovery_design(model: str, prediction: dict) -> np.ndarray:
    q4_shape = {
        row["id"]: _number(row["E4hat_over_i"])
        for row in prediction["moduli"]
    }
    complex_rows = []
    for n in SIZE_ORDER:
        x = n ** (-13 / 8)
        for modulus_index, _ in enumerate(MODULUS_ORDER):
            constants = [float(modulus_index == j) for j in range(3)]
            if model == "constant_by_modulus":
                features = constants
            elif model == "constant_by_modulus_plus_Q4_shape_tail":
                # Frozen Q4 modulus fingerprint from P275; only its overall
                # complex amplitude is fitted in this post-reveal diagnostic.
                q4 = q4_shape[MODULUS_ORDER[modulus_index]]
                features = constants + [x * q4]
            elif model == "constant_by_modulus_plus_free_tail":
                features = constants + [x * value for value in constants]
            else:  # pragma: no cover
                raise ValueError(model)
            complex_rows.append(features)
    p = len(complex_rows[0])
    design = np.zeros((18, 2 * p), dtype=float)
    for row, features in enumerate(complex_rows):
        design[2 * row, :p] = features
        design[2 * row + 1, p:] = features
    return design


def gls_fit(
    observation: np.ndarray, covariance: np.ndarray, model: str, prediction: dict
) -> dict:
    design = discovery_design(model, prediction)
    weight = np.linalg.pinv(covariance, rcond=1e-12, hermitian=True)
    normal_inv = np.linalg.pinv(design.T @ weight @ design, rcond=1e-12, hermitian=True)
    beta = normal_inv @ design.T @ weight @ observation
    residual = observation - design @ beta
    chi_square = float(residual @ weight @ residual)
    dof = len(observation) - design.shape[1]
    return {
        "chi_square": chi_square,
        "dof": dof,
        "survival_p": _survival(chi_square, dof),
        "coefficients_real_then_imag": beta.tolist(),
        "coefficient_standard_errors": np.sqrt(np.maximum(0, np.diag(normal_inv))).tolist(),
        "coefficient_covariance": normal_inv.tolist(),
        "fitted_observation": (design @ beta).tolist(),
    }


def zero_score(observation: np.ndarray, covariance: np.ndarray) -> dict:
    weight = np.linalg.pinv(covariance, rcond=1e-12, hermitian=True)
    chi_square = float(observation @ weight @ observation)
    return {
        "chi_square": chi_square,
        "dof": len(observation),
        "survival_p": _survival(chi_square, len(observation)),
    }


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def source_inventory(source_dir: Path) -> list[dict]:
    paths = [source_dir / "score.json"]
    paths.extend(sorted((source_dir / "raw").glob("N*/*.batches.csv")))
    paths.extend(sorted((source_dir / "raw").glob("N*/*.metadata.json")))
    return [
        {"path_relative_to_source": str(path.relative_to(source_dir)), "sha256": sha256(path)}
        for path in paths
    ]


def build_report(source_dir: Path, prediction_path: Path) -> dict:
    prediction = _load_prediction(prediction_path)
    runs = discover_runs(source_dir)
    provenance = validate_runs(runs, prediction)
    estimates = jackknife(runs, prediction)
    source_score = json.loads((source_dir / "score.json").read_text(encoding="utf-8"))
    source_y = np.asarray(source_score["estimates"]["observation_Y"], dtype=float)
    expected_gamma = source_y.copy()
    for size_index, n in enumerate(SIZE_ORDER):
        expected_gamma[6 * size_index:6 * (size_index + 1)] /= n ** (13 / 8)
    if not np.allclose(expected_gamma, estimates["gamma"], rtol=2e-12, atol=2e-12):
        raise ValueError("recomputed Gamma does not reproduce the frozen source score")
    if not np.allclose(
        estimates["gamma"], estimates["thermal"] + estimates["relative"], atol=5e-15
    ):
        raise ValueError("contact components do not reconstruct Gamma")
    sum_parts = np.concatenate((np.eye(18), np.eye(18)), axis=1)
    reconstructed_covariance = sum_parts @ estimates["part_covariance"] @ sum_parts.T
    if not np.allclose(reconstructed_covariance, estimates["covariance"], rtol=1e-9, atol=1e-18):
        raise ValueError("joint contact covariance does not reconstruct Gamma covariance")

    fits = {
        model: gls_fit(estimates["gamma"], estimates["covariance"], model, prediction)
        for model in (
            "constant_by_modulus",
            "constant_by_modulus_plus_Q4_shape_tail",
            "constant_by_modulus_plus_free_tail",
        )
    }
    thermal_covariance = estimates["part_covariance"][:18, :18]
    relative_covariance = estimates["part_covariance"][18:, 18:]
    thermal_fits = {
        model: gls_fit(estimates["thermal"], thermal_covariance, model, prediction)
        for model in (
            "constant_by_modulus",
            "constant_by_modulus_plus_Q4_shape_tail",
            "constant_by_modulus_plus_free_tail",
        )
    }
    base = fits["constant_by_modulus"]
    nested = {}
    for model in (
        "constant_by_modulus_plus_Q4_shape_tail",
        "constant_by_modulus_plus_free_tail",
    ):
        added = 2 if model.endswith("Q4_shape_tail") else 6
        delta = max(0.0, base["chi_square"] - fits[model]["chi_square"])
        nested[model] = {
            "delta_chi_square_from_constant": delta,
            "delta_dof": added,
            "survival_p": _survival(delta, added),
        }

    beta = np.asarray(base["coefficients_real_then_imag"])
    background = beta[:3] + 1j * beta[3:]
    q4_shape = np.asarray([
        _number(next(row for row in prediction["moduli"] if row["id"] == modulus)["E4hat_over_i"])
        for modulus in MODULUS_ORDER
    ], dtype=complex)
    projected_q4 = q4_shape - background * (
        np.vdot(background, q4_shape) / np.vdot(background, background)
    )
    normalization = np.vdot(projected_q4, q4_shape)
    weights = np.conjugate(projected_q4) / normalization
    if abs(np.dot(weights, background)) > 1e-10 or abs(np.dot(weights, q4_shape) - 1) > 1e-10:
        raise ValueError("background-annihilator construction failed")

    return {
        "schema": "matching-one/p275-cartan-contact-background/v1",
        "status": "post_reveal_discovery_analysis",
        "issues": [215, 258, 275],
        "source_boundary": {
            "frozen_score_is_unchanged": True,
            "source_archive_label": source_dir.name,
            "files": source_inventory(source_dir),
        },
        "provenance": provenance,
        "exact_transition_oracle": exact_transition_oracle(),
        "exact_identity": {
            "cartan_polynomial": "q^3=q for q in {-1,0,1}",
            "one_birth_identities": [
                "q_minus*D=(S-D)/2",
                "q_plus*D=(S+D)/2",
            ],
            "bernoulli_edge_balance": [
                "E_p[k*D_in]=p*T_D",
                "E_p[(N-k)*D_out]=(1-p)*T_D",
                "E_p[J_D]=T_D",
            ],
            "finite_N_connected_identity": "Cov_p(q,J_D)=E_p[J_S]/2+(p-1/2-E_p[q])E_p[J_D]",
            "provenance": "independently established in 8bef10b; rechecked here only as a numerical semantic control",
            "scope": "exact for the runner's global same-site line source; direct 0->2 jumps have null line character",
        },
        "estimates": {
            "coordinate_order": [
                f"N{n}:{modulus}:{part}"
                for n in SIZE_ORDER
                for modulus in MODULUS_ORDER
                for part in ("Re", "Im")
            ],
            "gamma": estimates["gamma"].tolist(),
            "ward_contact": estimates["contact"].tolist(),
            "ward_residual": estimates["residual"].tolist(),
            "thermal_contact": estimates["thermal"].tolist(),
            "relative_contact": estimates["relative"].tolist(),
            "max_abs_ward_residual": float(np.max(np.abs(estimates["residual"]))),
            "max_abs_delete_one_ward_residual": estimates["max_delete_one_ward_residual"],
            "covariance_gamma_18x18": estimates["covariance"].tolist(),
            "contact_part_order": "thermal_18_then_relative_18; each uses coordinate_order",
            "contact_part_covariance_36x36": estimates["part_covariance"].tolist(),
            "contact_point_and_covariance_reconstruction": "pass",
            "geometry_details": estimates["details"],
        },
        "post_reveal_discovery_fit": {
            "data": "unscaled Gamma with complete same-N, three-modulus jackknife covariance",
            "models": fits,
            "nested_gains": nested,
            "thermal_contact_models": thermal_fits,
            "relative_contact_zero_score": zero_score(
                estimates["relative"], relative_covariance
            ),
            "max_relative_to_total_complex_magnitude": float(max(
                abs(complex(*estimates["relative"][index:index + 2]))
                / abs(complex(*estimates["gamma"][index:index + 2]))
                for index in range(0, 18, 2)
            )),
            "interpretation_boundary": "descriptive only; these fits do not replace the frozen P275 selector",
        },
        "frozen_fixed_N_modulus_annihilator": {
            "training_data": "revealed N50/N130/N170 full-covariance block",
            "holdout": "new common N250 three-modulus block; no N250 target was used to form these coefficients",
            "background_C_hat_order_i_2i_5i_over_2": [
                [value.real, value.imag] for value in background
            ],
            "Q4_shape_F_order_i_2i_5i_over_2": [value.real for value in q4_shape],
            "complex_weights_w_order_i_2i_5i_over_2": [
                [value.real, value.imag] for value in weights
            ],
            "score": "Z250=sum_tau w_tau*Gamma(N250,tau)",
            "constraints": {
                "sum_w_C_hat": [
                    float(np.dot(weights, background).real),
                    float(np.dot(weights, background).imag),
                ],
                "sum_w_F_Q4": [
                    float(np.dot(weights, q4_shape).real),
                    float(np.dot(weights, q4_shape).imag),
                ],
            },
            "interpretation": "Z250 removes the trained scale-zero modulus profile at one fixed N and has unit response to the frozen Q4 modulus vector. Training uncertainty must be propagated; a nonzero value is not independent evidence unless N250 uses a new seed/counter.",
        },
        "scientific_layers": {
            "exact": "The observed Gamma is exactly the Cartan contact term at every finite N; contact subtraction leaves algebraic zero.",
            "continuum_shape": "A nonzero scale-zero limit by modulus is compatible with a conditional/projective-line polarization. Its limiting modulus function is inferred, not fixed by the Ward identity.",
            "N_minus_13_over_8": "No independent remainder is identifiable in this same-site statistic: the exact identity annihilates it before asymptotics.",
            "exploratory_conjecture": "Any genuine H4 field propagation must be measured with a separated or typed source/readout rather than a global same-site q-J_D covariance.",
        },
        "scientific_card": [
            "Question: why does revealed Gamma approach a nonzero modulus-dependent constant?",
            "Exact: three-state Cartan increment algebra plus Bernoulli edge balance makes Gamma a finite-N contact Ward term.",
            "Discovery: nine-geometry GLS separates modulus constants from optional N^-13/8 tails using the full covariance.",
            "Boundary: the same Gamma contains no separately identifiable H4 remainder after exact contact subtraction.",
            "Next test: freeze a background-annihilated held-out size and a separated typed-source observable; never reinterpret the frozen selector.",
        ],
    }


def render_markdown(report: dict) -> str:
    fits = report["post_reveal_discovery_fit"]["models"]
    details = report["estimates"]["geometry_details"]
    lines = [
        "# P275 Cartan contact background",
        "",
        "This is a post-reveal mechanism analysis. It does not modify or replace the frozen P275 score.",
        "",
        "## Exact regression constraint",
        "",
        "Commit `8bef10b` establishes the connected rank-gate identity",
        "",
        "`Cov_p(q,J_D) = E_p[J_S]/2 + (p-1/2-E_p[q])E_p[J_D]`.",
        "",
        "Here it is used as a numerical regression constraint, not reproved as a new result. "
        "The current Gamma is a contact term, not a contact term plus an identifiable small remainder.",
        "",
        "## Nine-geometry reproduction",
        "",
        "| geometry | Re Gamma | Re thermal | Re relative | max component residual |",
        "|---|---:|---:|---:|---:|",
    ]
    for n in SIZE_ORDER:
        for modulus in MODULUS_ORDER:
            row = details[f"N{n}:{modulus}"]
            gamma = row["gamma_canonical"]
            contact = row["contact_canonical"]
            thermal = row["thermal_contact_canonical"]
            relative = row["relative_contact_canonical"]
            residual = max(abs(v) for v in row["ward_residual_canonical"])
            lines.append(
                f"| N{n}/{modulus} | {gamma[0]:.10g} | {thermal[0]:.10g} | "
                f"{relative[0]:.10g} | {residual:.3g} |"
            )
    lines.extend([
        "",
        f"The relative-source term is statistically resolved but small: its largest complex-magnitude fraction is "
        f"`{report['post_reveal_discovery_fit']['max_relative_to_total_complex_magnitude']:.4%}`. "
        "The order-one background is therefore specifically thermal-contact dominated.",
        "",
        "## Post-reveal discovery GLS",
        "",
        "| model | chi2 | dof | survival p |",
        "|---|---:|---:|---:|",
    ])
    for name, row in fits.items():
        lines.append(f"| {name} | {row['chi_square']:.6g} | {row['dof']} | {row['survival_p']:.6g} |")
    lines.extend([
        "",
        "These fits describe the scale-zero contact shape and finite-size drift only. They are not a re-score of P275.",
        "Even the free per-modulus N^-13/8 tail fails, so the same statistic does not support a constant-plus-small-field decomposition.",
        "",
        "## Scientific layers",
        "",
    ])
    for key, value in report["scientific_layers"].items():
        lines.append(f"- **{key}:** {value}")
    lines.extend(["", "## Scientific card", ""])
    lines.extend(f"- {line}" for line in report["scientific_card"])
    lines.append("")
    return "\n".join(lines)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument(
        "--prediction",
        type=Path,
        default=ROOT / "predictions/p275_atop_q4_field_identity_20260829.yaml",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path)
    args = parser.parse_args(argv)
    report = build_report(args.source_dir.resolve(), args.prediction)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    if args.markdown:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(render_markdown(report), encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

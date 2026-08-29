#!/usr/bin/env python3
"""Freeze a normalization-free A_top field-identity selector.

The selected observable is gamma=Cov(A_top,J_D4)/B, where J_D4 is the
complex matching-odd rank-birth H4 source and B is the unmarked birth mass.
Three genuine moduli and three small common-N geometries turn ordinary Q4,
the inherited Q4 Jordan completion, and a generic allowed H4 completion into
different complex GLS subspaces.  No target data are read by this generator.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Sequence

import mpmath as mp

from derive_q4_jordan_log_slope_shape import e4hat


SIZES = (50, 130, 170)
ALPHA_N = Fraction(13, 8)
MODULI = (
    ("i", 1, 1),
    ("2i", 2, 1),
    ("5i_over_2", 5, 2),
)
REPRESENTATIVES = {
    50: {"i": (7, 1), "2i": (4, 3), "5i_over_2": (2, 1)},
    130: {"i": (11, 3), "2i": (8, 1), "5i_over_2": (3, 2)},
    170: {"i": (13, 1), "2i": (9, 2), "5i_over_2": (4, 1)},
}
INPUT_COMMITS = {
    "H4_selection": "fc14817bb8c0b2f6e7cbde41778e715dcb62bc64",
    "double_projector": "ddf41aa",
    "A_top_rank_birth_proxy": "078bd61",
    "eta_cocycle": "f5dcca6",
    "complex_character_archive": "cc1d43c",
}


def period_matrix(
    representation: tuple[int, int], numerator: int, denominator: int
) -> tuple[tuple[int, int], tuple[int, int]]:
    """Columns are denominator*z and numerator*i*z, giving tau=(p/q)i."""
    a, b = representation
    return (
        (denominator * a, -numerator * b),
        (denominator * b, numerator * a),
    )


def determinant(matrix: Sequence[Sequence[int]]) -> int:
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def smith_invariants(matrix: Sequence[Sequence[int]]) -> tuple[int, int]:
    n = abs(determinant(matrix))
    divisor = math.gcd(*(abs(value) for row in matrix for value in row))
    return divisor, n // divisor


def phase4(representation: tuple[int, int]) -> tuple[Fraction, Fraction]:
    """Return Re/Im[(a+ib)^4/(a^2+b^2)^2] exactly."""
    a, b = representation
    n = a * a + b * b
    return (
        Fraction(a**4 - 6 * a * a * b * b + b**4, n * n),
        Fraction(4 * a * b * (a * a - b * b), n * n),
    )


def fraction_text(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else str(value)


def eta_abs_imaginary(y: mp.mpf, dps: int = 90) -> mp.mpf:
    with mp.workdps(dps):
        q = mp.exp(-2 * mp.pi * y)
        log_value = -mp.pi * y / 12
        tolerance = mp.power(10, -(dps + 8))
        for n in range(1, 100_001):
            term = mp.log1p(-(q**n))
            log_value += term
            if abs(term) < tolerance:
                return +mp.exp(log_value)
    raise RuntimeError("eta product did not converge")


def exact_modulus_oracle(dps: int = 90) -> dict[str, object]:
    with mp.workdps(dps):
        values = {}
        reference_e4, _, _ = e4hat(mp.j, dps=dps)
        reference_eta = eta_abs_imaginary(mp.mpf(1), dps)
        for name, numerator, denominator in MODULI:
            y = mp.mpf(numerator) / denominator
            e4_value, terms, last = e4hat(mp.j * y, dps=dps)
            eta_value = eta_abs_imaginary(y, dps)
            e4_ratio = e4_value / reference_e4
            if abs(mp.im(e4_ratio)) > mp.power(10, -70):
                raise AssertionError("imaginary-axis E4hat ratio is not real")
            values[name] = {
                "tau": f"{fraction_text(Fraction(numerator, denominator))}i",
                "E4hat_over_i": mp.nstr(mp.re(e4_ratio), 70),
                "two_log_abs_eta_ratio_to_i": mp.nstr(
                    2 * mp.log(eta_value / reference_eta), 70
                ),
                "eta_abs": mp.nstr(eta_value, 70),
                "E4_terms": terms,
                "E4_last_term_abs": mp.nstr(last, 12),
            }
        ratio2 = mp.mpf(values["2i"]["E4hat_over_i"])
        if abs(ratio2 - mp.mpf(11) / 4) > mp.power(10, -60):
            raise AssertionError("E4hat(2i)/E4hat(i) lost the exact 11/4 ratio")
        eta2 = mp.mpf(values["2i"]["two_log_abs_eta_ratio_to_i"])
        eta25 = mp.mpf(values["5i_over_2"]["two_log_abs_eta_ratio_to_i"])
        return {
            "normalization": "E4hat(tau)=Im(tau)^2 E4(tau)",
            "values": values,
            "exact_E4hat_2i_over_i": "11/4",
            "eta_cocycle_ratio_2i_over_5i_over_2": mp.nstr(eta2 / eta25, 70),
            "eta_cocycle_residual": (
                "eta25*[R(2i)-R(i)]-eta2*[R(5i/2)-R(i)]=0"
            ),
        }


def geometry_payload() -> list[dict[str, object]]:
    rows = []
    for n in SIZES:
        for modulus, numerator, denominator in MODULI:
            representation = REPRESENTATIVES[n][modulus]
            matrix = period_matrix(representation, numerator, denominator)
            if determinant(matrix) != n:
                raise AssertionError(f"N={n}/{modulus}: determinant mismatch")
            smith = smith_invariants(matrix)
            if smith != (1, n):
                raise AssertionError(f"N={n}/{modulus}: selected quotient is not cyclic")
            real, imag = phase4(representation)
            if real * real + imag * imag != 1:
                raise AssertionError("spin-four frame phase is not on the unit circle")
            rows.append({
                "N": n,
                "modulus": modulus,
                "tau_imaginary": fraction_text(Fraction(numerator, denominator)),
                "representation_z": list(representation),
                "period_matrix": [list(row) for row in matrix],
                "smith_invariants": list(smith),
                "lab_to_canonical_H4_transport": {
                    "operation": "multiply gamma_lab by conjugate((z/abs(z))^4)",
                    "real": fraction_text(real),
                    "imag": fraction_text(-imag),
                },
            })
    return rows


def build_payload(dps: int = 90) -> dict[str, object]:
    modulus = exact_modulus_oracle(dps)
    f = {
        name: modulus["values"][name]["E4hat_over_i"]
        for name, _, _ in MODULI
    }
    eta = {
        name: modulus["values"][name]["two_log_abs_eta_ratio_to_i"]
        for name, _, _ in MODULI
    }
    return {
        "schema": "matching-one/p275-atop-q4-field-identity-design/v1",
        "status": "design_frozen_before_any_multi_modulus_gamma_target",
        "issues": [205, 215, 275],
        "input_commits": INPUT_COMMITS,
        "decision_boundary": (
            "H4 is treated as selected. No H8/H12 score is permitted in this design."
        ),
        "archive_audit": {
            "fc14817_quotient_prism": {
                "available": "A_top H4 at square modulus for N25/N50/N125 with full pair covariance",
                "result": "H4 chi2=1.087585/2 and A4=0.8039184543+/-0.0131464871",
                "post_reveal_log_slope_diagnostic": "-0.08258+/-0.08018; not a frozen field-identity score",
                "missing": "ell/H4 rank-birth source and same-sample A_top*J_D4 moments",
            },
            "Gaussian_Hecke_children": {
                "available": "scale and angular character at tau=i",
                "missing": "a genuine modulus axis and A_top*J_D4 moments",
            },
            "cc1d43c_complex_character": {
                "available": "same-parent charged local pivotal response with 4x4 covariance",
                "missing": "global q=A_top and q*J_D4; it is a different observable and cannot be relabelled",
            },
            "078bd61_rank_birth_proxy": {
                "available": "exact finite definition, parity controls, and tiny nonzero complex coupling",
                "missing": "multi-N multi-modulus common-field stream",
            },
            "conclusion": "no committed archive can reconstruct the field-identity score",
        },
        "observable": {
            "only_global_channel": "A_top=P2-P0=q=rank(H1 image)-1",
            "complex_source": "J_D4=sum_v chi4(ell_v)*(I12-I01)_v",
            "birth_mass": "B=E[sum_v(I01+I12)_v]=M_prime",
            "normalized_response": "gamma=Cov(A_top,J_D4)/B",
            "canonical_complex_response": "Gamma=conjugate((z/abs(z))^4)*gamma_lab",
            "scaled_response": "Y=N^(13/8)*Gamma",
            "excluded": [
                "E_top", "primitive_rank1_character", "charged_local_O_Schi_response",
                "H8_or_H12_refit", "free_radial_exponent",
            ],
        },
        "modulus_oracle": modulus,
        "selected_geometries": geometry_payload(),
        "common_field_blocks": {
            "sizes": list(SIZES),
            "within_each_N": (
                "all i/2i/5i_over_2 geometries use identical seed, counter, priority field, "
                "uniform-root schedule, and batch boundaries"
            ),
            "across_N": "independent seeds",
            "all_selected_quotients": "cyclic Z/N; Smith invariants (1,N)",
        },
        "frozen_model_subspaces": {
            "Q4_epsilon_ordinary": {
                "equation": "Y(N,tau)=c*F(tau)",
                "F_E4hat_over_i": f,
                "nuisance": "one common complex c",
                "real_observations_parameters_df": [18, 2, 16],
            },
            "Q4_energy_Jordan": {
                "equation": "Y(N,tau)=F(tau)*[c0+c1*log(N)+c2*e(tau)]",
                "F_E4hat_over_i": f,
                "e_two_log_eta_ratio_to_i": eta,
                "nuisance": "three common complex coefficients c0,c1,c2",
                "real_observations_parameters_df": [18, 6, 12],
                "parameter_free_relations": [
                    "log-slope vector is proportional to E4hat(tau)",
                    modulus["eta_cocycle_residual"],
                ],
            },
            "generic_allowed_H4_pure": {
                "equation": "Y(N,tau)=c_tau",
                "nuisance": "one unrelated complex intercept per modulus",
                "real_observations_parameters_df": [18, 6, 12],
            },
            "generic_allowed_H4_affine_log": {
                "equation": "Y(N,tau)=a_tau+b_tau*log(N)",
                "nuisance": "one unrelated complex intercept and slope per modulus",
                "real_observations_parameters_df": [18, 12, 6],
            },
        },
        "frozen_score_order": [
            "Q4_epsilon_ordinary",
            "Q4_energy_Jordan",
            "generic_allowed_H4_pure",
            "generic_allowed_H4_affine_log",
            "zero_response",
        ],
        "acquisition": {
            "status": "not_launched_runner_extension_required",
            "samples_per_geometry": 20_000_000,
            "batches": 100,
            "sizes": list(SIZES),
            "moduli": [name for name, _, _ in MODULI],
            "seeds_by_N": {"50": 2026275050, "130": 2026275130, "170": 2026275170},
            "replica_counter": [9_500_000_000, 9_520_000_000],
            "random_root": "counter-derived uniform label, identical across moduli within N",
            "production_authorized": False,
        },
        "minimal_sufficient_statistics": {
            "per_aligned_batch_per_geometry": [
                "samples and replica interval",
                "sum q and sum q^2",
                "sum I01, I12, and direct I02",
                "sum Re/Im J_S4 and J_D4",
                "sum q*Re/Im J_S4 and q*Re/Im J_D4",
                "unmarked birth mass sum N*(I01+I12)",
                "optional landing-H4 S/D values and q-products as a correlated diagnostic",
                "priority-field digest and random-root digest",
            ],
            "covariance": (
                "recompute every covariance, gamma ratio, canonical phase transport, N^(13/8) "
                "scaling, and GLS residual inside synchronized delete-one batches; preserve the "
                "full 18x18 real covariance in the frozen observation order"
            ),
            "not_sufficient": [
                "marginal threshold histograms", "unmarked birth mass alone",
                "marked H4 means without q*J", "separate marginal z scores",
            ],
        },
        "decision_map": {
            "ordinary_Q4_fits_and_no_log_gain": "thermal Q4 epsilon selected",
            "ordinary_fails_Jordan_fits_and_c1_nonzero": "thermal Q4 Jordan completion selected",
            "Q4_models_fail_generic_H4_fits": "H4 is real but its field completion is not the thermal Q4 module",
            "all_nonzero_models_fail": "multi-field flow or J_D4 does not isolate the global H4 coupling",
        },
        "scientific_card": [
            "Question: which field completes the already-selected global H4 channel?",
            "Observable: only gamma=Cov(A_top,J_D4)/birth_mass, transported as one complex character.",
            "Selector: three moduli times three cyclic sizes, scored with one full covariance block structure.",
            "Q4 fingerprint: E4hat modulus vector; Jordan adds the gauge-free eta cocycle, never a free exponent.",
            "Stop rule: do not collect more harmonic votes; failure promotes another H4 completion, not H8/H12.",
        ],
    }


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dps", type=int, default=90)
    parser.add_argument(
        "--output", type=Path,
        default=root / "results/p275-atop-q4-field-identity/latest.json",
    )
    args = parser.parse_args()
    payload = build_payload(args.dps)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

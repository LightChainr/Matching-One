#!/usr/bin/env python3
"""Clock-anchored N100 shear deformation from the existing 200 paired batches."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.stats import chi2


FIELDS = ["A_top", "E_top", "C", "W"]
REST = [0, 1, 3]
PARAMETERS = ["r_C", "D_A/D_C", "D_E/D_C", "D_W/D_C",
              "R_A", "R_E", "R_W", "R_A/D_C", "R_E/D_C", "R_W/D_C"]
R_E4 = (47-60*np.sqrt(2))/(47+60*np.sqrt(2))
SOURCE_COMMIT = "7b30648be558df0652a7ff22143cc87ed399d042"


def contrasts(mean):
    y = np.asarray(mean).reshape(3, 4)
    return y[1]-y[0], y[2]-y[0]


def parameters(mean):
    d, u = contrasts(mean)
    r = u[2]/d[2]
    residual = u-r*d
    return np.r_[r, d[REST]/d[2], residual[REST], residual[REST]/d[2]]


def jackknife_parameters(mean, batches):
    count = len(batches)
    rows = np.array([parameters((count*mean-row)/(count-1)) for row in batches])
    center = rows.mean(axis=0)
    covariance = (count-1)/count*(rows-center).T@(rows-center)
    return parameters(mean), covariance


def wald(mean, covariance):
    value = float(mean@np.linalg.solve(covariance, mean))
    return {"chi_square": value, "df": len(mean),
            "nominal_p_value": float(chi2.sf(value, len(mean))),
            "status": "post-reveal correlated mechanism diagnostic, not independent validation"}


def pinned_input(path):
    raw = path.read_bytes()
    data = json.loads(raw)
    if data.get("schema") == "matching-one/p267-clock-quotient-input/v1":
        return data
    return {
        "schema": "matching-one/p267-clock-quotient-input/v1",
        "original_score_sha256": hashlib.sha256(raw).hexdigest(),
        "original_score_path": "results/etop-n100-three-modulus/score.json",
        "original_score_state": "committed production score; supplied working-tree bytes match this source commit exactly",
        "original_score_commit": SOURCE_COMMIT,
        "prediction_freeze_commit": data["prediction_freeze_commit"],
        "contract": data["contract"], "shape_order": data["shape_order"],
        "field_order": data["field_order"], "mean": data["mean"],
        "covariance": data["covariance"], "batch_vectors": data["batch_vectors"],
        "raw_sha256": data["raw_sha256"],
    }


def calculate(data):
    mean = np.array(data["mean"])
    covariance = np.array(data["covariance"])
    batches = np.array(data["batch_vectors"])
    if data["field_order"] != FIELDS or batches.shape != (200, 12):
        raise ValueError("this analysis is only the declared 200x12 N100 block")
    d, u = contrasts(mean)
    point, theta_cov = jackknife_parameters(mean, batches)
    r, residual = point[0], point[4:7]
    vr = theta_cov[4:7, 4:7]
    q = wald(residual, vr)

    # Fixed-E4 residuals, conditioned on the C row being exactly zero under
    # that hypothesis. This is a Schur decomposition of the frozen joint
    # discrepancy, not another independent test or a declaration C is E4.
    L = np.c_[(R_E4-1)*np.eye(4), -R_E4*np.eye(4), np.eye(4)]
    z = L@mean
    vz = L@covariance@L.T
    beta = vz[REST, 2]/vz[2, 2]
    cond = z[REST]-beta*z[2]
    cond_cov = vz[np.ix_(REST, REST)]-np.outer(vz[REST, 2], vz[2, REST])/vz[2, 2]

    # Whitening is descriptive in the declared physical units. The matched
    # readout is trained here and is frozen only for future independent data.
    eig, vectors = np.linalg.eigh(vr)
    for k in range(3):
        if vectors[np.argmax(np.abs(vectors[:, k])), k] < 0:
            vectors[:, k] *= -1
    whitened = vectors.T@residual/np.sqrt(eig)
    dual = np.linalg.solve(vr, residual)
    dual /= dual[2]
    axes = {}
    for j, name in enumerate(["A_only", "E_only", "W_only"]):
        kept = [k for k in range(3) if k != j]
        axes[name] = wald(residual[kept], vr[np.ix_(kept, kept)])

    # A clock profile plus a single shear defect is the unique saturated
    # three-shape representation after this gauge choice; it is not a rank
    # or field-count discovery.
    y = mean.reshape(3, 4)
    full_r = np.zeros(4)
    full_r[REST] = residual
    reconstruction = y[0][None, :]+np.outer([0, 1, r], d)+np.outer([0, 0, 1], full_r)
    assert np.max(np.abs(reconstruction-y)) < 1e-15

    return {
        "schema": "matching-one/p267-clock-shape-quotient/v1",
        "dependency_group": "p267-N100-seed20260831125401-counter267100000000-fresh2M",
        "original_score_sha256": data["original_score_sha256"],
        "original_score_commit": data["original_score_commit"],
        "status": "post-reveal localization and future source freeze; zero new samples",
        "shape_order": data["shape_order"], "field_order": FIELDS,
        "contrasts": {"D_Y4_minus_Y2": d.tolist(), "U_Ys_minus_Y2": u.tolist()},
        "parameter_order": PARAMETERS, "parameters": point.tolist(),
        "parameter_se": np.sqrt(np.diag(theta_cov)).tolist(),
        "parameter_covariance": theta_cov.tolist(),
        "clock_quotient": {
            "formula": "r_C=U_C/D_C; R_j=U_j-r_C D_j",
            "r_C": float(r), "r_C_se": float(np.sqrt(theta_cov[0, 0])),
            "residual_field_order": [FIELDS[i] for i in REST],
            "residual": residual.tolist(), "covariance": vr.tolist(),
            "z_by_coordinate": (residual/np.sqrt(np.diag(vr))).tolist(),
            "joint_zero": q,
            "even_area_zero_R_W": wald(residual[2:], vr[2:, 2:]),
            "correlation": (vr/np.sqrt(np.outer(np.diag(vr), np.diag(vr)))).tolist(),
        },
        "exact_functional_crosswalk": {
            "input_identities": ["integral P4[A(p)] dp=-2 C", "integral P4[E(p)] dp=-W"],
            "clock_quotient_odd_area": "integral R_A(p) dp=0 exactly by the r_C definition",
            "clock_quotient_even_area": "integral R_E(p) dp=-R_W",
            "even_area_estimate": -float(residual[2]),
            "even_area_se": float(np.sqrt(vr[2, 2])),
            "boundary": "The zero odd area is an algebraic gauge condition, not an additional empirical pass; the full p-dependent curve is scored separately by the coordinator.",
        },
        "E4_schur_decomposition": {
            "fixed_r": R_E4, "four_residuals": z.tolist(),
            "four_covariance": vz.tolist(),
            "clock_z": float(z[2]/np.sqrt(vz[2, 2])),
            "clock_chi_square_1df": float(z[2]**2/vz[2, 2]),
            "noise_regression_beta_C_to_AEW": beta.tolist(),
            "conditional_AEW_mean": cond.tolist(),
            "conditional_AEW_covariance": cond_cov.tolist(),
            "conditional_AEW_zero": wald(cond, cond_cov),
            "boundary": "Conditions on the exact C-E4 mean-zero hypothesis. Different from fitting r_C; C compatibility is not its proof.",
        },
        "descriptive_deformation_axes": axes,
        "whitened_deformation": {
            "noise_eigenvalues": eig.tolist(), "noise_eigenvectors_columns": vectors.tolist(),
            "whitened_components": whitened.tolist(),
            "squared_norm": float(whitened@whitened),
            "future_matched_readout_AEW_weights_W_one": dual.tolist(),
            "source_readout_estimate": float(dual@residual),
            "source_readout_fixed_weight_se": float(np.sqrt(dual@vr@dual)),
            "boundary": "Weights and whitening are source-selected. Their source SNR is descriptive and does not have a new one-df significance interpretation.",
        },
        "minimal_two_profile_representation": {
            "formula": "Y_sj=Y_2j+D_j*[0,1,r_C]_s+R_j*[0,0,1]_s, R_C=0",
            "baseline": y[0].tolist(), "clock_profile": [0, 1, float(r)],
            "clock_loadings_D": d.tolist(), "shear_profile": [0, 0, 1],
            "shear_loadings_R": full_r.tolist(),
            "boundary": "Saturated after centering three shapes; gauge-fixed descriptive coordinates, not evidence of two fields or a unique continuum mechanism.",
        },
        "future_predictions_not_acquisition": {
            "scope": "A future independent block with the same three moduli, observable definitions, O map, and declared lineage/quotient change. No homogeneity or N400 recommendation is assumed.",
            "direction_transfer": "R_A,100 R_E,next - R_E,100 R_A,next=0 and R_A,100 R_W,next - R_W,100 R_A,next=0; profile source and target covariance, not source-fixed errors",
            "direction_training_vector_AEW": residual.tolist(),
            "direction_training_covariance": vr.tolist(),
            "stronger_clock_normalized_transfer": {
                "hypothesis": "R_next/D_C,next = R_100/D_C,100; an additional testable scale/geometry law, not inferred homogeneity",
                "mean": point[7:].tolist(), "covariance": theta_cov[7:, 7:].tolist(),
            },
            "pure_even_zero_area_alternative": "R_W,next=0; current R_W is 2.39 SE, so even-area zero is not rejected at the declared .01 boundary",
            "input_for_matched_readout": "Compute new R using its own measured r_C; apply the frozen three weights without reselecting the source direction.",
        },
    }


def report(result):
    q = result["clock_quotient"]
    s = result["E4_schur_decomposition"]
    lines = ["# N100 clock-anchored shear deformation", "",
             "Zero new samples; post-reveal reuse of the existing 200 paired batches.", "",
             f"Clock shape r_C = {q['r_C']:.9g} +/- {q['r_C_se']:.5g}.", "",
             "| residual | estimate | SE | marginal z |", "|---|---:|---:|---:|"]
    for i, name in enumerate(q["residual_field_order"]):
        lines.append(f"| {name} | {q['residual'][i]:.9g} | {np.sqrt(q['covariance'][i][i]):.5g} | {q['z_by_coordinate'][i]:.5g} |")
    z = q["joint_zero"]
    lines += ["", f"Joint clock-closure diagnostic: chi2={z['chi_square']:.6g}/3, nominal p={z['nominal_p_value']:.6g}.",
              "This localizes the coordinator's common-secant failure; it is not an independent test.", "",
              "The exact integral identities give integral R_A(p) dp=0 and integral R_E(p) dp=-R_W.",
              "The odd zero is imposed by the clock gauge. The even area estimate is nonzero in sign but only 2.39 SE; it does not independently reject zero at alpha .01.", "",
              f"Under fixed E4, C contributes {s['clock_chi_square_1df']:.6g} chi-square and conditional A/E/W contribute {s['conditional_AEW_zero']['chi_square']:.6g}/3. Conditioning on a true E4 clock is stronger than using an empirical clock.", "",
              "## Which deformation direction?", "",
              "Marginal A is negative, W positive and E unresolved. The complete covariance matters: E/W residual correlation is about -.763.",
              "The descriptive single-axis nuisance fits are:", ""]
    for name, value in result["descriptive_deformation_axes"].items():
        lines.append(f"- {name}: chi2={value['chi_square']:.6g}/2, nominal p={value['nominal_p_value']:.6g}.")
    dual = result["whitened_deformation"]["future_matched_readout_AEW_weights_W_one"]
    lines += ["", "A W-only deformation remains compatible at .01; a unique A/E/W loading is not identified. These source-selected diagnostics are not corrected independent model elections.",
              f"The source-selected covariance-matched readout is Psi={dual[0]:.8g} R_A + {dual[1]:.8g} R_E + R_W. Freeze it for future data; its in-sample SNR is not a new significance result.", "",
              "## Minimal representation and next falsifiers", "",
              "After subtracting Y2, use clock profile [0,1,r_C] and shear profile [0,0,1]. Their loadings are D and R. Three shapes make this saturated: it is a transparent coordinate choice, not two-field evidence.",
              "A next independent same-semantics block can test the source R direction via two cross-products, retaining source uncertainty. A stronger optional law keeps R/D_C fixed; its three source values and full covariance are in score.json. No area exponent, same-lineage homogeneity or new production is assumed.", "",
              "## Scientific card", "",
              "- Changed space: a readable C clock does not close the fixed-p A/E/W morphology; the remainder has an exact zero-area odd interpretation and an even lifetime-area coordinate.",
              "- Not proved: independent fields, physical Jordan identity, universal two-profile closure, or a separately significant nonzero even integral.",
              "- Observer/source/geometry: signed normalized A/E/C/W rank/clock contrasts at p=.59274605079, N100, three fixed moduli and the same cyclic/noncyclic O map.",
              "- Dependency: one common-random 2M/shape, 200-batch block; all new readouts are correlated post-reveal reuse.",
              "- Next upgrade: an independent geometry/scale block that preserves the frozen deformation direction or rejects it; whole-p curves separately determine how the zero-area redistribution is carried.", ""]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--input-copy", type=Path)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    data = pinned_input(args.source)
    result = calculate(data)
    for path in (args.json, args.report, args.input_copy):
        if path:
            path.parent.mkdir(parents=True, exist_ok=True)
    if args.input_copy:
        args.input_copy.write_text(json.dumps(data, indent=2)+"\n")
    args.json.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    args.report.write_text(report(result))


if __name__ == "__main__":
    main()

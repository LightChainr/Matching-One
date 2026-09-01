#!/usr/bin/env python3
"""Project the archived P334 covariance onto additive directional responses."""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import platform
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
NAMES = ("branch_success_Q", "total_gap", "common_gate", "between_checkpoints", "within_checkpoint")


def cos4(a: int, b: int) -> Fraction:
    return Fraction(a**4 - 6*a*a*b*b + b**4, (a*a+b*b)**2)


def rational(x: Fraction) -> dict:
    return {"numerator": x.numerator, "denominator": x.denominator, "decimal": float(x)}


def summaries(vector: np.ndarray, covariance: np.ndarray) -> dict:
    se = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    return {name: {"estimate": float(vector[i]), "checkpoint_cluster_se": float(se[i]),
                   "descriptive_standardized_contrast": float(vector[i]/se[i]) if se[i] else None}
            for i, name in enumerate(NAMES)}


def project(parent: dict, design: dict) -> dict:
    source_order = parent["vector_order"]
    source = np.asarray(parent["vector"], dtype=float)
    covariance = np.asarray(parent["checkpoint_cluster_covariance"], dtype=float)
    if source.shape != (22,) or covariance.shape != (22, 22):
        raise ValueError("Expected the full parent 22-coordinate checkpoint covariance")
    x = np.empty(10)
    jacobian = np.zeros((10, 22))
    for oi, orientation in enumerate(("first", "second")):
        idx = {key: source_order.index(f"{orientation}:{key}") for key in
               ("delta_coop_rb", "closure_baseline", "branch_success_rb", "s1", "s2")}
        within, baseline, q, s1, s2 = (source[idx[k]] for k in idx)
        j = 5*oi
        x[j:j+5] = q, q-s2*s2, (1/s1-1)*s2*s2, baseline-s2*s2/s1, within
        jacobian[j, idx["branch_success_rb"]] = 1
        jacobian[j+1, idx["branch_success_rb"]] = 1
        jacobian[j+1, idx["s2"]] = -2*s2
        jacobian[j+2, idx["s1"]] = -s2*s2/s1**2
        jacobian[j+2, idx["s2"]] = 2*s2*(1/s1-1)
        jacobian[j+3, idx["closure_baseline"]] = 1
        jacobian[j+3, idx["s1"]] = s2*s2/s1**2
        jacobian[j+3, idx["s2"]] = -2*s2/s1
        jacobian[j+4, idx["delta_coop_rb"]] = 1
    raw_covariance = jacobian @ covariance @ jacobian.T
    contrast = np.hstack((-np.eye(5), np.eye(5)))
    delta = contrast @ x
    delta_covariance = contrast @ raw_covariance @ contrast.T
    first = cos4(*design["first_rep"])
    second = cos4(*design["second_rep"])
    delta_cos4 = second-first
    if not delta_cos4:
        raise ValueError("The archived pair has no cos4 contrast")
    normalized = delta/float(delta_cos4)
    normalized_covariance = delta_covariance/float(delta_cos4)**2
    addition = np.array([0., 1., -1., -1., -1.])
    parent_indices = [1, 2, 3, 4, 6, 7, 8, 9]
    original = parent["secondary_pooled_gap_decomposition"]
    audit = {
        "first_addition_residual": float(addition @ x[:5]),
        "second_addition_residual": float(addition @ x[5:]),
        "directional_addition_residual": float(addition @ delta),
        "normalized_addition_residual": float(addition @ normalized),
        "directional_addition_covariance_residual_max": float(np.max(np.abs(addition @ delta_covariance))),
        "parent_decomposition_value_difference_max": float(np.max(np.abs(x[parent_indices]-original["vector"]))),
        "parent_decomposition_covariance_difference_max": float(np.max(np.abs(
            raw_covariance[np.ix_(parent_indices, parent_indices)]-original["checkpoint_cluster_covariance"]))),
    }
    if max(abs(v) for v in audit.values()) > 1e-12:
        raise ValueError(f"Additive projection does not match the archived decomposition: {audit}")
    return {
        "physical_first": design["first_rep"], "physical_second": design["second_rep"],
        "exact_cos4_first": rational(first), "exact_cos4_second": rational(second),
        "exact_delta_cos4": rational(delta_cos4),
        "source_vector_order": source_order, "source_vector": source.tolist(),
        "source_checkpoint_cluster_covariance": covariance.tolist(),
        "raw_vector_order": [f"{o}:{name}" for o in ("first", "second") for name in NAMES],
        "raw_vector": x.tolist(), "raw_checkpoint_cluster_covariance": raw_covariance.tolist(),
        "jacobian_from_parent": jacobian.tolist(),
        "first": summaries(x[:5], raw_covariance[:5, :5]),
        "second": summaries(x[5:], raw_covariance[5:, 5:]),
        "contrast_vector_order": list(NAMES),
        "directional_vector": delta.tolist(), "directional_covariance": delta_covariance.tolist(),
        "directional": summaries(delta, delta_covariance),
        "cos4_normalized_vector": normalized.tolist(),
        "cos4_normalized_covariance": normalized_covariance.tolist(),
        "cos4_normalized": summaries(normalized, normalized_covariance),
        "directional_component_fractions": None,
        "audit": audit,
    }


def report(payload: dict) -> str:
    lines = [
        "# P334 fork directional allocation: nonclosure is not yet an H4 mechanism",
        "",
        "**The within-checkpoint degree-dispersion excess has a resolved positive mean, but its directional contrast is unresolved in this design.** This zero-new-sample result projects the parent's complete checkpoint covariance; it does not rerun the already completed cooperative decomposition.",
        "",
        "The exact real-checkpoint counterexamples remain decisive against the specified scalar state `(N, orientation, k0, H2, b2, age, ell)`. That is a different claim from identifying which coordinate carries an orientation response. No new significance threshold, exponent fit, Monte Carlo block, or scalar-state test is introduced here.",
        "",
        "## Definitions and dependence",
        "",
        "Let `d=N-k0`, `b1=d-H2`, `c_v` be the safe second-site count after safe first site v, and `b2=sum(c_v)/2`. Per checkpoint, `s1=b1/d`, `s2=2*b2/[d(d-1)]`, and the one-common-update/two-clone probability is `sum(c_v²)/[d(d-1)²]`. Capital `S1,S2` denote checkpoint means, `Q` the mean exact branching probability, and `B=E[s2²/s1]`.",
        "",
        "`total_gap = Q-S2² = common_gate + between_checkpoints + within_checkpoint`, where the three parts are `(1/S1-1)S2²`, `B-S2²/S1`, and `Q-B=E[s1 Var(c_v/(d-1) | safe,C)]`. `Q` itself is separately reported and must not be confused with `total_gap`.",
        "",
        "The source has 20k base permutations per size. A permutation, both orientations and all clone/exact-count rows form one checkpoint cluster. The input 22×22 covariance is propagated through the parent's differentiable aggregate definitions and then the linear `second-first` contrast. The two sizes remain separate; all output rows reuse their original dependency groups, not additional independent evidence.",
        "",
        "## Direction-resolved outputs",
        "",
        "All ± quantities below are parent checkpoint-cluster/delta-method SEs. Standardized contrasts are descriptive post-reveal summaries, not multiplicity-adjusted claims. Cos4 normalization is a geometric comparison only: with two directions it cannot distinguish H4 from other angular contributions.",
    ]
    for size, item in payload["sizes"].items():
        lines += ["", f"### {size}", "",
                  f"Physical pair: `{item['physical_first']}` → `{item['physical_second']}`. Exact Δcos4 = `{item['exact_delta_cos4']['numerator']}/{item['exact_delta_cos4']['denominator']}`.", "",
                  "| observable | first ± SE | second ± SE | second−first ± SE | contrast / Δcos4 ± SE |",
                  "|---|---:|---:|---:|---:|"]
        for name in NAMES:
            cells = []
            for mode in ("first", "second", "directional", "cos4_normalized"):
                entry = item[mode][name]
                cells.append(f"{entry['estimate']:.8g} ± {entry['checkpoint_cluster_se']:.5g}")
            lines.append(f"| {name} | " + " | ".join(cells) + " |")
        z = item["directional"]["within_checkpoint"]["descriptive_standardized_contrast"]
        lines += ["", f"The within-checkpoint contrast / SE is `{z:.3f}`. No component/total-directional-gap fractions are formed: signed, uncertain directional denominators would turn additive responses into unstable mechanism shares."]
    lines += [
        "", "## What changes next", "",
        "1. **Keep the exact nonclosure result; do not resample its already saved inputs.** The safe-insertion graph's degree second moment (equivalently its 2-star overlap count) distinguishes real checkpoints with the same scalar tuple. It is a concrete microscopic coordinate, not an unexplained generic memory label.",
        "2. **Do not call successor-H2 one-step closure a new mechanism experiment.** In the new archive, on `branch_common_safe=1`, `H2_after=(N-k0-1)-branch_q_after_safe_count`. Given this value, the independent one-site clones have exact success `q=1-H2_after/(N-k0-1)` and product expectation `q²`; this is calibration. Absorbed rows use q=0 and do not define a rank-one successor H2.",
        "3. **Use existing states to target a genuinely richer question.** Frozen seed/counter/k0/period-matrix metadata reconstructs the current microscopic configuration. A deterministic replay can expose the safe-insertion degree distribution and ask what boundary organization carries its 2-star variation; a specifically chosen third-clone fan-out reads the cubic degree moment. These are possible follow-ups on existing states, not requests for another generic production block. The present direction contrast supplies no identified H4 carrier and does not justify enlarging the same experiment by default.",
        "", "## Provenance, reproducibility and narrow check", "",
        f"Source branch `{payload['provenance']['source_branch']}` at `{payload['provenance']['source_commit']}`; production result `{payload['provenance']['production_result_commit']}`. Input score SHA256 `{payload['provenance']['source_score_sha256']}`. Parent exact witnesses: `notes/p334-real-checkpoint-scalar-nonclosure.md`. This is a retrospective derived analysis, not independent confirmation.",
        "",
        "`score.json` preserves the parent 22×22 covariance, Jacobians, both-direction 10×10 covariance, 5×5 raw and normalized contrast covariances, exact rational geometry, environment versions and the additive check. Only one additive/parent-projection check is run; no test suite and no Monte Carlo are invoked.",
        "", "```bash",
        "/Users/lc/python-envs/research-py311/bin/python scripts/analyze_p334_fork_directional_allocation.py",
        "```", "",
        "Assessment: share with the stated finite-design, post-reveal and angular-identification caveats. No full-state temporal memory, scaling law, or continuum field identity follows from this projection.", "",
    ]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=ROOT/"analysis/p334_fork_directional_allocation_manifest.json")
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    source_path = ROOT/manifest["source_score"]
    source_bytes = source_path.read_bytes()
    if hashlib.sha256(source_bytes).hexdigest() != manifest["source_score_sha256"]:
        raise ValueError("Archived parent score differs from the pinned input")
    parent = json.loads(source_bytes)
    freeze = json.loads((ROOT/manifest["source_freeze"]).read_text())
    payload = {
        "schema": "matching-one/p334-fork-directional-allocation/v1",
        "status": manifest["status"], "new_samples": 0,
        "provenance": {key: manifest[key] for key in ("source_branch", "source_commit", "production_result_commit", "source_score", "source_score_sha256", "dependency_groups")},
        "environment": {"python": sys.version, "numpy": np.__version__, "platform": platform.platform(), "machine": platform.machine()},
        "definitions": {"contrast": "second-minus-first", "total_gap": "Q-S2^2", "claim_boundary": manifest["claim_boundary"], "normalization": manifest["normalization"], "directional_fraction_policy": manifest["directional_fraction_policy"]},
        "sizes": {n: project(parent["sizes"][n], freeze["runs"][n]) for n in ("N325", "N425")},
    }
    for key, content in (("json", json.dumps(payload, indent=2, sort_keys=True)+"\n"), ("markdown", report(payload))):
        path = ROOT/manifest["outputs"][key]
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)
    print(json.dumps({n: {"within": s["directional"]["within_checkpoint"], "audit": s["audit"]} for n, s in payload["sizes"].items()}, indent=2))


if __name__ == "__main__":
    main()

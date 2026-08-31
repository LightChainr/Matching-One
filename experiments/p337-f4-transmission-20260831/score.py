#!/usr/bin/env python3
"""Score the frozen fresh F4 transmission block, with paired root-complete LOO.

Usage: python score.py --raw-dir raw --contract CONTRACT.json --output results
The output directory receives score.json and REPORT.md. No simulation is run.
Each nN.hist.csv[.gz] must have its completed nN.run.json and nN.metadata.json.
Only these fresh streams provide the ordinary baseline and its matching root.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
from fractions import Fraction
from functools import lru_cache
import gzip
import hashlib
import json
from pathlib import Path
import platform
import sys

import numpy as np
import scipy
from scipy.optimize import brentq
from scipy.special import gammaln
from scipy.stats import norm


NS = (65, 85, 130, 170)
ORIENTATIONS = ("first", "second")
MODES = ("ordinary", "forced")
BIRTHS = ("first", "second")
DESIGNS = {
    65: ((8, 1), (7, 4)), 85: ((9, 2), (7, 6)),
    130: ((11, 3), (9, 7)), 170: ((13, 1), (11, 7)),
}
HEADER = ("n", "batch", "orientation", "mode", "degree", "replicas", "birth", "k", "count")
BATCHES = 100
SAMPLES_PER_BATCH = 200000
ROOT_BRACKET = (0.55, 0.65)
FAMILY_ALPHA = 0.05
PRACTICAL_HALFWIDTH = 0.5


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def cos4(a, b):
    return Fraction(a**4 - 6*a*a*b*b + b**4, (a*a + b*b)**2)


def load_contract(path):
    contract = json.loads(path.read_text())
    require(contract["Ns"] == list(NS), "contract Ns differ from the four frozen primary sizes")
    require(contract["batches"] == BATCHES, "contract batch budget changed")
    require(contract["samples_per_batch"] == SAMPLES_PER_BATCH, "contract batch size changed")
    require(contract["root_bracket"] == list(ROOT_BRACKET), "contract root bracket changed")
    require(contract["family_alpha"] == FAMILY_ALPHA and contract["family_size"] == 4,
            "contract primary family changed")
    require(contract["practical_halfwidth"] == PRACTICAL_HALFWIDTH,
            "contract material-resolution band changed")
    require(contract["old_data_pooled"] is False, "old-data pooling is prohibited")
    require(contract["primary_labels"] == [f"N{n}.V_F4" for n in NS],
            "contract primary labels changed")
    require(contract["canonical_degree"] == {"ordinary": "N", "forced": "N-4"},
            "ordinary and forced canonical degrees must remain different")
    require(len(set(contract["seeds"][str(n)] for n in NS)) == len(NS),
            "the four sizes must have independent seed domains")
    for n in NS:
        require(contract["samples_per_N"][str(n)] == BATCHES*SAMPLES_PER_BATCH,
                f"N{n}: sample budget changed")
        require(contract["geometries"][str(n)] == [list(x) for x in DESIGNS[n]],
                f"N{n}: frozen first/second geometry order changed")
    return contract


def load_size(raw_dir, n, contract, contract_digest):
    """Read the completed dense histogram, without discarding a single batch."""
    prefix = raw_dir / f"n{n}"
    csv_path = prefix.with_suffix(".hist.csv")
    gzip_path = prefix.with_suffix(".hist.csv.gz")
    path = gzip_path if gzip_path.is_file() else csv_path
    receipt_path = prefix.with_suffix(".run.json")
    metadata_path = prefix.with_suffix(".metadata.json")
    receipt = json.loads(receipt_path.read_text())
    require(receipt.get("status") == "completed" and receipt.get("exit_code") == 0,
            f"N{n}: production is not completed with exit 0")
    require(receipt.get("N") == n and receipt.get("samples") == BATCHES*SAMPLES_PER_BATCH,
            f"N{n}: wrong receipt size or sample count")
    require(receipt.get("batch_begin") == 0 and receipt.get("batch_end") == BATCHES,
            f"N{n}: incomplete original batch domain")
    require(receipt.get("old_data_pooled") is False and receipt.get("contract") == contract,
            f"N{n}: receipt does not identify this exact fresh contract")
    freeze = receipt.get("freeze_commit")
    require(isinstance(freeze, str) and len(freeze) == 40 and
            all(c in "0123456789abcdef" for c in freeze.lower()),
            f"N{n}: missing full freeze commit")
    source_hashes = receipt["source_sha256"]
    contract_key = "experiments/p337-f4-transmission-20260831/CONTRACT.json"
    require(source_hashes.get(contract_key) == contract_digest,
            f"N{n}: producer contract bytes differ from --contract")
    metadata_digest = sha256(metadata_path)
    require(metadata_digest == receipt.get("producer_metadata_sha256"),
            f"N{n}: metadata hash mismatch")
    metadata = json.loads(metadata_path.read_text())
    require(metadata.get("status") == "completed" and metadata.get("N") == n and
            metadata.get("freeze_commit") == freeze and
            metadata.get("master_seed") == contract["seeds"][str(n)],
            f"N{n}: producer metadata does not match the frozen fresh domain")
    require(metadata.get("batch_begin") == 0 and metadata.get("batch_end_exclusive") == BATCHES and
            metadata.get("samples_per_batch") == SAMPLES_PER_BATCH and
            metadata.get("paired_permutations") == BATCHES*SAMPLES_PER_BATCH,
            f"N{n}: incomplete producer batch/replica domain")
    require(metadata.get("representations") == contract["geometries"][str(n)] and
            metadata.get("ordinary_degree") == n and metadata.get("forced_degree") == n-4,
            f"N{n}: producer geometry or canonical degree changed")
    raw_digest = sha256(path)
    digest_key = "gzip_sha256" if path.suffix == ".gz" else "csv_sha256"
    require(raw_digest == receipt.get(digest_key), f"N{n}: raw hash mismatch")

    degrees = {"ordinary": n, "forced": n-4}
    counts = {mode: np.zeros((BATCHES, 2, 2, degree+1), dtype=np.int64)
              for mode, degree in degrees.items()}
    seen = {mode: np.zeros_like(value, dtype=bool) for mode, value in counts.items()}
    opener = gzip.open if path.suffix == ".gz" else open
    row_count = 0
    with opener(path, "rt", newline="") as stream:
        reader = csv.DictReader(stream)
        require(tuple(reader.fieldnames or ()) == HEADER, f"N{n}: unexpected dense CSV schema")
        for line, row in enumerate(reader, start=2):
            try:
                row_n, batch, degree, replicas, k, count = (
                    int(row[key]) for key in ("n", "batch", "degree", "replicas", "k", "count"))
                mode = row["mode"]
                orientation = ORIENTATIONS.index(row["orientation"])
                birth = BIRTHS.index(row["birth"])
            except (ValueError, KeyError, TypeError) as error:
                raise ValueError(f"{path.name}:{line}: malformed histogram row") from error
            require(row_n == n and mode in degrees, f"{path.name}:{line}: wrong N or mode")
            require(degree == degrees[mode] and 0 <= k <= degree,
                    f"{path.name}:{line}: degree/k mismatch (forced degree is N-4)")
            require(0 <= batch < BATCHES and replicas == SAMPLES_PER_BATCH,
                    f"{path.name}:{line}: wrong batch domain/denominator")
            require(0 <= count <= replicas, f"{path.name}:{line}: invalid histogram count")
            index = (batch, orientation, birth, k)
            require(not seen[mode][index], f"{path.name}:{line}: duplicate dense histogram bin")
            counts[mode][index] = count
            seen[mode][index] = True
            row_count += 1
    for mode in MODES:
        require(bool(seen[mode].all()), f"N{n}/{mode}: missing bins/batches (including k=0)")
        require(bool((counts[mode].sum(axis=-1) == SAMPLES_PER_BATCH).all()),
                f"N{n}/{mode}: birth histogram does not sum to the original replica denominator")
        cdf = counts[mode].cumsum(axis=-1)
        require(bool((cdf[:, :, 0] >= cdf[:, :, 1]).all()),
                f"N{n}/{mode}: first-birth CDF is below second-birth CDF")
    provenance = {
        "file": path.name, "raw_sha256": raw_digest, "rows": row_count,
        "receipt_file": receipt_path.name, "receipt_sha256": sha256(receipt_path),
        "metadata_file": metadata_path.name, "metadata_sha256": metadata_digest,
        "metadata": metadata, "freeze_commit": freeze,
        "source_sha256": source_hashes, "seed": contract["seeds"][str(n)],
        "samples": BATCHES*SAMPLES_PER_BATCH,
    }
    return counts, provenance


@lru_cache(maxsize=None)
def binomial_basis(n):
    k = np.arange(n+1, dtype=float)
    log_choose = gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1)
    return k, log_choose


def bernstein_weights(n, p):
    """Stable probability weights; p is always inside the frozen bracket."""
    k, log_choose = binomial_basis(n)
    log_weight = log_choose + k*np.log(p) + (n-k)*np.log1p(-p)
    weight = np.exp(log_weight - log_weight.max())
    return weight/weight.sum()


def birth_jets(cdf, p, max_order):
    """Evaluate CDF Bernstein polynomials and their exact forward-difference jets."""
    degree = cdf.shape[-1]-1
    jets = []
    coefficient = cdf
    multiplier = 1
    for order in range(max_order+1):
        if order:
            coefficient = np.diff(coefficient, axis=-1)
            multiplier *= degree-order+1
        jets.append(multiplier*np.einsum("...k,k->...", coefficient,
                                        bernstein_weights(degree-order, p)))
    return np.stack(jets, axis=-1)


def topology_jets(birth):
    q = birth[:, 0]+birth[:, 1]
    e = -birth[:, 0]+birth[:, 1]
    q[:, 0] -= 1
    e[:, 0] += 1
    return q, e


def score_retained(n, retained_counts, denominator, delta):
    cdf = {mode: value.cumsum(axis=-1)/denominator for mode, value in retained_counts.items()}

    def pooled_q(p):
        f = birth_jets(cdf["ordinary"], p, 0)
        return float(np.mean(-1+f[:, 0, 0]+f[:, 1, 0]))

    left, right = ROOT_BRACKET
    q_left, q_right = pooled_q(left), pooled_q(right)
    require(q_left <= 0 <= q_right,
            f"N{n}: retained ordinary root is outside frozen bracket {ROOT_BRACKET}")
    p0 = brentq(pooled_q, left, right, xtol=5e-15,
                rtol=4*np.finfo(float).eps, maxiter=100)
    ordinary = birth_jets(cdf["ordinary"], p0, 2)
    forced = birth_jets(cdf["forced"], p0, 1)
    q, e = topology_jets(ordinary)
    q_forced, e_forced = topology_jets(forced)
    q_delta, e_delta = q_forced-q[:, :2], e_forced-e[:, :2]
    jq = n*p0**4*q_delta[:, 0]
    jqp = n*(4*p0**3*q_delta[:, 0]+p0**4*q_delta[:, 1])
    je = n*p0**4*e_delta[:, 0]
    jep = n*(4*p0**3*e_delta[:, 0]+p0**4*e_delta[:, 1])
    pooled = q.mean(axis=0)
    y = (e[0]-e[1])/delta
    d = float(pooled[1])
    require(d > 0 and np.isfinite(d), f"N{n}: retained ordinary root has no positive slope")
    amplitude = n**(13/8)/2
    rootdot = -float(jq.mean())/d
    ddot = float(jqp.mean())+float(pooled[2])*rootdot
    p4_je, p4_jep = float((je[0]-je[1])/delta), float((jep[0]-jep[1])/delta)
    fixed_p = amplitude/d*p4_jep
    root_motion = amplitude/d*float(y[2])*rootdot
    denominator_term = -amplitude/d*(float(y[1])/d)*ddot
    value = {
        "p0": p0, "U": amplitude*float(y[1])/d,
        "V_F4": fixed_p+root_motion+denominator_term, "rootdot": rootdot,
        "Q": float(pooled[0]), "Qp": d, "Qpp": float(pooled[2]),
        "Y": float(y[0]), "Yp": float(y[1]), "Ypp": float(y[2]),
        "mean_Jq": float(jq.mean()), "mean_Jqp": float(jqp.mean()),
        "P4_JE": p4_je, "P4_JEp": p4_jep, "Ddot": ddot,
        "V_fixed_p_source_Ejet": fixed_p, "V_root_motion_Ecurvature": root_motion,
        "V_denominator": denominator_term,
    }
    for oi, orientation in enumerate(ORIENTATIONS):
        for mode, jets in (("ordinary", ordinary), ("forced", forced)):
            for bi, birth in enumerate(("F1", "F2")):
                for order, suffix in enumerate(("value", "p", "pp")[:jets.shape[-1]]):
                    value[f"{orientation}.{mode}.{birth}.{suffix}"] = float(jets[oi, bi, order])
        for name, vector in (("Jq", jq), ("Jqp", jqp), ("JE", je), ("JEp", jep)):
            value[f"{orientation}.{name}"] = float(vector[oi])
    require(bool(np.isfinite(list(value.values())).all()), f"N{n}: nonfinite score")
    return value


def score_size(n, counts, provenance, zcrit):
    first, second = DESIGNS[n]
    delta_exact = cos4(*first)-cos4(*second)
    require(delta_exact != 0, f"N{n}: degenerate P4 normalization")
    delta = float(delta_exact)
    total = {mode: value.sum(axis=0) for mode, value in counts.items()}
    point = score_retained(n, total, BATCHES*SAMPLES_PER_BATCH, delta)
    labels = list(point)
    omissions = []
    for batch in range(BATCHES):
        retained = {mode: total[mode]-counts[mode][batch] for mode in MODES}
        result = score_retained(n, retained, (BATCHES-1)*SAMPLES_PER_BATCH, delta)
        omissions.append([result[label] for label in labels])
    loo = np.asarray(omissions)
    factor = np.sqrt((BATCHES-1)/BATCHES)*(loo-loo.mean(axis=0))
    covariance = factor.T@factor
    se = np.sqrt(np.maximum(np.diag(covariance), 0))
    estimate = np.asarray([point[label] for label in labels])
    index = labels.index("V_F4")
    v, v_se = float(estimate[index]), float(se[index])
    lower, upper = v-zcrit*v_se, v+zcrit*v_se
    excludes_zero = lower > 0 or upper < 0
    inside_band = lower >= -PRACTICAL_HALFWIDTH and upper <= PRACTICAL_HALFWIDTH
    material = lower > PRACTICAL_HALFWIDTH or upper < -PRACTICAL_HALFWIDTH
    resolution_status = ("material_positive" if lower > PRACTICAL_HALFWIDTH else
                         "material_negative" if upper < -PRACTICAL_HALFWIDTH else
                         "within_practical_band" if inside_band else "inconclusive")
    return {
        "N": n, "dependency_group": f"fresh-F4-N{n}-seed-{provenance['seed']}",
        "batch_ids": list(range(BATCHES)), "samples_per_batch": SAMPLES_PER_BATCH,
        "canonical_degrees": {"ordinary": n, "forced": n-4},
        "geometries_first_second": [list(first), list(second)],
        "delta_cos4_exact": str(delta_exact), "delta_cos4": delta,
        "labels": labels, "estimate": estimate.tolist(), "se": se.tolist(),
        "point": point, "se_by_label": dict(zip(labels, se.tolist())),
        "leave_one_batch_out": loo.tolist(), "factor": factor.tolist(),
        "factor_convention": "sqrt((B-1)/B)*(LOO-mean_LOO); covariance=factor.T@factor",
        "covariance": covariance.tolist(), "provenance": provenance,
        "primary": {"label": f"N{n}.V_F4", "estimate": v, "se": v_se,
                    "family95_interval": [lower, upper], "excludes_zero": excludes_zero,
                    "wholly_within_practical_band": inside_band,
                    "wholly_outside_practical_band": material,
                    "resolution_status": resolution_status},
    }


def report_text(score):
    lines = [
        "# Fresh F4 → global U: fixed-budget transmission", "",
        f"Zero-projection family: **{score['family_decision']['zero_projection']}**.",
        f"Finite-resolution action: **{score['family_decision']['action']}**.", "",
        "| N | V_F4 | paired jackknife SE | simultaneous family95 interval | resolution |",
        "|---:|---:|---:|---:|---|",
    ]
    for n in NS:
        primary = score["sizes"][str(n)]["primary"]
        low, high = primary["family95_interval"]
        lines.append(f"| {n} | {primary['estimate']:.9g} | {primary['se']:.9g} | "
                     f"[{low:.9g}, {high:.9g}] | {primary['resolution_status']} |")
    lines.extend([
        "", "## Definition and inference", "",
        "The sole source is the unscaled number F4 of fully occupied elementary faces, "
        "with source measure exp(t F4). Ordinary degree is N; forced-face degree is N−4, "
        "including births at forced k=0. J=N p^4(forced−ordinary).",
        "Every central estimate and each of the 100 paired batch deletions refits its "
        "own fresh ordinary pooled root inside [0.55,0.65], then differentiates the root, "
        "numerator and denominator of U. No old anchor or old sample enters.",
        f"Four primary V coordinates use Bonferroni normal critical value "
        f"{score['inference']['normal_critical_value']:.12g}. Independent N seed domains "
        "form separate covariance blocks; matching batch numbers across N are not pairing.",
        "The ±0.5 practical band uses the same bulk-source/global-U units. A resolved "
        "small nonzero effect can reject zero while all intervals remain inside that band. "
        "A non-excluding interval does not establish zero; an overlapping interval is inconclusive.",
        "", "## Raw/source diagnostics", "",
        "score.json retains the full central vector, all original delete-one vectors and "
        "their covariance. Its three V terms are fixed-p source E jet, moving-root E "
        "curvature and denominator response; their errors are correlated and are not added independently.",
        "", "## Scope", "",
        "Prediction scope is zero projection at the four declared finite sizes. This is "
        "not field identification, an asymptotic model, a source winner, or a validation of "
        "another source. No sign was predicted. The fixed block ends here regardless of outcome: "
        "no sample top-up, source substitution or size substitution is authorized by this score.",
        "", f"Freeze commit: `{score['freeze_commit']}`; old_data_pooled: false.", "",
    ])
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True, help="new directory for score.json and REPORT.md")
    parser.add_argument("--contract", type=Path, required=True, help="frozen CONTRACT.json")
    args = parser.parse_args()
    require(not (args.output/"score.json").exists() and not (args.output/"REPORT.md").exists(),
            "refusing to overwrite an existing final score/report")
    contract = load_contract(args.contract)
    contract_digest = sha256(args.contract)
    zcrit = float(norm.ppf(1-FAMILY_ALPHA/(2*len(NS))))
    sizes = {}
    reference_freeze = reference_sources = None
    for n in NS:
        counts, provenance = load_size(args.raw_dir, n, contract, contract_digest)
        if reference_freeze is None:
            reference_freeze = provenance["freeze_commit"]
            reference_sources = provenance["source_sha256"]
        require(provenance["freeze_commit"] == reference_freeze and
                provenance["source_sha256"] == reference_sources,
                "the four completed streams do not share one source/contract freeze")
        sizes[str(n)] = score_size(n, counts, provenance, zcrit)
    primary = [sizes[str(n)]["primary"] for n in NS]
    zero_rejected = any(row["excludes_zero"] for row in primary)
    all_inside = all(row["wholly_within_practical_band"] for row in primary)
    material_labels = [row["label"] for row in primary if row["wholly_outside_practical_band"]]
    action = ("STOP_F4_AS_MAJOR_SOURCE_AT_DECLARED_FINITE_RESOLUTION" if all_inside else
              "MATERIAL_FINITE_SIZE_RESPONSE_RESOLVED" if material_labels else
              "INCONCLUSIVE_STOP_FIXED_BLOCK_WITHOUT_TOP_UP")
    score = {
        "schema": "matching-one.p337-f4-fresh-transmission-score.v1",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "status": "completed_fixed_budget_score", "old_data_pooled": False,
        "prediction_scope": "zero-projection at the four declared finite sizes only",
        "nonclaims": ["field identity", "asymptotic model", "source winner", "predicted sign"],
        "freeze_commit": reference_freeze, "contract": contract,
        "contract_sha256": contract_digest, "scorer_sha256": sha256(Path(__file__)),
        "runtime": {"python": sys.version, "executable": sys.executable,
                    "numpy": np.__version__, "scipy": scipy.__version__, "platform": platform.platform()},
        "inference": {"estimator": "full fresh block; paired root-complete delete-one-batch jackknife",
                      "family_alpha": FAMILY_ALPHA, "family_size": len(NS),
                      "normal_critical_value": zcrit, "practical_halfwidth": PRACTICAL_HALFWIDTH,
                      "cross_N_dependence": "independent seed groups; primary covariance block diagonal"},
        "primary_labels": [row["label"] for row in primary],
        "primary_estimate": [row["estimate"] for row in primary],
        "primary_covariance": np.diag([row["se"]**2 for row in primary]).tolist(),
        "family_decision": {
            "zero_projection": "REJECTED" if zero_rejected else "NOT_EXCLUDED",
            "zero_rejecting_labels": [row["label"] for row in primary if row["excludes_zero"]],
            "all_intervals_inside_practical_band": all_inside, "material_labels": material_labels,
            "action": action, "additional_sampling": "none; fixed block is terminal",
        },
        "sizes": sizes,
    }
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output/"score.json").write_text(json.dumps(score, indent=2, allow_nan=False)+"\n")
    (args.output/"REPORT.md").write_text(report_text(score))
    print(json.dumps({"score": str(args.output/"score.json"), "family_decision": score["family_decision"]}))


if __name__ == "__main__":
    main()

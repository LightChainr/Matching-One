#!/usr/bin/env python3
"""One-pass secondary clock-line score AFTER immutable official P154 delivery.

Read only unmarked sum_q/sum_e from completed raw shards. Read the source
response and its saved omissions from the official result; never rescore it.
"""
import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
import re
import subprocess
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.stats import norm

FREEZE = "0820b8d203e2dc534bb883d6fdb4d6d1e0acb11f"
RULE = "83f3eba88d7f1290704f82610c28669dc5e12f3c"
MAP_COMMIT = "c2828e3430fe1ac7e02fbe0e5ddc0e6a24c99847"
MAP_SHA256 = "4c1f9a513324a3fd3167be48ffdd544c0e4e0905a7defcdd07eea1d1d24f96b5"
PINNED = {
    "CONTRACT.json": "0201a6c20db366c25c9a3fbf10a74b6efb198cc47d8127eaa47c0ee21a352f6e",
    "producer.cpp": "1bddf8aa3947b3a38a624bac669133680225af2765b648c9454659cfba3dd8ba",
    "run_production.py": "598336939239e601d9d0bd7367ecf9e982aca4788f399568687bf0bb28097c54",
    "score_production.py": "8fcd12bdf2954a342f2242c7b158a0779a74d09c05d79092629ec7fc2827e063",
    "archive_channel_split.py": "1e9bba472c0e8f0b9d0bcae13e66b4476eda806642d43ca7ab21df64288e8c3d",
    "vendor/primitive.cpp": "7893216c66802b28eb67eb27ac61976835291c4ad734f94a0d255a3e6d7e179a",
    "vendor/integer_period.cpp": "7df0d9362b31111eab8fc73ab4032eea37458fde0ff2d720bded8e7b530fa94a",
}
NS, BATCHES, MODES = (85, 340), 200, ("M10", "M11")
PRIMARY = tuple(f"N{n}.{c}.v" for n in NS for c in ("entry", "completion", "total"))
SECONDARY = tuple(f"N{n}.{m}" for n in NS for m in MODES)


def require(condition, message):
    if not condition:
        raise ValueError("UNSCORABLE: " + message)


def sha(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fixed_core(package):
    for name, expected in PINNED.items():
        require(sha(package / name) == expected, "official frozen file changed: " + name)
    spec = importlib.util.spec_from_file_location("p154_frozen_moments", package / "archive_channel_split.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)  # imports definitions only; guarded main is NOT run.
    return module.moments


def official_gate(package, result_path):
    contract = json.loads((package / "CONTRACT.json").read_text())
    result = json.loads(result_path.read_text())
    require(result["status"] == "FRESH_INDEPENDENT_FIXED_BUDGET_COMPLETED", "official status not final")
    require(result["freeze_commit"] == FREEZE, "wrong official freeze")
    require(result["contract"] == contract and result["contract_sha256"] == PINNED["CONTRACT.json"], "contract mismatch")
    require(result["code_sha256"] == PINNED["score_production.py"], "wrong official scorer")
    require(result["old_data_pooled"] is False, "official result pooled old data")
    require(result["sample_totals"] == contract["samples_per_N"], "wrong sample totals")
    require(result["primary_labels"] == list(PRIMARY), "official primary ordering changed")
    labels = result["labels"]
    require(len(set(labels)) == len(labels) and all(x in labels for x in PRIMARY), "bad official labels")
    require(set(result["covariance_groups"]) == {f"fresh_N{n}" for n in NS}, "unexpected official dependency groups")
    for n in NS:
        group = result["covariance_groups"][f"fresh_N{n}"]
        require(group["Ns"] == [n] and group["delete_one_batch_ids"] == list(range(BATCHES)), "omission IDs are not the original200 batches")
        values = np.asarray(group["delete_one_vectors"], dtype=float)
        require(values.shape == (BATCHES, len(labels)) and np.isfinite(values).all(), "invalid official omission vectors")
    return contract, result


def load_unmarked(raw_dir, contract, official):
    """Check all nine final receipts first; convert only q/e and row identities."""
    expected = {"n85-b000-200.csv.gz": (85, 0, 200)}
    expected.update({f"n340-b{b:03d}-{b+25:03d}.csv.gz": (340, b, b+25) for b in range(0, 200, 25)})
    paths = sorted(raw_dir.glob("n*-b*.csv.gz"))
    require({p.name for p in paths} == set(expected) and len(paths) == 9, "require exactly nine frozen completed shards")
    official_receipts = {Path(r["path"]).name: r for r in official["receipts"]}
    require(len(official["receipts"]) == len(official_receipts) == 9, "official receipt list is incomplete")
    receipts = []
    for path in paths:
        n, begin, end = expected[path.name]
        rp = path.with_name(path.name.replace(".csv.gz", ".run.json"))
        receipt = json.loads(rp.read_text())
        expected_samples = contract["samples_per_N"][str(n)] * (end - begin) // BATCHES
        require(receipt["status"] == "completed" and receipt["freeze_commit"] == FREEZE, "incomplete or differently frozen shard")
        require((receipt["N"], receipt["batch_begin"], receipt["batch_end"], receipt["samples"]) == (n, begin, end, expected_samples), "shard domain/sample mismatch")
        require(all(receipt["frozen_sha256"].get(k) == v for k, v in PINNED.items()), "receipt frozen source hashes mismatch")
        digest = sha(path)
        require(receipt["gzip_sha256"] == digest, "altered raw gzip")
        saved = official_receipts.get(rp.name)
        require(saved is not None and saved["sha256"] == sha(rp) and saved["samples"] == expected_samples and saved["N"] == n,
                "raw receipt is not the one used by the official result")
        receipts.append({"raw_path": str(path), "raw_sha256": digest, "receipt_path": str(rp),
                         "receipt_sha256": saved["sha256"], "receipt": receipt})
    arrays = {n: np.zeros((BATCHES, 2, n+1, 2), dtype=np.int64) for n in NS}
    seen = {n: np.zeros((BATCHES, 2, n+1), dtype=bool) for n in NS}
    for path in paths:
        n, begin, end = expected[path.name]
        per_batch = contract["samples_per_N"][str(n)] // BATCHES
        with gzip.open(path, "rt") as stream:
            for row in csv.DictReader(stream):
                b, g, k = int(row["batch"]), ("first", "second").index(row["orientation"]), int(row["k"])
                require(int(row["n"]) == n and begin <= b < end and 0 <= k <= n and int(row["samples"]) == per_batch, "unplanned unmarked row")
                require(not seen[n][b, g, k], "duplicate unmarked row")
                seen[n][b, g, k] = True
                arrays[n][b, g, k] = int(row["sum_q"]), int(row["sum_e"])
    require(all(s.all() for s in seen.values()), "missing unmarked batch/geometry/K rows")
    return arrays, receipts


def flat_gains(n, q, e, qq, ee, delta):
    """The value-only M10/M11 subset of c2828e34 evaluate; no old gains read."""
    d, t, a = q.mean(), qq.mean(), n ** (13 / 8) / 2
    require(np.isfinite(d) and d > 0, "nonpositive/nonfinite pooled slope")
    p4 = lambda v: (v[0] - v[1]) / delta
    f, fp = ((q-e)/2, (q+e)/2), ((qq-ee)/2, (qq+ee)/2)
    gains = []
    for alpha in ((np.array([1., -1.]), np.array([1., -1.])),
                  (np.array([-1., 1.]), np.array([1., -1.]))):
        j = tuple(al * ff / n for al, ff in zip(alpha, f))
        jp = tuple(al * ffp / n for al, ffp in zip(alpha, fp))
        rootdot = -(j[0]+j[1]).mean() / d
        ddot = (jp[0]+jp[1]).mean() + rootdot*t
        gains.append([sign*a/d * (p4(jjp) + rootdot*p4(ffp) - p4(ff)*ddot/d)
                      for sign, ff, ffp, jjp in zip((-1, 1), f, fp, jp)])
    gains = np.asarray(gains)
    require(np.isfinite(gains).all() and np.all(np.hypot(gains[:, 0], gains[:, 1]) > 0), "nonfinite/zero gain line")
    return gains


def baseline_gains(sums, samples, n, contract, moments):
    bracket = contract["root_bracket"]
    q_at = lambda p: sum(moments(sums[g], samples, p, n)[0][0] for g in range(2)) / 2
    lo, hi = q_at(bracket[0]), q_at(bracket[1])
    require(np.isfinite([lo, hi]).all() and lo*hi <= 0, "fixed root bracket failed")
    p0 = brentq(q_at, *bracket, xtol=5e-14, rtol=5e-14)
    jets = np.array([moments(sums[g], samples, p0, n) for g in range(2)])
    require(np.isfinite(jets).all(), "nonfinite baseline jets")
    delta = float(Fraction(contract["normalization"]["delta_cos4"][str(n)]))
    gains = flat_gains(n, jets[:, 1, 0], jets[:, 1, 1], jets[:, 2, 0], jets[:, 2, 1], delta)
    return gains, {"p0": float(p0), "direction_order": ["first", "second"],
                   "moment_axis": ["value", "p_derivative", "p_second_derivative"],
                   "field_axis": ["q", "E"], "jets": jets.tolist(), "delta_cos4": delta}


def perpendicular(gains, response):
    ce, cc = gains[:, 0], gains[:, 1]
    return (ce*response[1] - cc*response[0]) / np.hypot(ce, cc)


def analyze(raw, contract, official, moments):
    labels = official["labels"]
    primary_idx = [labels.index(k) for k in PRIMARY]
    primary = np.array([official["estimates"][k]["value"] for k in PRIMARY], dtype=float)
    totals = {n: raw[n].sum(axis=0) for n in NS}
    central, gains, baselines = np.empty(4), {}, {}
    for i, n in enumerate(NS):
        gains[n], baselines[n] = baseline_gains(totals[n], contract["samples_per_N"][str(n)], n, contract, moments)
        central[2*i:2*i+2] = perpendicular(gains[n], primary[3*i:3*i+2])
    joint_center = np.r_[primary, central]
    require(np.isfinite(joint_center).all(), "nonfinite central result")
    covariance = np.zeros((10, 10))
    groups = {}
    for i, n in enumerate(NS):
        saved = np.asarray(official["covariance_groups"][f"fresh_N{n}"]["delete_one_vectors"], dtype=float)
        joint_loo = np.tile(joint_center, (BATCHES, 1))
        joint_loo[:, :6] = saved[:, primary_idx]
        other = slice(3*(1-i), 3*(1-i)+3)
        require(np.array_equal(joint_loo[:, other], np.tile(primary[other], (BATCHES, 1))), "official independent-N omission alignment failed")
        gain_loo, root_loo = [], []
        sample_loo = contract["samples_per_N"][str(n)] * (BATCHES-1) // BATCHES
        for b in range(BATCHES):
            c, baseline = baseline_gains(totals[n]-raw[n][b], sample_loo, n, contract, moments)
            gain_loo.append(c.tolist()); root_loo.append(baseline["p0"])
            joint_loo[b, 6+2*i:8+2*i] = perpendicular(c, joint_loo[b, 3*i:3*i+2])
        require(np.isfinite(joint_loo).all(), "nonfinite omission result")
        factor = np.sqrt((BATCHES-1)/BATCHES) * (joint_loo-joint_loo.mean(axis=0))
        covariance += factor.T @ factor
        groups[f"fresh_N{n}"] = {"Ns": [n], "delete_one_batch_ids": list(range(BATCHES)),
                                  "joint_delete_one_vectors": joint_loo.tolist(), "joint_factor": factor.tolist(),
                                  "gain_order": ["M10", "M11"], "gain_channel_order": ["entry", "completion"],
                                  "gain_delete_one": gain_loo, "baseline_root_delete_one": root_loo}
    require(np.allclose(covariance[:6, :6], official["primary_covariance"], rtol=5e-11, atol=5e-14), "joint factor does not retain official primary covariance")
    require(np.isfinite(covariance).all(), "nonfinite covariance")
    critical = float(norm.ppf(1-.05/8))
    se = np.sqrt(np.maximum(0, np.diag(covariance)[6:]))
    estimates = {label: {"value": float(v), "se": float(s),
                          "simultaneous_interval": [float(v-critical*s), float(v+critical*s)]}
                 for label, v, s in zip(SECONDARY, central, se)}
    models = {}
    for mode in MODES:
        failures = [f"N{n}.{mode}" for n in NS if estimates[f"N{n}.{mode}"]["simultaneous_interval"][0] > 0
                    or estimates[f"N{n}.{mode}"]["simultaneous_interval"][1] < 0]
        models[mode] = {"status": "rejected_pure_flat_jet_restriction" if failures else "not_excluded",
                        "contradicting_coordinates": failures}
    return {"secondary_labels": list(SECONDARY), "estimates": estimates,
            "secondary_covariance": covariance[6:, 6:].tolist(), "joint_labels": list(PRIMARY+SECONDARY),
            "joint_center": joint_center.tolist(), "joint_covariance": covariance.tolist(), "covariance_groups": groups,
            "baseline_by_N": {str(n): {**baselines[n], "flat_gain_entry_completion": gains[n].tolist()} for n in NS},
            "decision": {"family_alpha": .05, "family_size": 4, "normal_critical": critical, "models": models,
                         "equivalence_test": False, "forced_winner": False},
            "official_primary_decision_unchanged": official["decision"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--official-package", type=Path, required=True)
    parser.add_argument("--official-result", type=Path, required=True)
    parser.add_argument("--raw-dir", type=Path, help="default: OFFICIAL_PACKAGE/production")
    parser.add_argument("--raw-source-commit", required=True)
    parser.add_argument("--official-result-commit", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists(), "refuse to overwrite output")
    for value in (args.raw_source_commit, args.official_result_commit):
        require(re.fullmatch(r"[0-9a-f]{40}", value) is not None, "full source/result commit IDs required")
    package, result_path = args.official_package.resolve(), args.official_result.resolve()
    moments = fixed_core(package)
    contract, official = official_gate(package, result_path)
    raw, receipts = load_unmarked((args.raw_dir or package/"production").resolve(), contract, official)
    result = analyze(raw, contract, official, moments)
    repo = Path(__file__).resolve().parents[1]
    implementation = subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip()
    result.update(status="COMPLETED_SECONDARY_CLOCK_LINE_SCORE_SAME_FRESH_BLOCK", rule_commit=RULE,
                  production_freeze=FREEZE, implementation_checkout_commit=implementation, implementation_sha256=sha(Path(__file__)),
                  gain_formula_commit=MAP_COMMIT, gain_formula_script_sha256=MAP_SHA256,
                  raw_source_commit_declared=args.raw_source_commit, official_result_commit_declared=args.official_result_commit,
                  official_result_path=str(result_path), official_result_sha256=sha(result_path),
                  official_contract=contract, official_frozen_sha256=PINNED, input_receipts=receipts,
                  raw_numeric_fields_read=["sum_q", "sum_e"], official_source_score_recomputed=False,
                  source_alpha_fitted=False, old_data_pooled=False, new_samples=0,
                  dependence="Same N-wise200-batch experiment; joint baseline/source omissions. Four secondary residuals are not independent evidence from six official primary coordinates.",
                  boundary="Bonferroni normal marginal intervals, not finite-sample coverage. Reject only pure locally flat jets, not all clocks; no equivalence, winner, mixture fit or primary revision.")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x") as stream:
        json.dump(result, stream, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps(result["decision"], indent=2))


if __name__ == "__main__":
    main()

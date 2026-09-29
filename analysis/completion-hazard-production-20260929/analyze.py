#!/usr/bin/env python3
"""Score the independent completion-hazard block; Python standard library only.

Read contract.json beside --data-dir, and batches from --data-dir/run.json.
Each compressed batch is read once for SHA256 and CSV parsing. Jackknife
replicates reuse batch sufficient statistics, with stratum weights recomputed.
"""

import argparse
import csv
from dataclasses import dataclass, field
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import sys


REQUIRED_COLUMNS = ("J1", "J2", "dx_b", "dy_b", "nu_b", "tau", "nu_tau")
PRIMARY = (
    "delta_survival",
    "delta_integrated_hazard",
    "closure_residual",
    "delta_nu_instantaneous_scaled",
    "delta_overlap_D",
    "delta_overlap_D_nu",
    "closure_early",
    "closure_late",
)
LABELS = {
    "delta_survival": "Pooled survival: early minus late",
    "delta_integrated_hazard": "Integrated-hazard prediction: -I early + I late",
    "closure_residual": "Survival difference minus integrated prediction",
    "delta_nu_instantaneous_scaled": "Instantaneous nu difference, window-scaled (not integrated)",
    "delta_overlap_D": "Exact-D overlap-weighted survival difference",
    "delta_overlap_D_nu": "Exact-(D,nu_b) overlap-weighted survival difference",
    "closure_early": "Early map residual: mean Z + mean I - 1",
    "closure_late": "Late map residual: mean Z + mean I - 1",
}


def integer(value, label):
    if isinstance(value, bool):
        raise ValueError(f"{label}: Boolean is not an integer parameter")
    if isinstance(value, int):
        return value
    if isinstance(value, str):
        try:
            return int(value)
        except ValueError:
            pass
    raise ValueError(f"{label}: expected an integer, got {value!r}")


def named_parameter(document, name, required=False):
    """Allow top-level or nested named parameters, never conflicting values."""
    matches = []

    def visit(value, path):
        if isinstance(value, dict):
            for key, item in value.items():
                location = f"{path}.{key}" if path else key
                if key == name and not isinstance(item, (dict, list)):
                    matches.append((location, integer(item, location)))
                if isinstance(item, (dict, list)):
                    visit(item, location)
        elif isinstance(value, list):
            for index, item in enumerate(value):
                if isinstance(item, (dict, list)):
                    visit(item, f"{path}[{index}]")

    visit(document, "")
    values = {value for _, value in matches}
    if len(values) > 1:
        raise ValueError(f"Conflicting {name} fields: {matches}")
    if not values:
        if required:
            raise ValueError(f"contract.json must contain named integer {name}")
        return None
    return values.pop()


@dataclass
class Moments:
    n: int = 0
    z: int = 0
    integrated: float = 0.0
    nu: int = 0

    def observe(self, z, integrated, nu):
        self.n += 1
        self.z += z
        self.integrated += integrated
        self.nu += nu

    def merge(self, other):
        self.n += other.n
        self.z += other.z
        self.integrated += other.integrated
        self.nu += other.nu

    def report(self):
        return {
            "count": self.n,
            "survivors": self.z,
            "sum_I": self.integrated,
            "sum_nu_b": self.nu,
            "mean_Z": self.z/self.n if self.n else None,
            "mean_I": self.integrated/self.n if self.n else None,
            "mean_nu_b": self.nu/self.n if self.n else None,
        }


def pair():
    return [Moments(), Moments()]


@dataclass
class BatchStats:
    samples: int = 0
    sign_canonicalizations: int = 0
    groups: list = field(default_factory=pair)
    directions: dict = field(default_factory=dict)
    direction_nu: dict = field(default_factory=dict)

    def observe(self, group, direction, z, integrated, nu):
        self.groups[group].observe(z, integrated, nu)
        for cells, key in ((self.directions, direction),
                           (self.direction_nu, (*direction, nu))):
            if key not in cells:
                cells[key] = pair()
            cells[key][group].observe(z, integrated, nu)

    def merge(self, other):
        self.samples += other.samples
        self.sign_canonicalizations += other.sign_canonicalizations
        for mine, theirs in zip(self.groups, other.groups):
            mine.merge(theirs)
        for mine, theirs in ((self.directions, other.directions),
                             (self.direction_nu, other.direction_nu)):
            for key, values in theirs.items():
                if key not in mine:
                    mine[key] = pair()
                for target, source in zip(mine[key], values):
                    target.merge(source)


def aggregate(blocks, omit=None):
    result = BatchStats()
    for index, block in enumerate(blocks):
        if index != omit:
            result.merge(block)
    return result


def read_batch(data_dir, entry, parameters):
    if not isinstance(entry, dict):
        raise ValueError("run.json batches must be objects")
    for name in ("file", "samples", "seed", "sha256"):
        if name not in entry:
            raise ValueError(f"Batch entry is missing {name}")
    path = (data_dir / entry["file"]).resolve()
    if not path.is_relative_to(data_dir):
        raise ValueError(f"Batch file is outside data-dir: {entry['file']}")
    expected_samples = integer(entry["samples"], "batch.samples")
    if expected_samples <= 0:
        raise ValueError("Every listed batch must have positive samples")
    seed = str(integer(entry["seed"], "batch.seed"))
    expected_hash = str(entry["sha256"]).strip().lower()
    if len(expected_hash) != 64 or any(ch not in "0123456789abcdef" for ch in expected_hash):
        raise ValueError(f"Invalid batch SHA256: {entry['file']}")

    # Exactly one filesystem read of the batch; both checks use these bytes.
    raw = path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    if digest != expected_hash:
        raise ValueError(f"SHA256 mismatch for {entry['file']}")
    payload = gzip.decompress(raw) if path.suffix == ".gz" else raw
    reader = csv.DictReader(io.StringIO(payload.decode("utf-8-sig"), newline=""))
    if reader.fieldnames is None or len(set(reader.fieldnames)) != len(reader.fieldnames):
        raise ValueError(f"Missing or duplicate CSV header: {entry['file']}")
    missing = set(REQUIRED_COLUMNS)-set(reader.fieldnames)
    if missing:
        raise ValueError(f"Missing CSV columns {sorted(missing)} in {entry['file']}")

    a, b, c, n = (parameters[key] for key in ("a", "b", "c", "N"))
    stats = BatchStats()
    for line_number, row in enumerate(reader, start=2):
        try:
            if None in row:
                raise ValueError("Extra CSV fields")
            j1, j2, dx, dy, nu, tau, nu_tau = (
                integer(row[key], key) for key in REQUIRED_COLUMNS)
            if not 1 <= j1 <= j2 <= n:
                raise ValueError("Require 1 <= J1 <= J2 <= N")
            if not b <= tau < c:
                raise ValueError("tau is outside [b,c-1]")
            at_b = j1 <= b < j2
            at_tau = j1 <= tau < j2
            if at_b:
                if not 0 <= nu <= n-b:
                    raise ValueError("nu_b is invalid on rank-one risk")
                if math.gcd(abs(dx), abs(dy)) != 1:
                    raise ValueError("Rank-one direction must be nonzero and primitive")
            elif nu != -1:
                raise ValueError("nu_b must be -1 outside rank-one risk")
            if at_tau:
                if not 0 <= nu_tau <= n-tau:
                    raise ValueError("nu_tau is invalid on rank-one risk")
            elif nu_tau != -1:
                raise ValueError("nu_tau must be -1 outside rank-one risk")
            if tau == b and nu_tau != nu:
                raise ValueError("nu_tau must equal nu_b when tau=b")
        except (ValueError, KeyError, TypeError) as error:
            raise ValueError(f"{entry['file']} CSV line {line_number}: {error}") from error
        stats.samples += 1
        if not at_b:
            continue
        if dx < 0 or (dx == 0 and dy < 0):
            dx, dy = -dx, -dy
            stats.sign_canonicalizations += 1
        group = 0 if j1 <= a else 1
        integrated = (c-b)*nu_tau/(n-tau) if j2 > tau else 0.0
        stats.observe(group, (dx, dy), int(j2 > c), integrated, nu)
    if stats.samples != expected_samples:
        raise ValueError(f"Sample count mismatch in {entry['file']}: "
                         f"expected {expected_samples}, found {stats.samples}")
    provenance = {
        "file": str(path.relative_to(data_dir)), "samples": stats.samples,
        "seed": seed, "sha256": digest, "sha256_verified": True,
    }
    return stats, provenance


def overlap(cells, totals, detailed=False):
    weights, early, late = [], [], []
    supported = [0, 0]
    common = early_only = late_only = 0
    rows = []
    for key, (e, l) in sorted(cells.items()):
        weight = e.n*l.n/(e.n+l.n) if e.n and l.n else 0.0
        if weight:
            common += 1
            supported[0] += e.n
            supported[1] += l.n
            weights.append(weight)
            early.append(weight*e.z/e.n)
            late.append(weight*l.z/l.n)
        elif e.n:
            early_only += 1
        elif l.n:
            late_only += 1
        if detailed:
            record = {
                "direction": list(key[:2]), "early_count": e.n, "late_count": l.n,
                "early_survivors": e.z, "late_survivors": l.z,
                "overlap_weight": weight,
            }
            if len(key) == 3:
                record["nu_b"] = key[2]
            rows.append(record)
    total_weight = math.fsum(weights)
    early_mean = math.fsum(early)/total_weight if total_weight else None
    late_mean = math.fsum(late)/total_weight if total_weight else None
    risk = sum(group.n for group in totals)
    return {
        "status": "scoreable" if total_weight else "not_scoreable",
        "reason": None if total_weight else "No cell has both early and late risk",
        "difference": early_mean-late_mean if total_weight else None,
        "weighted_survival_early": early_mean,
        "weighted_survival_late": late_mean,
        "overlap_weight_sum": total_weight,
        "common_cells": common,
        "early_only_cells": early_only,
        "late_only_cells": late_only,
        "supported_early_count": supported[0],
        "supported_late_count": supported[1],
        "drop_fraction_early": 1-supported[0]/totals[0].n if totals[0].n else None,
        "drop_fraction_late": 1-supported[1]/totals[1].n if totals[1].n else None,
        "drop_fraction_all_risk": 1-sum(supported)/risk if risk else None,
        **({"cells": rows} if detailed else {}),
    }


def score(stats, parameters, detailed=False):
    e, l = (group.report() for group in stats.groups)
    by_d = overlap(stats.directions, stats.groups, detailed)
    by_dn = overlap(stats.direction_nu, stats.groups, detailed)
    values = {name: None for name in PRIMARY}
    reasons = {name: None for name in PRIMARY}
    if e["count"]:
        values["closure_early"] = e["mean_Z"]+e["mean_I"]-1
    else:
        reasons["closure_early"] = "Early cohort is empty"
    if l["count"]:
        values["closure_late"] = l["mean_Z"]+l["mean_I"]-1
    else:
        reasons["closure_late"] = "Late cohort is empty"
    if e["count"] and l["count"]:
        values["delta_survival"] = e["mean_Z"]-l["mean_Z"]
        values["delta_integrated_hazard"] = -e["mean_I"]+l["mean_I"]
        values["closure_residual"] = values["delta_survival"]-values["delta_integrated_hazard"]
        values["delta_nu_instantaneous_scaled"] = (
            (parameters["c"]-parameters["b"])/(parameters["N"]-parameters["b"])
            * (e["mean_nu_b"]-l["mean_nu_b"]))
    else:
        for name in PRIMARY[:4]:
            reasons[name] = "Both initial cohorts must be nonempty"
    for name, table in (("delta_overlap_D", by_d), ("delta_overlap_D_nu", by_dn)):
        values[name] = table["difference"]
        reasons[name] = table["reason"]
    result = {
        "samples": stats.samples,
        "rank1_at_b_count": e["count"]+l["count"],
        "excluded_non_rank1_at_b": stats.samples-e["count"]-l["count"],
        "cohorts": {"early": e, "late": l},
        "scalars": {
            name: {"estimate": values[name],
                   "status": "scoreable" if values[name] is not None else "not_scoreable",
                   "reason": reasons[name]}
            for name in PRIMARY
        },
        "conditional_D": by_d,
        "conditional_D_nu": by_dn,
    }
    if detailed:
        result["per_direction"] = [
            {"direction": list(direction), "early": groups[0].report(),
             "late": groups[1].report()}
            for direction, groups in sorted(stats.directions.items())
        ]
    return result


def jackknife(blocks, parameters, point):
    b_count = len(blocks)
    replicates = [score(aggregate(blocks, omit=i), parameters) for i in range(b_count)]
    matrix = [[replicate["scalars"][name]["estimate"] for name in PRIMARY]
              for replicate in replicates]
    covariance = [[None for _ in PRIMARY] for _ in PRIMARY]
    means, eligible, reasons = {}, {}, {}
    global_problem = ("At least two batches are required" if b_count < 2 else
                      "Unequal batch sizes: equal-size delete-one covariance is not reported"
                      if len({batch.samples for batch in blocks}) != 1 else None)
    for j, name in enumerate(PRIMARY):
        missing = [i for i, row in enumerate(matrix) if row[j] is None]
        problem = global_problem
        if point["scalars"][name]["estimate"] is None:
            problem = "Full-data point estimate is not scoreable"
        elif missing:
            problem = f"Undefined delete-one estimate after omitting batch indices {missing}"
        eligible[name] = problem is None
        reasons[name] = problem
        means[name] = math.fsum(row[j] for row in matrix)/b_count if problem is None else None
    for j, left in enumerate(PRIMARY):
        for k, right in enumerate(PRIMARY):
            if eligible[left] and eligible[right]:
                covariance[j][k] = (b_count-1)/b_count * math.fsum(
                    (row[j]-means[left])*(row[k]-means[right]) for row in matrix)
    uncertainty = {}
    for j, name in enumerate(PRIMARY):
        uncertainty[name] = {
            "status": "scoreable" if eligible[name] else "not_scoreable",
            "reason": reasons[name],
            "standard_error": math.sqrt(max(0.0, covariance[j][j])) if eligible[name] else None,
        }
    return {
        "method": "aligned equal-size whole-batch delete-one; weights and support recomputed",
        "batch_count": b_count,
        "point_estimator": "pooled plug-in, not jackknife-bias-corrected",
        "scalar_order": list(PRIMARY),
        "leave_one_out_estimates": matrix,
        "leave_one_out_mean": means,
        "covariance": covariance,
        "uncertainty": uncertainty,
        "no_partial_batch_subsets_used_for_covariance": True,
        "exact_linear_relations_can_make_covariance_singular": True,
    }


def fmt(value, digits=7):
    return "not_scoreable" if value is None else f"{value:.{digits}g}"


def markdown(result):
    p, point, jk = result["parameters"], result["point"], result["jackknife"]
    lines = [
        "# Completion-hazard production: fixed-cohort readout", "",
        f"Square L={p['L']}; N={p['N']}; fixed insertion cutoffs "
        f"a={p['a']}, b={p['b']}, c={p['c']}.", "",
        f"Scored {point['samples']:,} filtrations in {len(result['batches'])} batches. "
        "This is one shared random block, not independent evidence for each column.", "",
        "## Primary estimates", "",
        "| Estimand | Estimate | Aligned batch SE | Status |",
        "|---|---:|---:|---|",
    ]
    for name in PRIMARY:
        value, uncertainty = point["scalars"][name], jk["uncertainty"][name]
        status = value["status"] if value["estimate"] is None else uncertainty["status"]
        reason = value["reason"] or uncertainty["reason"]
        if reason:
            status += ": " + reason
        lines.append(f"| {LABELS[name]} | {fmt(value['estimate'])} | "
                     f"{fmt(uncertainty['standard_error'])} | {status} |")
    lines += ["", "The nu contrast has the hazard sign (early minus late) and is only "
              "instantaneous. It is not an integrated survival prediction.", "",
              "## Initial cohort means", "",
              "| Cohort | Count | mean Z | mean I | mean nu_b |",
              "|---|---:|---:|---:|---:|"]
    for name, group in point["cohorts"].items():
        lines.append(f"| {name} | {group['count']} | {fmt(group['mean_Z'])} | "
                     f"{fmt(group['mean_I'])} | {fmt(group['mean_nu_b'])} |")
    lines += ["", "## Exact common-support conditional estimands", "",
              "| Stratification | Common cells | Supported E / L | Dropped E | Dropped L | Dropped total |",
              "|---|---:|---:|---:|---:|---:|"]
    for label, key in (("D", "conditional_D"), ("(D, integer nu_b)", "conditional_D_nu")):
        table = point[key]
        lines.append(f"| {label} | {table['common_cells']} | "
                     f"{table['supported_early_count']} / {table['supported_late_count']} | "
                     f"{fmt(table['drop_fraction_early'])} | {fmt(table['drop_fraction_late'])} | "
                     f"{fmt(table['drop_fraction_all_risk'])} |")
    lines += ["", "Weights are e*l/(e+l) in each exact shared cell. These targets differ "
              "from the pooled difference and from each other. No mediated/explained "
              "fraction is defined. Missing support is reported, not imputed or binned.", "",
              "## Direction-resolved descriptive means", "",
              "| Primitive unoriented D | E n | L n | E Z | L Z | E I | L I | E nu_b | L nu_b |",
              "|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for row in point["per_direction"]:
        e, l = row["early"], row["late"]
        lines.append(f"| {tuple(row['direction'])} | {e['count']} | {l['count']} | "
                     + " | ".join(fmt(group[field]) for field in ("mean_Z", "mean_I", "mean_nu_b")
                                  for group in (e, l)) + " |")
    lines += ["", "## Interpretation and uncertainty", "",
              "- E[Z+I | initial cohort]=1 is an expectation identity, not a pathwise equality. "
              "Closure assesses the implemented geometry-to-exit map; it does not establish nu_b closure.",
              "- The (D,nu_b) finite-lag contrast is a distinct necessary prediction of a "
              "state-sufficiency model. A zero weighted contrast can hide opposite cell effects "
              "and does not prove recursive Markovness.",
              "- All delete-one replicates remove the same whole batch for every metric. "
              "All stratum weights and shared-cell support are recomputed. Full covariance and "
              "replicate vectors are in result.json; no pairwise subset covariance is used.",
              "- Per-direction rows are descriptive. No claim about causal mediation, the iid-label "
              "clock, or a large-size non-Markov limit is made.", "",
              "## Input provenance", "",
              "Each listed batch was read once, checked against its SHA256 and expected sample count, "
              "then retained only as sufficient statistics for resampling. Checksums do not prove "
              "tau independence or independent RNG design; those are generator-contract assumptions.", "",
              "| Batch | Samples | Seed | SHA256 |", "|---|---:|---|---|"]
    for row in result["batches"]:
        lines.append(f"| {row['file']} | {row['samples']} | {row['seed']} | {row['sha256']} |")
    lines += ["", "No old archive was pooled into these estimates.", ""]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    data_dir = args.data_dir.resolve()
    contract_path = data_dir.parent / "contract.json"
    contract = json.loads(contract_path.read_text())
    run = json.loads((data_dir / "run.json").read_text())
    if run.get("status") != "completed":
        raise ValueError("Only run.status='completed' may be scored; partial production blocks are not accepted")
    expected_batches = integer(contract["batches"], "contract.batches")
    expected_samples = integer(contract["samples_per_batch"], "contract.samples_per_batch")
    parameters = {name: named_parameter(contract, name, required=True) for name in ("a", "b", "c", "L")}
    parameters["N"] = named_parameter(contract, "N")
    if parameters["N"] is None:
        parameters["N"] = parameters["L"]**2
    a, b, c, n, size = (parameters[name] for name in ("a", "b", "c", "N", "L"))
    if size < 2 or n != size**2 or not 0 <= a < b < c <= n:
        raise ValueError("Require square N=L^2 and 0 <= a < b < c <= N")
    metadata_checks = {}
    for name, expected in parameters.items():
        observed = named_parameter(run, name)
        if observed is not None and observed != expected:
            raise ValueError(f"run.json {name}={observed} differs from contract {expected}")
        metadata_checks[name] = "matches" if observed is not None else "not_recorded_by_name"
    entries = run.get("batches")
    if not isinstance(entries, list) or not entries:
        raise ValueError("run.json must contain a nonempty batches list")
    if len(entries) != expected_batches:
        raise ValueError(f"Incomplete fixed block: contract requires {expected_batches} batches, "
                         f"run lists {len(entries)}")
    for index, entry in enumerate(entries):
        if integer(entry["samples"], f"batch[{index}].samples") != expected_samples:
            raise ValueError(f"Batch {index} does not have contract.samples_per_batch={expected_samples}; "
                             "partial batches are not accepted")
    blocks, provenance = [], []
    seen_files, seen_seeds, seen_hashes = set(), set(), set()
    for entry in entries:
        block, info = read_batch(data_dir, entry, parameters)
        for name, seen in (("file", seen_files), ("seed", seen_seeds), ("sha256", seen_hashes)):
            if info[name] in seen:
                raise ValueError(f"Repeated batch {name}; cannot count as an independent batch: {info[name]}")
            seen.add(info[name])
        blocks.append(block)
        provenance.append(info)
    total = aggregate(blocks)
    point = score(total, parameters, detailed=True)
    result = {
        "schema": "matching-one.completion-hazard-production.v1",
        "parameters": parameters,
        "fixed_block": {
            "run_status": run["status"], "required_batches": expected_batches,
            "required_samples_per_batch": expected_samples,
            "required_total_samples": expected_batches*expected_samples,
        },
        "contract_file": str(contract_path),
        "run_metadata_parameter_checks": metadata_checks,
        "clock": "finite uniform-permutation insertion count",
        "data_policy": "only listed new-block files; no old-archive pooling",
        "tau_design_assumption": "independent uniform integer in [b,c-1] for each full filtration",
        "direction_convention": "primitive, first nonzero coordinate positive; sign normalized only",
        "direction_sign_canonicalizations": total.sign_canonicalizations,
        "batches": provenance,
        "point": point,
        "jackknife": jackknife(blocks, parameters, point),
        "batch_cohort_sufficient_statistics": [
            {"batch_index": i, "early": block.groups[0].report(), "late": block.groups[1].report()}
            for i, block in enumerate(blocks)
        ],
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_json = args.output_dir / "result.json"
    output_md = args.output_dir / "RESULT.md"
    output_json.write_text(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True,
                                      allow_nan=False)+"\n")
    output_md.write_text(markdown(result))
    print(json.dumps({"result_json": str(output_json), "report": str(output_md),
                      "batches": len(blocks), "samples": total.samples}))


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, KeyError, TypeError, csv.Error, EOFError) as error:
        print(f"analysis failed: {error}", file=sys.stderr)
        sys.exit(2)

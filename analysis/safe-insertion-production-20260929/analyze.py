#!/usr/bin/env python3
"""Score only the fixed prospective square-L512 safe-insertion block (stdlib).

No acquisition, historical pooling, cutoff selection, or fitted models. Each
prefix contributes its four-probe mean; all uncertainty uses aligned batches.
"""

import argparse
import csv
from dataclasses import dataclass, field, fields
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import re
import sys


ROOT = Path(__file__).resolve().parent
FIXED = dict(N=262144, a=154646, b=155385, probes=4,
             batches=14, samples_per_batch=10000)
COLUMNS = ("J1", "rank_b", "dx_b", "dy_b", "nu_b", "Y0", "Y1", "Y2", "Y3")
METRICS = ("deltaY", "deltaE", "deltaZ2", "secondary_deltaNu",
           "weighted_earlyY", "weighted_lateY")
GROUPS = ("early", "late")


def integer(value, label):
    if type(value) is int:
        return value
    if isinstance(value, str) and re.fullmatch(r"-?[0-9]+", value):
        return int(value)
    raise ValueError(f"{label}: expected an integer, got {value!r}")


def sha256(raw):
    return hashlib.sha256(raw).hexdigest()


def read_json(path):
    raw = path.read_bytes()
    doc = json.loads(raw)
    if not isinstance(doc, dict):
        raise ValueError(f"{path}: expected a JSON object")
    return doc, {"path": str(path), "bytes": len(raw), "sha256": sha256(raw)}


def check_contract(contract, label="contract.json"):
    for key, expected in FIXED.items():
        if type(contract.get(key)) is not int or contract[key] != expected:
            raise ValueError(f"{label}.{key} must be the fixed integer {expected}")
    for key, expected in (("lattice", "square"), ("L", 512)):
        if key in contract and contract[key] != expected:
            raise ValueError(f"{label}.{key} must be {expected!r}")


@dataclass
class Moments:
    n: int = 0
    sum_nu: int = 0
    sum_nu_squared: int = 0
    safe_prefixes: int = 0
    sum_probe_Y: int = 0
    sum_probe_Y_squared: int = 0
    sum_squared_prefix_probe_sum: int = 0
    nonzero_probes: int = 0
    nonzero_prefixes: int = 0

    def observe(self, nu, ys):
        self.n += 1
        self.sum_nu += nu
        self.sum_nu_squared += nu * nu
        if ys is not None:
            total = sum(ys)
            self.safe_prefixes += 1
            self.sum_probe_Y += total
            self.sum_probe_Y_squared += sum(y * y for y in ys)
            self.sum_squared_prefix_probe_sum += total * total
            self.nonzero_probes += sum(y > 0 for y in ys)
            self.nonzero_prefixes += int(total > 0)

    def merge(self, other):
        for item in fields(self):
            setattr(self, item.name, getattr(self, item.name) + getattr(other, item.name))

    def report(self, q):
        return {
            "prefixes": self.n, "safe_prefixes": self.safe_prefixes,
            "dead_prefixes": self.n - self.safe_prefixes,
            "sum_nu_b": self.sum_nu, "sum_nu_b_squared": self.sum_nu_squared,
            "mean_nu_b": self.sum_nu / self.n if self.n else None,
            "valid_probes": q * self.safe_prefixes,
            "nonzero_probes": self.nonzero_probes,
            "nonzero_prefixes": self.nonzero_prefixes,
            "sum_probe_Y": self.sum_probe_Y,
            "sum_probe_Y_squared": self.sum_probe_Y_squared,
            "sum_squared_prefix_probe_sum": self.sum_squared_prefix_probe_sum,
            "sum_Ybar": self.sum_probe_Y / q,
            "sum_Ybar_squared": self.sum_squared_prefix_probe_sum / (q * q),
            "mean_Ybar": (self.sum_probe_Y / (q * self.safe_prefixes)
                          if self.safe_prefixes else None),
        }


def pair():
    return [Moments(), Moments()]


@dataclass
class BatchStats:
    samples: int = 0
    rank_counts: list = field(default_factory=lambda: [0, 0, 0])
    sign_canonicalizations: int = 0
    groups: list = field(default_factory=pair)
    primary_cells: dict = field(default_factory=dict)
    secondary_cells: dict = field(default_factory=dict)

    def observe(self, values, parameters):
        j1, rank, dx, dy, nu, *ys = values
        a, b, m = parameters["a"], parameters["b"], parameters["N"] - parameters["b"]
        if rank not in (0, 1, 2):
            raise ValueError("rank_b must be 0, 1, or 2")
        if rank == 0 and j1 != 0:
            raise ValueError("rank zero requires J1=0")
        if rank in (1, 2) and not 1 <= j1 <= b:
            raise ValueError("rank one/two requires 1<=J1<=b")
        if rank != 1:
            if (dx, dy) != (0, 0) or nu != -1 or ys != [-1] * parameters["probes"]:
                raise ValueError("outside rank one require direction=(0,0), nu_b=Y0=...=Y3=-1")
        else:
            if math.gcd(abs(dx), abs(dy)) != 1:
                raise ValueError("rank-one direction must be nonzero and primitive")
            if not 0 <= nu <= m:
                raise ValueError("rank one requires 0<=nu_b<=m")
            if nu == m:
                if ys != [-1] * parameters["probes"]:
                    raise ValueError("nu_b=m requires all Y=-1; safe continuation is undefined")
            elif any(y < 0 or y > m - 1 - nu for y in ys):
                raise ValueError("safe probes require 0<=Y<=m-1-nu_b")
        self.samples += 1
        self.rank_counts[rank] += 1
        if rank != 1:
            return
        if dx < 0 or (dx == 0 and dy < 0):
            dx, dy = -dx, -dy
            self.sign_canonicalizations += 1
        group = int(j1 > a)
        probes = ys if nu < m else None
        self.groups[group].observe(nu, probes)
        self.secondary_cells.setdefault((dx, dy), pair())[group].observe(nu, probes)
        if probes is not None:
            self.primary_cells.setdefault((dx, dy, nu), pair())[group].observe(nu, probes)

    def merge(self, other):
        self.samples += other.samples
        self.rank_counts = [x + y for x, y in zip(self.rank_counts, other.rank_counts)]
        self.sign_canonicalizations += other.sign_canonicalizations
        for target, source in zip(self.groups, other.groups):
            target.merge(source)
        for name in ("primary_cells", "secondary_cells"):
            target = getattr(self, name)
            for key, source in getattr(other, name).items():
                for mine, theirs in zip(target.setdefault(key, pair()), source):
                    mine.merge(theirs)

    def totals(self, q):
        risk = Moments()
        for group in self.groups:
            risk.merge(group)
        return {
            "rows": self.samples,
            "rank_counts": dict(zip(("rank0", "rank1", "rank2"), self.rank_counts)),
            "outside_rank1": self.rank_counts[0] + self.rank_counts[2],
            "sign_canonicalizations": self.sign_canonicalizations,
            "rank1": risk.report(q),
            "cohorts": {key: group.report(q) for key, group in zip(GROUPS, self.groups)},
        }


def aggregate(blocks, omit=None):
    total = BatchStats()
    for index, block in enumerate(blocks):
        if index != omit:
            total.merge(block)
    return total


def batch_entries(run, data_dir, parameters):
    if run.get("status") != "completed":
        raise ValueError("run.status must be completed; partial production is not scoreable")
    entries = run.get("batches")
    if not isinstance(entries, list) or len(entries) != parameters["batches"]:
        raise ValueError("run must contain all 14 contracted batches")
    seen_ids, seen_paths, seen_hashes = set(), set(), set()
    normalized = []
    for entry in entries:
        if not isinstance(entry, dict) or not {"file", "samples", "batch", "sha256"} <= entry.keys():
            raise ValueError("each batch needs file, samples, batch, sha256")
        identity = integer(entry["batch"], "batch ID")
        if identity < 0 or identity in seen_ids:
            raise ValueError(f"invalid or duplicate batch ID {identity}")
        seen_ids.add(identity)
        if integer(entry["samples"], "batch.samples") != parameters["samples_per_batch"]:
            raise ValueError(f"batch {identity} must declare 10000 samples")
        if not isinstance(entry["file"], str) or not entry["file"]:
            raise ValueError(f"batch {identity}: file must be a nonempty relative path")
        relative = Path(entry["file"])
        path = (data_dir / relative).resolve()
        if relative.is_absolute() or data_dir not in path.parents or path.suffix != ".gz":
            raise ValueError(f"batch {identity}: expected a .gz file inside data-dir")
        digest = entry["sha256"]
        if not isinstance(digest, str) or not re.fullmatch(r"[0-9a-fA-F]{64}", digest):
            raise ValueError(f"batch {identity}: invalid SHA256")
        digest = digest.lower()
        if path in seen_paths or digest in seen_hashes:
            raise ValueError(f"batch {identity}: duplicate batch file or compressed content")
        seen_paths.add(path)
        seen_hashes.add(digest)
        normalized.append((identity, path, digest, entry))
    return sorted(normalized, key=lambda item: item[0])


def read_batch(item, parameters):
    identity, path, expected_hash, entry = item
    raw = path.read_bytes()
    if sha256(raw) != expected_hash:
        raise ValueError(f"batch {identity}: SHA256 mismatch for {path}")
    payload = gzip.decompress(raw).decode("utf-8-sig")
    reader = csv.reader(io.StringIO(payload, newline=""), strict=True)
    if next(reader, None) != list(COLUMNS):
        raise ValueError(f"batch {identity}: require exact CSV header {','.join(COLUMNS)}")
    stats = BatchStats()
    for row in reader:
        try:
            if len(row) != len(COLUMNS):
                raise ValueError(f"expected {len(COLUMNS)} fields, found {len(row)}")
            stats.observe([integer(v, key) for key, v in zip(COLUMNS, row)], parameters)
        except ValueError as error:
            raise ValueError(f"batch {identity}, CSV line {reader.line_num}: {error}") from error
    if stats.samples != parameters["samples_per_batch"]:
        raise ValueError(f"batch {identity}: row count {stats.samples}, expected 10000")
    manifest = {"batch": identity, "path": str(path), "bytes": len(raw),
                "sha256": expected_hash, "sha256_verified": True,
                "rows_verified": stats.samples, "declared": entry}
    return stats, manifest


def overlap(cells, groups, parameters, primary, detailed=False):
    """Every call reselects common support and recomputes overlap weights."""
    q, m = parameters["probes"], parameters["N"] - parameters["b"]
    eligible = [g.safe_prefixes if primary else g.n for g in groups]
    supported = [0, 0]
    cell_counts = {"common": 0, "early_only": 0, "late_only": 0}
    weights, weighted_early, weighted_late = [], [], []
    weighted_delta, weighted_edges, weighted_survival = [], [], []
    rows = []
    for key, (early, late) in sorted(cells.items()):
        common = bool(early.n and late.n)
        weight = early.n * late.n / (early.n + late.n) if common else 0.0
        mean_e = ((early.sum_probe_Y / q if primary else early.sum_nu) / early.n
                  if early.n else None)
        mean_l = ((late.sum_probe_Y / q if primary else late.sum_nu) / late.n
                  if late.n else None)
        delta = mean_e - mean_l if common else None
        if common:
            cell_counts["common"] += 1
            supported[0] += early.n
            supported[1] += late.n
            weights.append(weight)
            weighted_early.append(weight * mean_e)
            weighted_late.append(weight * mean_l)
            weighted_delta.append(weight * delta)
            if primary:
                s = m - key[2]
                weighted_edges.append(weight * (s / 2) * delta)
                weighted_survival.append(weight * (-s / (m * (m - 1))) * delta)
        else:
            cell_counts["early_only" if early.n else "late_only"] += 1
        if detailed:
            record = {"direction": list(key[:2]), "common_support": common,
                      "overlap_weight": weight, "early": early.report(q),
                      "late": late.report(q)}
            if primary:
                s = m - key[2]
                record.update(nu_b=key[2], safe_sites=s, deltaY=delta,
                              deltaE=s / 2 * delta if common else None,
                              deltaZ2=-s / (m * (m - 1)) * delta if common else None)
            else:
                record["deltaNu"] = delta
            rows.append(record)
    weight_sum = math.fsum(weights)
    available = weight_sum > 0
    coverage = {}
    for i, group in enumerate(GROUPS):
        coverage[group] = {
            "rank1_prefixes": groups[i].n, "eligible_prefixes": eligible[i],
            "supported_prefixes": supported[i],
            "unsupported_eligible_prefixes": eligible[i] - supported[i],
            "supported_fraction_of_eligible": supported[i] / eligible[i] if eligible[i] else None,
            "supported_fraction_of_rank1": supported[i] / groups[i].n if groups[i].n else None,
            "dead_prefixes": groups[i].n - groups[i].safe_prefixes,
        }
    coverage["combined"] = {
        "rank1_prefixes": sum(g.n for g in groups), "eligible_prefixes": sum(eligible),
        "supported_prefixes": sum(supported),
        "unsupported_eligible_prefixes": sum(eligible) - sum(supported),
        "supported_fraction_of_eligible": sum(supported) / sum(eligible) if sum(eligible) else None,
        "supported_fraction_of_rank1": (sum(supported) / sum(g.n for g in groups)
                                        if sum(g.n for g in groups) else None),
    }
    result = {
        "status": "scoreable" if available else "not_scoreable",
        "reason": None if available else "no common early/late support in eligible cells",
        "stratification": "exact (D,nu_b), nu_b<m" if primary else "exact D, all rank-one prefixes",
        "cell_counts": cell_counts, "overlap_weight_sum": weight_sum,
        "coverage": coverage,
        "weighted_early_mean": math.fsum(weighted_early) / weight_sum if available else None,
        "weighted_late_mean": math.fsum(weighted_late) / weight_sum if available else None,
        "delta": math.fsum(weighted_delta) / weight_sum if available else None,
    }
    if primary:
        result["deltaE"] = math.fsum(weighted_edges) / weight_sum if available else None
        result["deltaZ2"] = math.fsum(weighted_survival) / weight_sum if available else None
    if detailed:
        for row in rows:
            row["normalized_weight"] = row["overlap_weight"] / weight_sum if available else None
        result["cells"] = rows
    return result


def estimate(stats, parameters, detailed=False):
    primary = overlap(stats.primary_cells, stats.groups, parameters, True, detailed)
    secondary = overlap(stats.secondary_cells, stats.groups, parameters, False, detailed)
    vector = [primary["delta"], primary["deltaE"], primary["deltaZ2"],
              secondary["delta"], primary["weighted_early_mean"], primary["weighted_late_mean"]]
    return vector, primary, secondary


def score(blocks, batch_ids, parameters):
    total = aggregate(blocks)
    point, primary, secondary = estimate(total, parameters, detailed=True)
    stops = []
    if any(v is None for v in point):
        stops.append("full block: absent common support; see primary/secondary diagnostics")
    replicates = []
    for i, identity in enumerate(batch_ids):
        vector, p, s = estimate(aggregate(blocks, omit=i), parameters)
        if any(v is None for v in vector):
            stops.append(f"delete batch {identity}: absent common support")
        replicates.append({"deleted_batch": identity, "vector": vector,
                           "primary_support": p, "secondary_support": s})
    B, dimension = len(blocks), len(METRICS)
    if B != parameters["batches"] or len(batch_ids) != B or len(set(batch_ids)) != B:
        raise ValueError("scoring requires all distinct contracted batch IDs")
    if any(block.samples != parameters["samples_per_batch"] for block in blocks):
        raise ValueError("scoring requires equal contracted batch sizes")
    covariance = [[None] * dimension for _ in METRICS]
    jackknife_mean = [None] * dimension
    standard_errors = [None] * dimension
    if not stops:
        vectors = [row["vector"] for row in replicates]
        jackknife_mean = [math.fsum(v[k] for v in vectors) / B for k in range(dimension)]
        covariance = [[(B - 1) / B * math.fsum(
            (v[j] - jackknife_mean[j]) * (v[k] - jackknife_mean[k]) for v in vectors)
            for k in range(dimension)] for j in range(dimension)]
        standard_errors = [math.sqrt(covariance[k][k]) for k in range(dimension)]
    return {
        "schema": "matching-one.safe-insertion-score.v1",
        "status": "not_scoreable" if stops else "scored",
        "stop_reasons": stops, "parameters": dict(parameters),
        "m": parameters["N"] - parameters["b"], "metric_order": list(METRICS),
        "point_vector": point, "point_estimates": dict(zip(METRICS, point)),
        "standard_errors": dict(zip(METRICS, standard_errors)),
        "totals": total.totals(parameters["probes"]),
        "batch_totals": [{"batch": identity, **block.totals(parameters["probes"])}
                         for identity, block in zip(batch_ids, blocks)],
        "primary": primary, "secondary": secondary,
        "joint_uncertainty": {
            "method": "aligned whole-batch delete-one jackknife; all weights/support recomputed",
            "covariance_formula": "(B-1)/B * sum_j (theta_(-j)-mean)(theta_(-j)-mean)^T",
            "metric_order": list(METRICS), "batch_count": B,
            "status": "stopped_absent_support" if stops else "available",
            "delete_one_mean_vector": jackknife_mean,
            "covariance": covariance, "delete_one": replicates,
            "point_estimator": "pooled plug-in; no jackknife bias correction",
            "dependence": "four probes share one prefix; all six metrics share one new block",
            "singularity": "deltaZ2=-2*deltaE/[m*(m-1)] and deltaY=weighted_earlyY-weighted_lateY",
        },
        "interpretation": [
            "The primary tests one safe-step completion-creation mean at fixed exact (D,nu_b).",
            "Derived edge and two-step-survival contrasts use cell-specific s before aggregation.",
            "Secondary exact-D nu replication is independent of the old block, not of this primary.",
            "A null weighted moment can hide cell cancellation and does not establish recursive closure.",
            "Finite square L512 insertion clock only; no limiting, label-clock, or universality claim.",
            "No historical data pooling, subgroup search, fitted exponent, or acquisition is performed.",
        ],
    }


def source_manifest(run, contract_info):
    local = []
    for name in ("analyze.py", "ESTIMANDS.md", "engine.cpp", "run.py"):
        path = ROOT / name
        if path.is_file():
            raw = path.read_bytes()
            local.append({"path": str(path), "bytes": len(raw), "sha256": sha256(raw)})
    declared = run.get("source_sha256", {})
    if not isinstance(declared, dict):
        raise ValueError("run.source_sha256 must be an object if provided")
    if "contract.json" in declared and declared["contract.json"] != contract_info["sha256"]:
        raise ValueError("run source contract.json hash disagrees with input contract")
    return {
        "local_files_at_scoring": local,
        "acquisition_declared_source_commit": run.get("source_commit"),
        "acquisition_declared_source_sha256": declared,
        "note": "Local scoring hashes and acquisition-declared hashes are distinct provenance; "
                "hash/range checks do not prove topology, uniform safe sampling, or seed independence.",
    }


def markdown(result):
    def number(value):
        return "undefined" if value is None else f"{value:.10g}"

    lines = ["# Prospective square-L512 safe-insertion score", "",
             f"Status: **{result['status']}**. Fixed 14 batches x 10000 prefixes; q=4.", "",
             "All differences are early minus late. The point estimates are pooled plug-in estimates.", "",
             "| Metric | Estimate | Aligned batch SE |", "|---|---:|---:|"]
    for key in METRICS:
        lines.append(f"| {key} | {number(result['point_estimates'][key])} | "
                     f"{number(result['standard_errors'][key])} |")
    if result["stop_reasons"]:
        lines += ["", "Joint inference STOPPED; no deletions were omitted or replaced.", ""]
        lines += [f"- {reason}" for reason in result["stop_reasons"]]
    totals = result["totals"]
    risk = totals["rank1"]
    lines += ["", f"Rows: {totals['rows']}; rank counts: {totals['rank_counts']}.",
              f"Rank-one prefixes: {risk['prefixes']}; safe: {risk['safe_prefixes']}; "
              f"dead (nu=m): {risk['dead_prefixes']}.",
              f"Valid probes: {risk['valid_probes']}; nonzero probes: {risk['nonzero_probes']}; "
              f"prefixes with any nonzero probe: {risk['nonzero_prefixes']}.", "",
              "| Comparison / cohort | Eligible | Supported | Coverage of eligible |",
              "|---|---:|---:|---:|"]
    for name in ("primary", "secondary"):
        for group in GROUPS:
            coverage = result[name]["coverage"][group]
            lines.append(f"| {name} / {group} | {coverage['eligible_prefixes']} | "
                         f"{coverage['supported_prefixes']} | "
                         f"{number(coverage['supported_fraction_of_eligible'])} |")
    lines += ["", "Primary excludes dead prefixes because a safe-step mean is undefined there; "
              "secondary includes every rank-one prefix, including dead prefixes.", "",
              "Within each shared exact (D,nu) cell: deltaE=(m-nu)*deltaY/2 and "
              "deltaZ2=-(m-nu)*deltaY/[m(m-1)]. The same normalized overlap weights "
              "combine all primary quantities. The exact-D secondary recomputes its own weights.", "",
              "The JSON records exact integer sufficient statistics, all cells, support coverage, "
              "all 14 aligned deletion vectors, the full 6x6 covariance and source/input manifests. "
              "Probes are averaged within prefix; metrics are not independent and covariance is "
              "necessarily singular. No covariance inverse or independent-metric test is used.", ""]
    lines += [f"- {item}" for item in result["interpretation"]]
    return "\n".join(lines) + "\n"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=ROOT / "data")
    parser.add_argument("--contract", type=Path,
                        help="default: contract.json beside data-dir")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results")
    parser.add_argument("--output", "--output-json", dest="output", type=Path,
                        help="JSON path (default: output-dir/result.json)")
    parser.add_argument("--md", "--output-md", dest="md", type=Path,
                        help="Markdown path (default: output-dir/RESULT.md)")
    args = parser.parse_args(argv)
    data_dir = args.data_dir.resolve()
    output = (args.output or args.output_dir / "result.json").resolve()
    report = (args.md or args.output_dir / "RESULT.md").resolve()
    try:
        if output == report:
            raise ValueError("JSON and Markdown output paths must differ")
        for path in (output, report):
            if path.exists():
                raise ValueError(f"refusing to overwrite {path}; choose new output paths")
        contract_path = (args.contract or data_dir.parent / "contract.json").resolve()
        contract, contract_info = read_json(contract_path)
        check_contract(contract)
        run, run_info = read_json(data_dir / "run.json")
        if "contract" in run:
            if not isinstance(run["contract"], dict):
                raise ValueError("run.contract must be an object")
            check_contract(run["contract"], "run.contract")
        sources = source_manifest(run, contract_info)
        entries = batch_entries(run, data_dir, FIXED)
        blocks, manifests = [], []
        for item in entries:
            block, manifest = read_batch(item, FIXED)
            blocks.append(block)
            manifests.append(manifest)
        result = score(blocks, [item[0] for item in entries], FIXED)
        result["created_utc"] = datetime.now(timezone.utc).isoformat()
        result["input_manifest"] = {"contract": contract_info, "run": run_info,
                                    "contract_document": contract, "run_document": run,
                                    "batches": manifests}
        result["source_manifest"] = sources
        result["outputs"] = {"json": str(output), "markdown": str(report)}
        # Serialize before creating either output. Exclusive creation prevents overwrites.
        encoded = json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n"
        prose = markdown(result)
        for path in (output, report):
            path.parent.mkdir(parents=True, exist_ok=True)
        with output.open("x", encoding="utf-8") as handle:
            handle.write(encoded)
        with report.open("x", encoding="utf-8") as handle:
            handle.write(prose)
        print(json.dumps({"status": result["status"], "json": str(output), "md": str(report)}))
        return 2 if result["stop_reasons"] else 0
    except (OSError, EOFError, ValueError, csv.Error) as error:
        print(f"safe-insertion scoring error: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())

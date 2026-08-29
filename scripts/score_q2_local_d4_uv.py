#!/usr/bin/env python3
"""Score the preregistered q2 local-landing D4 UV annihilator.

The primary is a mean source, never Cov(A_top,source).  R=4 minus R=2
removes an R-independent contact term.  The two N=65 orientations train the
thermal-shell nuisance and one H4 amplitude; the two N=130 q2-child
orientations are held out with the fixed -2^(-13/8) transfer.  q-products are
retained solely to verify the exact rank-gate contact identities.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Any, Sequence

import mpmath as mp
import yaml


RADII = (2, 4)
ORIENTATIONS = ("first", "second")
CHANNELS = ("S", "D", "qS", "qD", "I02")
VALUE_COLUMNS = (
    "sum_q", "sum_q2", "sum_gate01", "sum_gate12",
    "sum_inactive_gate01", "sum_inactive_gate12", "sum_site_S", "sum_site_D",
    "sum_local_R2_S", "sum_local_R2_D", "sum_q_local_R2_S", "sum_q_local_R2_D", "sum_local_R2_I02",
    "sum_local_R4_S", "sum_local_R4_D", "sum_q_local_R4_S", "sum_q_local_R4_D", "sum_local_R4_I02",
)


@dataclass
class Row:
    n: int
    a: int
    b: int
    orientation: str
    batch: int
    samples: int
    k: int
    values: dict[str, mp.mpf]


@dataclass
class Run:
    prefix: Path
    metadata: dict[str, Any]
    groups: dict[tuple[str, int], list[Row]]
    audit: dict[str, int]


def cos4(a: int, b: int) -> mp.mpf:
    n = a*a+b*b
    return mp.mpf(a**4-6*a*a*b*b+b**4)/(n*n)


def load_run(prefix: Path) -> Run:
    metadata = json.loads(Path(str(prefix)+".metadata.json").read_text(encoding="utf-8"))
    groups: dict[tuple[str, int], list[Row]] = {}
    with Path(str(prefix)+".path.csv").open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = set(VALUE_COLUMNS)-set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"missing local multiradius fields: {sorted(missing)}")
        for raw in reader:
            row = Row(
                n=int(raw["n"]), a=int(raw["a"]), b=int(raw["b"]),
                orientation=raw["orientation"], batch=int(raw["batch"]),
                samples=int(raw["samples"]), k=int(raw["k"]),
                values={name: mp.mpf(raw[name]) for name in VALUE_COLUMNS},
            )
            groups.setdefault((row.orientation, row.batch), []).append(row)
    for key, rows in groups.items():
        rows.sort(key=lambda row: row.k)
        n = rows[0].n
        if len(rows) != n or [row.k for row in rows] != list(range(n)):
            raise ValueError(f"incomplete microcanonical grid for {key}")
        if len({row.samples for row in rows}) != 1:
            raise ValueError(f"inconsistent samples for {key}")
    audit = {
        name: 0 for name in (
            "endpoint_failures", "site_failures", "line_failures",
            "local_mark_failures", "index_mismatches",
        )
    }
    with Path(str(prefix)+".complement_audit.csv").open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            for name in audit:
                audit[name] += int(raw[name])
    return Run(prefix, metadata, groups, audit)


def validate_run(run: Run, size_contract: dict[str, Any], phase: dict[str, Any]) -> None:
    metadata = run.metadata
    expected_n = int(size_contract["N"])
    expected_batches = int(phase["batches"])
    checks = {
        "schema": metadata.get("schema") == "matching-one/threshold-rank-marked-multiradius/v2",
        "N": metadata.get("designs", [{}])[0].get("N") == expected_n,
        "samples": metadata.get("samples_per_pair") == int(phase["samples_per_size"]),
        "batches": metadata.get("batches") == expected_batches,
        "seed": metadata.get("seed") == int(size_contract["seed"]),
        "counter_first": metadata.get("replica_counter_first") == int(size_contract["replica_counter"][0]),
        "counter_last": metadata.get("replica_counter_last_exclusive") == int(size_contract["replica_counter"][1]),
        "runner": metadata.get("git_commit") == phase["runner_commit"],
        "binary": metadata.get("binary_sha256") == phase["binary_sha256_arm64"],
        "multiradius": "per-k S,D,qS,qD,I02 Horvitz sums" in metadata.get("local_multiradius_schema", ""),
        "primary": metadata.get("local_multiradius_primary") == "mean D4 shell; q-products are exact contact controls only",
        "audit": not any(run.audit.values()),
    }
    design = metadata.get("designs", [{}])[0]
    checks.update({
        "id": design.get("id") == size_contract["id"],
        "first_rep": design.get("first") == size_contract["first_rep"],
        "second_rep": design.get("second") == size_contract["second_rep"],
        "first_matrix": design.get("first_period_matrix") == size_contract["first_matrix"],
        "second_matrix": design.get("second_period_matrix") == size_contract["second_matrix"],
    })
    if not all(checks.values()):
        raise ValueError(f"run contract failed for N{expected_n}: {checks}")
    expected_groups = {(orientation, batch) for orientation in ORIENTATIONS for batch in range(expected_batches)}
    if set(run.groups) != expected_groups:
        raise ValueError(f"batch/orientation set failed for N{expected_n}")


def combine(groups: Sequence[list[Row]]) -> list[Row]:
    if not groups:
        raise ValueError("cannot combine zero batches")
    result = []
    for k in range(len(groups[0])):
        rows = [group[k] for group in groups]
        first = rows[0]
        result.append(Row(
            n=first.n, a=first.a, b=first.b, orientation=first.orientation,
            batch=-1, samples=sum(row.samples for row in rows), k=k,
            values={name: mp.fsum(row.values[name] for row in rows) for name in VALUE_COLUMNS},
        ))
    return result


def binomial_weights(n: int, p: mp.mpf) -> list[mp.mpf]:
    if p <= 0:
        return [mp.mpf(1)]+[mp.mpf(0)]*n
    if p >= 1:
        return [mp.mpf(0)]*n+[mp.mpf(1)]
    q = 1-p
    values = [q**n]
    for k in range(n):
        values.append(values[-1]*(n-k)*p/((k+1)*q))
    return values


def matching_curve(rows: Sequence[Row], p: mp.mpf) -> mp.mpf:
    q = [row.values["sum_q"]/row.samples for row in rows]+[mp.mpf(1)]
    return mp.fsum(w*v for w, v in zip(binomial_weights(len(rows), p), q))


def intrinsic_center(first: Sequence[Row], second: Sequence[Row]) -> mp.mpf:
    lo, hi = mp.mpf(0), mp.mpf(1)
    for _ in range(100):
        p = (lo+hi)/2
        value = (matching_curve(first, p)+matching_curve(second, p))/2
        if value < 0:
            lo = p
        else:
            hi = p
    return (lo+hi)/2


def canonical(rows: Sequence[Row], column: str, p: mp.mpf) -> mp.mpf:
    n = len(rows)
    weights = binomial_weights(n-1, p)
    return mp.fsum(
        weights[k]*mp.mpf(n)/(n-k)*rows[k].values[column]/rows[k].samples
        for k in range(n)
    )


def q_pre(rows: Sequence[Row], p: mp.mpf) -> mp.mpf:
    weights = binomial_weights(len(rows)-1, p)
    return mp.fsum(
        weights[k]*rows[k].values["sum_q"]/rows[k].samples
        for k in range(len(rows))
    )


def orientation_sources(rows: Sequence[Row], p: mp.mpf) -> dict[str, mp.mpf]:
    result: dict[str, mp.mpf] = {}
    qmean = matching_curve(rows, p)
    for radius in RADII:
        for channel, suffix in (("S", "S"), ("D", "D"), ("qS", "S"), ("qD", "D")):
            prefix = "sum_q_local" if channel.startswith("q") else "sum_local"
            result[f"R{radius}_{channel}"] = canonical(
                rows, f"{prefix}_R{radius}_{suffix}", p)
        result[f"R{radius}_I02"] = canonical(rows, f"sum_local_R{radius}_I02", p)
        coefficient = p-mp.mpf("0.5")-qmean
        connected_d = (
            result[f"R{radius}_qD"]+p*result[f"R{radius}_D"]
            - qmean*result[f"R{radius}_D"]
        )
        connected_s = (
            result[f"R{radius}_qS"]+p*result[f"R{radius}_S"]
            - qmean*result[f"R{radius}_S"]
        )
        result[f"R{radius}_contact_D"] = connected_d-(
            result[f"R{radius}_S"]/2+coefficient*result[f"R{radius}_D"]
            - result[f"R{radius}_I02"]
        )
        result[f"R{radius}_contact_S"] = connected_s-(
            result[f"R{radius}_D"]/2+coefficient*result[f"R{radius}_S"]
            - result[f"R{radius}_I02"]
        )
    result["shell_D"] = (result["R4_D"]-result["R2_D"])/mp.log(2)
    result["shell_S"] = (result["R4_S"]-result["R2_S"])/mp.log(2)
    result["mean_q"] = qmean
    return result


def size_estimate(run: Run, omitted_batch: int | None = None) -> dict[str, Any]:
    batches = sorted({batch for _, batch in run.groups})
    selected = [batch for batch in batches if batch != omitted_batch]
    first = combine([run.groups[("first", batch)] for batch in selected])
    second = combine([run.groups[("second", batch)] for batch in selected])
    p = intrinsic_center(first, second)
    orientations = {
        "first": orientation_sources(first, p),
        "second": orientation_sources(second, p),
    }
    c_first, c_second = cos4(first[0].a, first[0].b), cos4(second[0].a, second[0].b)
    delta = c_first-c_second
    if delta == 0:
        raise ValueError("zero paired-orientation H4 leverage")
    p4_shell_d = (orientations["first"]["shell_D"]-orientations["second"]["shell_D"])/delta
    p4_shell_s = (orientations["first"]["shell_S"]-orientations["second"]["shell_S"])/delta
    base = []
    coordinate_order = []
    for orientation in ORIENTATIONS:
        for radius in RADII:
            for channel in CHANNELS:
                base.append(orientations[orientation][f"R{radius}_{channel}"])
                coordinate_order.append(f"N{first[0].n}:{orientation}:R{radius}:{channel}")
    contact = max(
        abs(orientations[o][f"R{r}_contact_{d}"])
        for o in ORIENTATIONS for r in RADII for d in ("D", "S")
    )
    return {
        "N": first[0].n, "p": p, "orientations": orientations,
        "cos4": {"first": c_first, "second": c_second, "delta": delta},
        "P4_shell_D": p4_shell_d, "P4_shell_S": p4_shell_s,
        "base_vector": base, "coordinate_order": coordinate_order,
        "contact_max_residual": contact,
    }


def solve_training(parent: dict[str, Any], child: dict[str, Any]) -> dict[str, Any]:
    s1 = parent["orientations"]["first"]["shell_S"]
    s2 = parent["orientations"]["second"]["shell_S"]
    c1 = parent["cos4"]["first"]
    c2 = parent["cos4"]["second"]
    d1 = parent["orientations"]["first"]["shell_D"]
    d2 = parent["orientations"]["second"]["shell_D"]
    determinant = s1*c2-s2*c1
    if abs(determinant) < mp.mpf("1e-20"):
        raise ValueError("N65 thermal/H4 training matrix is singular")
    beta = (d1*c2-d2*c1)/determinant
    amplitude = (s1*d2-s2*d1)/determinant
    scale = mp.power(2, -mp.mpf(13)/8)
    residual = []
    predictions = []
    for orientation in ORIENTATIONS:
        thermal = child["orientations"][orientation]["shell_S"]
        x = -scale*child["cos4"][orientation]
        prediction = beta*thermal+amplitude*x
        predictions.append(prediction)
        residual.append(child["orientations"][orientation]["shell_D"]-prediction)
    return {
        "thermal_beta": beta, "H4_amplitude_N65": amplitude,
        "training_determinant": determinant, "target_scale": -scale,
        "child_prediction": predictions, "child_residual": residual,
    }


def jackknife_covariance(rows: Sequence[Sequence[mp.mpf]]) -> list[list[mp.mpf]]:
    count = len(rows)
    width = len(rows[0])
    means = [mp.fsum(row[j] for row in rows)/count for j in range(width)]
    factor = mp.mpf(count-1)/count
    return [[
        factor*mp.fsum((row[i]-means[i])*(row[j]-means[j]) for row in rows)
        for j in range(width)
    ] for i in range(width)]


def add_matrix(left: list[list[mp.mpf]], right: list[list[mp.mpf]]) -> list[list[mp.mpf]]:
    return [[a+b for a, b in zip(lrow, rrow)] for lrow, rrow in zip(left, right)]


def quadratic_2(residual: Sequence[mp.mpf], covariance: list[list[mp.mpf]]) -> mp.mpf:
    a, b = covariance[0]
    _, d = covariance[1]
    determinant = a*d-b*b
    if determinant <= 0:
        raise ValueError("heldout residual covariance is not positive definite")
    x, y = residual
    return (d*x*x-2*b*x*y+a*y*y)/determinant


def block_covariance(parent_run: Run, child_run: Run, parent: dict[str, Any], child: dict[str, Any]) -> list[list[mp.mpf]]:
    blocks = []
    for run, full in ((parent_run, parent), (child_run, child)):
        batches = sorted({batch for _, batch in run.groups})
        deletes = [size_estimate(run, batch)["base_vector"] for batch in batches]
        blocks.append(jackknife_covariance(deletes))
    width = len(blocks[0])+len(blocks[1])
    result = [[mp.mpf(0) for _ in range(width)] for _ in range(width)]
    offset = 0
    for block in blocks:
        for i, row in enumerate(block):
            for j, value in enumerate(row):
                result[offset+i][offset+j] = value
        offset += len(block)
    return result


def _text(value: mp.mpf) -> str:
    return mp.nstr(value, 25)


def build_report(parent_run: Run, child_run: Run, prediction: dict[str, Any]) -> dict[str, Any]:
    phase = prediction["phase"]
    if not phase.get("production_authorized"):
        raise ValueError("production phase was not prereveal-authorized")
    validate_run(parent_run, phase["sizes"]["N65"], phase)
    validate_run(child_run, phase["sizes"]["N130"], phase)
    parent = size_estimate(parent_run)
    child = size_estimate(child_run)
    point = solve_training(parent, child)
    parent_batches = sorted({batch for _, batch in parent_run.groups})
    child_batches = sorted({batch for _, batch in child_run.groups})
    parent_delete = [
        solve_training(size_estimate(parent_run, batch), child)["child_residual"]
        for batch in parent_batches
    ]
    child_delete = [
        solve_training(parent, size_estimate(child_run, batch))["child_residual"]
        for batch in child_batches
    ]
    residual_covariance = add_matrix(
        jackknife_covariance(parent_delete), jackknife_covariance(child_delete))
    chi_square = quadratic_2(point["child_residual"], residual_covariance)
    base_covariance = block_covariance(parent_run, child_run, parent, child)
    transfer = child["P4_shell_D"]/parent["P4_shell_D"]
    return {
        "schema": "matching-one/q2-local-d4-uv-annihilator-score/v1",
        "status": "frozen_preregistered_reveal",
        "issues": [215, 225, 275],
        "runner_commit": phase["runner_commit"],
        "binary_sha256": phase["binary_sha256_arm64"],
        "primary": "heldout N130 two-orientation residual after N65 thermal/H4 training",
        "q_observer_role": "exact contact control only; never a primary field matrix element",
        "sizes": {
            "N65": _serialize_size(parent), "N130": _serialize_size(child),
        },
        "fit": {
            "thermal_beta": _text(point["thermal_beta"]),
            "H4_amplitude_N65": _text(point["H4_amplitude_N65"]),
            "training_determinant": _text(point["training_determinant"]),
            "target_child_over_parent": _text(point["target_scale"]),
            "observed_P4_shell_child_over_parent": _text(transfer),
            "child_prediction": [_text(v) for v in point["child_prediction"]],
            "child_residual": [_text(v) for v in point["child_residual"]],
            "child_residual_covariance_2x2": [[_text(v) for v in row] for row in residual_covariance],
            "chi_square": _text(chi_square), "dof": 2,
            "survival_p": _text(mp.exp(-chi_square/2)),
        },
        "full_sufficient_statistics": {
            "coordinate_order": parent["coordinate_order"]+child["coordinate_order"],
            "point_vector": [_text(v) for v in parent["base_vector"]+child["base_vector"]],
            "delete_one_covariance_40x40_block_diagonal_across_independent_sizes": [
                [_text(v) for v in row] for row in base_covariance
            ],
        },
        "claim_boundary": [
            "R4-R2 cancels only an R-independent local contact term under the frozen cutoff model.",
            "The N65 thermal-shell nuisance and H4 amplitude are trained before the N130 child is scored.",
            "The q-product contact identities are exact controls, not evidence for a field coupling.",
            "A surviving two-coordinate heldout score nominates the local source; it does not prove Q4 epsilon.",
        ],
    }


def _serialize_size(value: dict[str, Any]) -> dict[str, Any]:
    return {
        "N": value["N"], "intrinsic_center": _text(value["p"]),
        "cos4": {name: _text(v) for name, v in value["cos4"].items()},
        "P4_shell_D": _text(value["P4_shell_D"]),
        "P4_shell_S": _text(value["P4_shell_S"]),
        "contact_identity_max_residual": _text(value["contact_max_residual"]),
        "orientations": {
            orientation: {name: _text(v) for name, v in values.items()}
            for orientation, values in value["orientations"].items()
        },
    }


def render_markdown(report: dict[str, Any]) -> str:
    fit = report["fit"]
    sizes = report["sizes"]
    return "\n".join([
        "# q2 local landing-D4 UV-annihilator", "",
        "Primary: paired-orientation mean local-D4. `q=A_top` is an exact contact control only.", "",
        "| size | intrinsic p | P4 D shell | P4 thermal S shell | contact max residual |",
        "|---:|---:|---:|---:|---:|",
        f"| 65 | {sizes['N65']['intrinsic_center']} | {sizes['N65']['P4_shell_D']} | {sizes['N65']['P4_shell_S']} | {sizes['N65']['contact_identity_max_residual']} |",
        f"| 130 | {sizes['N130']['intrinsic_center']} | {sizes['N130']['P4_shell_D']} | {sizes['N130']['P4_shell_S']} | {sizes['N130']['contact_identity_max_residual']} |",
        "", "## Heldout q2 score", "",
        f"- target child/parent: `{fit['target_child_over_parent']}`",
        f"- raw P4-shell child/parent: `{fit['observed_P4_shell_child_over_parent']}`",
        f"- N65-trained thermal beta: `{fit['thermal_beta']}`",
        f"- N130 residual: `{fit['child_residual']}`",
        f"- chi-square: `{fit['chi_square']}/2`, survival p `{fit['survival_p']}`",
        "", "## Boundary", "",
        *[f"- {line}" for line in report["claim_boundary"]], "",
    ])


def main(argv: Sequence[str] | None = None) -> int:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prediction", type=Path, default=root/"experiments/p275_q2_local_d4_uv_20260829.yaml")
    parser.add_argument("--parent-prefix", type=Path, required=True)
    parser.add_argument("--child-prefix", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path)
    parser.add_argument("--dps", type=int, default=60)
    args = parser.parse_args(argv)
    mp.mp.dps = args.dps
    prediction = yaml.safe_load(args.prediction.read_text(encoding="utf-8"))
    report = build_report(load_run(args.parent_prefix), load_run(args.child_prefix), prediction)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    if args.markdown:
        args.markdown.write_text(render_markdown(report), encoding="utf-8")
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Phase-first score of the preregistered N65->N130->N260 mean-JD4 chain."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any, Sequence

import mpmath as mp

from score_marked_birth_path import combine, cos4, projected, read_path


SIZES = (65, 130, 260)
NAMES = {65: "N65_5m", 130: "N130_5m", 260: "N260_5m"}
RUNNER_COMMIT = "6899b119db5b16e9918db53abf5280d990eb6653"
PREREG_COMMIT = "0b3cebf9cd1b536c859d3fad1591345725a068f7"
BINARY_SHA256 = "46d8a2690b9a3b1899b3fe61e9a2c16019cb39487d493998c477ca302eaa1223"
EXPECTED_METADATA = {
    "git_commit": RUNNER_COMMIT,
    "samples_per_pair": 5_000_000,
    "batches": 100,
    "seed": 202608290315,
    "replica_counter_first": 9_300_000_000,
    "replica_counter_last_exclusive": 9_305_000_000,
    "threads_requested": 16,
}


def target_edge() -> mp.mpf:
    return -mp.power(2, -mp.mpf(13) / 8)


def target_two_step() -> mp.mpf:
    return mp.power(4, -mp.mpf(13) / 8)


def _complex_point(first, second) -> tuple[mp.mpf, mp.mpc]:
    p, _, orientations = projected(first, second)
    delta = cos4(first[0].a, first[0].b) - cos4(second[0].a, second[0].b)
    value = mp.mpc(
        (orientations["first"]["J_D_re"] - orientations["second"]["J_D_re"]) / delta,
        (orientations["first"]["J_D_im"] - orientations["second"]["J_D_im"]) / delta,
    )
    return p, value


def load_size(prefix: Path) -> tuple[mp.mpf, mp.mpc, list[mp.mpc], list[mp.mpf]]:
    groups = read_path(Path(str(prefix) + ".path.csv"))
    n = next(iter(groups))[0]
    batches = sorted(
        set(key[2] for key in groups if key[1] == "first")
        & set(key[2] for key in groups if key[1] == "second")
    )
    if len(batches) != 100:
        raise ValueError(f"expected 100 batches for N={n}, got {len(batches)}")
    first = combine([groups[(n, "first", batch)] for batch in batches])
    second = combine([groups[(n, "second", batch)] for batch in batches])
    center, point = _complex_point(first, second)
    deletes: list[mp.mpc] = []
    delete_centers: list[mp.mpf] = []
    for omitted in batches:
        first_delete = combine(
            [groups[(n, "first", batch)] for batch in batches if batch != omitted]
        )
        second_delete = combine(
            [groups[(n, "second", batch)] for batch in batches if batch != omitted]
        )
        p_delete, value_delete = _complex_point(first_delete, second_delete)
        delete_centers.append(p_delete)
        deletes.append(value_delete)
    return center, point, deletes, delete_centers


def jackknife_covariance(rows: Sequence[Sequence[mp.mpf]]) -> list[list[mp.mpf]]:
    n = len(rows)
    means = [mp.fsum(row[index] for row in rows) / n for index in range(len(rows[0]))]
    factor = mp.mpf(n - 1) / n
    return [
        [
            factor
            * mp.fsum(
                (row[i] - means[i]) * (row[j] - means[j]) for row in rows
            )
            for j in range(len(means))
        ]
        for i in range(len(means))
    ]


def _jackknife_se(values: Sequence[mp.mpf]) -> mp.mpf:
    covariance = jackknife_covariance([[value] for value in values])
    return mp.sqrt(max(mp.mpf(0), covariance[0][0]))


def _phase(value: mp.mpc) -> mp.mpf:
    return mp.atan2(mp.im(value), mp.re(value))


def _text(value: mp.mpf | mp.mpc) -> str:
    return mp.nstr(value, 24)


def _complex_json(value: mp.mpc) -> dict[str, str]:
    return {"re": _text(mp.re(value)), "im": _text(mp.im(value)), "abs": _text(abs(value))}


def _inverse2(matrix: Sequence[Sequence[mp.mpf]]) -> list[list[mp.mpf]]:
    determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
    if determinant <= 0:
        raise ValueError("ratio covariance is not positive definite")
    return [
        [matrix[1][1] / determinant, -matrix[0][1] / determinant],
        [-matrix[1][0] / determinant, matrix[0][0] / determinant],
    ]


def score_ratio(
    numerator: mp.mpc,
    denominator: mp.mpc,
    delete_numerator: Sequence[mp.mpc],
    delete_denominator: Sequence[mp.mpc],
    target: mp.mpf,
) -> dict[str, Any]:
    ratio = numerator / denominator
    deletes = [left / right for left, right in zip(delete_numerator, delete_denominator)]
    residual_rows = [[mp.re(value - target), mp.im(value)] for value in deletes]
    covariance = jackknife_covariance(residual_rows)
    inverse = _inverse2(covariance)
    residual = [mp.re(ratio - target), mp.im(ratio)]
    chi2 = mp.fsum(residual[i] * inverse[i][j] * residual[j] for i in range(2) for j in range(2))
    phase_values = [_phase(value / target) for value in deletes]
    magnitude_values = [abs(value) for value in deletes]
    exponent_values = [-mp.log(abs(value), 2) for value in deletes]
    phase_se = _jackknife_se(phase_values)
    return {
        "ratio": _complex_json(ratio),
        "target": _text(target),
        "phase_residual_radians": _text(_phase(ratio / target)),
        "phase_residual_se": _text(phase_se),
        "phase_residual_z": _text(_phase(ratio / target) / phase_se),
        "magnitude": _text(abs(ratio)),
        "magnitude_se": _text(_jackknife_se(magnitude_values)),
        "effective_area_exponent": _text(-mp.log(abs(ratio), 2)),
        "effective_area_exponent_se": _text(_jackknife_se(exponent_values)),
        "target_residual": {"re": _text(residual[0]), "im": _text(residual[1])},
        "target_residual_covariance": [[_text(value) for value in row] for row in covariance],
        "target_chi2_2dof": _text(chi2),
        "target_chi2_p_2dof": _text(mp.exp(-chi2 / 2)),
    }


def audit_file(path: Path) -> dict[str, int]:
    totals = {
        name: 0
        for name in (
            "endpoint_failures",
            "site_failures",
            "line_failures",
            "local_mark_failures",
            "index_mismatches",
        )
    }
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            for name in totals:
                totals[name] += int(row[name])
    return totals


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def provenance(raw: Path, provenance_dir: Path, n: int) -> dict[str, Any]:
    name = NAMES[n]
    metadata_path = raw / f"{name}.metadata.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    sums_path = provenance_dir / f"N{n}.SHA256SUMS.remote"
    remote_sums = {}
    for line in sums_path.read_text(encoding="utf-8").splitlines():
        digest, filename = line.split(maxsplit=1)
        remote_sums[Path(filename).name] = digest
    sufficient = {}
    for suffix in ("path.csv", "hist.csv", "moments.csv", "metadata.json", "complement_audit.csv"):
        path = raw / f"{name}.{suffix}"
        actual = _sha256(path)
        expected = remote_sums[path.name]
        if actual != expected:
            raise ValueError(f"SHA mismatch for {path}")
        sufficient[path.name] = actual
    for key, expected in EXPECTED_METADATA.items():
        if metadata.get(key) != expected:
            raise ValueError(
                f"metadata mismatch for N={n}, {key}: {metadata.get(key)!r} != {expected!r}"
            )
    if remote_sums.get("threshold_rank_integer_period_mc") != BINARY_SHA256:
        raise ValueError(f"binary SHA mismatch for N={n}")
    audit = audit_file(raw / f"{name}.complement_audit.csv")
    if any(audit.values()):
        raise ValueError(f"nonzero complement audit for N={n}: {audit}")
    return {
        "metadata": {
            key: metadata[key]
            for key in (
                "git_commit",
                "samples_per_pair",
                "batches",
                "seed",
                "replica_counter_first",
                "replica_counter_last_exclusive",
                "threads_requested",
                "elapsed_seconds",
            )
        },
        "binary_sha256": remote_sums["threshold_rank_integer_period_mc"],
        "remote_sparse_marked_births_sha256": remote_sums[f"{name}.marked_births.csv"],
        "local_sufficient_statistic_sha256": sufficient,
        "complement_audit": audit,
    }


def build_report(raw: Path, provenance_dir: Path) -> dict[str, Any]:
    centers: dict[int, mp.mpf] = {}
    points: dict[int, mp.mpc] = {}
    deletes: dict[int, list[mp.mpc]] = {}
    for n in SIZES:
        center, point, delete, _ = load_size(raw / NAMES[n])
        centers[n], points[n], deletes[n] = center, point, delete

    vectors = [
        [coordinate for n in SIZES for coordinate in (mp.re(deletes[n][index]), mp.im(deletes[n][index]))]
        for index in range(100)
    ]
    cross_covariance = jackknife_covariance(vectors)
    edge1 = score_ratio(points[130], points[65], deletes[130], deletes[65], target_edge())
    edge2 = score_ratio(points[260], points[130], deletes[260], deletes[130], target_edge())
    two_step = score_ratio(points[260], points[65], deletes[260], deletes[65], target_two_step())
    joint_exponents = [
        (
            -mp.log(abs(deletes[130][index] / deletes[65][index]), 2)
            -mp.log(abs(deletes[260][index] / deletes[130][index]), 2)
        )
        / 2
        for index in range(100)
    ]
    joint_exponent = (
        -mp.log(abs(points[130] / points[65]), 2)
        -mp.log(abs(points[260] / points[130]), 2)
    ) / 2
    return {
        "schema": "matching-one/mean-jd4-q2-chain-score/v1",
        "preregistration_commit": PREREG_COMMIT,
        "runner_commit": RUNNER_COMMIT,
        "metric_order": ["N65_re", "N65_im", "N130_re", "N130_im", "N260_re", "N260_im"],
        "points": {
            str(n): {"intrinsic_center": _text(centers[n]), "P4_mean_J_D4": _complex_json(points[n])}
            for n in SIZES
        },
        "primary_phase_first": {
            "N65_to_N130": edge1,
            "N130_to_N260": edge2,
            "N65_to_N260": two_step,
            "joint_edge_exponent": _text(joint_exponent),
            "joint_edge_exponent_se": _text(_jackknife_se(joint_exponents)),
            "joint_residual_from_13_over_8": _text(joint_exponent - mp.mpf(13) / 8),
        },
        "full_cross_size_delete_one_covariance": [[_text(value) for value in row] for row in cross_covariance],
        "provenance": {str(n): provenance(raw, provenance_dir, n) for n in SIZES},
        "failed_setup_excluded": {
            "server": "DevEnvC_ZyTrST",
            "path": "/workspace/Matching-One-mean-jd4-q2-chain",
            "event": "initial clone/fetch stalled (PIDs 24929/24937); first N65 process disappeared without metadata or accepted output",
            "action": "stalled fetch stopped; accepted N65 rerun used explicit /workspace/mean-jd4-N65-5m with verified binary",
        },
        "remote_raw": {
            "65": "/workspace/mean-jd4-N65-5m/raw/N65_5m.* on DevEnvC_ZyTrST",
            "130": "/workspace/Matching-One-mean-jd4-q2-chain/results/server-20260829/mean-jd4-q2-chain/raw/N130_5m.* on DevEnvC_XPk2PZ",
            "260": "/workspace/Matching-One-mean-jd4-q2-chain/results/server-20260829/mean-jd4-q2-chain/raw/N260_5m.* on DevEnvC_HZsCM6",
        },
    }


def render_markdown(report: dict[str, Any]) -> str:
    lines = [
        "# Three-generation mean-JD4 q2 score",
        "",
        "The score reads only the preregistered complex mean source; the contact-closed connected response is excluded.",
        "",
        "| N | intrinsic p | Re P4 J_D4 | Im P4 J_D4 | |J| |",
        "|---:|---:|---:|---:|---:|",
    ]
    for n in SIZES:
        point = report["points"][str(n)]
        z = point["P4_mean_J_D4"]
        lines.append(f"| {n} | {point['intrinsic_center']} | {z['re']} | {z['im']} | {z['abs']} |")
    lines += ["", "## Preregistered phase-first score", ""]
    for name, label in (("N65_to_N130", "65->130"), ("N130_to_N260", "130->260"), ("N65_to_N260", "65->260")):
        score = report["primary_phase_first"][name]
        lines.append(
            f"- `{label}` ratio `{score['ratio']['re']}+({score['ratio']['im']})i`; "
            f"phase residual `{score['phase_residual_radians']} +/- {score['phase_residual_se']}` rad "
            f"(`z={score['phase_residual_z']}`); effective exponent "
            f"`{score['effective_area_exponent']} +/- {score['effective_area_exponent_se']}`; "
            f"full-target `chi2={score['target_chi2_2dof']}` (2 dof, "
            f"`p={score['target_chi2_p_2dof']}`)."
        )
    lines += [
        "",
        f"Joint edge exponent: `{report['primary_phase_first']['joint_edge_exponent']} +/- {report['primary_phase_first']['joint_edge_exponent_se']}`.",
        "",
        "**Decision.** Both one-step q2 transfers reject their preregistered negative-real "
        "`-2^(-13/8)` phase target. The N65->N130 rejection is already decisive; the "
        "held-out N130->N260 edge independently rejects it. Thus the 20k-pilot phase "
        "coincidence did not reproduce at 5m. The nearly positive two-step phase is only "
        "a cancellation of two incompatible one-step phases and cannot rescue the frozen mechanism.",
        "",
        "This rejects the **primitive-line mean `J_D4` source as the proposed q2 H4 carrier**. "
        "It does not reject the established global H4 response, nor a different external/local "
        "observer coupling to the rank-birth stream.",
        "",
        "## Scientific card",
        "",
        "1. **Mechanism-space change:** remove the bare primitive-line mean-`J_D4` q2 staircase from the live H4 mechanisms.",
        "2. **Not proved:** no claim is made against global H4, the thermal `Q4 epsilon` channel, or externally marked rank births.",
        "3. **Observer / sector / source / geometry:** complex mean `chi4(ell)(I12-I01)` / rank-one Alexander-odd / topology-only primitive line / q2 square-torus chain.",
        "4. **Dependency group:** one frozen runner, one frozen counter range, paired orientations within each size, three independent Huawei size jobs; full aligned-batch covariance retained.",
        "5. **Upweighted observation:** a local landing, bulk Betti, seam-charged, or other non-topology-only observer is now required before reconnecting rank births to the H4 field.",
        "",
        "## Provenance boundary",
        "",
        "All downloaded sufficient statistics match the remote SHA manifests. The GiB-scale sparse birth tables remain on their declared Huawei paths.",
        "The failed first Zy setup is explicitly excluded; only the metadata-bearing rerun is primary.",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw", type=Path, required=True)
    parser.add_argument("--provenance", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--markdown", type=Path)
    parser.add_argument("--dps", type=int, default=50)
    args = parser.parse_args()
    mp.mp.dps = args.dps
    report = build_report(args.raw, args.provenance)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if args.markdown:
        args.markdown.parent.mkdir(parents=True, exist_ok=True)
        args.markdown.write_text(render_markdown(report), encoding="utf-8")


if __name__ == "__main__":
    main()

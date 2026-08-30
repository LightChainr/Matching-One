#!/usr/bin/env python3
"""Exact mark-fiber decomposition of the P334 N=9 Y=0 obstruction."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
from typing import Any, Iterable

from p334_n9_reservoir_obstruction import (
    _coarse_classes,
    candidate_rows,
    two_output_mark_targets as build_two_output_mark_targets,
)
from p334_tm_translation_orbit_hall import (
    descriptor,
    normalize_target,
    transverse_reservoir_targets,
)


SCHEMA = "matching-one/p334-n9-mark-fiber-theorem/v1"
PROTOCOL = Path("analysis/p334_n9_mark_fiber_protocol.json")


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def payload_sha256(payload: Any) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
    ).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def count_histogram(values: Iterable[int]) -> dict[str, int]:
    return {
        str(value): count
        for value, count in sorted(Counter(values).items())
    }


def active_side(marks, face) -> str:
    base, left, right = face
    left_rank = marks[base | (1 << left)][0]
    right_rank = marks[base | (1 << right)][0]
    if (left_rank == 2) == (right_rank == 2):
        raise AssertionError("an M face must have exactly one rank-2-producing mark")
    return "L" if left_rank == 2 else "R"


def image_for_row(row):
    n, _matrix, _geometry, _carrier, marks, line, _layer, faces = row
    _keys, representatives, permutations, inverses = _coarse_classes(row)
    old_targets = set()
    two_mark_targets = set()
    for source in representatives:
        old_targets.update(
            normalize_target(target, inverses)
            for target in transverse_reservoir_targets(
                marks, line, source, permutations, n, transport=True
            )
        )
        two_mark_targets.update(
            build_two_output_mark_targets(
                marks, line, source, permutations, inverses, n
            )
        )
    full_targets = {
        normalize_target(("MM", first, second), inverses)
        for first in faces["M"]
        for second in faces["M"]
    }
    return old_targets, two_mark_targets, full_targets


def fiber_signature(row, old_targets, two_mark_targets, full_targets) -> dict[str, Any]:
    _n, _matrix, _geometry, _carrier, marks, _line, _layer, _faces = row
    anchors = sorted({target[1] for target in full_targets})

    def anchor_histogram(targets):
        counts = Counter(target[1] for target in targets)
        return count_histogram(counts.get(anchor, 0) for anchor in anchors)

    def activity_counts(targets):
        return {
            activity: count
            for activity, count in sorted(
                Counter(
                    active_side(marks, target[1])
                    + active_side(marks, target[2])
                    for target in targets
                ).items()
            )
        }

    def base_pair_histogram(targets):
        counts = Counter((target[1][0], target[2][0]) for target in targets)
        return {
            "base_pair_count": len(counts),
            "targets_per_base_pair_histogram": count_histogram(counts.values()),
        }

    def activity_base_pair_histograms(targets):
        counts: dict[str, Counter] = defaultdict(Counter)
        for target in targets:
            activity = active_side(marks, target[1]) + active_side(marks, target[2])
            counts[activity][(target[1][0], target[2][0])] += 1
        return {
            activity: {
                "base_pair_count": len(rows),
                "targets_per_base_pair_histogram": count_histogram(rows.values()),
            }
            for activity, rows in sorted(counts.items())
        }

    signature = {
        "anchor_count": len(anchors),
        "old_anchor_partner_count_histogram": anchor_histogram(old_targets),
        "two_mark_anchor_partner_count_histogram": anchor_histogram(two_mark_targets),
        "full_anchor_partner_count_histogram": anchor_histogram(full_targets),
        "old_activity_counts": activity_counts(old_targets),
        "two_mark_activity_counts": activity_counts(two_mark_targets),
        "full_activity_counts": activity_counts(full_targets),
        "old_base_pair_fibers": base_pair_histogram(old_targets),
        "two_mark_base_pair_fibers": base_pair_histogram(two_mark_targets),
        "full_base_pair_fibers": base_pair_histogram(full_targets),
        "old_activity_base_pair_fibers": activity_base_pair_histograms(old_targets),
        "two_mark_activity_base_pair_fibers": activity_base_pair_histograms(
            two_mark_targets
        ),
    }
    signature["sha256"] = payload_sha256(signature)
    return signature


def row_result(index: int, row) -> dict[str, Any]:
    n, matrix, _geometry, carrier, _marks, line, layer, faces = row
    old_targets, two_mark_targets, full_targets = image_for_row(row)
    signature = fiber_signature(row, old_targets, two_mark_targets, full_targets)
    if not old_targets <= two_mark_targets:
        raise AssertionError("the old image must be contained in the two-mark image")
    if two_mark_targets != full_targets:
        raise AssertionError("the two-mark image must equal the full MM orbit reservoir")
    return {
        "row_index": index,
        **descriptor(n, matrix, carrier, line, layer, faces),
        "image_counts": {
            "old": len(old_targets),
            "two_mark": len(two_mark_targets),
            "full_MM": len(full_targets),
            "two_mark_gain": len(two_mark_targets - old_targets),
        },
        "set_relations": {
            "old_subset_two_mark": True,
            "two_mark_equals_full_MM": True,
        },
        "fiber_signature": signature,
    }


def validate_protocol(repo: Path) -> dict[str, Any]:
    protocol_path = repo / PROTOCOL
    protocol = json.loads(protocol_path.read_text(encoding="utf-8"))
    if protocol["status"] != "frozen_before_fiber_enumeration":
        raise AssertionError("protocol is not the frozen pre-enumeration contract")
    for relative, expected in protocol["frozen_inputs_sha256"].items():
        actual = sha256_path(repo / relative)
        if actual != expected:
            raise AssertionError(f"frozen input changed: {relative}")
    return {
        "path": str(PROTOCOL),
        "sha256": sha256_path(protocol_path),
        "parent_commit": protocol["parent_commit"],
        "input_hashes_verified": True,
    }


def build_result(repo: Path) -> dict[str, Any]:
    protocol = validate_protocol(repo)
    rows = [row_result(index, row) for index, row in candidate_rows()]
    signature_hashes = {row["fiber_signature"]["sha256"] for row in rows}
    if len(signature_hashes) != 1:
        raise AssertionError("the six N9 rows do not share one frozen fiber signature")

    first = rows[0]
    counts = first["image_counts"]
    motifs = first["motifs"]
    n = first["N"]
    m_count = motifs["M"]
    demand = 4 * motifs["D"] * motifs["F"] // n
    if (demand, counts["old"], counts["two_mark"], counts["two_mark_gain"]) != (
        6912,
        4752,
        20736,
        15984,
    ):
        raise AssertionError("the frozen N9 image arithmetic changed")

    signature = first["fiber_signature"]
    if signature["old_anchor_partner_count_histogram"] != {
        "0": 12,
        "67": 8,
        "70": 4,
        "164": 24,
    }:
        raise AssertionError("the primary old anchor-fiber theorem changed")
    if signature["two_mark_anchor_partner_count_histogram"] != {"432": 48}:
        raise AssertionError("the two-mark complete-fiber theorem changed")

    derived = {
        "M": m_count,
        "normalized_first_M_anchors": m_count // n,
        "source_demand": demand,
        "source_demand_in_M_units": demand // m_count,
        "old_image": counts["old"],
        "old_image_in_M_units": counts["old"] // m_count,
        "deficiency": demand - counts["old"],
        "deficiency_in_M_units": (demand - counts["old"]) // m_count,
        "old_coverage_fraction": "11/16",
        "deficiency_fraction": "5/16",
        "two_mark_image": counts["two_mark"],
        "two_mark_image_in_M_units": counts["two_mark"] // m_count,
        "two_mark_gain": counts["two_mark_gain"],
        "two_mark_gain_in_M_units": counts["two_mark_gain"] // m_count,
        "full_MM_formula": "M^2/N = 432^2/9 = 20736 = 48M",
        "old_anchor_sum": "12*0 + 8*67 + 4*70 + 24*164 = 4752 = 11M",
    }
    if (
        derived["source_demand_in_M_units"],
        derived["old_image_in_M_units"],
        derived["deficiency_in_M_units"],
        derived["two_mark_image_in_M_units"],
        derived["two_mark_gain_in_M_units"],
    ) != (16, 11, 5, 48, 37):
        raise AssertionError("the M-unit theorem changed")

    return {
        "schema": SCHEMA,
        "status": "exact_bounded_mark_fiber_theorem",
        "issue": 334,
        "protocol": protocol,
        "rows": rows,
        "common_fiber_signature": signature,
        "derived_identities": derived,
        "conclusion": {
            "positive": (
                "All six N9 Y=0 rows have the same exact nonuniform first-M anchor "
                "fiber histogram. Its sum gives 11M old targets against 16M demand, "
                "hence the 5M = 5/16 deficit without a max-flow calculation. Releasing "
                "two output marks makes all 48 anchors complete with M partners, giving 48M."
            ),
            "negative": (
                "The old 11M image is not a union of eleven complete M-sized first-anchor "
                "fibers. The 11/16 fraction is an exact global sum over nonuniform fibers, "
                "not a literal eleven-of-sixteen local-slot rule."
            ),
        },
        "crosswalk_to_birth_age_review": {
            "commit": "fee33287cf4830e07ccef6177f43034add02256e",
            "relationship": (
                "The birth-age review concerns path-time memory of (tau1,ell,tau2) "
                "and the direct 0-to-2 collision mass D_N. This certificate concerns "
                "the image capacity of a bounded N9 D x F switching reservoir. Neither "
                "identity implies the other, and the exact N10 1/57 witness is not rerun here."
            ),
            "next_existing_archive_test": {
                "inputs": "same-batch (tau1, primitive ell, tau2, DIRECT_RANK2) records already emitted by the projective-birth archive",
                "birth_age_estimand": (
                    "At each occupied count k and line ell, compare the discrete exit "
                    "hazards P(tau2=k+1 | tau2>k,tau1=j,ell) across entry times j; "
                    "freeze pooling bins before scoring."
                ),
                "collision_estimand": "D_N = mean[DIRECT_RANK2] = P(tau1=tau2)",
                "covariance_contract": (
                    "Compute both estimands in every original batch and retain their full "
                    "joint covariance. A covariance-whitened common-age score must count "
                    "DIRECT_RANK2 paths as a separate no-plateau channel rather than silently "
                    "dropping them from the risk set."
                ),
                "primary_null": (
                    "Conditional on (k,ell,tau2>k), exit hazard is independent of tau1; "
                    "D_N is reported jointly but is not combined as an independent vote."
                ),
                "compute": "no new simulation; one pass over existing sparse joint records plus batch-level scoring",
            },
        },
        "claim_boundary": [
            "Exact for the six frozen N9 matching/layer4/Y=0 rows.",
            "The target-image deficit follows without max flow; source injection still uses the separately certified Hall result.",
            "No arbitrary-HNF formula or saturation theorem is claimed.",
            "No slot, source, phase or provenance label contributes target capacity.",
        ],
    }


def render_markdown(result: dict[str, Any]) -> str:
    derived = result["derived_identities"]
    signature = result["common_fiber_signature"]
    lines = [
        "# P334 N9 mark-fiber theorem",
        "",
        "## Exact result",
        "",
        "All six frozen matching/layer-4/Y=0 rows have the same normalized target-fiber signature.",
        "The old one-carrier/one-mark image, grouped by the normalized first ordered M face, is",
        "",
        "| partners in anchor fiber | number of anchors |",
        "|---:|---:|",
    ]
    for partners, anchors in signature["old_anchor_partner_count_histogram"].items():
        lines.append(f"| {partners} | {anchors} |")
    lines.extend(
        [
            "",
            f"Therefore `{derived['old_anchor_sum']}`. The coarse source demand is "
            f"`{derived['source_demand']}=16M`, so the exact image deficit is "
            f"`{derived['deficiency']}=5M`, or `{derived['deficiency_fraction']}` of demand.",
            "",
            "This is a target-image counting identity; it does not invoke maximum flow.",
            "",
            "## Two-output-mark closure",
            "",
            "After keeping both bases fixed and releasing exactly two output marks, every one of the 48 normalized first-M anchors has all 432 partners. Thus",
            "",
            f"`{derived['full_MM_formula']}`.",
            "",
            f"The new image is `{derived['two_mark_gain']}=37M` beyond the old image.",
            "",
            "## Structural correction",
            "",
            result["conclusion"]["negative"],
            "The finer nonuniform histogram is the theorem: 12 anchors are absent, 8 have 67 partners, 4 have 70, and 24 have 164.",
            "",
            "## Crosswalk to the birth-age result",
            "",
            result["crosswalk_to_birth_age_review"]["relationship"],
            "The next production statistic should reuse the existing same-batch `(tau1,ell,tau2,DIRECT_RANK2)` archive: test equality of exit hazards across `tau1` within each `(k,ell)` risk set, estimate `D_N=mean(DIRECT_RANK2)` in the same batches, and retain their full covariance. Direct births are a separate no-plateau channel, not missing risk-set rows.",
            "",
            "## Boundary",
            "",
        ]
    )
    lines.extend(f"- {item}" for item in result["claim_boundary"])
    lines.append("")
    return "\n".join(lines)


def render_card(result: dict[str, Any]) -> str:
    derived = result["derived_identities"]
    return "\n".join(
        [
            "# Scientific card: P334 N9 mark-fiber theorem",
            "",
            "- **Mechanism space changed:** the common N9 Hall failure is resolved into an exact output-mark fiber anatomy rather than six unexplained equal max-flow totals.",
            f"- **Result:** the nonuniform anchor sum is `{derived['old_anchor_sum']}` against `16M` demand, hence `5M=5/16` deficit; two-mark release gives all `48M` targets.",
            "- **Unexpected correction:** `11M` is not eleven complete first-anchor fibers; it is a global sum over the frozen `0/67/70/164` partner-count classes.",
            "- **Not proved:** no arbitrary-HNF theorem, no universal Y=0 formula, and no new injection beyond the separately certified Hall result.",
            "- **Observer/sector/source/geometry:** ordinary untagged MM orbit targets | fixed projective line | corrected D x F reservoir | six N9 matching/layer4/Y=0 HNF rows.",
            "- **Dependency group:** one exact bounded N9 certificate built from the same six rows as `a5de4d6`; it is not six independent votes.",
            "- **Crosswalk:** `fee33287` proves N10 birth-age non-Markovianity and a separate direct-birth mass. This N9 target-capacity theorem neither repeats nor implies that path-time result.",
            "- **Next production test:** from the existing same-batch `(tau1,ell,tau2,DIRECT_RANK2)` archive, score entry-time dependence of the exit hazard within each `(k,ell)` risk set jointly with `D_N`, preserving their full batch covariance and treating direct births as a separate no-plateau channel.",
            "- **Next upweight:** derive the `0/67/70/164` classes from a symbolic local incidence rule and test that rule on the next Y=0 hard row before running flow.",
            "",
        ]
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, default=Path("."))
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    parser.add_argument("--card", type=Path, required=True)
    args = parser.parse_args()
    repo = args.repo.resolve()
    result = build_result(repo)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.markdown.parent.mkdir(parents=True, exist_ok=True)
    args.card.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    args.markdown.write_text(render_markdown(result), encoding="utf-8")
    args.card.write_text(render_card(result), encoding="utf-8")
    print(json.dumps(result["derived_identities"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

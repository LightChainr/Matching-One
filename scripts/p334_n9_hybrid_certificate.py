#!/usr/bin/env python3
"""Close the N=9 corrected-reservoir shell with two exact proof paths.

Twenty-four rows are immutable direct translation-orbit matching artifacts.
The four omitted heavy rows are recomputed by the proved twin-class capacity
compression.  This is deliberately a one-shell certificate generator, not a
new general scan interface.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

from p334_tm_coarse_reservoir_hall import coarse_row_audit
from p334_tm_corrected_reservoir_scan import rows_for_order
from p334_tm_translation_orbit_hall import descriptor


SCHEMA = "p334-n9-hybrid-certificate-v1"
SOURCE_SCHEMA = "p334-tm-corrected-reservoir-scan-v1"
SOURCE_COMMIT = "4bb75176c56558084c8397917995026e54420b9f"
ALL_INDICES = tuple(range(28))
HEAVY_INDICES = (11, 17, 20, 26)
DIRECT_INDICES = tuple(index for index in ALL_INDICES if index not in HEAVY_INDICES)
EXPECTED_FAILURES = (1, 3, 6, 9, 15, 24)
EXPECTED_DIRECT_SHA256 = {
    0: "f7e75032eb485d9b85cdd0966a8b0af76357941651146fc923481eba975b88f5",
    1: "95d83b543e87ded76d17601956502529a8d4808f1c243bebbe84741e04733f20",
    2: "4e476905499e3dd829c78556e6a19f6834541612ea671168fce0199b0bb619e5",
    3: "13d87c5406026a8ea60381c9749155507f04973956a9244d30d72ea5d83624a8",
    4: "f13ed8e4d1dadef683198c206e56ea9c406b476a2f4e70007696c980f319b46b",
    5: "53770f6720520dd6d0cf3838b199e814661a7971ce12e6cc478111941dbe6c64",
    6: "a1b794df578e8525b5cdca1464dd6a407a0c89760cadcc39f917f35089f1247d",
    7: "8b0d748277ca426348af45bc2e436c48e34120aec2ad50ada439b2fb61fc7c11",
    8: "53580a9d38ac0fc505be7520c62802b438679f031322fec45b07ad2bb2ef0016",
    9: "6f401a92323ac13d6d93440912570a76b661f48fbaf766eb999b86995d32a198",
    10: "e905393d233837347ff326a5bb6796ca50051f67ca6c453c668538aab993f652",
    12: "632a7457b69ba46ceaf4f410e07fc88a4fde97c77357d659dfaf37500b5732da",
    13: "481260ea308eff2bbdcfe2d287cbc6bd53136eaaf05e21f35ad846e6f0b48aa1",
    14: "dcfd7109daa76bf29625fffe2f5cc8c48ee09f483814c00f155031c573492b11",
    15: "01686d89c19fbc20d0ba33c58671b61ede1edb4b2b5968b907a4acd8717fd0fc",
    16: "2d2d6b223677bf5dd9dc8e8956823c04cda1b27deab9bf861b3dd80a6976a153",
    18: "6e3a8a66d01cc0f588d2cce2b4c1a164a098bfe1a05374dbd325ea14d67a47e6",
    19: "f2d51be0549ad5798247d7757e4274eae4c87a46f33d3479e44b9fa8bbb8f404",
    21: "193310a0a89a46084140308a9df92e7366cb87c7c4bcef772d4bdecd3f791125",
    22: "a922568066b8118e064b24b35d4d1fc156ab39561d3797c9afb42525c48d70a0",
    23: "72778848d3c7b3bb1326dec05189f3650233221aa359b1c1ded707b004460c96",
    24: "b18051d9e85010a555caf154f638d45222de48095b7d9cb84c30699c01d46140",
    25: "2871e3051bff0117bbd98d52b1dac5fb6a70a709eb2de24f1a1c50ed38093ab9",
    27: "5cd658a0aedd91a40524a6cd324656c85cafd81573508e95ceb99b59f146a478",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _jsonable(value: Any) -> Any:
    return json.loads(json.dumps(value))


def validate_direct_artifact(path: Path, index: int, row: tuple) -> dict[str, Any]:
    if path.name != f"row{index}.json":
        raise AssertionError(f"row {index}: unexpected filename {path.name}")
    actual_hash = sha256(path)
    if actual_hash != EXPECTED_DIRECT_SHA256[index]:
        raise AssertionError(f"row {index}: immutable artifact SHA256 mismatch")
    payload = json.loads(path.read_text())
    if payload.get("schema") != SOURCE_SCHEMA:
        raise AssertionError(f"row {index}: source schema mismatch")
    required = {
        "schema", "order", "shard", "complete_order_row_count",
        "selected_row_count", "rows", "summary", "scientific_boundary",
    }
    if set(payload) != required:
        raise AssertionError(f"row {index}: top-level schema keys changed")
    if payload["order"] != 9 or payload["complete_order_row_count"] != 28:
        raise AssertionError(f"row {index}: shell contract mismatch")
    if payload["shard"] != {"index": index, "count": 28}:
        raise AssertionError(f"row {index}: shard identity mismatch")
    if payload["selected_row_count"] != 1 or len(payload["rows"]) != 1:
        raise AssertionError(f"row {index}: artifact is not a singleton shard")
    report = payload["rows"][0]
    if report.get("row_index") != index:
        raise AssertionError(f"row {index}: row index mismatch")
    expected_description = _jsonable(descriptor(*(
        row[0], row[1], row[3], row[5], row[6], row[7]
    )))
    observed_description = {
        key: report[key]
        for key in expected_description
    }
    if observed_description != expected_description:
        raise AssertionError(f"row {index}: deterministic row descriptor mismatch")
    graph = report.get("corrected_reservoir_orbit_graph")
    graph_keys = {
        "Hall_deficiency", "compression_factor", "maximum_degree",
        "maximum_matching", "minimum_degree", "raw_source_tokens",
        "reachable_cover_tokens", "saturates", "source_action_free",
        "source_tokens", "target_action_free", "zero_degree_sources",
    }
    if not isinstance(graph, dict) or set(graph) != graph_keys:
        raise AssertionError(f"row {index}: direct graph schema mismatch")
    if graph["raw_source_tokens"] != 9 * graph["source_tokens"]:
        raise AssertionError(f"row {index}: translation compression mismatch")
    deficiency = graph["source_tokens"] - graph["maximum_matching"]
    if graph["Hall_deficiency"] != deficiency:
        raise AssertionError(f"row {index}: flow deficiency mismatch")
    if graph["saturates"] != (deficiency == 0):
        raise AssertionError(f"row {index}: saturation flag mismatch")
    summary = payload["summary"]
    failed = int(not graph["saturates"])
    if summary["saturated_rows"] != 1 - failed or summary["failed_rows"] != failed:
        raise AssertionError(f"row {index}: shard summary mismatch")
    return {
        "row_index": index,
        "proof_path": "direct_translation_orbit_matching",
        "source_artifact": f"direct/{path.name}",
        "source_sha256": actual_hash,
        "descriptor": expected_description,
        "flow": {
            "source_tokens": graph["source_tokens"],
            "maximum_flow": graph["maximum_matching"],
            "Hall_deficiency": graph["Hall_deficiency"],
            "saturates": graph["saturates"],
        },
        "direct_audit": graph,
    }


def validate_direct_directory(direct_dir: Path) -> list[dict[str, Any]]:
    actual_names = {path.name for path in direct_dir.glob("*.json")}
    expected_names = {f"row{index}.json" for index in DIRECT_INDICES}
    if actual_names != expected_names:
        raise AssertionError(
            f"direct coverage mismatch: missing={sorted(expected_names-actual_names)}, "
            f"extra={sorted(actual_names-expected_names)}"
        )
    rows = rows_for_order(9)
    if len(rows) != 28:
        raise AssertionError("deterministic N9 shell no longer contains 28 rows")
    return [
        validate_direct_artifact(direct_dir / f"row{index}.json", index, rows[index])
        for index in DIRECT_INDICES
    ]


def recompute_heavy_row(index: int, row: tuple) -> dict[str, Any]:
    audit = coarse_row_audit(row, verify_all_twins=False)
    combined = audit["channel_flows"]["combined"]
    compression = audit["source_compression"]
    if (
        combined["total_demand"] != 45360
        or combined["maximum_flow"] != 45360
        or combined["Hall_deficiency"] != 0
        or not combined["saturates"]
    ):
        raise AssertionError(f"heavy row {index}: expected exact 45360/45360 saturation")
    if (
        compression["raw_sources"] != 408240
        or compression["translation_orbit_sources"] != 45360
        or compression["coarse_twin_classes"] != 5040
        or compression["twin_class_size"] != 9
    ):
        raise AssertionError(f"heavy row {index}: coarse compression contract changed")
    return {
        "row_index": index,
        "proof_path": "coarse_twin_capacitated_hall",
        "source_artifact": None,
        "source_sha256": None,
        "descriptor": _jsonable(descriptor(*(
            row[0], row[1], row[3], row[5], row[6], row[7]
        ))),
        "flow": {
            "source_tokens": combined["total_demand"],
            "maximum_flow": combined["maximum_flow"],
            "Hall_deficiency": combined["Hall_deficiency"],
            "saturates": combined["saturates"],
        },
        "coarse_capacity_audit": audit,
    }


def build_result(direct_dir: Path) -> dict[str, Any]:
    rows = rows_for_order(9)
    reports = validate_direct_directory(direct_dir)
    reports.extend(recompute_heavy_row(index, rows[index]) for index in HEAVY_INDICES)
    reports.sort(key=lambda report: report["row_index"])
    indices = [report["row_index"] for report in reports]
    if indices != list(ALL_INDICES):
        raise AssertionError("hybrid proof paths do not cover every N9 row exactly once")
    failures = [
        report["row_index"] for report in reports if not report["flow"]["saturates"]
    ]
    if failures != list(EXPECTED_FAILURES):
        raise AssertionError(f"unexpected N9 failure set: {failures}")
    if any(reports[index]["flow"]["Hall_deficiency"] != 2160 for index in failures):
        raise AssertionError("each frozen N9 failure must retain deficiency 2160")
    return {
        "schema": SCHEMA,
        "source_commit": SOURCE_COMMIT,
        "order": 9,
        "proof_paths": {
            "direct_translation_orbit_matching": {
                "row_indices": list(DIRECT_INDICES),
                "count": 24,
                "contract": "byte-identified singleton scan shards",
            },
            "coarse_twin_capacitated_hall": {
                "row_indices": list(HEAVY_INDICES),
                "count": 4,
                "contract": "exact demand-9 twin compression and integral capacity flow",
            },
        },
        "coverage": {"expected": list(ALL_INDICES), "observed": indices, "complete": True},
        "rows": reports,
        "summary": {
            "row_count": 28,
            "saturated_rows": 22,
            "failed_rows": 6,
            "failed_row_indices": failures,
            "failure_Hall_deficiency": 2160,
            "heavy_row_indices": list(HEAVY_INDICES),
            "heavy_exact_flow": "45360/45360",
            "status": "complete_N9_hybrid_exact_certificate",
        },
        "scientific_card": {
            "question": "Does the corrected two-carrier reservoir saturate every deterministic N9 HNF/line row?",
            "design": "24 immutable direct orbit matchings plus four exact twin-class capacity flows",
            "result": "22/28 saturate; failures are exactly rows 1,3,6,9,15,24, each deficient by 2160",
            "heavy_gate": "rows 11,17,20,26 each saturate exactly at 45360/45360 after proved compression",
            "meaning": "The N9 obstruction is a six-row topology class, not a missing-heavy-row artifact; the repaired MM output-mark channel remains the minimal known local repair for those six rows.",
        },
    }


def render_markdown(result: dict[str, Any]) -> str:
    lines = [
        "# P334 N9 hybrid corrected-reservoir certificate",
        "",
        f"Source commit: `{result['source_commit']}`.",
        "",
        "This closes all deterministic N9 rows through two non-interchangeable exact paths: 24 byte-identified direct translation-orbit matchings and four recomputed coarse twin-class capacity flows.",
        "",
        "| row | path | HNF | line | flow | deficiency | status |",
        "|---:|---|---|---|---:|---:|---|",
    ]
    for row in result["rows"]:
        flow = row["flow"]
        lines.append(
            f"| {row['row_index']} | {row['proof_path']} | "
            f"`{row['descriptor']['matrix']}` | `{row['descriptor']['line']}` | "
            f"{flow['maximum_flow']}/{flow['source_tokens']} | "
            f"{flow['Hall_deficiency']} | {'saturated' if flow['saturates'] else 'FAIL'} |"
        )
    card = result["scientific_card"]
    lines.extend([
        "", "## Scientific card", "",
        f"- Question: {card['question']}",
        f"- Design: {card['design']}",
        f"- Result: {card['result']}",
        f"- Heavy gate: {card['heavy_gate']}",
        f"- Meaning: {card['meaning']}",
        "",
        "The certificate is exhaustive only for the frozen deterministic N9 shell. It does not claim the arbitrary-HNF reservoir theorem.",
        "",
    ])
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--direct-dir", type=Path, required=True)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--markdown", type=Path, required=True)
    args = parser.parse_args()
    result = build_result(args.direct_dir)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.markdown.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    args.markdown.write_text(render_markdown(result))


if __name__ == "__main__":
    main()

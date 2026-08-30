#!/usr/bin/env python3
"""Fail-closed verifier for the dedicated P250 radius-six Level-S certificate."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Mapping

from scipy.stats import f

from score_z5_projective_leg_bivariate_state import read_batches as read_old
from score_z5_projective_leg_radius5_morphism import read_new as read_radius5
from score_z5_projective_leg_radius6_flat import read_new as read_radius6, score


SCHEMA = "matching-one/p250-radius6-level-s-elimination-certificate/v1"
EXPECTED_ELIMINATED_RANKS = (5, 6, 7)
COMPATIBLE_RANK = 8
FLOAT_TOLERANCE = 1e-12


class VerificationError(RuntimeError):
    pass


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerificationError(message)


def close(first: float, second: float) -> bool:
    return abs(first - second) <= FLOAT_TOLERANCE * max(1.0, abs(first), abs(second))


def verify_hashes(root: Path, certificate: Mapping[str, object]) -> dict[str, str]:
    checked = {}
    for name, row in certificate["immutable_inputs"].items():
        path = root / row["path"]
        require(path.is_file(), f"missing immutable input: {name}")
        observed = sha256(path)
        require(observed == row["sha256"], f"immutable input hash mismatch: {name}")
        checked[name] = observed
    return checked


def verify_one_rank(
    certificate_row: Mapping[str, object], score_row: Mapping[str, object],
    rank: int, alpha: float, expected_status: str,
) -> dict:
    frozen = score_row["rank_nulls"][str(rank)]["score"]
    modes = int(frozen["resolved_covariance_modes"])
    denominator = 400 - modes
    critical = float(f.ppf(1.0 - alpha, modes, denominator))
    expected = {
        "T_squared": float(frozen["asymptotic_chi_square"]),
        "Hotelling_F": float(frozen["finite_batch_Hotelling_F"]),
        "df_numerator": modes,
        "df_denominator": denominator,
        "critical_F_alpha_0_01": critical,
        "p": float(frozen["finite_batch_survival_p"]),
    }
    for key in ("T_squared", "Hotelling_F", "critical_F_alpha_0_01", "p"):
        require(close(float(certificate_row[key]), expected[key]), f"rank {rank} statistic mismatch: {key}")
    for key in ("df_numerator", "df_denominator"):
        require(int(certificate_row[key]) == expected[key], f"rank {rank} degree-of-freedom mismatch: {key}")
    observed_status = "eliminated" if expected["p"] < alpha else "compatible_not_eliminated"
    require(observed_status == expected_status, f"rank {rank} frozen decision differs from certificate class")
    require(certificate_row["status"] == expected_status, f"rank {rank} certificate status was edited")
    if expected_status == "eliminated":
        require(expected["Hotelling_F"] > critical, f"rank {rank} did not cross the frozen F threshold")
    else:
        require(expected["Hotelling_F"] <= critical, f"rank {rank} crosses the frozen F threshold")
    return expected


def verify_certificate(root: Path, certificate_path: Path) -> dict:
    certificate = json.loads(certificate_path.read_text())
    require(certificate.get("schema") == SCHEMA, "wrong certificate schema")
    require(certificate.get("claim_level") == "Level-S statistical elimination certificate", "claim level changed")
    require(certificate.get("status") == "verified_by_dedicated_recomputation_required", "certificate status changed")
    require(certificate["invariant_model_class"]["rank_null"] == "rank(H3)<=r", "rank model class changed")
    require(certificate["typed_observer"]["matrix_shape_per_hand"] == [20, 10], "H3 shape changed")
    require(certificate["decision_rule"]["alpha"] == 0.01, "decision alpha changed")
    require(certificate["decision_rule"]["batches_per_dependency_block"] == 400, "batch calibration changed")
    checked_hashes = verify_hashes(root, certificate)

    raw_lock = json.loads((root / certificate["immutable_inputs"]["raw_lock_manifest"]["path"]).read_text())
    stored_score = json.loads((root / certificate["immutable_inputs"]["frozen_score"]["path"]).read_text())
    recomputed = score(
        read_old(root / raw_lock["old4_batches"]),
        read_radius5(root / raw_lock["old5_batches"]),
        read_radius6(root / raw_lock["new6_batches"]),
        raw_lock,
    )
    require(recomputed == stored_score, "frozen score does not reproduce exactly from locked artifacts")
    require(stored_score["decision"] == "rank5_flat_extension_rejected", "stored decision changed")
    require(stored_score["R2_conjugate_kernel_projector_bridge"] == "LOCKED_RANK5_FLAT_EXTENSION_FAILED", "R2 bridge lock changed")

    summaries = {}
    for certificate_hand, score_hand in (("plus", "plus"), ("minus_R2_gauge", "minus")):
        summaries[certificate_hand] = {}
        certificate_rows = certificate["hand_certificates"][certificate_hand]
        score_rows = stored_score["hand_rank_ladders"][score_hand]
        for rank in EXPECTED_ELIMINATED_RANKS:
            summaries[certificate_hand][str(rank)] = verify_one_rank(
                certificate_rows[f"rank_le_{rank}"], score_rows, rank, 0.01, "eliminated"
            )
        summaries[certificate_hand][str(COMPATIBLE_RANK)] = verify_one_rank(
            certificate_rows["rank_le_8"], score_rows, COMPATIBLE_RANK, 0.01, "compatible_not_eliminated"
        )
        require(certificate_rows["certified_rank_lower_bound"] == 8, f"{certificate_hand} lower bound changed")
        require(score_rows["rank_lower_bound_at_alpha"] == 8, f"{certificate_hand} scorer lower bound changed")

    expected_eliminated = [
        "plus rank<=5", "plus rank<=6", "plus rank<=7",
        "minus-R2 rank<=5", "minus-R2 rank<=6", "minus-R2 rank<=7",
    ]
    require(certificate["classification"]["eliminated"] == expected_eliminated, "eliminated class list changed")
    require(certificate["classification"]["survives_as_compatible_only"] == ["plus rank<=8", "minus-R2 rank<=8"], "compatible class list changed")
    require(certificate["bridge_gate"]["status"] == "not_reached", "R2 bridge was reclassified")
    require(certificate["bridge_gate"]["observed_support"] is False, "R2 bridge support was reclassified")

    return {
        "schema": "matching-one/p250-radius6-level-s-verification/v1",
        "verified": True,
        "certificate_id": certificate["certificate_id"],
        "claim_level": certificate["claim_level"],
        "checked_hashes": checked_hashes,
        "recomputed_decision": recomputed["decision"],
        "certified_rank_lower_bounds": {"plus": 8, "minus_R2_gauge": 8},
        "eliminated_classes": expected_eliminated,
        "compatible_only": certificate["classification"]["survives_as_compatible_only"],
        "R2_bridge": "not_reached",
        "statistics": summaries,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    try:
        result = verify_certificate(args.root.resolve(), args.certificate.resolve())
    except (VerificationError, KeyError, TypeError, ValueError, OSError, json.JSONDecodeError) as error:
        print(json.dumps({"verified": False, "error": str(error)}, indent=2))
        return 1
    payload = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(payload)
    print(payload, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from score_norm5_chiral_phase import (  # noqa: E402
    gls_model,
    score_payload,
    sha256,
    synthetic_payload,
)


MANIFEST_PATH = ROOT / "experiments" / "p226_norm5_chiral_fixedp_production_20260829.json"
CONTRACT_PATH = ROOT / "predictions" / "norm5_chiral_phase_scorer_20260829.json"


def fixtures() -> tuple[dict, dict]:
    return json.loads(MANIFEST_PATH.read_text()), json.loads(CONTRACT_PATH.read_text())


def test_production_manifest_hash_is_frozen() -> None:
    _manifest, contract = fixtures()
    assert sha256(MANIFEST_PATH) == contract["production_manifest_sha256"]
    assert contract["production_manifest_sha256"] == (
        "9b381f7ecf651d482bc4cbb2a63d2217893a44feb2007f06302bc45bc32ade9f"
    )


def test_each_exact_synthetic_target_is_recovered_jointly() -> None:
    manifest, contract = fixtures()
    for target in ("H4", "H8", "H12"):
        payload = synthetic_payload(manifest, contract, target)
        score = score_payload(payload, manifest, contract)
        assert score["best_model"] == target
        assert score["models"][target]["chi_square"] < 1e-20
        assert score["models"][target]["degrees_of_freedom"] == 2
        assert score["reflection_conjugacy_null"]["status"] == "exact_configurationwise_null"
        assert score["evidence_accounting"].startswith("one joint 4D")


def test_offdiagonal_covariance_enters_the_gls_fit() -> None:
    mean = [0.8, -0.2, 1.1, 0.3]
    correlated = [
        [0.04, 0.012, 0.009, -0.006],
        [0.012, 0.03, 0.004, 0.008],
        [0.009, 0.004, 0.025, -0.007],
        [-0.006, 0.008, -0.007, 0.035],
    ]
    diagonal = [
        [correlated[i][i] if i == j else 0.0 for j in range(4)]
        for i in range(4)
    ]
    joint = gls_model(mean, correlated, complex(-527 / 625, -336 / 625))
    fake_independent = gls_model(mean, diagonal, complex(-527 / 625, -336 / 625))
    assert abs(joint["chi_square"] - fake_independent["chi_square"]) > 1.0
    assert joint["fitted_common_complex_normalization_re_im"] != (
        fake_independent["fitted_common_complex_normalization_re_im"]
    )


def test_run_metadata_tampering_is_rejected() -> None:
    manifest, contract = fixtures()
    payload = synthetic_payload(manifest, contract, "H8")
    payload["run"]["seed"] += 1
    try:
        score_payload(payload, manifest, contract)
    except ValueError as error:
        assert "seed differs" in str(error)
    else:
        raise AssertionError("tampered production metadata was accepted")

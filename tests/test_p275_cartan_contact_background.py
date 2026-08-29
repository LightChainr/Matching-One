from __future__ import annotations

import importlib.util
from fractions import Fraction
from pathlib import Path

import numpy as np
import yaml


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "p275_cartan_contact",
    ROOT / "scripts/analyze_p275_cartan_contact_background.py",
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def test_transition_oracle_is_exact() -> None:
    rows = MODULE.exact_transition_oracle()
    assert {row["transition"] for row in rows} == {"01", "12", "02_line_null"}
    assert all(row["pass"] for row in rows)


def test_frozen_modulus_weights_annihilate_training_profile() -> None:
    result = yaml.safe_load(
        (ROOT / "analysis/p275_cartan_background_holdout_20260829.yaml").read_text()
    )
    weights = np.asarray([
        complex(real, imag)
        for real, imag in result["secondary_negative_control"][
            "weights_order_i_2i_5i_over_2"
        ]
    ])
    background = np.asarray([
        0.33820875354274704 - 7.12221571533348e-06j,
        0.4323381206314479 - 0.00017807975668895294j,
        0.4633403652728997 - 0.0003041720626276482j,
    ])
    q4 = np.asarray([1.0, 2.75, 4.293436854374923])
    assert abs(weights @ background) < 1e-12
    assert abs(weights @ q4 - 1) < 1e-12


def test_n250_geometry_and_phases() -> None:
    manifest = yaml.safe_load(
        (ROOT / "analysis/p275_cartan_background_holdout_20260829.yaml").read_text()
    )
    for geometry in manifest["heldout_acquisition"]["moduli"].values():
        matrix = np.asarray(geometry["matrix"], dtype=int)
        assert round(np.linalg.det(matrix)) == 250
        phase = complex(
            float(Fraction(geometry["transport"]["real"])),
            float(Fraction(geometry["transport"]["imag"])),
        )
        assert abs(abs(phase) - 1) < 1e-12


def test_discovery_design_dimensions() -> None:
    prediction = MODULE._load_prediction(
        ROOT / "predictions/p275_atop_q4_field_identity_20260829.yaml"
    )
    assert MODULE.discovery_design("constant_by_modulus", prediction).shape == (18, 6)
    assert MODULE.discovery_design(
        "constant_by_modulus_plus_Q4_shape_tail", prediction
    ).shape == (18, 8)
    assert MODULE.discovery_design(
        "constant_by_modulus_plus_free_tail", prediction
    ).shape == (18, 12)

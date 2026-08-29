import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "p118", ROOT / "scripts" / "p118_cross_microscopic_invariants.py"
)
assert SPEC and SPEC.loader
MOD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)


def test_raw_single_modulus_has_no_nonconstant_gauge_invariant():
    assert [MOD.raw_charge(p, q) for p in range(-3, 4) for q in range(-3, 4) if MOD.raw_charge(p, q) == (0, 0)] == [(0, 0)]


def test_hex_child_cross_microscopic_double_ratio_is_one():
    a4, a6 = MOD.hex_child_forms(scale4=complex(7, -2), scale6=complex(-3, 5))
    for j in range(3):
        for k in range(j + 1, 3):
            assert abs(MOD.double_ratio(a4[j], a4[k], a6[j], a6[k]) - 1) < 1e-14
            assert abs(MOD.polynomial_null(a4[j], a4[k], a6[j], a6[k])) < 1e-10


def test_double_ratio_cancels_independent_microscopic_scales():
    a4, a6 = MOD.hex_child_forms()
    base = MOD.double_ratio(a4[0], a4[1], a6[0], a6[1])
    scaled = MOD.double_ratio(11j * a4[0], 11j * a4[1], -7 * a6[0], -7 * a6[1])
    assert abs(base - scaled) < 1e-14


def test_projective_jordan_coordinate_rescaling_and_shear_boundary():
    gram = (0.2, 0.7, 1.3)
    scaled = MOD.rescale_gram(*gram, local_scale=-3, top_scale=5)
    assert abs(MOD.jordan_k(*gram) - MOD.jordan_k(*scaled)) < 1e-14

    # At finite cutoff a Jordan shear is a real gauge ambiguity.
    sheared = MOD.shear_top(*gram, shear=2)
    assert abs(MOD.jordan_k(*gram) - MOD.jordan_k(*sheared)) > 1e-3

    # On the null-bottom LCFT boundary it disappears exactly.
    boundary = (0.0, 0.7, 1.3)
    for alpha in (-10.0, 0.0, 4.5):
        assert MOD.jordan_k(*MOD.shear_top(*boundary, shear=alpha)) == 0.0


def test_artifact_freezes_all_three_child_pairs():
    artifact = MOD.build_artifact()
    values = artifact["cross_microscopic_double_ratio"]["hex_degree2_child_pair_values"]
    assert set(values) == {"0:1", "0:2", "1:2"}
    assert all(v["exact"] == "1" for v in values.values())
    assert all(v["character_exponent_mod_3"] == 0 for v in values.values())
    assert all(v["floating_oracle_absolute_error"] < 1e-14 for v in values.values())

import importlib.util
import sys
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
SPEC = importlib.util.spec_from_file_location(
    "p180", ROOT / "scripts" / "p180_affine_clock_hankel.py"
)
assert SPEC and SPEC.loader
MOD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)


def test_cayley_hamilton_ordinary_q2_clock():
    a = [[Fraction(1), Fraction(0)], [Fraction(0), Fraction(1, 2)]]
    alpha, beta = MOD.square_clock_coefficients(Fraction(3, 2), Fraction(1, 2))
    assert (alpha, beta) == (Fraction(-1, 2), Fraction(3, 2))
    assert MOD.matmul2(a, a) == MOD.affine_matrix2(alpha, beta, a)


def test_cayley_hamilton_jordan_clock():
    a = [[Fraction(1), Fraction(1)], [Fraction(0), Fraction(1)]]
    alpha, beta = MOD.square_clock_coefficients(Fraction(2), Fraction(1))
    assert (alpha, beta) == (Fraction(-1), Fraction(2))
    assert MOD.matmul2(a, a) == MOD.affine_matrix2(alpha, beta, a)


def test_enriched_context_breaks_endpoint_affine_pencil():
    counterexample = MOD.exact_oracle()["enriched_sector_counterexample"]
    assert counterexample["endpoint_affine_prediction"] == "-5"
    assert counterexample["U_charged"] == "0"
    assert counterexample["defect"] == "5"


def test_affine_prediction_and_independent_covariance():
    parent = [1.0, 2.0]
    child = [3.0, 5.0]
    assert MOD.affine_predict(parent, child, -1.0, 2.0) == [5.0, 8.0]
    cov0 = [[1.0, 0.25], [0.25, 4.0]]
    cov1 = [[9.0, -0.5], [-0.5, 16.0]]
    assert MOD.affine_covariance(cov0, cov1, -1.0, 2.0) == [
        [37.0, -1.75],
        [-1.75, 68.0],
    ]


def test_state_coordinates():
    stat = {
        "mean_slope": 2.0,
        "P4_S": 3.0,
        "P4_D_prime": 5.0,
        "P4_D": 7.0,
        "P4_S_prime": 11.0,
    }
    state = MOD.state_from_statistics(stat, 4)
    assert state[0] == 12.0
    assert state[1] == 10.0
    assert state[2] == 4 ** (13 / 8) * 7
    assert state[3] == 4 ** (13 / 8) * 11 / 2

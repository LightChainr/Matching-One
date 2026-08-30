import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "src" / "c4_general_period_thermal_null_mc.cpp"


def test_exact_oracle_and_r8_injectivity(tmp_path):
    binary = tmp_path / "p155_mc"
    subprocess.run(
        ["g++", "-O1", "-std=c++17", str(SOURCE), "-o", str(binary)],
        cwd=ROOT,
        check=True,
    )
    exact = subprocess.run([str(binary), "--self-test"], check=True, text=True,
                           capture_output=True)
    assert "alpha*=3/64" in exact.stdout
    validation = subprocess.run(
        [str(binary), "--validate-only", "--radii", "2,4,8"],
        check=True,
        text=True,
        capture_output=True,
    )
    assert "injective Euclidean R=2,4,8" in validation.stdout

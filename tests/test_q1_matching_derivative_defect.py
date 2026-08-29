import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from q1_matching_derivative_defect import build_oracle  # noqa: E402


class Q1MatchingDerivativeDefectTests(unittest.TestCase):
    def test_planar_dual_defect_vanishes_configurationwise(self) -> None:
        control = build_oracle()["edge_FK_planar_dual_control"]
        self.assertTrue(control["configurationwise_zero"])
        self.assertEqual(control["first_log_Q_tangent_power_coefficients"], [0, 0, 0, 0])

    def test_nondual_pair_has_exact_first_tangent(self) -> None:
        obstruction = build_oracle()["edge_FK_nondual_obstruction"]
        self.assertFalse(obstruction["configurationwise_zero"])
        self.assertEqual(
            obstruction["first_log_Q_tangent_bernstein_coefficients"], [0, -3, -6, -2]
        )
        self.assertEqual(
            obstruction["first_log_Q_tangent_power_coefficients"], [0, -3, 0, 1]
        )

    def test_site_matching_tangent_is_not_a_formal_parity_proof(self) -> None:
        site = build_oracle()["site_matching_C4_to_K4"]
        self.assertEqual(
            site["first_log_Q_tangent_bernstein_coefficients"], [-1, 0, 2, 0, 1]
        )
        self.assertEqual(
            site["first_log_Q_tangent_power_coefficients"], [-1, 4, -4, 0, 2]
        )
        self.assertEqual(site["p_half"]["mean_defect"], "1/8")

    def test_committed_oracle_is_reproducible(self) -> None:
        committed = json.loads(
            (ROOT / "results" / "exact-q1-matching-derivative-defect" / "tiny_oracle.json").read_text()
        )
        self.assertEqual(committed, build_oracle())

    def test_cli_round_trip(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / "oracle.json"
            subprocess.run(
                [sys.executable, str(ROOT / "scripts" / "q1_matching_derivative_defect.py"), "--json", str(target)],
                check=True,
            )
            self.assertEqual(json.loads(target.read_text()), build_oracle())


if __name__ == "__main__":
    unittest.main()

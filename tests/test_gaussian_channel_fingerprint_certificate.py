from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path
import tempfile
import unittest

from scripts.gaussian_channel_fingerprint_certificate import (
    EXPONENT_EIGHTHS,
    NORM5_H12_OVER_H4,
    build_contract,
    fingerprints_at,
    h12_fingerprint,
    validate_contract,
)


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT / "analysis" / "gaussian_channel_fingerprint_contract.json"


class GaussianChannelFingerprintCertificateTests(unittest.TestCase):
    def test_checked_in_contract_closes_exactly(self) -> None:
        self.assertEqual(validate_contract(), build_contract())

    def test_four_channels_do_not_alias_at_either_norm(self) -> None:
        for norm in (2, 5):
            values = fingerprints_at(norm)
            self.assertEqual(len(set(values.values())), 4)

    def test_derivative_pairs_differ_by_q_cubed(self) -> None:
        for norm in (2, 5):
            values = fingerprints_at(norm)
            self.assertEqual(values["P4_Dprime"] / values["P4_S"], norm**3)
            self.assertEqual(values["P4_Sprime"] / values["P4_D"], norm**3)

    def test_norm5_h12_sign_and_eighth_power_are_exact(self) -> None:
        self.assertLess(NORM5_H12_OVER_H4, 0)
        for channel in ("P4_D", "P4_Sprime"):
            exponent = EXPONENT_EIGHTHS[channel]
            expected = Fraction(1679**8, 625**8 * 5**exponent)
            self.assertEqual(h12_fingerprint(5, exponent), expected)

    def test_contract_drift_fails_closed(self) -> None:
        frozen = build_contract()
        frozen["exponent_eighths"]["P4_D"] = 12
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "drift.json"
            path.write_text(json.dumps(frozen), encoding="utf-8")
            with self.assertRaisesRegex(AssertionError, "contract drifted"):
                validate_contract(path)


if __name__ == "__main__":
    unittest.main()

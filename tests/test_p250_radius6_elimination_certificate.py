from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from verify_p250_radius6_elimination_certificate import VerificationError, verify_certificate  # noqa: E402


CERTIFICATE = ROOT / "analysis/p250_radius6_level_s_elimination_certificate.json"


class P250Radius6EliminationCertificateTests(unittest.TestCase):
    def test_repository_certificate_recomputes_and_verifies(self) -> None:
        result = verify_certificate(ROOT, CERTIFICATE)
        self.assertTrue(result["verified"])
        self.assertEqual(result["certified_rank_lower_bounds"], {"plus": 8, "minus_R2_gauge": 8})
        self.assertEqual(result["R2_bridge"], "not_reached")

    def test_tampered_input_hash_fails_closed_before_claim(self) -> None:
        certificate = json.loads(CERTIFICATE.read_text())
        certificate["immutable_inputs"]["fresh6_batches"]["sha256"] = "0" * 64
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "certificate.json"
            path.write_text(json.dumps(certificate))
            with self.assertRaisesRegex(VerificationError, "hash mismatch"):
                verify_certificate(ROOT, path)

    def test_tampered_statistical_status_fails_closed(self) -> None:
        certificate = json.loads(CERTIFICATE.read_text())
        certificate["hand_certificates"]["plus"]["rank_le_7"]["status"] = "compatible_not_eliminated"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "certificate.json"
            path.write_text(json.dumps(certificate))
            with self.assertRaisesRegex(VerificationError, "certificate status was edited"):
                verify_certificate(ROOT, path)


if __name__ == "__main__":
    unittest.main()

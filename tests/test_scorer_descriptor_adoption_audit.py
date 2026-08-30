from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from scorer_descriptor_adoption_audit import audit, git_blob_sha  # noqa: E402


TYPED = """from wrapping_channels import ObservableDescriptor, map_observable\n"""
KERNEL = """def score():\n    return 1\n"""


class ScorerDescriptorAdoptionAuditTests(unittest.TestCase):
    def fixture(self, root: Path) -> dict:
        scripts = root / "scripts"
        scripts.mkdir()
        (scripts / "score_direct.py").write_text(TYPED, encoding="utf-8")
        (scripts / "score_wrapper.py").write_text(TYPED, encoding="utf-8")
        (scripts / "score_kernel.py").write_text(KERNEL, encoding="utf-8")
        (scripts / "score_generic.py").write_text(KERNEL, encoding="utf-8")
        (scripts / "score_migration.py").write_text(KERNEL, encoding="utf-8")
        (scripts / "score_unclassified.py").write_text(KERNEL, encoding="utf-8")
        (scripts / "scorer_descriptor_adoption_audit.py").write_text(
            KERNEL, encoding="utf-8"
        )
        return {
            "schema": "fixture",
            "corpus_glob": "scripts/*score*.py",
            "corpus_exclusions": ["scripts/scorer_descriptor_adoption_audit.py"],
            "typed_wrapper_kernels": {
                "scripts/score_wrapper.py": "scripts/score_kernel.py"
            },
            "direct_typed_standalone": ["scripts/score_direct.py"],
            "descriptor_not_applicable_generic_utilities": {
                "scripts/score_generic.py": "opaque generic score helper"
            },
            "channel_bearing_migration_required": {
                "scripts/score_migration.py": "fixture migration"
            },
            "boundary": "fixture boundary",
        }

    def test_complete_partition_and_git_blob_identity(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result = audit(root, self.fixture(root))
            self.assertEqual(
                result["counts"],
                {
                    "total": 6,
                    "direct_typed_entrypoint": 2,
                    "covered_frozen_kernel": 1,
                    "descriptor_not_applicable_generic_utility": 1,
                    "channel_bearing_migration_required": 1,
                    "outside_registered_typed_path": 1,
                },
            )
            rows = {row["path"]: row for row in result["rows"]}
            self.assertEqual(
                rows["scripts/score_kernel.py"]["git_blob_sha"],
                git_blob_sha(KERNEL.encode()),
            )

    def test_unregistered_direct_import_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = self.fixture(root)
            manifest["direct_typed_standalone"] = []
            with self.assertRaisesRegex(ValueError, "do not equal"):
                audit(root, manifest)

    def test_stale_corpus_exclusion_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = self.fixture(root)
            manifest["corpus_exclusions"] = ["scripts/missing_score_audit.py"]
            with self.assertRaisesRegex(ValueError, "exclusions do not match"):
                audit(root, manifest)

    def test_missing_wrapped_kernel_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = self.fixture(root)
            manifest["typed_wrapper_kernels"]["scripts/score_wrapper.py"] = "scripts/score_missing.py"
            with self.assertRaisesRegex(ValueError, "kernels missing"):
                audit(root, manifest)

    def test_not_applicable_cannot_overlap_typed_path(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = self.fixture(root)
            manifest["descriptor_not_applicable_generic_utilities"] = {
                "scripts/score_direct.py": "invalid overlap"
            }
            with self.assertRaisesRegex(ValueError, "overlap typed paths"):
                audit(root, manifest)

    def test_migration_cannot_overlap_resolved_class(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = self.fixture(root)
            manifest["channel_bearing_migration_required"] = {
                "scripts/score_generic.py": "invalid overlap"
            }
            with self.assertRaisesRegex(ValueError, "overlap resolved classes"):
                audit(root, manifest)

    def test_checked_result_has_current_partition(self) -> None:
        manifest = json.loads(
            (ROOT / "analysis/scorer_descriptor_adoption_manifest.json").read_text(
                encoding="utf-8"
            )
        )
        result = json.loads(
            (ROOT / "results/scorer-descriptor-adoption/latest.json").read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(result, audit(ROOT, manifest))
        self.assertEqual(result["counts"]["total"], 39)
        self.assertEqual(result["counts"]["direct_typed_entrypoint"], 8)
        self.assertEqual(result["counts"]["covered_frozen_kernel"], 7)
        self.assertEqual(
            result["counts"]["descriptor_not_applicable_generic_utility"], 1
        )
        self.assertEqual(result["counts"]["channel_bearing_migration_required"], 23)
        self.assertEqual(result["counts"]["outside_registered_typed_path"], 0)
        self.assertEqual(len(result["rows"]), 39)
        statuses = {row["path"]: row["status"] for row in result["rows"]}
        for path in (
            "scripts/score_axis_pair_annihilator.py",
            "scripts/score_axis_pair_annihilator_stable.py",
            "scripts/score_c4_tangent_orthogonal_holdout.py",
            "scripts/score_intrinsic_quantile_center_n145_n290.py",
            "scripts/score_issue43_full_curve.py",
            "scripts/score_issue43_full_curve_locked.py",
            "scripts/score_issue43_secondary.py",
            "scripts/score_norm4_production.py",
            "scripts/score_norm4_thermal_jet.py",
            "scripts/score_norm5_thermal_jet.py",
            "scripts/score_p159_pell_hex_filter.py",
            "scripts/score_p231_vacuum_kdv_sector.py",
            "scripts/score_p48_new_geometry_channels.py",
            "scripts/score_p49_fullcurve_doubling.py",
            "scripts/score_p50_fullcurve_n290.py",
            "scripts/score_p50_sprime_n290.py",
            "scripts/score_rank_gap_boundary_targets.py",
            "scripts/score_threshold_rank_root_doubling.py",
            "scripts/threshold_score_modes.py",
            "scripts/score_v14_fixedp_scalar_projector.py",
            "scripts/score_v14_scalar_root_projector.py",
        ):
            self.assertEqual(
                statuses[path], "channel_bearing_migration_required"
            )
        self.assertEqual(
            statuses["scripts/score_c4_self_matching_n26.py"],
            "covered_frozen_kernel",
        )
        self.assertEqual(
            statuses["scripts/score_c4_self_matching_n26_typed.py"],
            "direct_typed_entrypoint",
        )
        self.assertEqual(
            statuses["scripts/score_issue50_n290.py"],
            "covered_frozen_kernel",
        )
        self.assertEqual(
            statuses["scripts/score_issue50_n290_typed.py"],
            "direct_typed_entrypoint",
        )
        self.assertEqual(
            statuses["scripts/score_matching_odd_synthesis.py"],
            "covered_frozen_kernel",
        )
        self.assertEqual(
            statuses["scripts/score_matching_odd_synthesis_typed.py"],
            "direct_typed_entrypoint",
        )
        self.assertEqual(
            statuses["scripts/score_p48_sprime_frozen.py"],
            "covered_frozen_kernel",
        )
        self.assertEqual(
            statuses["scripts/score_p48_sprime_frozen_typed.py"],
            "direct_typed_entrypoint",
        )


if __name__ == "__main__":
    unittest.main()

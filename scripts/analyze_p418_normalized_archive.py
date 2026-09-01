#!/usr/bin/env python3
"""Re-score the pinned P418 archives in per-sample, not per-batch-sum units.

No Monte Carlo configurations are generated.  The original reader, finite
geometry, CRT masks, covariance whitening, NNLS, and parametric bootstrap are
loaded from immutable Git objects.  Only the archive exposure normalization
changes; historical code and results are not modified.
"""

from __future__ import annotations

import argparse
import ast
import csv
from datetime import datetime, timezone
import hashlib
import importlib
import io
import json
from pathlib import Path
import platform
import subprocess
import sys
import tempfile
import time

import numpy as np
import scipy
from threadpoolctl import threadpool_info, threadpool_limits


ROOT = Path(__file__).resolve().parents[1]
BASE_COMMIT = "8704eee790403e14e5ad75d3465ee1496eaa9c0e"
FAMILY_COMMIT = "588ca452dedd47213a424d79fc119ad67f8f77df"
OLD_SCORE_PATH = "results/huawei-20260830/P418-radius-resolved-elimination/score.json"
BOOTSTRAP = 250
SEED = 40610120260830
RADIUS_NAMES = ("radius4", "radius5", "radius6")
SCHEMA = "matching-one/p418-normalized-archive/v1"


def git(*args: str) -> bytes:
    return subprocess.check_output(["git", "-C", str(ROOT), *args])


def git_bytes(commit: str, path: str) -> bytes:
    return git("show", f"{commit}:{path}")


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def materialize_modules(directory: Path) -> list[dict]:
    """Extract only the recursively imported pinned script modules, never run MC."""
    paths = git("ls-tree", "-r", "--name-only", BASE_COMMIT, "scripts").decode().splitlines()
    available = {
        Path(path).stem: path for path in paths
        if path.endswith(".py") and Path(path).parent == Path("scripts")
    }
    family_module = "score_p418_radius_resolved_elimination"
    available[family_module] = f"scripts/{family_module}.py"
    pending = ["score_p406_spatial_fourier_cone", "score_p418_crt_degauging", family_module]
    seen = set()
    provenance = []
    while pending:
        name = pending.pop()
        if name in seen:
            continue
        seen.add(name)
        commit = FAMILY_COMMIT if name == family_module else BASE_COMMIT
        path = available[name]
        content = git_bytes(commit, path)
        (directory / f"{name}.py").write_bytes(content)
        provenance.append({"module": name, "path": path, "commit": commit, "sha256": digest(content)})
        for node in ast.walk(ast.parse(content, filename=path)):
            candidates = []
            if isinstance(node, ast.Import):
                candidates = [alias.name.split(".")[0] for alias in node.names]
            elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
                candidates = [node.module.split(".")[0]]
            pending.extend(candidate for candidate in candidates if candidate in available and candidate not in seen)
    return sorted(provenance, key=lambda row: row["path"])


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def fit_family_with_spectra(p406, family, blocks, design_builder, seed):
    """Call the unchanged historical family fit; retain its four point spectra."""
    original_fit = p406.fit_nonnegative
    point_fits = []

    def recording_fit(matrix, vector):
        weights, statistic = original_fit(matrix, vector)
        if len(point_fits) < 4:
            point_fits.append({"weights": weights.tolist(), "total_mass": float(weights.sum())})
        return weights, statistic

    p406.fit_nonnegative = recording_fit
    try:
        result = family.fit_family(blocks, design_builder, BOOTSTRAP, seed)
    finally:
        p406.fit_nonnegative = original_fit
    result["shared_spectrum"].update(point_fits[0])
    for row, spectrum in zip(result["per_radius_separate_spectra"], point_fits[1:]):
        row.update(spectrum)
    return result


def report(result: dict) -> str:
    lines = [
        "# P418: normalized production-archive reanalysis",
        "",
        "This is a per-sample-unit correction of the existing radius4/radius5/radius6 production archives, not a new simulation.",
        "",
        (
            "All four joint raw/masked spectra are compatible at the original alpha=0.01. The historical radius-flow interpretation loses the support of this score. Compatibility does not identify a unique physical spectrum."
            if not result["decision"]["shared_raw_rejections"] and not result["decision"]["shared_mask_rejections"] else
            "The corrected channel decisions are reported below without assuming that normalization must remove every tension."
        ),
        "",
        "## Corrected shared-spectrum result",
        "",
        "| channel | raw shared d2 (p) | masked shared d2 (p) | inherited masked sharing penalty (p; unreliable) |",
        "|---|---:|---:|---:|",
    ]
    for key, channel in result["channels"].items():
        raw = channel["raw_family"]["shared_spectrum"]
        masked = channel["masked_family"]["shared_spectrum"]
        penalty = channel["masked_family"]["cross_radius_sharing_penalty"]
        lines.append(
            f"| {key} | {raw['distance_squared']:.6g} ({raw['bootstrap_p']:.6g}) | "
            f"{masked['distance_squared']:.6g} ({masked['bootstrap_p']:.6g}) | "
            f"{penalty['distance_squared']:.6g} ({penalty['bootstrap_p_under_shared_spectrum']:.6g}) |"
        )
    lines.extend([
        "",
        f"At the unchanged alpha=0.01, corrected shared raw rejections: `{result['decision']['shared_raw_rejections']}`; "
        f"corrected shared masked rejections: `{result['decision']['shared_mask_rejections']}`. "
        "Inherited sharing decisions are not scoreable because the radius5 fit is numerically unreliable.",
        "",
        "## What changed from the historical result",
        "",
        "| channel | historical masked shared d2 | corrected masked shared d2 | historical masked sharing penalty | corrected penalty |",
        "|---|---:|---:|---:|---:|",
    ])
    for key, channel in result["channels"].items():
        old = result["historical_comparison"][key]
        masked = channel["masked_family"]
        lines.append(
            f"| {key} | {old['masked_shared_distance_squared']:.6g} | "
            f"{masked['shared_spectrum']['distance_squared']:.6g} | "
            f"{old['masked_sharing_penalty']:.6g} | {masked['cross_radius_sharing_penalty']['distance_squared']:.6g} |"
        )
    lines.extend([
        "",
        "## Separate-radius compatibility",
        "",
        "The radius5 rows and sharing penalties below are preserved outputs of the inherited solver, not reliable mechanism statistics: the saved-point diagnostic detects amplification of numerical null directions. Joint fits do not have this defect.",
        "",
        "| channel | radius | raw d2 (p) | masked d2 (p) | resolved covariance modes |",
        "|---|---|---:|---:|---:|",
    ])
    for key, channel in result["channels"].items():
        for raw, masked in zip(channel["raw_family"]["per_radius_separate_spectra"], channel["masked_family"]["per_radius_separate_spectra"]):
            lines.append(
                f"| {key} | {raw['radius']} | {raw['distance_squared']:.6g} ({raw['bootstrap_p']:.6g}) | "
                f"{masked['distance_squared']:.6g} ({masked['bootstrap_p']:.6g}) | {masked['resolved_modes']} |"
            )
    lines.extend([
        "",
        "## Method and finite scope",
        "",
        "The CSV coordinates are sums over each batch. Every radius4 batch contains 200 samples, whereas every radius5/radius6 batch contains 3000. The original reader used these sums as observations without rescaling the common Fourier design. This wrapper divides each row by its own `samples` before calling the unchanged historical whitening and family fitter.",
        "",
        "For block exposure n, replacing batch sums by per-sample means sends mean to mean/n and covariance to covariance/n². The whitened response is unchanged (up to the eigensystem basis), but the correctly whitened design acquires the factor n. An independent spectrum for each radius can absorb this factor; one common spectrum cannot absorb three inconsistent units. Thus the old common-spectrum rejection cannot itself be interpreted as physical radius flow.",
        "",
        "The method keeps the original covariance-of-the-batch-mean formula, relative eigenvalue cutoff 1e-10, nonnegative 101-frequency Fourier cone, exact CRT mask, 250 Gaussian parametric-bootstrap draws, original channel seeds and alpha=0.01. The sharing penalty is calibrated under the common fitted center for both nested fits. `score.json` retains every original family-fit output plus the point spectra; `normalized-inputs.json` retains the normalized means and full within-channel/radius covariance matrices.",
        "",
        "Different radius streams retain the original block-diagonal covariance convention. Hands/charges in one radius share configurations; their separate p values are not independent pieces of evidence and are not combined. The bootstrap draws are uncertainty calculations on archived means, not new Monte Carlo lattice configurations. Basis/sign choices in covariance eigendecomposition and NNLS nonuniqueness can change finite-bootstrap details without changing the statistical convention.",
        "",
        "Non-rejection is compatibility of this finite spatial-spectrum model, not a unique recovered spectrum, physical state count, local field identification, or proof of a continuum mechanism. A saturated radius6 fit is not affirmative model identification. If corrected tension remains, it belongs to this normalized finite observation contract; no new mechanism is assigned automatically.",
        "",
        "The radius4 block observes zero displacement, while radius5 and radius6 do not. Adding a constant to every spectral weight changes only zero displacement, so positivity alone cannot restrict an isolated nonzero-displacement shell beyond its signed Fourier span. This argument does not trivialize the common fit, whose total spectral mass is constrained by radius4's zero-lag row. Equal common raw/masked distances here are not a proof that their cones are always equal.",
        "",
        "A small saved-point diagnostic (solver-note.json; no NNLS or bootstrap rerun) finds all eight joint fits consistent with their reconstructed residuals within 7e-13, with the resolved-rank-69 least-squares floor within 7e-11 and scaled KKT violations below 1e-14. The radius5 fits instead amplify floating-point null directions: their saved-weight residuals fall about 10–24 below the exact/resolved-rank-20 Fourier-span floor and scaled negative-gradient violations reach 0.066. The missing zero-lag row permits arbitrarily large uniform null mass in that isolated shell. Those radius5 distances and the derived sharing penalties are numerically unreliable, not physical changes caused by rescaling. The robust result is the disappearance of the large joint masked penalty; the inherited sharing-penalty p values are not used for mechanism inference.",
        "",
        "The exact CRT/root-translation certificates, correctly normalized P250 Hankel/radius5/radius6 scorers and paired-anchor pilot's own statistics are not modified. Historical sum-unit outputs remain available; they are included for comparison, not silently overwritten.",
        "",
        "## Provenance and reproduction",
        "",
        f"- Data, geometry and mask source: `{BASE_COMMIT}`.",
        f"- Radius-family fitter and historical comparison: `{FAMILY_COMMIT}`.",
        "- Input byte SHA256 values must match the historical score's frozen inputs. No replacement archive is accepted.",
        "- Exact imported module commits and SHA256 values, Python/NumPy/SciPy/BLAS information, elapsed time and command are saved in `score.json` and `manifest.json`.",
        "",
        "```bash",
        "/Users/lc/python-envs/research-py311/bin/python scripts/analyze_p418_normalized_archive.py",
        "/Users/lc/python-envs/research-py311/bin/python scripts/analyze_p418_normalized_archive.py --diagnose-existing",
        "```",
        "",
        f"Elapsed wall time: {result['execution']['elapsed_seconds']:.3f} seconds.",
        "",
    ])
    return "\n".join(lines)


def diagnose_saved_point_fits(output: Path) -> None:
    """Read saved point weights: compute objectives/KKT and LS bounds, no NNLS/bootstrap."""
    started = time.perf_counter()
    score_path = output / "score.json"
    score = json.loads(score_path.read_text())
    inputs = json.loads((output / "normalized-inputs.json").read_text())
    diagnostics = {
        "schema": SCHEMA + "/saved-point-solver-diagnostic",
        "operation": "reconstruct saved whitened designs; evaluate saved weights and one unconstrained least-squares bound per fit; no NNLS or bootstrap rerun",
        "input_score_sha256_before_annotation": digest(score_path.read_bytes()),
        "wrapper_sha256": digest(Path(__file__).read_bytes()),
        "rank_contract": "LS floor uses rcond=1e-10, the exact/resolved Fourier rank; spurious floating-point null directions with huge coefficients are not additional physical modes",
        "fits": {},
    }
    for key, channel in score["channels"].items():
        hand, charge = key.split("_")
        masks = score["exact_CRT_section_and_masks"]["hands"][hand]["masks"][charge]["values"]
        channel_rows = {}
        for kind in ("raw_family", "masked_family"):
            matrices, vectors = [], []
            for block in inputs["blocks"][key]:
                rows = []
                for a, b in block["coordinate_order"]:
                    residue = (a - 10 * b) % 101
                    phase = 2.0 * np.pi * np.arange(101) * residue / 101
                    real, imag = np.cos(phase), -np.sin(phase)
                    if kind == "masked_family":
                        mask = masks[residue]
                        real, imag = mask["real"] * real - mask["imag"] * imag, mask["imag"] * real + mask["real"] * imag
                    rows.extend((real, imag))
                design = np.asarray(rows)
                covariance = np.asarray(block["covariance_of_mean"])
                eigenvalues, eigenvectors = np.linalg.eigh((covariance + covariance.T) / 2)
                keep = eigenvalues > max(float(eigenvalues[-1]) * 1e-10, 0.0)
                transform = eigenvectors[:, keep].T / np.sqrt(eigenvalues[keep])[:, None]
                matrices.append(transform @ design)
                vectors.append(transform @ np.asarray(block["mean"]))
            for scope, matrix, vector, saved in (
                ("joint", np.vstack(matrices), np.concatenate(vectors), channel[kind]["shared_spectrum"]),
                ("radius5", matrices[1], vectors[1], channel[kind]["per_radius_separate_spectra"][1]),
            ):
                weights = np.asarray(saved["weights"])
                residual = matrix @ weights - vector
                objective = float(residual @ residual)
                gradient = matrix.T @ residual
                solution, _, rank, singular = np.linalg.lstsq(matrix, vector, rcond=1e-10)
                ls_residual = matrix @ solution - vector
                lower_bound = float(ls_residual @ ls_residual)
                scale = max(float(singular[0]) * float(np.linalg.norm(residual)), 1.0)
                active = weights > max(float(weights.max()) * 1e-12, 1e-15)
                channel_rows[f"{kind}:{scope}"] = {
                    "rows": len(vector), "columns": matrix.shape[1], "least_squares_rank_rcond_1e10": int(rank),
                    "reported_nnls_distance_squared": saved["distance_squared"],
                    "reconstructed_saved_weight_distance_squared": objective,
                    "reported_minus_reconstructed_distance_squared": saved["distance_squared"] - objective,
                    "unconstrained_ls_lower_bound": lower_bound,
                    "saved_weight_objective_minus_ls_lower_bound": objective - lower_bound,
                    "minimum_saved_weight": float(weights.min()),
                    "minimum_KKT_gradient": float(gradient.min()),
                    "negative_gradient_violation_scaled": max(0.0, -float(gradient.min())) / scale,
                    "active_gradient_max_abs_scaled": float(np.max(np.abs(gradient[active]))) / scale if np.any(active) else 0.0,
                    "maximum_complementarity_abs": float(np.max(np.abs(weights * gradient))),
                    "KKT_scale": scale,
                }
        diagnostics["fits"][key] = channel_rows
    diagnostics["elapsed_seconds"] = time.perf_counter() - started
    diagnostics["finished_utc"] = datetime.now(timezone.utc).isoformat()
    write_json(output / "solver-note.json", diagnostics)
    score["numerical_limitations"] = {
        "radius5_scale_invariance": "In exact arithmetic per-shell cone distance is unchanged by constant exposure rescaling. Historical versus corrected radius5 NNLS point distances drift; this is not a physical effect of units.",
        "diagnostic": "solver-note.json evaluates saved point weights and LS bounds without refitting NNLS or rerunning bootstrap. Interpret returned residual norms and sharing penalties in light of that diagnostic.",
        "joint_saved_weight_fit": "stable: all raw/masked joint fits attain the rank-69 LS floor within 7e-11; scaled KKT violations below 1e-14",
        "radius5_and_derived_sharing_penalties": "numerically_unreliable: large null-direction weights amplify floating-point Fourier design errors; retained for provenance, not mechanism inference",
        "claim_scope": "The large common-mask rejection is not reproduced and the corrected joint point fits are numerically supported; this is not an exact mechanism certification or unique-spectrum identification.",
    }
    write_json(score_path, score)
    manifest_path = output / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["post_run_changes"] = {
        "scope": "reporting text, numerical limitation annotation, and saved-point residual/LS/KKT diagnostic only; no archive-score rerun",
        "wrapper_sha256_after_reporting_update": digest(Path(__file__).read_bytes()),
        "solver_diagnostic_command": [sys.executable, str(Path(__file__).relative_to(ROOT)), "--diagnose-existing"],
    }
    manifest["outputs"] = [{"path": name, "sha256": digest((output / name).read_bytes())} for name in ("score.json", "normalized-inputs.json", "REPORT.md", "solver-note.json")]
    write_json(manifest_path, manifest)
    print(json.dumps(diagnostics, indent=2), flush=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results/p418-normalized-archive")
    parser.add_argument("--diagnose-existing", action="store_true", help="read saved point fits and compute residual/LS/KKT diagnostics only; no scorer run")
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if args.diagnose_existing:
        with threadpool_limits(limits=1):
            diagnose_saved_point_fits(output)
        return 0
    started = time.perf_counter()
    started_at = datetime.now(timezone.utc).isoformat()
    old_bytes = git_bytes(FAMILY_COMMIT, OLD_SCORE_PATH)
    historical = json.loads(old_bytes)
    execution = {
        "started_utc": started_at,
        "command": [sys.executable, *sys.argv],
        "python": sys.version,
        "executable": sys.executable,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "git_head_at_run": git("rev-parse", "HEAD").decode().strip(),
        "wrapper_sha256": digest(Path(__file__).read_bytes()),
        "blas_threads_limit": 1,
    }
    result = {
        "schema": SCHEMA,
        "status": "normalized_archived_production_analysis",
        "new_monte_carlo": False,
        "issues": [418, 406, 250],
        "only_estimand_input_change": "batch_sum / row.samples, before unchanged covariance whitening and model fit",
        "source_commits": {"archive_geometry_mask": BASE_COMMIT, "radius_family_fitter": FAMILY_COMMIT},
        "bootstrap": {"replicates": BOOTSTRAP, "base_seed": SEED, "channel_seed": "base_seed + 1000*hand_index + charge", "decision_alpha": 0.01},
        "historical_score": {"commit": FAMILY_COMMIT, "path": OLD_SCORE_PATH, "sha256": digest(old_bytes)},
        "execution": execution,
        "inputs": [],
        "channels": {},
        "historical_comparison": {},
    }
    normalized_inputs = {"schema": SCHEMA + "/normalized-inputs", "covariance": "covariance of the batch-mean estimator", "blocks": {}}
    with tempfile.TemporaryDirectory(prefix="matching-p418-normalized-") as temporary:
        isolated = Path(temporary)
        module_directory = isolated / "scripts"
        module_directory.mkdir()
        result["pinned_module_sources"] = materialize_modules(module_directory)
        sys.path.insert(0, str(module_directory))
        p406 = importlib.import_module("score_p406_spatial_fourier_cone")
        p418 = importlib.import_module("score_p418_crt_degauging")
        family = importlib.import_module("score_p418_radius_resolved_elimination")
        original_read = p406.read_block

        def normalized_read(path, expected_hash):
            blocks = original_read(path, expected_hash)
            with path.open(newline="") as handle:
                rows = list(csv.DictReader(handle))
            exposures = np.asarray([int(row["samples"]) for row in rows], dtype=float)
            if np.any(exposures <= 0) or len(set(exposures)) != 1:
                raise ValueError("This correction retains the historical equal-exposure within-radius contract")
            for block in blocks.values():
                block["values"] = block["values"] / exposures[:, None]
            return blocks

        p406.read_block = normalized_read
        archives = []
        for radius, source in zip(RADIUS_NAMES, historical["inputs"]):
            content = git_bytes(BASE_COMMIT, source["path"])
            if digest(content) != source["sha256"]:
                raise ValueError(f"Pinned archive no longer matches the historical bytes: {source['path']}")
            rows = list(csv.DictReader(io.StringIO(content.decode())))
            exposures = sorted({int(row["samples"]) for row in rows})
            blob = git("rev-parse", f"{BASE_COMMIT}:{source['path']}").decode().strip()
            source_path = isolated / f"{radius}.csv"
            source_path.write_bytes(content)
            archives.append(p406.read_block(source_path, source["sha256"]))
            result["inputs"].append({**source, "commit": BASE_COMMIT, "git_blob": blob, "radius": radius, "batches": len(rows), "samples_per_batch": exposures, "total_samples": sum(int(row["samples"]) for row in rows), "stored_unit": "batch_sum", "analysis_unit": "per_sample_mean"})
        with threadpool_limits(limits=1):
            execution["threadpools"] = threadpool_info()
            exact = p418.exact_section_and_masks()
            result["exact_CRT_section_and_masks"] = exact
            for hand_index, hand in enumerate(p406.HANDS):
                for charge in p406.CHARGES:
                    key = f"{hand}_r{charge}"
                    blocks = [archive[(hand, charge)] for archive in archives]
                    channel_seed = SEED + 1000 * hand_index + charge
                    normalized_inputs["blocks"][key] = []
                    for radius, block in zip(RADIUS_NAMES, blocks):
                        values = block["values"]
                        mean = values.mean(axis=0)
                        centered = values - mean
                        covariance = centered.T @ centered / (len(values) * (len(values) - 1))
                        normalized_inputs["blocks"][key].append({"radius": radius, "coordinate_order": [list(value) for value in block["coordinates"]], "component_order": "re,im for each coordinate", "mean": mean.tolist(), "covariance_of_mean": covariance.tolist()})
                    channel = {}
                    for name, builder in (
                        ("raw_family", lambda _index, coordinates: p406.design(coordinates)),
                        ("masked_family", lambda _index, coordinates: p418.masked_design(exact, hand, charge, coordinates)),
                    ):
                        step_start = time.perf_counter()
                        print(f"{key} {name}: fitting 250-draw archive bootstrap", flush=True)
                        channel[name] = fit_family_with_spectra(p406, family, blocks, builder, channel_seed)
                        channel[name]["elapsed_seconds"] = time.perf_counter() - step_start
                        print(f"{key} {name}: d2={channel[name]['shared_spectrum']['distance_squared']:.9g}, p={channel[name]['shared_spectrum']['bootstrap_p']:.9g}", flush=True)
                    result["channels"][key] = channel
                    old = historical["channels"][key]
                    result["historical_comparison"][key] = {
                        "raw_shared_distance_squared": old["raw_family"]["shared_spectrum"]["distance_squared"],
                        "masked_shared_distance_squared": old["masked_family"]["shared_spectrum"]["distance_squared"],
                        "raw_sharing_penalty": old["raw_family"]["cross_radius_sharing_penalty"]["distance_squared"],
                        "masked_sharing_penalty": old["masked_family"]["cross_radius_sharing_penalty"]["distance_squared"],
                    }
    result["decision"] = {
        "shared_raw_rejections": [key for key, row in result["channels"].items() if row["raw_family"]["shared_spectrum"]["bootstrap_p"] < .01],
        "shared_mask_rejections": [key for key, row in result["channels"].items() if row["masked_family"]["shared_spectrum"]["bootstrap_p"] < .01],
        "sharing_penalty_scoreability": "not_scoreable_radius5_numerical_null_mode_leakage",
        "single_radius5_scoreability": "not_scoreable_numerical_null_mode_leakage",
    }
    result["raw_solver_decision"] = {
        "masked_sharing_penalty_rejections": [key for key, row in result["channels"].items() if row["masked_family"]["cross_radius_sharing_penalty"]["bootstrap_p_under_shared_spectrum"] < .01],
        "single_radius_mask_rejections": [f"{key}:{radius['radius']}" for key, row in result["channels"].items() for radius in row["masked_family"]["per_radius_separate_spectra"] if radius["bootstrap_p"] < .01],
        "status": "uninterpretable_trace_only_for_radius5_and_derived_sharing_penalties",
    }
    execution["elapsed_seconds"] = time.perf_counter() - started
    execution["finished_utc"] = datetime.now(timezone.utc).isoformat()
    write_json(output / "score.json", result)
    write_json(output / "normalized-inputs.json", normalized_inputs)
    (output / "REPORT.md").write_text(report(result), encoding="utf-8")
    write_json(output / "manifest.json", {
        "schema": SCHEMA + "/manifest", "inputs": result["inputs"],
        "pinned_module_sources": result["pinned_module_sources"], "execution": execution,
        "outputs": [{"path": name, "sha256": digest((output / name).read_bytes())} for name in ("score.json", "normalized-inputs.json", "REPORT.md")],
    })
    print(json.dumps({"decision": result["decision"], "elapsed_seconds": execution["elapsed_seconds"]}, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

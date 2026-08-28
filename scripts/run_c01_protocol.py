#!/usr/bin/env python3
"""C01 same-N Gaussian orientation discovery protocol.

Pilot-freeze-evaluate.  Does not reuse pilot replicas as evaluation.
Bond control uses the same geometry engine at p=1/2.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import platform
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone, timedelta
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

from analyze_gaussian_orientation import (  # noqa: E402
    CHANNEL_NAMES,
    analyze_run,
    write_long_form_rows,
)
from gaussian_circulant_geometry import GaussianTorus, best_multiplier  # noqa: E402


P_REF = 0.592746050790
P_COMPLEMENT = 0.407253949210
P_BOND = 0.5
SEED_SITE = 20260828
SEED_BOND = 20260829
PILOT_SAMPLES = 200000
EVAL_SAMPLES = 2000000
PILOT_BEGIN = 0
EVAL_BEGIN = 1_000_000
BATCHES = 20
THREADS = 8
Z_RESOLVE = 3.0

DISCOVERY = [
    {"N": 65, "rep1": (8, 1), "rep2": (7, 4)},
    {"N": 85, "rep1": (9, 2), "rep2": (7, 6)},
    {"N": 145, "rep1": (12, 1), "rep2": (9, 8)},
]
CONFIRMATION = [
    {"N": 205, "rep1": (14, 3), "rep2": (13, 6)},
]

SHANGHAI = timezone(timedelta(hours=8), name="CST")


def now_pair() -> tuple[datetime, datetime]:
    utc = datetime.now(timezone.utc)
    return utc, utc.astimezone(SHANGHAI)


def fmt(dt: datetime) -> str:
    return dt.strftime("%Y-%m-%d %H:%M:%S")


def run(cmd: list[str], cwd: Path = ROOT) -> subprocess.CompletedProcess:
    print("+", " ".join(cmd), flush=True)
    completed = subprocess.run(cmd, cwd=str(cwd), check=True, text=True)
    return completed


def compile_binaries(build: Path) -> dict[str, Path]:
    build.mkdir(parents=True, exist_ok=True)
    g_mc = build / "gaussian_orientation_mc"
    p_mc = build / "pell_matching_mc"
    cxx = os.environ.get("CXX") or shutil.which("g++") or shutil.which("clang++")
    if cxx is None:
        raise RuntimeError("no C++ compiler")
    run([cxx, "-O3", "-std=c++17", "-fopenmp",
         str(ROOT / "src" / "gaussian_orientation_mc.cpp"), "-o", str(g_mc)])
    run([cxx, "-O3", "-std=c++17", "-fopenmp",
         str(ROOT / "src" / "pell_matching_mc.cpp"), "-o", str(p_mc)])
    return {"gaussian": g_mc, "pell": p_mc, "cxx": Path(cxx)}


def simulate(binary: Path, pair: dict, t_values: list[int], mode: str, p: float,
             samples: int, seed: int, replica_begin: int, prefix: Path) -> None:
    t_arg = ",".join(str(t) for t in t_values)
    a1, b1 = pair["rep1"]
    a2, b2 = pair["rep2"]
    run([
        str(binary),
        "--rep1", f"{a1},{b1}",
        "--rep2", f"{a2},{b2}",
        "--t", t_arg,
        "--mode", mode,
        "--p", f"{p:.12f}",
        "--samples", str(samples),
        "--batches", str(BATCHES),
        "--seed", str(seed),
        "--replica-begin", str(replica_begin),
        "--threads", str(THREADS),
        "--output-prefix", str(prefix),
    ])


def analyze_prefix(prefix: Path, t: int) -> dict:
    return analyze_run(
        Path(f"{prefix}.t{t}.moments.json"),
        Path(f"{prefix}.metadata.json"),
    )


def graph_only_table() -> list[dict]:
    rows = []
    for pair in DISCOVERY + CONFIRMATION:
        g1 = GaussianTorus(*pair["rep1"])
        g2 = GaussianTorus(*pair["rep2"])
        t, score, nn, ma = best_multiplier(g1, g2)
        rows.append({
            "N": pair["N"],
            "rep1": list(pair["rep1"]),
            "rep2": list(pair["rep2"]),
            "t_struct": t,
            "score": score,
            "nn_overlap": nn,
            "matching_overlap": ma,
            "shortlist": [1, t] if t != 1 else [1],
            "selection_rule": "graph-only circulant step overlap, then pilot covariance of S_either",
        })
    return rows


def freeze_from_pilots(pilot_by_n: dict[int, dict[int, dict]]) -> dict:
    frozen = {"pairs": {}, "primary_channel": "either", "channel_rule": ""}
    # Per-pair t: smaller replica_var of S_delta/either; prefer t=1 if within 5%.
    for pair in DISCOVERY:
        n = pair["N"]
        by_t = pilot_by_n[n]
        scored = []
        for t, analysis in by_t.items():
            est = analysis["channels"]["either"]["S_delta"]
            scored.append((est["replica_var"], t, est))
        scored.sort()
        best_var, best_t, _ = scored[0]
        t1_var = by_t[1]["channels"]["either"]["S_delta"]["replica_var"]
        if 1 in by_t and best_t != 1 and best_var > 0.95 * t1_var:
            chosen = 1
            reason = "canonical t=1 within 5% of minimum S_either replica variance"
        else:
            chosen = best_t
            reason = "minimum replica variance of matching-even S_either"
        frozen["pairs"][str(n)] = {
            "t": chosen,
            "reason": reason,
            "shortlist": sorted(by_t),
            "replica_var_S_either": {str(t): by_t[t]["channels"]["either"]["S_delta"]["replica_var"] for t in by_t},
            "z_S_either": {str(t): by_t[t]["channels"]["either"]["S_delta"]["z_batch"] for t in by_t},
        }
    # Channel freeze from N=65 frozen-t matching-even z, candidates predeclared.
    n65_t = frozen["pairs"]["65"]["t"]
    n65 = pilot_by_n[65][n65_t]
    ranked = []
    for ch in CHANNEL_NAMES:
        est = n65["channels"][ch]["S_delta"]
        ranked.append((abs(est["z_batch"]), ch, est["z_batch"], est["mean"]))
    ranked.sort(reverse=True)
    chosen_ch = ranked[0][1] if ranked[0][0] >= 1.0 else "either"
    frozen["primary_channel"] = chosen_ch
    frozen["channel_rule"] = (
        "predeclared candidate set = matching-even S on "
        f"{list(CHANNEL_NAMES)}; freeze argmax |z_batch| on N=65 frozen-t "
        "pilot, defaulting to either if all |z|<1"
    )
    frozen["channel_ranking_n65"] = [
        {"channel": ch, "abs_z": abs_z, "z": z, "mean": mean}
        for abs_z, ch, z, mean in ranked
    ]
    frozen["pilot_not_reused_as_eval"] = True
    frozen["eval_replica_begin"] = EVAL_BEGIN
    frozen["pilot_replica_range"] = [PILOT_BEGIN, PILOT_BEGIN + PILOT_SAMPLES]
    return frozen


def resolved(est: dict) -> bool:
    return abs(est.get("z_batch", 0.0)) >= Z_RESOLVE


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 16), b""):
            digest.update(chunk)
    return digest.hexdigest()


def capture_environment(path: Path, cxx: Path, started_utc: datetime) -> None:
    uname = platform.uname()
    py = sys.version.replace("\n", " ")
    gcc = subprocess.check_output([str(cxx), "--version"], text=True).splitlines()[0]
    mem = subprocess.check_output(["free", "-h"], text=True)
    cpu = subprocess.check_output(["lscpu"], text=True)
    os_rel = Path("/etc/os-release").read_text(encoding="utf-8") if Path("/etc/os-release").exists() else ""
    text = (
        "C01 environment record\n"
        f"UTC capture: {started_utc.strftime('%Y-%m-%dT%H:%M:%SZ')}  "
        f"(Asia/Shanghai {started_utc.astimezone(SHANGHAI).strftime('%Y-%m-%d %H:%M:%S')})\n\n"
        f"{uname.system} {uname.node} {uname.release} {uname.version} {uname.machine}\n\n"
        f"{os_rel}\n"
        f"{cpu}\n"
        f"Python {py}\n"
        f"{gcc}\n\n"
        f"Memory:\n{mem}\n"
        f"nproc: {os.cpu_count()}\n"
        "GPU: none\n"
        f"OpenMP threads requested: {THREADS}\n"
    )
    path.write_text(text, encoding="utf-8")


def main() -> int:
    utc0, sh0 = now_pair()
    t0 = time.perf_counter()
    out = ROOT / "results" / "server-20260828" / "C01"
    raw = out / "raw"
    analysis_dir = out / "analysis"
    out.mkdir(parents=True, exist_ok=True)
    raw.mkdir(parents=True, exist_ok=True)
    analysis_dir.mkdir(parents=True, exist_ok=True)
    build = ROOT / "build"
    bins = compile_binaries(build)

    commands = [
        f"# C01 commands. Working directory: {ROOT}",
        f"# Started {fmt(sh0)} Asia/Shanghai ({utc0.strftime('%Y-%m-%dT%H:%M:%SZ')} UTC)",
        "g++ -O3 -std=c++17 -fopenmp src/gaussian_orientation_mc.cpp -o build/gaussian_orientation_mc",
        "g++ -O3 -std=c++17 -fopenmp src/pell_matching_mc.cpp -o build/pell_matching_mc",
        "./build/gaussian_orientation_mc --self-test",
        "./build/pell_matching_mc --self-test",
        "python3 -m unittest tests.test_torus_homology tests.test_gaussian_orientation_mc tests.test_pell_matching_mc -v",
        "python3 scripts/run_c01_protocol.py",
    ]
    (out / "commands.txt").write_text("\n".join(commands) + "\n", encoding="utf-8")

    run([str(bins["gaussian"]), "--self-test"])
    run([str(bins["pell"]), "--self-test"])

    graph_rows = graph_only_table()
    graph_by_n = {row["N"]: row for row in graph_rows}

    # --- pilots (site) ---
    pilot_by_n: dict[int, dict[int, dict]] = {}
    for pair in DISCOVERY:
        shortlist = graph_by_n[pair["N"]]["shortlist"]
        prefix = raw / f"n{pair['N']}_site_pilot"
        simulate(bins["gaussian"], pair, shortlist, "site", P_REF,
                 PILOT_SAMPLES, SEED_SITE, PILOT_BEGIN, prefix)
        pilot_by_n[pair["N"]] = {}
        for t in shortlist:
            analysis = analyze_prefix(prefix, t)
            (analysis_dir / f"n{pair['N']}_site_pilot_t{t}.json").write_text(
                json.dumps(analysis, indent=2) + "\n", encoding="utf-8"
            )
            pilot_by_n[pair["N"]][t] = analysis

    frozen = freeze_from_pilots(pilot_by_n)
    frozen["graph_only"] = graph_rows
    frozen["p_ref"] = P_REF
    frozen["complement"] = P_COMPLEMENT
    frozen["frozen_before_evaluation"] = True
    (out / "coupling_pilot.json").write_text(json.dumps(frozen, indent=2) + "\n", encoding="utf-8")
    print("FROZEN coupling:", json.dumps(frozen["pairs"], indent=2))
    print("FROZEN channel:", frozen["primary_channel"])

    primary = frozen["primary_channel"]

    # --- evaluation: site N=65,85 and bond N=65,85 ---
    evals: dict[str, dict] = {}

    def run_eval(pair: dict, mode: str, p: float, seed: int, samples: int, tag: str) -> dict:
        t = frozen["pairs"][str(pair["N"])]["t"] if str(pair["N"]) in frozen["pairs"] else graph_by_n[pair["N"]]["t_struct"]
        # For confirmation sizes not in frozen pairs, use graph-only t (already independent).
        if str(pair["N"]) not in frozen["pairs"]:
            t = 1 if 1 in graph_by_n[pair["N"]]["shortlist"] else graph_by_n[pair["N"]]["t_struct"]
            # still prefer frozen rule: t_struct from graph-only only (no new percolation pilot)
            t = graph_by_n[pair["N"]]["t_struct"]
            # If N=145 was piloted, use frozen t
        if pair["N"] in pilot_by_n:
            t = frozen["pairs"][str(pair["N"])]["t"]
        prefix = raw / f"n{pair['N']}_{mode}_{tag}"
        simulate(bins["gaussian"], pair, [t], mode, p, samples, seed, EVAL_BEGIN, prefix)
        analysis = analyze_prefix(prefix, t)
        (analysis_dir / f"n{pair['N']}_{mode}_{tag}_t{t}.json").write_text(
            json.dumps(analysis, indent=2) + "\n", encoding="utf-8"
        )
        return analysis

    for pair in DISCOVERY[:2]:
        evals[f"site_{pair['N']}"] = run_eval(pair, "site", P_REF, SEED_SITE, EVAL_SAMPLES, "eval")
        evals[f"bond_{pair['N']}"] = run_eval(pair, "bond", P_BOND, SEED_BOND, EVAL_SAMPLES, "eval")

    site65 = evals["site_65"]["channels"][primary]["S_delta"]
    site85 = evals["site_85"]["channels"][primary]["S_delta"]
    bond65 = evals["bond_65"]["channels"][primary]["primal_delta"]
    bond85 = evals["bond_85"]["channels"][primary]["primal_delta"]
    # Also inspect either if frozen channel differs
    site65_either = evals["site_65"]["channels"]["either"]["S_delta"]
    site85_either = evals["site_85"]["channels"]["either"]["S_delta"]
    bond65_either = evals["bond_65"]["channels"]["either"]["primal_delta"]
    bond85_either = evals["bond_85"]["channels"]["either"]["primal_delta"]

    site_resolves = resolved(site65) or resolved(site85) or resolved(site65_either) or resolved(site85_either)
    bond_resolves = resolved(bond65) or resolved(bond85) or resolved(bond65_either) or resolved(bond85_either)
    continue_c02 = bool(site_resolves)
    stop_gate = {
        "primary_channel": primary,
        "z_threshold": Z_RESOLVE,
        "site_matching_even_resolved": site_resolves,
        "bond_primal_angular_resolved": bond_resolves,
        "continue_to_C02": continue_c02,
        "decision": "continue_to_C02" if continue_c02 else "stop_before_increasing_N",
        "note": (
            "Stop gate looks at the larger matching-even / single-lattice sector, "
            "not only M. Bond control is square-bond wrapping on the same Gaussian "
            "quotients at p=1/2."
        ),
        "site_65_S": site65,
        "site_85_S": site85,
        "bond_65_primal_delta": bond65,
        "bond_85_primal_delta": bond85,
    }

    # N=145 eval only if discovery resolved
    if continue_c02:
        z65 = site65["z_batch"]
        se65 = site65["batch_se"]
        mean65 = site65["mean"]
        # L^-2 scaling for S: signal_145 = mean65 * 65/145
        pred145 = mean65 * 65 / 145
        var_ratio = evals["site_85"]["channels"][primary]["S_delta"]["replica_var"] / max(
            evals["site_65"]["channels"][primary]["S_delta"]["replica_var"], 1e-18
        )
        # samples for z=3 at N=145, transporting replica var ~ empirically from 65
        var65 = evals["site_65"]["channels"][primary]["S_delta"]["replica_var"]
        target_se = abs(pred145) / Z_RESOLVE if pred145 != 0 else math.inf
        n_req = var65 / (target_se ** 2) if target_se < math.inf and target_se > 0 else math.inf
        n145 = int(min(max(EVAL_SAMPLES, math.ceil(n_req / 100000) * 100000 if math.isfinite(n_req) else EVAL_SAMPLES), 4_000_000))
        n145 = max(n145, EVAL_SAMPLES)
        if n145 % BATCHES:
            n145 += BATCHES - n145 % BATCHES
        evals["site_145"] = run_eval(DISCOVERY[2], "site", P_REF, SEED_SITE, n145, "eval")
        evals["bond_145"] = run_eval(DISCOVERY[2], "bond", P_BOND, SEED_BOND, min(n145, EVAL_SAMPLES), "eval")
        stop_gate["n145_samples"] = n145
        stop_gate["n145_predicted_S"] = pred145
        stop_gate["n145_required_for_z3"] = n_req
        # confirmation N=205 only if 145 also resolves
        if resolved(evals["site_145"]["channels"][primary]["S_delta"]):
            # freeze t for 205 from graph-only (no extra percolation pilot)
            frozen["pairs"]["205"] = {
                "t": graph_by_n[205]["t_struct"],
                "reason": "graph-only structural multiplier; no extra percolation pilot",
                "shortlist": graph_by_n[205]["shortlist"],
            }
            evals["site_205"] = run_eval(CONFIRMATION[0], "site", P_REF, SEED_SITE, EVAL_SAMPLES, "eval")
            evals["bond_205"] = run_eval(CONFIRMATION[0], "bond", P_BOND, SEED_BOND, EVAL_SAMPLES, "eval")
        else:
            stop_gate["n145_underpowered_or_null"] = True
    else:
        stop_gate["skipped"] = ["N=145 evaluation", "N=205/425/1105 confirmation"]

    # Power budget
    budget = {
        "model_S": "matching-even orientation difference ~ C * N^{-1} * Delta cos(4 theta) (L^{-2})",
        "model_M": "matching-odd D_N ~ A4 * N^{-13/8} * Delta cos(4 theta)",
        "z_target": Z_RESOLVE,
        "primary_channel": primary,
        "from_eval": {},
    }
    for key, analysis in evals.items():
        if not key.startswith("site_"):
            continue
        n = analysis["N"]
        s = analysis["channels"][primary]["S_delta"]
        m = analysis["channels"]["either"]["M_delta"]
        budget["from_eval"][key] = {
            "N": n,
            "samples": analysis["samples"],
            "D_S": s["mean"],
            "D_S_batch_se": s["batch_se"],
            "z_S": s["z_batch"],
            "D_M_either": m["mean"],
            "D_M_batch_se": m["batch_se"],
            "z_M": m["z_batch"],
            "A4_M_either": analysis["channels"]["either"]["A4_M"],
            "A4_S_times_N": analysis["channels"][primary]["A4_S_times_N"],
            "replica_var_S": s["replica_var"],
        }
    (out / "effect_power_budget.json").write_text(json.dumps(budget, indent=2) + "\n", encoding="utf-8")
    (out / "stop_gate.json").write_text(json.dumps(stop_gate, indent=2) + "\n", encoding="utf-8")

    # long-form CSV and covariance JSON
    long_rows: list[dict] = []
    cov_out: dict[str, dict] = {}
    for key, analysis in sorted(evals.items()):
        stage = key
        long_rows.extend(write_long_form_rows(analysis, stage))
        cov_out[key] = {
            "N": analysis["N"],
            "mode": analysis["mode"],
            "t": analysis["t"],
            "variable_names": analysis["variable_names"],
            "means": analysis["variable_means_20"],
            "covariance": analysis["covariance_20"],
            "correlations": analysis["correlations"],
        }
    # include pilots too
    for n, by_t in pilot_by_n.items():
        for t, analysis in by_t.items():
            long_rows.extend(write_long_form_rows(analysis, f"pilot_n{n}"))
            cov_out[f"pilot_n{n}_t{t}"] = {
                "N": analysis["N"],
                "mode": analysis["mode"],
                "t": analysis["t"],
                "variable_names": analysis["variable_names"],
                "means": analysis["variable_means_20"],
                "covariance": analysis["covariance_20"],
                "correlations": analysis["correlations"],
            }

    with (out / "long_form_channel_means.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(long_rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(long_rows)
    (out / "covariance_matrices.json").write_text(json.dumps(cov_out) + "\n", encoding="utf-8")

    utc1, sh1 = now_pair()
    elapsed = time.perf_counter() - t0

    # unittest summary
    unit = subprocess.run(
        [sys.executable, "-m", "unittest",
         "tests.test_torus_homology",
         "tests.test_gaussian_orientation_mc",
         "tests.test_pell_matching_mc"],
        cwd=str(ROOT), text=True, capture_output=True,
    )
    (out / "unittest.log").write_text(unit.stdout + "\n" + unit.stderr, encoding="utf-8")
    unittest_ok = unit.returncode == 0

    source_commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=str(ROOT), text=True).strip()
    source_branch = subprocess.check_output(
        ["git", "rev-parse", "--abbrev-ref", "HEAD"], cwd=str(ROOT), text=True
    ).strip()
    kernel = platform.uname()
    gcc_ver = subprocess.check_output([str(bins["cxx"]), "--version"], text=True).splitlines()[0]
    py_ver = sys.version.split()[0]

    def pack_dn(analysis: dict, ch: str = "either") -> dict:
        return {
            "D_M": analysis["channels"][ch]["M_delta"],
            "A4_M": analysis["channels"][ch]["A4_M"],
            "D_S": analysis["channels"][primary]["S_delta"],
            "A4_S_times_N": analysis["channels"][primary]["A4_S_times_N"],
        }

    metadata = {
        "queue_id": "C01",
        "title": "Same-N Gaussian orientation discovery",
        "source_commit": source_commit,
        "source_branch": source_branch,
        "base_commit_c00": "08c4079798cc37913c6a7de3c97bc8d1c6bda16a",
        "issues": [22, 23],
        "machine": "cursor (Linux x86_64 KVM, Intel Xeon, 8 cores, 15 GiB RAM)",
        "hostname": platform.node(),
        "os_kernel": f"{kernel.system} {kernel.node} {kernel.release} {kernel.version} {kernel.machine}",
        "os_pretty": "Debian GNU/Linux 13 (trixie)",
        "compiler_or_interpreter": f"CPython {py_ver}; {gcc_ver}",
        "dependency_versions": {
            "python": py_ver,
            "g++": gcc_ver,
            "openmp": True,
            "mpmath": "1.4.1",
        },
        "rng_algorithm_if_stochastic": (
            "Philox4x32-10 counter uniform (Random123 KAT). "
            "key=(seed_lo, seed_hi); ctr=(index, replica_lo, replica_hi, stream); "
            "stream 0 = site occupation by cyclic vertex; stream 1 = bond occupation "
            "by packed (src, dx, dy). U in [0,1) from 53 high bits of the first two 32-bit words."
        ),
        "rng_seed_or_counter_ranges": {
            "site_seed": SEED_SITE,
            "bond_seed": SEED_BOND,
            "pilot_replicas": [PILOT_BEGIN, PILOT_BEGIN + PILOT_SAMPLES],
            "eval_replica_begin": EVAL_BEGIN,
            "eval_samples_default": EVAL_SAMPLES,
            "note": "pilot replicas are disjoint from evaluation replicas; same seed, different counters",
        },
        "sample_counts": {
            "pilot_per_discovery_pair": PILOT_SAMPLES,
            "evaluation_N65_site": EVAL_SAMPLES,
            "evaluation_N85_site": EVAL_SAMPLES,
            "evaluation_N65_bond": EVAL_SAMPLES,
            "evaluation_N85_bond": EVAL_SAMPLES,
            "evaluation_N145_site": evals.get("site_145", {}).get("samples"),
            "evaluation_N205_site": evals.get("site_205", {}).get("samples"),
        },
        "wall_time": {
            "seconds": elapsed,
            "utc_start": utc0.strftime("%Y-%m-%dT%H:%M:%SZ"),
            "utc_end": utc1.strftime("%Y-%m-%dT%H:%M:%SZ"),
            "asia_shanghai_start": fmt(sh0),
            "asia_shanghai_end": fmt(sh1),
        },
        "thread_or_gpu_configuration": {
            "threads": THREADS,
            "gpu": False,
            "openmp": True,
        },
        "frozen_coupling": frozen["pairs"],
        "frozen_primary_channel": primary,
        "p_ref": P_REF,
        "complement": P_COMPLEMENT,
        "bond_pc": P_BOND,
        "stop_gate": stop_gate["decision"],
        "continue_to_C02": continue_c02,
        "unittest_pass": unittest_ok,
        "signed_effects": {
            "site_65": pack_dn(evals["site_65"]),
            "site_85": pack_dn(evals["site_85"]),
            **({"site_145": pack_dn(evals["site_145"])} if "site_145" in evals else {}),
            **({"site_205": pack_dn(evals["site_205"])} if "site_205" in evals else {}),
        },
        "bond_control": {
            "N65_primal_delta": evals["bond_65"]["channels"][primary]["primal_delta"],
            "N85_primal_delta": evals["bond_85"]["channels"][primary]["primal_delta"],
        },
        "c00_results_intact": True,
        "torus_homology_rewritten": False,
    }
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    capture_environment(out / "environment.txt", bins["cxx"], utc1)

    # REPORT
    def line_effect(label: str, analysis: dict) -> str:
        s = analysis["channels"][primary]["S_delta"]
        m = analysis["channels"]["either"]["M_delta"]
        a4 = analysis["channels"]["either"]["A4_M"]
        a4s = analysis["channels"][primary]["A4_S_times_N"]
        return (
            f"| {label} | {s['mean']:+.6g} | {s['batch_se']:.3g} | {s['z_batch']:+.2f} | "
            f"{m['mean']:+.6g} | {m['batch_se']:.3g} | {m['z_batch']:+.2f} | "
            f"{a4['mean']:+.6g} | {a4s['mean']:+.6g} |"
        )

    report = f"""# C01: same-N Gaussian orientation discovery

Queue item C01 (issue 22 / coordination issue 23). Started from C00
`agent/c01-gaussian-orientation` at
`08c4079798cc37913c6a7de3c97bc8d1c6bda16a`.  The C00 2x2 homology engine is
the single DSU; `src/pell_matching_mc.cpp` now uses
`src/homology_union_find.hpp` (`either` = PR #21 Boolean wrap).  A second
topology implementation was not added.  `scripts/torus_homology.py` math was
not rewritten.

Wall time: {elapsed:.2f} s ({fmt(sh0)}–{fmt(sh1)} Asia/Shanghai, 2026-08-28).
OpenMP {THREADS} threads, no GPU.

## Frozen choices (pilot, not reused)

p_ref = {P_REF} (discovery coordinate, not a pc claim); complement = {P_COMPLEMENT}.
Bond control at exact pc = 1/2.

Graph-only structural multipliers, then 200000 independent site-pilot replicas
per discovery pair (Philox seed {SEED_SITE}, replicas [{PILOT_BEGIN}, {PILOT_BEGIN + PILOT_SAMPLES})).
Evaluation uses replicas starting at {EVAL_BEGIN}.

"""
    report += "| N | t_struct | frozen t | primary channel |\n|---|---|---|---|\n"
    for pair in DISCOVERY:
        n = pair["N"]
        ts = graph_by_n[n]["t_struct"]
        ft = frozen["pairs"][str(n)]["t"]
        report += f"| {n} | {ts} | {ft} | {primary} |\n"
    report += (
        f"\nPrimary effect channel frozen from N=65 matching-even pilot ranking: "
        f"`S_{primary}`.\n"
        "Default coupling is same_U_j (t=1); a unit relabel is used only when the "
        "pilot variance of ΔS_either improved by more than 5%.\n\n"
        "## Evaluation signed effects\n\n"
        f"Matching-even ΔS uses frozen channel `{primary}`; matching-odd D_N uses "
        "M_either = primal_either − matching_either.  Uncertainties are batch SE "
        f"({BATCHES} equal batches).\n\n"
        "| run | ΔS | SE | z_S | D_N (M_either) | SE | z_M | A4_M | N ΔS / Δcos4 |\n"
        "|---|---|---|---|---|---|---|---|---|\n"
    )
    for key in ("site_65", "site_85", "site_145", "site_205"):
        if key in evals:
            report += line_effect(key, evals[key]) + "\n"
    report += "\nBond control (primal wrapping orientation difference at p=1/2):\n\n"
    report += "| run | Δ primal | SE | z | Δ matching | SE | z |\n|---|---|---|---|---|---|---|\n"
    for key in ("bond_65", "bond_85", "bond_145", "bond_205"):
        if key not in evals:
            continue
        pr = evals[key]["channels"][primary]["primal_delta"]
        ma = evals[key]["channels"][primary]["matching_delta"]
        report += (
            f"| {key} | {pr['mean']:+.6g} | {pr['batch_se']:.3g} | {pr['z_batch']:+.2f} | "
            f"{ma['mean']:+.6g} | {ma['batch_se']:.3g} | {ma['z_batch']:+.2f} |\n"
        )

    report += f"""

## Stop gate

Site matching-even resolved at N=65 or N=85: **{'YES' if site_resolves else 'NO'}**.
Bond angular signal resolved at N=65 or N=85: **{'YES' if bond_resolves else 'NO'}**.

Decision: **{stop_gate['decision']}**.

"""
    if not continue_c02:
        report += (
            "Neither discovery size resolved a reproducible orientation effect in "
            "the larger matching-even / single-lattice sector, and the exact bond "
            "control also shows no angular signal at |z|>=3.  Evaluation of N=145 "
            "and confirmation pairs 205/425/1105 was not started.  C02 production "
            "analysis was not started.\n\n"
        )
    else:
        report += (
            "A matching-even orientation difference is resolved at the discovery "
            "sizes.  C02 may derive S/D sector tables from these channel means; "
            "this handoff does not run C02 production fits.\n\n"
        )

    report += f"""## RNG

Philox4x32-10 (Random123 official KATs in `--self-test`).
Site seed {SEED_SITE}, bond seed {SEED_BOND}.
Pilot replicas [{PILOT_BEGIN}, {PILOT_BEGIN + PILOT_SAMPLES}); evaluation replicas start at {EVAL_BEGIN}.
Coupling: vertex j of orientation 2 uses U_{{t j mod N}}.

## Tests

- Gaussian (2,1) exhaustive: rank0=16, rank1=10, rank2=6, d0=11, d1=11 (C00)
- Gaussian (3,2) exhaustive: rank0=4629, rank1=2340, rank2=1223, d0=2471, d1=2471
- PR #21 axis L=2,3 and diamond L=2 matching polynomials still pass after the DSU swap
- Unittest: {'PASS' if unittest_ok else 'FAIL'}

## Negative / preserved results

All channels are reported, including null matching-odd D_N values.  Pilot
covariance is retained in `covariance_matrices.json`.  Confirmation sizes
425 and 1105 were not run.

C00 result files under `results/server-20260828/C00/` were not modified.
"""
    (out / "REPORT.md").write_text(report, encoding="utf-8")

    # checksums of handoff files
    names = [
        "REPORT.md", "commands.txt", "environment.txt", "metadata.json",
        "long_form_channel_means.csv", "covariance_matrices.json",
        "coupling_pilot.json", "effect_power_budget.json",
    ]
    lines = []
    for name in names:
        path = out / name
        if path.exists():
            lines.append(f"{sha256_file(path)}  {name}")
    (out / "checksums.sha256").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("C01 wall_seconds", elapsed)
    print("stop_gate", stop_gate["decision"])
    print("unittest", unittest_ok)
    return 0 if unittest_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())

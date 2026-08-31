#!/usr/bin/env python3
"""Compile and enumerate the two frozen single-A-vacancy parents once; no score."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
CONTRACT = ROOT/"analysis/p337_endpoint_defect_contract.json"
CPP = ROOT/"scripts/p337_endpoint_defect_exact.cpp"


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    start_utc = datetime.now(timezone.utc).isoformat()
    contract = json.loads(CONTRACT.read_text())
    if contract["parent_geometries"] != [[5, 5], [1, 7]] or contract["coefficient_degree"] != 25:
        raise ValueError("this producer is restricted to the frozen parent pair and degree")
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    binary = Path(tempfile.mkdtemp(prefix="p337-endpoint-defect-"))/"enumerate"
    compiler = subprocess.check_output(["/usr/bin/clang++", "--version"], text=True).splitlines()[0]
    compile_command = ["/usr/bin/clang++", "-O3", "-std=c++17", str(CPP), "-o", str(binary)]
    subprocess.run(compile_command, check=True)
    binary_hash = sha(binary)

    def enumerate_one(item):
        geometry, label = item
        path = output/f"{label}.csv"
        command = [str(binary), *(str(x) for x in geometry), str(path)]
        completed = subprocess.run(command, check=True, capture_output=True, text=True)
        receipt = json.loads(completed.stdout)
        receipt.update(geometry=geometry, file=path.name, sha256=sha(path),
                       command=command, exit_code=completed.returncode)
        print(json.dumps(receipt), flush=True)
        return receipt

    with ThreadPoolExecutor(max_workers=contract["workers"]) as pool:
        receipts = list(pool.map(enumerate_one, zip(contract["parent_geometries"], ("first", "second"))))
    code_commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    run = {
        "schema": "matching-one.p337-single-A-vacancy.counts.v1", "status": "completed",
        "started_utc": start_utc, "finished_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": time.perf_counter()-started, "code_commit": code_commit,
        "contract_freeze_commit": "6c65157f8211565a374148e087dde542f2502c06",
        "contract": contract, "enumerations": receipts,
        "compile_command": compile_command, "compiler": compiler, "binary_sha256": binary_hash,
        "source_sha256": {"contract": sha(CONTRACT), "cpp": sha(CPP), "runner": sha(Path(__file__))},
        "python": sys.version, "machine": platform.machine(), "command": sys.argv,
        "saturated_baseline_enumerated": False, "response_scored": False,
        "new_random_samples": 0, "cloud_jobs": 0, "tests_run": 0,
    }
    (output/"run.json").write_text(json.dumps(run, indent=2)+"\n")
    (output/"README.md").write_text("""# Exact single-A-vacancy endpoint coefficients

`first.csv` is parent (5,5), corresponding to child (5,0).
`second.csv` is parent (1,7), corresponding to child (4,3).
In each parent, origin A is vacant, the other 24 A vertices occupied, and
all 2^25 free-B configurations are enumerated once on the actual N50 graph.

Rows k=0..25 store count and sums of q, E, Sstar, q*Sstar and E*Sstar.
The free-coordinate degree is **25**, Ktotal=24+k; integer sums already
include multiplicity. Weight them by p^k*(1-p)^(25-k), not an extra binomial.

Sstar=CB_NN+CW_matching+F+Bvac, Bvac=100−4*Ktotal+Bocc.
All fixed A neighbors are active independently of vertex IDs. Faces close
at their last free B activation. No conventional child graph is substituted
for the defective parent. Only the origin defect was computed; saturated
coefficients must come from the preexisting N25 complement dictionary.

The declared finite endpoint question is H_gain: U(s,t)=g(s)U_child(t),
with zero determinant R=U*U_st−U_s*U_t at s=1,t=0. This producer does not
score R or any response; root-owned rational scoring remains separate.

Contract, code/source hashes, compiler and exact enumeration receipts are
in run.json. No Monte Carlo, old-script replay, tests or cloud jobs ran.
""")
    print(json.dumps({"status": "completed", "code_commit": code_commit,
                      "output": str(output), "elapsed_seconds": run["elapsed_seconds"]}), flush=True)


if __name__ == "__main__":
    main()

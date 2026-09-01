#!/usr/bin/env python3
"""Export only frozen Git blobs; never modify the source checkout."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

FORK = "e32a85939279b8574278024d647b56d2d1485247"
CONTACT = "959a7fa26677c416b874d272f1ba66523fb38f73"
COMMON = "4db356e1b026853468f94d59d938895a2367ceb7"
RANK = "73608ba9d3eef34c6980cb5a049f726cfebdd72d"

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path(__file__).parent / "inputs")
    args = parser.parse_args()
    entries = []
    for n in (325, 425):
        for batch in range(20):
            entries.extend([
                (FORK, f"results/p334-nested-next-label-forks/N{n}/N{n}.batch{batch:02}.csv.gz", f"forks/N{n}.batch{batch:02}.csv.gz"),
                (CONTACT, f"results/p334-next-label-contact-coordinates/N{n}/N{n}.batch{batch:02}.csv.gz", f"contact/N{n}.batch{batch:02}.csv.gz"),
            ])
    entries.extend([
        (FORK, "src/p334_nested_next_label_forks.cpp", "provenance/p334_nested_next_label_forks.cpp"),
        (COMMON, "results/p334-common-label-euler-tangent/score.json", "common_label_score.json"),
        (RANK, "results/p334-common-label-response-rank/score.json", "ensemble_rank_score.json"),
    ])
    manifest = {"source_repository": "https://github.com/LightChainr/Matching-One", "files": []}
    args.output.mkdir(parents=True, exist_ok=True)
    for commit, path, target in entries:
        data = subprocess.check_output(["git", "show", f"{commit}:{path}"], cwd=args.repo)
        destination = args.output / target
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists() and destination.read_bytes() != data:
            raise RuntimeError(f"Refuse to replace differing input: {destination}")
        destination.write_bytes(data)
        manifest["files"].append({"commit": commit, "git_path": path, "local_path": target,
                                  "sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data)})
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"files": len(entries), "bytes": sum(x["bytes"] for x in manifest["files"])}))

if __name__ == "__main__":
    main()

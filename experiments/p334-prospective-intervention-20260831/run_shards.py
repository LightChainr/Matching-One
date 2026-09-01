#!/usr/bin/env python3
"""One fixed shard pipeline: prefix/census -> saved predictions -> independent tails.

No global F1 approval is required. Never execute until the supplied commit
has frozen this package; model/file hashes are checked against FREEZE.json.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def predictions(n, shard, output, model):
    stem = f"N{n}.shard{shard:03d}"
    path = output / (stem + ".prediction.csv.gz")
    if path.exists():
        raise ValueError("refuse to overwrite predictions")
    with gzip.open(output / (stem + ".prefix.csv.gz"), "rt") as source, gzip.open(path, "wt", newline="") as dest:
        reader = csv.DictReader(source)
        writer = csv.writer(dest)
        writer.writerow(["index", "batch", "counter", "rank_first", "rank_second"] + [f"{o}_{s}_{f}" for o in ("first", "second") for s in ("first", "second") for f in ("C", "W")])
        count = 0
        for r in reader:
            pred = []
            for ori in ("first", "second"):
                if r["rank_first"] != "0" or r["rank_second"] != "0":
                    pred.extend([0.]*4)
                    continue
                m = model["sizes"][str(n)]["point"][ori]
                features = [float(r[f"{ori}_{name}"]) for name in ("mass", "energy", "degree", "loop")]
                pred.extend([m["mean_responses"][j] + sum((features[i]-m["mean_features"][i])*m["beta"][i][j] for i in range(4)) for j in range(4)])
            writer.writerow([r[k] for k in ("index", "batch", "counter", "rank_first", "rank_second")] + pred)
            count += 1
    if count != 5000:
        raise ValueError("incomplete all-prefix prediction shard")
    digest = sha(path)
    (output / (stem + ".prediction.sha256")).write_text(digest + "\n")
    return digest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--freeze-commit", required=True)
    p.add_argument("--sizes", default="325,425")
    p.add_argument("--shards", required=True, help="fixed list/range, e.g. 0:30 or 30:60")
    p.add_argument("--workers", type=int, default=14)
    p.add_argument("--output", type=Path, default=ROOT/"production")
    args = p.parse_args()
    if len(args.freeze_commit)!=40 or any(c not in "0123456789abcdef" for c in args.freeze_commit):
        p.error("a full committed F0 sha is required")
    freeze = json.loads((ROOT/"FREEZE.json").read_text())
    if freeze.get("status") != "ready_for_root_freeze" or freeze["prefixes_per_size"]!=300000 or freeze["shards_per_size"]!=60:
        raise ValueError("not a completed fixed-budget freeze record")
    for path,digest in freeze["file_sha256"].items():
        if sha(ROOT/path)!=digest:
            raise ValueError(f"frozen file changed: {path}")
    if not 1 <= args.workers <= 14:
        raise ValueError("workers must respect 14.5-CPU cgroup")
    if ":" in args.shards:
        lo,hi=map(int,args.shards.split(":")); shards=list(range(lo,hi))
    else:
        shards=list(map(int,args.shards.split(",")))
    if len(set(shards))!=len(shards) or any(not 0<=s<60 for s in shards):
        raise ValueError("duplicate or invalid shard")
    sizes=list(map(int,args.sizes.split(",")))
    if len(set(sizes))!=len(sizes) or any(n not in (325,425) for n in sizes):
        raise ValueError("invalid sizes")
    model=json.loads((ROOT/"existing_model.json").read_text())
    args.output.mkdir(parents=True,exist_ok=True)
    def job(item):
        n,shard=item; stem=f"N{n}.shard{shard:03d}"
        spec=freeze["sizes"][str(n)]
        common=[str(n),str(shard),"5000",str(spec["counter_begin"]),str(spec["prefix_seed"]),str(args.output),args.freeze_commit]
        receipt={"N":n,"shard":shard,"freeze_commit":args.freeze_commit,"features_started":datetime.now(timezone.utc).isoformat()}
        log=args.output/(stem+".log")
        if log.exists():
            raise ValueError("shard already started; inspect before any retry")
        with log.open("w") as out:
            subprocess.run([str(ROOT/"prospective"),"features",*common],check=True,stdout=out,stderr=subprocess.STDOUT)
            digest=predictions(n,shard,args.output,model)
            receipt.update(predictions_sealed=datetime.now(timezone.utc).isoformat(),prediction_sha256=digest,
                           prefix_sha256=sha(args.output/(stem+".prefix.csv.gz")),census_sha256=sha(args.output/(stem+".census.csv.gz")))
            (args.output/(stem+".run.json")).write_text(json.dumps(receipt,indent=2)+"\n")
            subprocess.run([str(ROOT/"prospective"),"tails",*common,digest],check=True,stdout=out,stderr=subprocess.STDOUT)
            receipt.update(tails_completed=datetime.now(timezone.utc).isoformat(),tail_sha256=sha(args.output/(stem+".tails.csv.gz")))
            subprocess.run([sys.executable,str(ROOT/"score_prospective.py"),"--summarize-shard",stem,"--input",str(args.output)],check=True,stdout=out,stderr=subprocess.STDOUT)
            receipt.update(sufficient_sha256=sha(args.output/(stem+".sufficient.json")))
            (args.output/(stem+".run.json")).write_text(json.dumps(receipt,indent=2)+"\n")
        print(stem+" complete",flush=True)
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        list(pool.map(job,[(n,s) for n in sizes for s in shards]))


if __name__=="__main__":
    main()

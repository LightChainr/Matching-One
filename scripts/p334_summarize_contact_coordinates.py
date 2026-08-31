#!/usr/bin/env python3
"""Descriptive contact architecture of the immutable sampled-label readout."""
from collections import Counter
import csv
import gzip
from hashlib import sha256
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/"results/p334-next-label-contact-coordinates"


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def main():
    sizes = {}
    for n in (325, 425):
        metadata = json.loads((OUT/f"N{n}/metadata.json").read_text())
        counts = Counter()
        unique_double, prefix_double = set(), set()
        double_arch = Counter()
        for path in sorted((OUT/f"N{n}").glob("*.csv.gz")):
            with gzip.open(path, "rt") as stream:
                for row in csv.DictReader(stream):
                    cell = row["first_oldrank"]+row["second_oldrank"]
                    for orientation in ("first", "second"):
                        if row[f"{orientation}_oldrank"] != "0":
                            continue
                        after = int(row[f"{orientation}_rank_after"])
                        e, c = (int(row[f"{orientation}_{name}"]) for name in ("e", "c"))
                        counts[f"{orientation}.cell{cell}.R0"] += 1
                        counts[f"{orientation}.cell{cell}.R0_to_{after}"] += 1
                        if after == 0:
                            counts[f"{orientation}.safe.e{e}.c{c}"] += 1
                        if after == 2:
                            architecture = int(row[f"{orientation}_r0_rank2_arch"])
                            double_arch[f"{orientation}.arch{architecture}.e{e}.c{c}"] += 1
                            unique_double.add((row["counter"],row["next_label"],orientation))
                            prefix_double.add((row["counter"],orientation))
        sizes[str(n)] = {"metadata": metadata, "R0_rank_cell_and_safe_contact_counts": dict(sorted(counts.items())),
                         "rank2_architecture_contacts": dict(sorted(double_arch.items())),
                         "unique_prefix_label_orientation_rank2_events": len(unique_double),
                         "prefix_orientation_pairs_with_sampled_rank2": len(prefix_double)}
    source_paths = [ROOT/f"results/p334-full-birth-archive/N{n}.csv" for n in (325,425)]
    source_paths += sorted((ROOT/"results/p334-nested-next-label-forks").rglob("*.csv.gz"))
    source_paths += sorted((ROOT/"results/p334-nested-next-label-forks").rglob("metadata.json"))
    result_paths = sorted(OUT.rglob("*.csv.gz"))+sorted(OUT.rglob("metadata.json"))
    manifest = {"code_commit":"b044e6452d3342f496d77a67f13319caa06c92fe",
                "fork_source_commit":"e32a85939279b8574278024d647b56d2d1485247",
                "original_prefix_source_commit":"9c495ab13e65f2bc93dc0849ee3b73f88724c4b1",
                "source_sha256":{str(p.relative_to(ROOT)):digest(p) for p in source_paths},
                "output_sha256":{str(p.relative_to(ROOT)):digest(p) for p in result_paths},
                "code_sha256":{path:sha256(subprocess.check_output(["git","show",f"b044e645:{path}"],cwd=ROOT)).hexdigest()
                               for path in ("src/p334_next_label_contact_coordinates.cpp","src/threshold_rank_integer_period_mc.cpp")}}
    (OUT/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    (OUT/"summary.json").write_text(json.dumps({"new_samples":0,"new_tail_replay":0,"DP_calls":0,"sizes":sizes},indent=2)+"\n")
    census = ROOT/"results/p334-exact-next-label-safe-census"
    census_paths = sorted(census.rglob("*.csv.gz"))+sorted(census.rglob("metadata.json"))
    census_manifest = {"data_lock_commit":"e9dc7a1078b2c64b319f4a36ffc1c844e8426aa0",
                       "code_commit":"0e4db1b8ccae26f2953522ff7428162b12e9e8fa",
                       "source_prefix_commit":"9c495ab13e65f2bc93dc0849ee3b73f88724c4b1",
                       "output_sha256":{str(p.relative_to(ROOT)):digest(p) for p in census_paths},
                       "code_sha256":{path:sha256(subprocess.check_output(["git","show",f"0e4db1b8:{path}"],cwd=ROOT)).hexdigest()
                                      for path in ("src/p334_exact_next_label_safe_census.cpp","src/p334_next_label_contact_coordinates.cpp","src/threshold_rank_integer_period_mc.cpp")}}
    (census/"manifest.json").write_text(json.dumps(census_manifest,indent=2)+"\n")
    for n, row in sizes.items():
        print(n, "rank2 contact architectures",row["rank2_architecture_contacts"],
              "unique labelled events",row["unique_prefix_label_orientation_rank2_events"])
        for c in row["metadata"]["R0_counts"]:
            classes=c["safe_loop_merger_00_01_10_11"]
            print(c["orientation"],"safe loops",classes[2]+classes[3],"/",c["safe_R0_draws"],
                  "fraction",(classes[2]+classes[3])/c["safe_R0_draws"],"classes",classes)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Export the two frozen pair/triple constraints for the bounded C++ oracle."""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results/p334-quartic-clock"


def main():
    triple_rows = json.loads((ROOT / "results/p334-pair-triple-clock/full_triples.json").read_text())["checkpoints"]
    prior = {row["counter"]: row for row in json.loads((ROOT / "results/p334-pair-triple-clock/pair_triple_survival.json").read_text())["records"]}
    lines = []
    for row in triple_rows:
        counter = row["replica_counter"]
        graph = json.loads((ROOT / prior[counter]["source_graph"]).read_text())
        pairs, triples = graph["minimal_trigger_pairs"], row["all_minimal_nonfaces"]
        lines.append(f"{counter} {len(pairs)} {len(triples)} {prior[counter]['independent_counts'][4]}")
        lines.extend(" ".join(map(str, e)) for e in pairs + triples)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    (OUTPUT / "frozen_constraints.txt").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()

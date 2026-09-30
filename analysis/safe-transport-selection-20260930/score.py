#!/usr/bin/env python3
"""Same-cell natural/reference/selection contrast, with one batch covariance."""
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
from itertools import zip_longest
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
OLD = HERE.parent / "safe-insertion-production-20260929/data"
A, B, M = 154646, 155385, 262144-155385
NAMES = ["natural_contrast", "safe_reference_contrast", "selection_contrast"]


def evaluate(batches, omit=None, detail=False):
    # Cell vectors: [early safe count,sum,natural count,sum, late ...].
    cells = defaultdict(lambda: [0., 0., 0., 0., 0., 0., 0., 0.])
    eligible = [0, 0]
    for batch in batches:
        if batch["batch"] == omit:
            continue
        for key, row in batch["cells"].items():
            for j, value in enumerate(row):
                cells[key][j] += value
            eligible[0] += int(row[2])
            eligible[1] += int(row[6])
    common = []
    for key, row in sorted(cells.items()):
        se, ye, ne, ze, sl, yl, nl, zl = row
        if ne and nl:
            weight = ne*nl/(ne+nl)
            common.append(dict(cell=list(key), weight=weight, safe_counts=[int(se), int(sl)],
                               natural_counts=[int(ne), int(nl)],
                               safe_means=[ye/se, yl/sl], natural_means=[ze/ne, zl/nl]))
    weight = math.fsum(c["weight"] for c in common)
    if not weight:
        raise ValueError("not_scoreable: no natural-cohort overlap")
    natural = math.fsum(c["weight"]*(c["natural_means"][0]-c["natural_means"][1]) for c in common)/weight
    safe = math.fsum(c["weight"]*(c["safe_means"][0]-c["safe_means"][1]) for c in common)/weight
    result = [natural, safe, natural-safe]
    if not detail:
        return result
    return result, dict(common_cells=len(common), weight_sum=weight,
                        eligible_natural_counts=eligible,
                        common_natural_counts=[sum(c["natural_counts"][g] for c in common) for g in range(2)],
                        common_safe_counts=[sum(c["safe_counts"][g] for c in common) for g in range(2)],
                        weighted_natural_means=[math.fsum(c["weight"]*c["natural_means"][g] for c in common)/weight for g in range(2)],
                        weighted_safe_means=[math.fsum(c["weight"]*c["safe_means"][g] for c in common)/weight for g in range(2)],
                        cells=common)


def main():
    manifest = json.loads((HERE / "data/run.json").read_text())
    assert manifest["status"] == "completed"
    batches = []
    for batch in manifest["batches"]:
        path = HERE / "data" / batch["file"]
        oldpath = OLD / batch["original_file"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == batch["sha256"]
        assert hashlib.sha256(oldpath.read_bytes()).hexdigest() == batch["original_sha256"]
        cells = defaultdict(lambda: [0., 0., 0., 0., 0., 0., 0., 0.])
        counts = Counter()
        with gzip.open(path, "rt") as f, gzip.open(oldpath, "rt") as g:
            for new, old in zip_longest(csv.DictReader(f), csv.DictReader(g)):
                assert new is not None and old is not None and new["J1"] == old["J1"]
                counts["prefixes"] += 1
                j1, entry, ref, alive, c, e = (int(new[k]) for k in ["J1", "entry_rank", "reference_rank_b", "natural_alive", "nu_b", "synergy_edges"])
                if entry == 0:
                    assert j1 == 0 and ref == 0 and alive == 0
                    counts["no_entrance_by_b"] += 1
                    continue
                if entry == 2:
                    assert 1 <= j1 <= B and ref == 2 and alive == 0
                    counts["direct_rank_two_entrance"] += 1
                    continue
                assert entry == 1 and 1 <= j1 <= B and alive in (0, 1)
                group = 0 if j1 <= A else 1
                label = "early" if group == 0 else "late"
                counts["rank_one_entrance_"+label] += 1
                if ref == -1:
                    assert alive == 0
                    counts["reference_cemetery_"+label] += 1
                    continue
                assert ref == 1 and 0 <= c <= M and e >= 0
                counts["reference_endpoint_"+label] += 1
                counts["natural_endpoint_"+label] += alive
                if c == M:
                    counts["no_safe_endpoint_"+label] += 1
                    continue
                d = (int(new["dx_entry"]), int(new["dy_entry"]))
                assert d != (0, 0)
                y = 2*e/(M-c)
                row = cells[d+(c,)]
                off = 4*group
                row[off] += 1
                row[off+1] += y
                row[off+2] += alive
                row[off+3] += alive*y
        assert counts["prefixes"] == batch["samples"] == 10000
        batches.append(dict(batch=batch["batch"], dependency_group=batch["dependency_group"],
                            counts=dict(counts), cells=dict(cells)))
    assert sorted(b["batch"] for b in batches) == list(range(14))
    estimate, support = evaluate(batches, detail=True)
    deletions = [evaluate(batches, omit=i) for i in range(14)]
    center = [math.fsum(d[j] for d in deletions)/14 for j in range(3)]
    covariance = [[13/14*math.fsum((d[i]-center[i])*(d[j]-center[j]) for d in deletions) for j in range(3)] for i in range(3)]
    se = [math.sqrt(covariance[i][i]) for i in range(3)]
    counts = sum((Counter(b["counts"]) for b in batches), Counter())
    result = dict(schema="matching-one.safe-transport-selection.v1", names=NAMES, estimate=estimate,
                  batch_delete_one_se=se, covariance=covariance, deletions=deletions, support=support,
                  counts=dict(counts), batch_counts=[dict(batch=b["batch"], counts=b["counts"], dependency_group=b["dependency_group"]) for b in batches],
                  all_original_births_match=True, new_independent_prefixes=0,
                  dependency_group="safe-insertion-independent-20260929/square/L512",
                  contract=json.loads((HERE / "contract.json").read_text()))
    target = HERE / "results"
    target.mkdir(exist_ok=False)
    (target / "result.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    lines = ["# Safe transport versus natural survival selection", "",
             "Same archived entrances, new conditional continuations; not an independent prefix block.", "",
             "| Contrast | Estimate | Whole-batch SE |", "|---|---:|---:|"]
    lines += [f"| {name} | {value:.12g} | {error:.12g} |" for name, value, error in zip(NAMES, estimate, se)]
    lines += ["", f"{support['common_cells']} common natural-cohort cells; one set of endpoint weights for all components.",
              f"Natural support {support['common_natural_counts']} / eligible {support['eligible_natural_counts']}; reference support {support['common_safe_counts']}.",
              "Natural = safe reference + selection exactly. Full singular3x3 covariance and14 aligned deletions in result.json.",
              "This decomposes a specified reference process; it does not distinguish imported birth geometry from subsequent safe growth, or identify original U/continuum fields.", ""]
    (target / "RESULT.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()

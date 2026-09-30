#!/usr/bin/env python3
"""One fixed finite-window mechanism readout with aligned whole-batch covariance."""
import csv
import gzip
import hashlib
from itertools import zip_longest
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
OLD = HERE.parent / "safe-insertion-production-20260929/data"
M, H = 262144-155385, 735
NAMES = ["finite_window_susceptibility", "baseline_survival", "two_step_susceptibility"]


def main():
    manifest = json.loads((HERE / "data/run.json").read_text())
    assert manifest["status"] == "completed"
    batches = []
    for batch in manifest["batches"]:
        path = HERE / "data" / batch["file"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == batch["sha256"]
        rows, values = 0, []
        with gzip.open(path, "rt") as f, gzip.open(OLD / batch["original_file"], "rt") as g:
            for new, old in zip_longest(csv.DictReader(f), csv.DictReader(g)):
                assert new is not None and old is not None
                assert all(new[k] == old[k] for k in ["J1", "rank_b", "dx_b", "dy_b", "nu_b"])
                rows += 1
                if int(new["rank_b"]) != 1:
                    continue
                c, e, d2, survives, total = (int(new[k]) for k in ["nu_b", "synergy_edges", "sum_degree_squared", "survival_h", "window_degree_sum"])
                s = M-c
                assert 0 <= c <= M and e >= 0 and survives in (0, 1)
                if not s:
                    assert survives == 0 and e == 0 and d2 == 0
                    values.append([0., 0., 0.])
                else:
                    variance_numerator = s*d2-4*e*e
                    assert variance_numerator >= 0
                    values.append([survives*(total/H-2*e/s), float(survives),
                                   -variance_numerator/(s*M*(M-1))])
        assert rows == batch["samples"] == 10000
        batches.append(dict(batch=batch["batch"], prefix_rows=rows, risk=len(values),
                            sums=[math.fsum(v[j] for v in values) for j in range(3)]))
    assert sorted(b["batch"] for b in batches) == list(range(14))
    def evaluate(omit=None):
        selected = [b for b in batches if b["batch"] != omit]
        n = sum(b["risk"] for b in selected)
        return [math.fsum(b["sums"][j] for b in selected)/n for j in range(3)]
    estimate = evaluate()
    deletions = [evaluate(i) for i in range(14)]
    center = [math.fsum(v[j] for v in deletions)/14 for j in range(3)]
    covariance = [[13/14*math.fsum((v[i]-center[i])*(v[j]-center[j]) for v in deletions)
                   for j in range(3)] for i in range(3)]
    se = [math.sqrt(covariance[i][i]) for i in range(3)]
    result = dict(schema="matching-one.geometric-source-window.v1", names=NAMES,
                  estimate=estimate, batch_delete_one_se=se, covariance=covariance,
                  deletions=deletions, batch_sufficient_statistics=batches,
                  prefix_rows=sum(b["prefix_rows"] for b in batches), risk_prefixes=sum(b["risk"] for b in batches),
                  all_first5_match_old=True, lag=H, clock="insertion count", source="one normalized safe-degree kick",
                  dependency_group="safe_insertion_independent_L512", new_independent_prefixes=0,
                  contract=json.loads((HERE / "contract.json").read_text()))
    target = HERE / "results"
    target.mkdir(exist_ok=False)
    (target / "result.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    lines = ["# Finite-window geometric-source response", "",
             "Same existing prefixes; new mechanism readout, not an independent block.", "",
             "| Metric | Estimate | Whole-batch SE |", "|---|---:|---:|"]
    lines += [f"| {name} | {value:.12g} | {error:.12g} |" for name, value, error in zip(NAMES, estimate, se)]
    lines += ["", f"{result['risk_prefixes']} rank-one prefixes / {result['prefix_rows']} total, 14 batches, lag={H}.",
              "Full3x3 covariance, aligned deletions and unchanged prefix correspondence are in result.json.",
              "An infinitesimal explicitly modified-source response; not a finite-theta effect, continuum limit, original-U claim or explained fraction of birth memory.", ""]
    (target / "RESULT.md").write_text("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()

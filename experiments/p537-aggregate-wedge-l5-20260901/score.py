#!/usr/bin/env python3
"""Combine exact shards and score the aggregate landing wedge at M(p)=0."""
from __future__ import annotations
import argparse, csv, json
from collections import defaultdict
from pathlib import Path
import mpmath as mp


def read_table(path: Path):
    meta = {}
    lines = []
    for line in path.read_text().splitlines():
        if line.startswith("# ") and "=" in line:
            key, value = line[2:].split("=", 1)
            meta[key] = value
        elif not line.startswith("#"):
            lines.append(line)
    rows = list(csv.DictReader(lines, delimiter="\t"))
    return meta, rows


def combine(paths):
    global_rows = defaultdict(lambda: [0, 0, 0])
    landing = defaultdict(lambda: [0, 0, 0])
    metas = []
    for path in paths:
        meta, rows = read_table(path); metas.append(meta)
        for row in rows:
            k = int(row["k"])
            if row["kind"] == "global":
                target = global_rows[k]
                target[0] += int(row["count"]); target[1] += int(row["sum_q"])
                target[2] += int(row["sum_source16"])
            else:
                target = landing[(row["transition"], k)]
                target[0] += int(row["signed_count"])
                target[1] += int(row["signed_source_mid16"])
                target[2] += int(row["unsigned_count"])
    Ls = {int(meta["L"]) for meta in metas}
    if len(Ls) != 1: raise ValueError("mixed L")
    shard_counts = {int(meta["shard_count"]) for meta in metas}
    shard_indices = {int(meta["shard_index"]) for meta in metas}
    if len(shard_counts) != 1 or shard_indices != set(range(shard_counts.pop())):
        raise ValueError("incomplete shard set")
    return Ls.pop(), global_rows, landing


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tables", nargs="+", type=Path)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--dps", type=int, default=80)
    args = ap.parse_args()
    if args.output.exists(): raise SystemExit(f"refusing to overwrite {args.output}")
    mp.mp.dps = args.dps
    L, global_rows, landing = combine(args.tables)
    N = L * L
    expected = 1 << N
    if sum(row[0] for row in global_rows.values()) != expected:
        raise ValueError("global mass mismatch")

    def bern(k, n, p): return p**k * (1-p)**(n-k)
    def matching(p): return mp.fsum(mp.mpf(row[1]) * bern(k,N,p) for k,row in global_rows.items())
    root = mp.findroot(matching, (mp.mpf("0.58"),mp.mpf("0.61")))
    mean_a = mp.fsum(mp.mpf(row[2]) * bern(k,N,root) for k,row in global_rows.items()) / (16*N*N)

    matrix = {}
    for tr in ("01","12"):
        B = mp.fsum(mp.mpf(landing[(tr,k)][0]) * bern(k,N-1,root) for k in range(N))
        T = mp.fsum(mp.mpf(landing[(tr,k)][0]) * (mp.mpf(k)+mp.mpf("0.5")-N*root) * bern(k,N-1,root) for k in range(N))
        raw_A = mp.fsum(mp.mpf(landing[(tr,k)][1]) * bern(k,N-1,root) for k in range(N)) / (32*N*N)
        A = raw_A - mean_a * B
        matrix[tr] = dict(B=B,T=T,raw_A=raw_A,A=A,unsigned_count=sum(landing[(tr,k)][2] for k in range(N)))
    psi4 = matrix["01"]["T"]*matrix["12"]["A"] - matrix["12"]["T"]*matrix["01"]["A"]
    thermal_norm2 = matrix["01"]["T"]**2 + matrix["12"]["T"]**2
    chi_perp = psi4 / thermal_norm2
    source_perp = psi4 / mp.sqrt(thermal_norm2)

    def s(x): return mp.nstr(x,50)
    payload = {
        "schema":"matching-one/p537-aggregate-wedge-score/v1",
        "L":L,"N":N,"matching_root":s(root),"matching_residual":s(matching(root)),
        "mean_canonical_source":s(mean_a),
        "matrix":{tr:{k:(v if isinstance(v,int) else s(v)) for k,v in row.items()} for tr,row in matrix.items()},
        "Psi4":s(psi4),"Psi4_sign":int(mp.sign(psi4)),
        "thermal_norm2":s(thermal_norm2),
        "chi_perp":s(chi_perp),
        "source_perp":s(source_perp),
        "scaled":{
            "L4_chi_perp":s(L**4*chi_perp),
            "L6_source_perp":s(L**6*source_perp),
            "L8_Psi4":s(L**8*psi4),
            "L2_thermal_norm":s(L**2*mp.sqrt(thermal_norm2)),
        },
        "interpretation":"complete finite-population radius-one ell4 aggregate; signed thermal-gauge-invariant birth/completion wedge",
        "boundary":"one finite axis torus; no asymptotic exponent or continuum-field identity",
    }
    args.output.write_text(json.dumps(payload,indent=2)+"\n")
    print(json.dumps(payload,indent=2))


if __name__ == "__main__": main()

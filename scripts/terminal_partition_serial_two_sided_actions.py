#!/usr/bin/env python3
"""Exact two-sided multiplication operators x -> axb for serial partitions."""
from __future__ import annotations
import argparse,json
from collections import Counter
from pathlib import Path
from typing import Any,Mapping,Optional,Sequence
try:
 from scripts.terminal_partition_canonical import enumerate_rgs
 from scripts.terminal_partition_serial_category import serial_compose
except ModuleNotFoundError:
 from terminal_partition_canonical import enumerate_rgs
 from terminal_partition_serial_category import serial_compose
ROOT=Path(__file__).resolve().parents[1];SCHEMA="matching-one/terminal-partition-serial-two-sided-actions/v1"
def table():
 states=enumerate_rgs(4);idx={s:i for i,s in enumerate(states)}
 return tuple(tuple(idx[serial_compose(a,b)] for b in states) for a in states)
def operator(a,b,p):
 if not all(type(x) is int and 0<=x<15 for x in (a,b)):raise ValueError("multiplier outside table")
 return tuple(p[p[a][x]][b] for x in range(15))
def compose(f,g):return tuple(f[g[i]] for i in range(len(f)))
def build_artifact()->dict[str,Any]:
 p=table();fibers={}
 for a in range(15):
  for b in range(15):fibers.setdefault(operator(a,b,p),[]).append((a,b))
 ordered=sorted(fibers.items(),key=lambda z:(len(z[1]),z[1]))
 ranks=Counter(len(set(f)) for f in fibers);sizes=Counter(len(v) for v in fibers.values())
 failures=sum(compose(operator(a,b,p),operator(c,d,p))!=operator(p[a][c],p[d][b],p) for a in range(15) for b in range(15) for c in range(15) for d in range(15))
 return {"schema":SCHEMA,"issue":13,"status":"complete_two_sided_operator_census","distinct_operator_count":len(fibers),
  "rank_histogram":{str(k):v for k,v in sorted(ranks.items())},"pair_fiber_size_histogram":{str(k):v for k,v in sorted(sizes.items())},
  "operator_fibers":[{"operator":list(f),"multiplier_pairs":[list(x) for x in pairs]} for f,pairs in ordered],
  "exact_checks":{"225_pairs_collapse_to_133_operators":len(fibers)==133,"composition_law_has_zero_failures":failures==0,
   "rank_profile_is_exact":ranks=={1:4,2:89,5:36,15:4},"fiber_profile_is_exact":sizes=={1:121,2:8,16:1,22:2,28:1}},
  "claim_boundary":{"included":"all maps x->axb, pair fibers, transformation ranks, and exact composition law","excluded":"endomorphism classification, physical transfer operators, probability kernels, reliability, or thresholds","parent_issue":"remain open"}}
def validate_artifact(v:Mapping[str,Any]):
 e=build_artifact()
 if v!=e:raise ValueError("two-sided artifact does not reproduce")
 if set(e["exact_checks"].values())!={True}:raise ValueError("checks failed")
 return {"schema":SCHEMA,"status":"valid","operators":133}
def main(argv:Optional[Sequence[str]]=None)->int:
 q=argparse.ArgumentParser(description=__doc__);q.add_argument("--output",type=Path);q.add_argument("--validate",type=Path);a=q.parse_args(argv)
 if a.validate:print(json.dumps(validate_artifact(json.loads(a.validate.read_text())),indent=2,sort_keys=True));return 0
 s=json.dumps(build_artifact(),indent=2,sort_keys=True)+"\n"
 if a.output:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(s,encoding="utf-8")
 else:print(s,end="")
 return 0
if __name__=="__main__":raise SystemExit(main())

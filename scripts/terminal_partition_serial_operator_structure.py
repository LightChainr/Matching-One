#!/usr/bin/env python3
"""Exact identity, units, idempotents, and ranks of two-sided operators."""
from __future__ import annotations
import argparse,json
from collections import Counter
from pathlib import Path
from typing import Any,Mapping,Optional,Sequence
try:
 from scripts.terminal_partition_serial_two_sided_actions import table,operator,compose
except ModuleNotFoundError:
 from terminal_partition_serial_two_sided_actions import table,operator,compose
SCHEMA="matching-one/terminal-partition-serial-operator-structure/v1"
def data():
 p=table();ops=sorted({operator(a,b,p) for a in range(15) for b in range(15)});ix={f:i for i,f in enumerate(ops)};m=tuple(tuple(ix[compose(f,g)] for g in ops) for f in ops);return ops,m
def build_artifact()->dict[str,Any]:
 ops,m=data();n=len(ops);ident=ops.index(tuple(range(15)));idem=[i for i in range(n) if m[i][i]==i];units=[]
 for i in range(n):
  if any(m[i][j]==m[j][i]==ident for j in range(n)):units.append(i)
 ranks=Counter(len(set(f)) for f in ops)
 return {"schema":SCHEMA,"issue":13,"status":"complete_two_sided_operator_structure","operator_count":n,"identity_index":ident,"unit_indices":units,"idempotent_count":len(idem),"idempotent_indices":idem,"transformation_rank_histogram":{str(k):v for k,v in sorted(ranks.items())},
 "exact_checks":{"operator_count_is_133":n==133,"identity_index_is_thirteen":ident==13,"units_form_klein_four":len(units)==4 and all(m[u][u]==ident for u in units),"idempotent_count_is_76":len(idem)==76,"rank_profile_is_exact":ranks=={1:4,2:89,5:36,15:4}},
 "claim_boundary":{"included":"identity, unit group, idempotents, and transformation ranks of the 133 two-sided operators","excluded":"operator Green relations, power profiles, physical transfer operators, reliability, or thresholds","parent_issue":"remain open"}}
def validate_artifact(v:Mapping[str,Any]):
 e=build_artifact()
 if v!=e:raise ValueError("operator structure artifact mismatch")
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

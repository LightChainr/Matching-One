#!/usr/bin/env python3
"""Exact index-period census of the two-sided operator semigroup."""
from __future__ import annotations
import argparse,json
from collections import Counter
from pathlib import Path
from typing import Any,Mapping,Optional,Sequence
try:
 from scripts.terminal_partition_serial_two_sided_actions import table,operator,compose
except ModuleNotFoundError:
 from terminal_partition_serial_two_sided_actions import table,operator,compose
SCHEMA="matching-one/terminal-partition-serial-operator-powers/v1"
def data():
 p=table();ops=sorted({operator(a,b,p) for a in range(15) for b in range(15)});ix={f:i for i,f in enumerate(ops)};m=tuple(tuple(ix[compose(f,g)] for g in ops) for f in ops);return ops,m
def profile(x,m):
 seen={};y=x;k=1;seq=[]
 while y not in seen:seen[y]=k;seq.append(y);y=m[y][x];k+=1
 return seen[y],k-seen[y],seq
def build_artifact()->dict[str,Any]:
 _,m=data();rows=[];hist=Counter()
 for x in range(len(m)):
  i,p,s=profile(x,m);hist[(i,p)]+=1;rows.append({"operator":x,"index":i,"period":p,"distinct_power_indices":s})
 return {"schema":SCHEMA,"issue":13,"status":"complete_two_sided_operator_power_census","profiles":rows,"index_period_histogram":{f"{i},{p}":v for (i,p),v in sorted(hist.items())},"maximum_index":max(r["index"] for r in rows),"maximum_period":max(r["period"] for r in rows),
 "exact_checks":{"all_133_profiles_present":len(rows)==133,"profile_is_exact":hist=={(1,1):76,(2,1):40,(1,2):17},"maximum_index_is_two":max(r["index"] for r in rows)==2,"maximum_period_is_two":max(r["period"] for r in rows)==2},
 "claim_boundary":{"included":"all power sequences and exact index-period profiles of the 133 operators","excluded":"element powers already certified, operator Green relations, physical transfer operators, reliability, or thresholds","parent_issue":"remain open"}}
def validate_artifact(v:Mapping[str,Any]):
 e=build_artifact()
 if v!=e:raise ValueError("operator power artifact mismatch")
 if set(e["exact_checks"].values())!={True}:raise ValueError("checks failed")
 return {"schema":SCHEMA,"status":"valid","profiles":133}
def main(argv:Optional[Sequence[str]]=None)->int:
 q=argparse.ArgumentParser(description=__doc__);q.add_argument("--output",type=Path);q.add_argument("--validate",type=Path);a=q.parse_args(argv)
 if a.validate:print(json.dumps(validate_artifact(json.loads(a.validate.read_text())),indent=2,sort_keys=True));return 0
 s=json.dumps(build_artifact(),indent=2,sort_keys=True)+"\n"
 if a.output:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(s,encoding="utf-8")
 else:print(s,end="")
 return 0
if __name__=="__main__":raise SystemExit(main())

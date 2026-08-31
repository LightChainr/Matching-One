#!/usr/bin/env python3
"""Exact generalized-inverse census of the two-sided operator semigroup."""
from __future__ import annotations
import argparse,json
from collections import Counter
from pathlib import Path
from typing import Any,Mapping,Optional,Sequence
try:
 from scripts.terminal_partition_serial_two_sided_actions import table,operator,compose
except ModuleNotFoundError:
 from terminal_partition_serial_two_sided_actions import table,operator,compose
SCHEMA="matching-one/terminal-partition-serial-operator-inverses/v1"
def data():
 p=table();ops=sorted({operator(a,b,p) for a in range(15) for b in range(15)});ix={f:i for i,f in enumerate(ops)};m=tuple(tuple(ix[compose(f,g)] for g in ops) for f in ops);return ops,m
def build_artifact()->dict[str,Any]:
 _,m=data();rows=[]
 for a in range(len(m)):
  inv=[b for b in range(len(m)) if m[m[a][b]][a]==a and m[m[b][a]][b]==b];rows.append({"operator":a,"generalized_inverses":inv})
 hist=Counter(len(r["generalized_inverses"]) for r in rows)
 return {"schema":SCHEMA,"issue":13,"status":"complete_two_sided_operator_generalized_inverse_census","regular_operator_count":sum(bool(r["generalized_inverses"]) for r in rows),"inverse_count_histogram":{str(k):v for k,v in sorted(hist.items())},"operators":rows,
 "exact_checks":{"all_133_operators_are_regular":all(r["generalized_inverses"] for r in rows),"inverse_count_profile_is_exact":hist=={1:4,4:28,6:16,9:4,16:16,24:32,36:24,54:8,81:1},"inverse_relation_is_symmetric":all(a in rows[b]["generalized_inverses"] for a,r in enumerate(rows) for b in r["generalized_inverses"])},
 "claim_boundary":{"included":"all mutual generalized inverses and regularity of the 133 operators","excluded":"base-monoid generalized inverses already certified, Green relations, physical transfer operators, reliability, or thresholds","parent_issue":"remain open"}}
def validate_artifact(v:Mapping[str,Any]):
 e=build_artifact()
 if v!=e:raise ValueError("operator inverse artifact mismatch")
 if set(e["exact_checks"].values())!={True}:raise ValueError("checks failed")
 return {"schema":SCHEMA,"status":"valid","regular":133}
def main(argv:Optional[Sequence[str]]=None)->int:
 q=argparse.ArgumentParser(description=__doc__);q.add_argument("--output",type=Path);q.add_argument("--validate",type=Path);a=q.parse_args(argv)
 if a.validate:print(json.dumps(validate_artifact(json.loads(a.validate.read_text())),indent=2,sort_keys=True));return 0
 s=json.dumps(build_artifact(),indent=2,sort_keys=True)+"\n"
 if a.output:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(s,encoding="utf-8")
 else:print(s,end="")
 return 0
if __name__=="__main__":raise SystemExit(main())

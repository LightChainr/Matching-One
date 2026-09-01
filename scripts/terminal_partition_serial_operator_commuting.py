#!/usr/bin/env python3
"""Exact center and commuting census of the two-sided operator semigroup."""
from __future__ import annotations
import argparse,json
from collections import Counter
from pathlib import Path
from typing import Any,Mapping,Optional,Sequence
try:
 from scripts.terminal_partition_serial_two_sided_actions import table,operator,compose
except ModuleNotFoundError:
 from terminal_partition_serial_two_sided_actions import table,operator,compose
SCHEMA="matching-one/terminal-partition-serial-operator-commuting/v1"
def data():
 p=table();ops=sorted({operator(a,b,p) for a in range(15) for b in range(15)});ix={f:i for i,f in enumerate(ops)};m=tuple(tuple(ix[compose(f,g)] for g in ops) for f in ops);return ops,m
def build_artifact()->dict[str,Any]:
 _,m=data();n=len(m);cs=[sorted(j for j in range(n) if m[i][j]==m[j][i]) for i in range(n)];hist=Counter(map(len,cs));center=[i for i,c in enumerate(cs) if len(c)==n]
 return {"schema":SCHEMA,"issue":13,"status":"complete_two_sided_operator_commuting_census","center":center,"ordered_commuting_pair_count":sum(map(len,cs)),"centralizer_size_histogram":{str(k):v for k,v in sorted(hist.items())},"centralizers":[{"operator":i,"commuting_operators":c} for i,c in enumerate(cs)],
 "exact_checks":{"all_133_centralizers_cataloged":len(cs)==133,"center_is_identity_only":center==[13],"ordered_commuting_pair_count_is_3239":sum(map(len,cs))==3239,"commutation_is_symmetric":all((j in cs[i])==(i in cs[j]) for i in range(n) for j in range(n))},
 "claim_boundary":{"included":"center and all elementwise centralizers of the 133 two-sided operators","excluded":"base-monoid centralizers, congruences, physical transfer operators, reliability, or thresholds","parent_issue":"remain open"}}
def validate_artifact(v:Mapping[str,Any]):
 e=build_artifact()
 if v!=e:raise ValueError("operator commuting artifact mismatch")
 if set(e["exact_checks"].values())!={True}:raise ValueError("checks failed")
 return {"schema":SCHEMA,"status":"valid","commuting_pairs":3239}
def main(argv:Optional[Sequence[str]]=None)->int:
 q=argparse.ArgumentParser(description=__doc__);q.add_argument("--output",type=Path);q.add_argument("--validate",type=Path);a=q.parse_args(argv)
 if a.validate:print(json.dumps(validate_artifact(json.loads(a.validate.read_text())),indent=2,sort_keys=True));return 0
 s=json.dumps(build_artifact(),indent=2,sort_keys=True)+"\n"
 if a.output:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(s,encoding="utf-8")
 else:print(s,end="")
 return 0
if __name__=="__main__":raise SystemExit(main())

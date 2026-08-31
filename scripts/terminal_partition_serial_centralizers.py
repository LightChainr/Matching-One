#!/usr/bin/env python3
"""Exact elementwise centralizers of the typed serial partition monoid."""
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
SCHEMA="matching-one/terminal-partition-serial-centralizers/v1"
def table():
 st=enumerate_rgs(4);idx={s:i for i,s in enumerate(st)}
 return tuple(tuple(idx[serial_compose(a,b)] for b in st) for a in st)
def build_artifact()->dict[str,Any]:
 p=table();cs=[sorted(y for y in range(15) if p[x][y]==p[y][x]) for x in range(15)];hist=Counter(map(len,cs));center=[x for x,c in enumerate(cs) if len(c)==15]
 return {"schema":SCHEMA,"issue":13,"status":"complete_elementwise_centralizer_catalog","centralizers":[{"element":i,"commuting_elements":c} for i,c in enumerate(cs)],"centralizer_size_histogram":{str(k):v for k,v in sorted(hist.items())},"center":center,"ordered_commuting_pair_count":sum(map(len,cs)),
 "exact_checks":{"all_fifteen_centralizers_cataloged":len(cs)==15,"size_profile_is_exact":hist=={3:6,4:4,5:2,7:2,15:1},"center_is_wire_identity_only":center==[6],"commutation_is_symmetric":all((y in cs[x])==(x in cs[y]) for x in range(15) for y in range(15))},
 "claim_boundary":{"included":"all elementwise centralizers and the center of the 15-state serial monoid","excluded":"operator-semigroup centralizers, congruences, planarity, reliability, or thresholds","parent_issue":"remain open"}}
def validate_artifact(v:Mapping[str,Any]):
 e=build_artifact()
 if v!=e:raise ValueError("centralizer artifact mismatch")
 if set(e["exact_checks"].values())!={True}:raise ValueError("checks failed")
 return {"schema":SCHEMA,"status":"valid","center_size":1}
def main(argv:Optional[Sequence[str]]=None)->int:
 q=argparse.ArgumentParser(description=__doc__);q.add_argument("--output",type=Path);q.add_argument("--validate",type=Path);a=q.parse_args(argv)
 if a.validate:print(json.dumps(validate_artifact(json.loads(a.validate.read_text())),indent=2,sort_keys=True));return 0
 s=json.dumps(build_artifact(),indent=2,sort_keys=True)+"\n"
 if a.output:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(s,encoding="utf-8")
 else:print(s,end="")
 return 0
if __name__=="__main__":raise SystemExit(main())

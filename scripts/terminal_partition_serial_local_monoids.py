#!/usr/bin/env python3
"""Exact local monoids eSe at every idempotent serial partition."""
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
ROOT=Path(__file__).resolve().parents[1];SCHEMA="matching-one/terminal-partition-serial-local-monoids/v1"
def table():
 states=enumerate_rgs(4);index={s:i for i,s in enumerate(states)}
 return tuple(tuple(index[serial_compose(a,b)] for b in states) for a in states)
def local_record(e:int,product):
 if type(e) is not int or not 0<=e<15 or product[e][e]!=e:raise ValueError("local identity must be an idempotent index")
 carrier=tuple(sorted({product[product[e][x]][e] for x in range(15)}))
 idempotents=tuple(x for x in carrier if product[x][x]==x)
 units=tuple(x for x in carrier if any(product[x][y]==e==product[y][x] for y in carrier))
 return {"identity":e,"carrier":list(carrier),"size":len(carrier),"idempotents":list(idempotents),"units":list(units),"commutative":all(product[a][b]==product[b][a] for a in carrier for b in carrier)}
def build_artifact()->dict[str,Any]:
 product=table();idem=[i for i in range(15) if product[i][i]==i];records=[local_record(e,product) for e in idem]
 return {"schema":SCHEMA,"issue":13,"status":"exact_local_monoid_census","idempotent_indices":idem,"local_monoids":records,
  "local_size_histogram":{str(k):v for k,v in sorted(Counter(x["size"] for x in records).items())},
  "local_unit_count_histogram":{str(k):v for k,v in sorted(Counter(len(x["units"]) for x in records).items())},
  "exact_checks":{"local_size_profile_is_four_singletons_seven_pairs_one_global":Counter(x["size"] for x in records)=={1:4,2:7,15:1},
   "all_proper_local_monoids_are_commutative":all(x["commutative"] for x in records if x["size"]<15),
   "only_global_identity_local_has_two_units":[x["identity"] for x in records if len(x["units"])==2]==[6],
   "all_other_local_unit_groups_are_trivial":all(len(x["units"])==1 for x in records if x["identity"]!=6),
   "global_local_monoid_recovers_all_states":next(x for x in records if x["identity"]==6)["carrier"]==list(range(15))},
  "claim_boundary":{"included":"all local monoids eSe, carriers, local units, local idempotents, and commutativity","excluded":"Schutzenberger representations, endomorphism or congruence lattices, planarity, reliability, or thresholds","parent_issue":"remain open"}}
def validate_artifact(value:Mapping[str,Any]):
 expected=build_artifact()
 if value!=expected:raise ValueError("local-monoid artifact does not exactly reproduce")
 if set(expected["exact_checks"].values())!={True}:raise ValueError("all exact checks must pass")
 return {"schema":SCHEMA,"status":"valid","local_monoids":12}
def main(argv:Optional[Sequence[str]]=None)->int:
 p=argparse.ArgumentParser(description=__doc__);p.add_argument("--output",type=Path);p.add_argument("--validate",type=Path);a=p.parse_args(argv)
 if a.validate:print(json.dumps(validate_artifact(json.loads(a.validate.read_text())),indent=2,sort_keys=True));return 0
 rendered=json.dumps(build_artifact(),indent=2,sort_keys=True)+"\n"
 if a.output:a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(rendered,encoding="utf-8")
 else:print(rendered,end="")
 return 0
if __name__=="__main__":raise SystemExit(main())

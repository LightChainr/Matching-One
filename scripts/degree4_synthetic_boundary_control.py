#!/usr/bin/env python3
"""Sensitivity control for the quartic census where it returned a null.

The degree-4 census returns zero survivors on the two narrowest method
intervals.  Those nulls are only informative if the pipeline would have found
a height-100 quartic root had one been present at that width.  This control
plants a root that is known to exist -- a quartic root witness already
committed by the census itself -- inside a synthetic interval of each such
width, and runs the unmodified census path
(``degree4_interval_exclusion.run_search``: certified fixed-point screen,
exact rational endpoints, Sturm isolation) on it.

The control deliberately covers only the widths where the census reported
zero survivors, and reads that set from the census artifacts rather than
naming it.  On the two wider intervals the census itself found 1 and 15
roots, so its sensitivity there is already demonstrated and a planted root
would add nothing.

Each planted root is tested at both polarities:

* **positive** -- the interval contains the planted root; the census must
  report that quartic among its root witnesses;
* **negative** -- the interval is shifted one full width away, so the planted
  root is certified outside; the census must not report that quartic.

The control therefore certifies detection sensitivity and boundary
specificity exactly where the paper rests on a null, using no new mathematics
and no modification of the census code.
"""

from __future__ import annotations
import argparse, hashlib, json
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

try:
    from scripts.exact_polynomial_root_certificate import isolate_roots
    from scripts.degree4_interval_exclusion import run_search
except ModuleNotFoundError:
    from exact_polynomial_root_certificate import isolate_roots
    from degree4_interval_exclusion import run_search

ROOT=Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT=ROOT/"analysis"/"pslq_search_contract.json"
DEFAULT_OUTPUT=ROOT/"results"/"pslq-degree4-synthetic-boundary-control"/"latest.json"
SCHEMA="matching-one/degree4-synthetic-boundary-control/v1"
ISSUE=551
ISOLATION_BITS=160
SEARCH_WINDOW=(Fraction(55,100),Fraction(65,100))

# Planted quartics are committed census root witnesses, not new constructions.
# Each is named with the artifact that already certifies its root, and is
# placed at a different offset inside the synthetic interval so the control
# does not depend on one symmetric placement.
PLANTED=(
    {"coefficients_ascending":(-84,99,-7,99,58),
     "committed_by":"results/pslq-degree4-mertens-2022-p-med/latest.json",
     "placement_fraction":Fraction(1,3)},
    {"coefficients_ascending":(-48,31,43,62,13),
     "committed_by":"results/pslq-degree4-mertens-2022-p-cell/latest.json",
     "placement_fraction":Fraction(2,3)},
)


def _require(condition:bool,message:str)->None:
    if not condition:raise ValueError(message)


def _text(value:Fraction)->str:
    return f"{value.numerator}/{value.denominator}"


def _digest(path:Path)->str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _planted_root(coefficients:Sequence[int])->tuple[Fraction,Fraction]:
    polynomial=[Fraction(value) for value in coefficients]
    roots=isolate_roots(polynomial,SEARCH_WINDOW[0],SEARCH_WINDOW[1],bits=ISOLATION_BITS)
    _require(len(roots)==1,"planted quartic must have exactly one root in the search window")
    return roots[0]


def _committed_witness(relative:str,coefficients:Sequence[int])->None:
    payload=json.loads((ROOT/relative).read_text(encoding="utf-8"))
    committed={tuple(row["coefficients_ascending"]) for row in payload["interval_result"]["root_witnesses"]}
    _require(tuple(coefficients) in committed,f"planted quartic is not a committed witness of {relative}")


def _run(interval_id:str,lower:Fraction,upper:Fraction)->dict[str,Any]:
    return run_search({"id":interval_id,"source_id":"synthetic-boundary-control",
                       "lower":_text(lower),"upper":_text(upper)})


def _trial(planted:Mapping[str,Any],bracket:tuple[Fraction,Fraction],width_id:str,width:Fraction,
           polarity:str)->dict[str,Any]:
    low,high=bracket
    if polarity=="positive":
        lower=low-width*planted["placement_fraction"]
        upper=lower+width
        _require(lower<low and high<upper,"planted root must lie strictly inside the positive interval")
    else:
        # shift a full width above the root, so the root is outside by >= width.
        lower=high+width
        upper=lower+width
        _require(high<lower,"planted root must lie strictly outside the negative interval")
    result=_run(f"synthetic-{width_id}-{polarity}",lower,upper)
    witnesses=[tuple(row["coefficients_ascending"]) for row in result["root_witnesses"]]
    detected=tuple(planted["coefficients_ascending"]) in witnesses
    return {
        "planted_coefficients_ascending":list(planted["coefficients_ascending"]),
        "width_id":width_id,"width_text":_text(width),"polarity":polarity,
        "lower":_text(lower),"upper":_text(upper),
        "planted_root_bracket":[_text(low),_text(high)],
        "placement_fraction":_text(planted["placement_fraction"]),
        "near_candidates_exactly_checked":result["near_candidates_exactly_checked"],
        "root_filter_candidates":result["root_filter_candidates"],
        "root_containing_polynomials":result["root_containing_polynomials"],
        "witness_coefficients":[list(row) for row in witnesses],
        "planted_quartic_detected":detected,
        "expected_detection":polarity=="positive",
        "passed":detected==(polarity=="positive"),
    }


def build_result(contract_path:Path=DEFAULT_CONTRACT)->dict[str,Any]:
    return _build_result(Path(contract_path).resolve())


@lru_cache(maxsize=2)
def _build_result(contract_path:Path)->dict[str,Any]:
    raw=contract_path.read_bytes()
    contract=json.loads(raw)
    provenance=contract["provenance"]
    provenance_digest=_digest(ROOT/provenance["path"])
    _require(provenance_digest==provenance["sha256"],"provenance digest drift")
    widths={row["id"]:Fraction(row["upper"])-Fraction(row["lower"]) for row in contract["intervals"]}
    # Cover exactly the intervals whose census result was a null.
    null_widths={}
    for interval_id in sorted(widths):
        census=json.loads((ROOT/f"results/pslq-degree4-{interval_id}/latest.json").read_text(encoding="utf-8"))
        if census["interval_result"]["root_containing_polynomials"]==0:
            null_widths[interval_id]=widths[interval_id]
    _require(bool(null_widths),"no census interval returned a null to certify")

    targets,trials=[],[]
    for planted in PLANTED:
        _committed_witness(planted["committed_by"],planted["coefficients_ascending"])
        bracket=_planted_root(planted["coefficients_ascending"])
        targets.append({
            "coefficients_ascending":list(planted["coefficients_ascending"]),
            "height":max(abs(value) for value in planted["coefficients_ascending"]),
            "committed_by":planted["committed_by"],
            "committed_by_sha256":_digest(ROOT/planted["committed_by"]),
            "root_bracket":[_text(value) for value in bracket],
            "isolation_bits":ISOLATION_BITS,
        })
        for width_id in sorted(null_widths):
            for polarity in ("positive","negative"):
                trials.append(_trial(planted,bracket,width_id,null_widths[width_id],polarity))

    positive=[row for row in trials if row["polarity"]=="positive"]
    negative=[row for row in trials if row["polarity"]=="negative"]
    return {
        "schema":SCHEMA,"issue":ISSUE,"status":"degree4_synthetic_boundary_control_complete",
        "contract_sha256":hashlib.sha256(raw).hexdigest(),"provenance_sha256":provenance_digest,
        "census_path":{
            "driver":"scripts/degree4_interval_exclusion.py::run_search",
            "screen":"scripts/degree4_fixed_point_screen.cpp",
            "modified_for_this_control":False,
        },
        "planted_targets":targets,
        "widths_tested":{width_id:_text(width) for width_id,width in sorted(null_widths.items())},
        "width_selection":{
            "rule":"every method interval whose committed degree-4 census reported zero root-containing polynomials",
            "covered":sorted(null_widths),
            "not_covered":sorted(set(widths)-set(null_widths)),
            "reason":"the census itself found roots on the uncovered intervals, so its sensitivity there is "
                     "already demonstrated",
        },
        "trials":trials,
        "conclusion":{
            "positive_trials":len(positive),
            "negative_trials":len(negative),
            "all_positive_trials_detected_the_planted_quartic":all(row["planted_quartic_detected"] for row in positive),
            "all_negative_trials_excluded_the_planted_quartic":not any(row["planted_quartic_detected"] for row in negative),
            "all_trials_passed":all(row["passed"] for row in trials),
            "meaning":"at every width where the census reported a null, the census path detects a height-100 "
                      "quartic root that is present and does not report one that is absent, so those nulls are "
                      "sensitivity-certified rather than blind spots",
        },
        "claim_boundary":{
            "included":"detection sensitivity and boundary specificity of the committed quartic census path at "
                       "the frozen method widths whose census result was a null",
            "excluded":"new census results, degree/height expansion, a null distribution, a p-value, a claim "
                       "about the planted quartics as candidate formulas, closed forms, or transcendence",
            "parent_issue":"remain open",
        },
    }


def validate_result(result:Mapping[str,Any],contract_path:Path=DEFAULT_CONTRACT)->Mapping[str,Any]:
    expected=build_result(contract_path)
    if result!=expected:raise ValueError("degree-4 boundary control does not exactly reproduce")
    return {"schema":SCHEMA,"status":"valid","all_trials_passed":expected["conclusion"]["all_trials_passed"]}


def main(argv:Optional[Sequence[str]]=None)->int:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract",type=Path,default=DEFAULT_CONTRACT)
    parser.add_argument("--output",type=Path)
    parser.add_argument("--validate",type=Path)
    args=parser.parse_args(argv)
    if args.validate:
        print(json.dumps(validate_result(json.loads(args.validate.read_text(encoding="utf-8")),args.contract),indent=2,sort_keys=True))
        return 0
    rendered=json.dumps(build_result(args.contract),indent=2,sort_keys=True)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(rendered,encoding="utf-8")
    else:
        print(rendered,end="")
    return 0


if __name__=="__main__":raise SystemExit(main())

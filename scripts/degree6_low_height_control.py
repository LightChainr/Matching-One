#!/usr/bin/env python3
"""Sensitivity control for the degree-6 height-3 exclusion.

The exclusion in ``degree6_low_height_exclusion`` returns zero screened
candidates on every frozen interval: the certified derivative screen discards
the whole class before any Sturm decision runs.  A null produced that way is
only informative if the same path would have reported a root had one been
present at that width, so this control plants one.

The planted polynomial is not a construction.  It is ``x^6 - 3x^4 + 1``, the
minimal polynomial of the (3,12^2) site threshold
``sqrt(1 - 2 sin(pi/18))``, taken from the committed lattice-native candidate
artifact -- the exact form whose absence from the degree-4 height-100 class
motivated the exclusion in the first place.  It is planted inside a synthetic
interval of each frozen method width, at both polarities:

* **positive** -- the interval contains the planted root; the scan must report
  that polynomial among its root witnesses;
* **negative** -- the interval is shifted one full width away, so the root is
  certified outside; the scan must not report it.
"""

from __future__ import annotations
import argparse, hashlib, json
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

try:
    from scripts.exact_polynomial_root_certificate import isolate_roots
    from scripts.degree6_low_height_exclusion import HEIGHT, _scan_degree, derivative_bound
except ModuleNotFoundError:
    from exact_polynomial_root_certificate import isolate_roots
    from degree6_low_height_exclusion import HEIGHT, _scan_degree, derivative_bound

ROOT=Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT=ROOT/"analysis"/"pslq_search_contract.json"
DEFAULT_OUTPUT=ROOT/"results"/"pslq-degree6-low-height-control"/"latest.json"
LATTICE_NATIVE="results/pslq-lattice-native-candidates/latest.json"
SCHEMA="matching-one/degree6-low-height-control/v1"
ISSUE=559
ISOLATION_BITS=160
PLANTED_CANDIDATE_ID="three-twelve-site"
PLACEMENT=Fraction(1,3)
SEARCH_WINDOW=(Fraction(1,2),Fraction(9,10))


def _require(condition:bool,message:str)->None:
    if not condition:raise ValueError(message)


def _text(value:Fraction)->str:
    return f"{value.numerator}/{value.denominator}"


def _digest(path:Path)->str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def planted_polynomial()->tuple[int,...]:
    """The (3,12^2) minimal polynomial, read from the certified artifact."""
    payload=json.loads((ROOT/LATTICE_NATIVE).read_text(encoding="utf-8"))
    rows=[row for row in payload["candidates"] if row["candidate_id"]==PLANTED_CANDIDATE_ID]
    _require(len(rows)==1,f"{PLANTED_CANDIDATE_ID} is not uniquely present in the lattice-native artifact")
    coefficients=tuple(rows[0]["minimal_polynomial_coefficients_ascending"])
    _require(len(coefficients)-1==6,"planted polynomial must have degree 6")
    _require(max(abs(value) for value in coefficients)<=HEIGHT,
             "planted polynomial must lie inside the height bound being certified")
    return coefficients


def planted_root(coefficients:Sequence[int])->tuple[Fraction,Fraction]:
    polynomial=[Fraction(value) for value in coefficients]
    roots=isolate_roots(polynomial,SEARCH_WINDOW[0],SEARCH_WINDOW[1],bits=ISOLATION_BITS)
    _require(len(roots)==1,"planted polynomial must have exactly one root in the search window")
    return roots[0]


def _trial(coefficients:Sequence[int],bracket:tuple[Fraction,Fraction],width_id:str,
           width:Fraction,polarity:str)->dict[str,Any]:
    low,high=bracket
    if polarity=="positive":
        lower=low-width*PLACEMENT
        upper=lower+width
        _require(lower<low and high<upper,"planted root must lie strictly inside the positive interval")
    else:
        lower=high+width
        upper=lower+width
        _require(high<lower,"planted root must lie strictly outside the negative interval")
    scan=_scan_degree(6,lower,upper)
    witnesses=[tuple(row["coefficients_ascending"]) for row in scan["root_witnesses"]]
    detected=tuple(coefficients) in witnesses
    return {
        "width_id":width_id,"width_text":_text(width),"polarity":polarity,
        "lower":_text(lower),"upper":_text(upper),
        "screen_candidates_exactly_decided":scan["screen_candidates_exactly_decided"],
        "root_containing_polynomials":scan["root_containing_polynomials"],
        "witness_coefficients":[list(row) for row in witnesses],
        "planted_polynomial_detected":detected,
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
    coefficients=planted_polynomial()
    bracket=planted_root(coefficients)
    widths={row["id"]:Fraction(row["upper"])-Fraction(row["lower"]) for row in contract["intervals"]}

    trials=[]
    for width_id in sorted(widths):
        for polarity in ("positive","negative"):
            trials.append(_trial(coefficients,bracket,width_id,widths[width_id],polarity))
    positive=[row for row in trials if row["polarity"]=="positive"]
    negative=[row for row in trials if row["polarity"]=="negative"]
    return {
        "schema":SCHEMA,"issue":ISSUE,"status":"degree6_low_height_control_complete",
        "contract_sha256":hashlib.sha256(raw).hexdigest(),"provenance_sha256":provenance_digest,
        "scan_path":{"function":"scripts/degree6_low_height_exclusion.py::_scan_degree",
                     "modified_for_this_control":False},
        "planted":{
            "candidate_id":PLANTED_CANDIDATE_ID,
            "closed_form":"sqrt(1-2*sin(pi/18))",
            "coefficients_ascending":list(coefficients),
            "degree":len(coefficients)-1,
            "height":max(abs(value) for value in coefficients),
            "derivative_bound_on_unit_interval":derivative_bound(len(coefficients)-1),
            "source":LATTICE_NATIVE,
            "source_sha256":_digest(ROOT/LATTICE_NATIVE),
            "root_bracket":[_text(bracket[0]),_text(bracket[1])],
            "isolation_bits":ISOLATION_BITS,
            "placement_fraction":_text(PLACEMENT),
        },
        "widths_tested":{width_id:_text(width) for width_id,width in sorted(widths.items())},
        "trials":trials,
        "conclusion":{
            "positive_trials":len(positive),
            "negative_trials":len(negative),
            "all_positive_trials_detected_the_planted_polynomial":all(row["planted_polynomial_detected"] for row in positive),
            "all_negative_trials_excluded_the_planted_polynomial":not any(row["planted_polynomial_detected"] for row in negative),
            "all_trials_passed":all(row["passed"] for row in trials),
            "meaning":"at every frozen method width the degree-6 scan reports a height-3 degree-6 root that is "
                      "present and does not report one that is absent, so the zero-candidate exclusion is a "
                      "sensitivity-certified null rather than a blind screen",
        },
        "claim_boundary":{
            "included":"detection sensitivity and boundary specificity of the degree-6 height-3 scan at the "
                       "four frozen method widths",
            "excluded":"new exclusion results, degree/height expansion, a null distribution, a p-value, any "
                       "claim that the planted (3,12^2) form is a candidate for square-site p_c, closed forms, "
                       "or transcendence",
            "parent_issue":"remain open",
        },
    }


def validate_result(result:Mapping[str,Any],contract_path:Path=DEFAULT_CONTRACT)->Mapping[str,Any]:
    expected=build_result(contract_path)
    if result!=expected:raise ValueError("degree-6 low-height control does not exactly reproduce")
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

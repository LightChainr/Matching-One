#!/usr/bin/env python3
"""Cross-check two independent implementations of the C(1..6, 3) census.

``degree6_low_height_exclusion`` and ``degree6_independent_replication`` decide
the same question over the same frozen class on the same frozen intervals, and
were written separately.  They differ where a shared bug would have to hide:

* **enumeration** -- each generates and sign-normalises ``C(d,3)`` in its own
  code, and each checks its count against the committed
  ``primitive_polynomial_count``;
* **screen** -- one rules a polynomial out when ``|P(l)| > D*(u-l)`` at the left
  endpoint, the other when ``|P(m)| > D*(u-l)/2`` at the midpoint.  Both are
  certified consequences of ``|P'| <= D`` on ``[0,1]``, but they evaluate
  different points and compare against different bounds.

Both import ``exact_polynomial_root_certificate`` unchanged, so the Sturm path
is shared rather than replicated.  It contributes nothing to the null: the
screens retain zero candidates everywhere, so root isolation never runs during
either census, and the agreement below is between two independently written
certified screens.

What is compared:

1. every one of the ``4 x 6`` interval-by-degree cells, on class size, screen
   survivors, root-containing polynomials and distinct roots;
2. the per-interval exclusion verdict;
3. the closest member of the whole class, coefficient by coefficient;
4. its two residuals against each other.  One implementation reports
   ``min(|P(l)|, |P(u)|)`` and the other ``|P(m)|``; for a single polynomial the
   mean-value theorem forces these to agree within ``D*(u-l)/2``, where ``D`` is
   that polynomial's own ``sum k|a_k|``.  This is the one check that would fail
   if the two agreed on the counts by coincidence rather than by computing the
   same thing.

All arithmetic is exact; the residual comparison is over ``Fraction``.
"""

from __future__ import annotations
import argparse, hashlib, json
from fractions import Fraction
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

ROOT=Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT=ROOT/"analysis"/"pslq_search_contract.json"
OUTPUT=ROOT/"results"/"pslq-degree6-implementation-agreement"/"latest.json"
SCHEMA="matching-one/degree6-implementation-agreement/v1"
ISSUE=559
PRIMARY="scripts/degree6_low_height_exclusion.py"
REPLICATION="scripts/degree6_independent_replication.py"


def _require(condition:bool,message:str)->None:
    if not condition:raise ValueError(message)


def _text(value:Fraction)->str:
    return f"{value.numerator}/{value.denominator}"


def _digest(path:Path)->str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load(directory:str)->dict[str,Any]:
    path=ROOT/"results"/directory/"latest.json"
    _require(path.is_file(),f"missing artifact {directory}")
    return json.loads(path.read_text(encoding="utf-8"))


def derivative_bound(coefficients:Sequence[int])->int:
    """sum k|a_k|, which bounds |P'| on [0,1] for this polynomial."""
    return sum(index*abs(value) for index,value in enumerate(coefficients))


def _cells(primary:Mapping[str,Any],replication:Mapping[str,Any])->list[dict[str,Any]]:
    left={row["degree"]:row for row in primary["by_degree"]}
    right={row["degree"]:row for row in replication["by_degree"]}
    _require(sorted(left)==sorted(right)==list(range(1,7)),"degree coverage differs")
    rows=[]
    for degree in range(1,7):
        a,b=left[degree],right[degree]
        pairs={
            "polynomials_in_class":(a["polynomials_in_class"],b["class_size"]),
            "screen_survivors":(a["screen_candidates_exactly_decided"],b["screen_survivors"]),
            "root_containing_polynomials":(a["root_containing_polynomials"],b["root_containing_polynomials"]),
            "distinct_roots_in_interval":(a["distinct_roots_in_interval"],b["distinct_roots_in_interval"]),
        }
        disagreements=sorted(name for name,(x,y) in pairs.items() if x!=y)
        rows.append({"degree":degree,
                     "values":{name:x for name,(x,_) in pairs.items()},
                     "fields_compared":sorted(pairs),
                     "agree":not disagreements,
                     "disagreements":disagreements})
    return rows


def _closest(primary:Mapping[str,Any],replication:Mapping[str,Any],
             lower:Fraction,upper:Fraction)->dict[str,Any]:
    best=min(primary["by_degree"],
             key=lambda row:Fraction(row["closest_polynomial"]["minimum_absolute_endpoint_residual"]))
    mine=best["closest_polynomial"]
    theirs=replication["closest_polynomial_at_midpoint"]
    coefficients=list(mine["coefficients_ascending"])
    same=coefficients==list(theirs["coefficients_ascending"])
    endpoint=Fraction(mine["minimum_absolute_endpoint_residual"])
    midpoint=Fraction(theirs["minimum_absolute_residual_at_midpoint"])
    bound=derivative_bound(coefficients)
    allowed=Fraction(bound)*(upper-lower)/2
    gap=abs(endpoint-midpoint)
    return {
        "coefficients_ascending":coefficients,
        "coefficients_agree":same,
        "endpoint_residual_text":_text(endpoint),
        "midpoint_residual_text":_text(midpoint),
        "residual_gap_text":_text(gap),
        "mean_value_allowance_text":_text(allowed),
        "polynomial_derivative_bound":bound,
        "within_mean_value_bound":gap<=allowed,
    }


def build_result(contract_path:Path=DEFAULT_CONTRACT)->dict[str,Any]:
    raw=Path(contract_path).read_bytes()
    contract=json.loads(raw)
    provenance=contract["provenance"]
    provenance_digest=_digest(ROOT/provenance["path"])
    _require(provenance_digest==provenance["sha256"],"provenance digest drift")

    intervals=[]
    for row in contract["intervals"]:
        interval_id=row["id"]
        lower,upper=Fraction(row["lower"]),Fraction(row["upper"])
        primary=_load(f"pslq-degree6-low-height-{interval_id}")
        replication=_load(f"pslq-degree6-low-height-replication-{interval_id}")
        left,right=primary["interval_result"],replication["interval_result"]
        _require(left["interval_id"]==right["interval_id"]==interval_id,"artifact interval mismatch")
        for artifact in (primary,replication):
            _require(artifact["contract_sha256"]==hashlib.sha256(raw).hexdigest(),
                     f"{interval_id}: artifact was built against a different contract")
            _require(artifact["provenance_sha256"]==provenance_digest,
                     f"{interval_id}: artifact was built against a different provenance manifest")
        cells=_cells(left,right)
        intervals.append({
            "interval_id":interval_id,
            "width_text":_text(upper-lower),
            "by_degree":cells,
            "cells_compared":len(cells),
            "cells_in_agreement":sum(1 for cell in cells if cell["agree"]),
            "exclusion_verdict":left["excluded"],
            "exclusion_verdicts_agree":left["excluded"]==right["excluded"],
            "closest_member":_closest(left,right,lower,upper),
        })

    cells_compared=sum(row["cells_compared"] for row in intervals)
    cells_agreed=sum(row["cells_in_agreement"] for row in intervals)
    return {
        "schema":SCHEMA,"issue":ISSUE,"status":"degree6_implementation_agreement_complete",
        "contract_sha256":hashlib.sha256(raw).hexdigest(),"provenance_sha256":provenance_digest,
        "implementations":[
            {"path":PRIMARY,"sha256":_digest(ROOT/PRIMARY),
             "screen":"exact integer evaluation at both endpoints; excludes when |P(l)| > D*(u-l)"},
            {"path":REPLICATION,"sha256":_digest(ROOT/REPLICATION),
             "screen":"exact rational evaluation at the midpoint; excludes when |P(m)| > D*(u-l)/2"},
        ],
        "shared_code":{
            "path":"scripts/exact_polynomial_root_certificate.py",
            "sha256":_digest(ROOT/"scripts"/"exact_polynomial_root_certificate.py"),
            "role":"exact Sturm isolation, imported unchanged by both",
            "contribution_to_this_result":"none: both screens retain zero candidates at every "
                                          "degree on every interval, so root isolation never runs "
                                          "during either census",
        },
        "fields_compared_per_cell":["polynomials_in_class","screen_survivors",
                                    "root_containing_polynomials","distinct_roots_in_interval"],
        "cells_compared":cells_compared,
        "cells_in_agreement":cells_agreed,
        "intervals":intervals,
        "implementations_agree":(cells_agreed==cells_compared
                                 and all(row["exclusion_verdicts_agree"] for row in intervals)
                                 and all(row["closest_member"]["coefficients_agree"] for row in intervals)
                                 and all(row["closest_member"]["within_mean_value_bound"] for row in intervals)),
        "claim_boundary":{
            "included":"the degree-1..6 height-3 census of section 6.5 is reproduced by a second "
                       "independently written enumeration and certified screen",
            "excluded":"the degree-4 height-100 census, whose C++ meet-in-the-middle screen and "
                       "retained-candidate Sturm decisions have no second implementation; also "
                       "higher degree or height, closed forms, or transcendence",
            "parent_issue":"remain open",
        },
    }


def validate_result(result:Mapping[str,Any],contract_path:Path=DEFAULT_CONTRACT)->Mapping[str,Any]:
    expected=build_result(contract_path)
    if result!=expected:raise ValueError("implementation agreement does not exactly reproduce")
    return {"schema":SCHEMA,"status":"valid",
            "cells_compared":expected["cells_compared"],
            "cells_in_agreement":expected["cells_in_agreement"],
            "implementations_agree":expected["implementations_agree"]}


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
    destination=args.output or OUTPUT
    destination.parent.mkdir(parents=True,exist_ok=True)
    destination.write_text(rendered,encoding="utf-8")
    print(f"wrote {destination.relative_to(ROOT)}")
    return 0


if __name__=="__main__":raise SystemExit(main())

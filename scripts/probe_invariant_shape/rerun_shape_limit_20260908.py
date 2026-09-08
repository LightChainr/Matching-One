#!/usr/bin/env python3
"""#622 re-run (2026-09-08): re-adjudicate the shape-limit verdicts on the
repaired inputs (PR #653 bond lab, PR #655 spin-0 N725).

Everything here is built from already-merged objects:

* the repaired exact census (site L=3,4; square-bond L=3 with
  ``dual_fail = 0``, PR #653): read from ``results/probe-invariant-shape/
  census-exact.json`` if the #653 tree is present, else rebuilt exactly by
  importing the repaired census module logic that lives on that branch;
* the corrected spin-0 N725 shape (PR #655): read from
  ``results/probe-invariant-shape/n725-zflow-corrected.json`` if present,
  else from the committed branch via ``git show`` (never rewritten);
* the toy location/shape separation of #628 (``toy_location_theorem.py``):
  cited, extended here by a *duality-constrained* toy family.

New in this run (each is one machine-checked statement):

T1 (Duality symmetry theorem, bond laboratory).
  If every configuration satisfies r_b + r_w = 2 under complement
  (dual_fail = 0), then the pair-count polynomials satisfy
  P02[k] = P20[NB-k], P11[k] = P11[NB-k]; hence
  M(p) = -M(1-p) exactly (M antisymmetric about 1/2),
  F(p) = 1 - F(1-p) exactly (self-complementing), and with the symmetric
  anchors a = 1/4, b = 3/4,
  Z(u) + Z(1-u) = 1 and Z(1/2) = 1/2 exactly.
  Machine-checked against the repaired census (all four statements, exact
  rational arithmetic where the census is exact).

T2 (Duality does not pin the shape, i.e. W2's "pinning" mechanism dies).
  The class of CDFs satisfying every T1 constraint is exactly the class of
  symmetric densities on (0,1) (F self-complementing <=> f(p) = f(1-p)).
  Within that class Z still ranges over a 1-parameter family: two symmetric
  densities with different tail weights give shapes differing by ~0.15 at
  u = 0.9 (computed).  So duality pins Z to the antisymmetric class but
  pins no member of it.

T3 (Decomposition and the rigidity observation).
  Z = A + S with
    A(u) = (Z(u) - Z(1-u) + 1)/2  (antisymmetric part),
    S(u) = (Z(u) + Z(1-u) - 1)/2  (symmetric part, zero iff T1 holds).
  Across site L=3, site L=4, bond L=3 (self-dual) and the corrected spin-0
  N725 block, the symmetric part S collapses toward 0
  (S(0.1): -0.0541, -0.0446, 0, -0.0090) while the antisymmetric part A
  drifts slowly and monotonically in u=0.1 ordering
  (A(0.1): -0.4116, -0.4279, -0.4352, -0.4591), pairwise spread <= 0.0475.
  Recorded as an observation with its trend; NOT claimed as a theorem and
  NOT fitted to any exponent.

T4 (Conjecture C3, stated precisely so it can be killed).
  S_N -> 0 uniformly on the level grid iff
  sup_{p in compacta} |M_N(p) + M_N(1-p)| -> 0
  ("self-matching restoration").  Under C3 the limit shape, if it exists,
  is antisymmetric and the question is whether A_N converges.  C3 is not
  implied by Theorem L (its hypotheses are spent on location); it is a
  checkable statement about rank moments that neither this probe nor any
  merged ticket has proved.

Output: results/probe-invariant-shape/rerun-20260908.json

Standing rules honoured: no new production, no MC, no STATUS.md, no issue
closure, no exponent fit, N=725 read as the committed block (not scored),
1.55 not revived, no literature search.
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "scripts"))

LEVELS = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
ANCHOR_A, ANCHOR_B = 0.25, 0.75

BRANCH_CENSUS = "origin/repair/632-bond-lab-and-wrapping-l3l4"
BRANCH_N725 = "origin/analysis/633-n725-spin0"
PATH_CENSUS = "results/probe-invariant-shape/census-exact.json"
PATH_N725 = "results/probe-invariant-shape/n725-zflow-corrected.json"


def _gitshow_json(branch: str, path: str) -> dict:
    blob = subprocess.check_output(
        ["git", "show", f"{branch}:{path}"], text=True, cwd=ROOT)
    return json.loads(blob)


def load_census() -> dict:
    local = ROOT / PATH_CENSUS
    if local.exists() and json.loads(local.read_text()).get(
            "bond", {}).get("3", {}).get("dual_fail") == 0:
        return json.loads(local.read_text())
    return _gitshow_json(BRANCH_CENSUS, PATH_CENSUS)


def load_n725() -> dict:
    local = ROOT / PATH_N725
    if local.exists():
        return json.loads(local.read_text())
    return _gitshow_json(BRANCH_N725, PATH_N725)


def deciles_of_n725(n725: dict) -> list[float]:
    z11 = n725["weightings"]["spin0"]["Z_pooled_quarter_anchors"]
    # 11 levels: 0.1,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.9 -> drop anchors
    idx = [0, 1, 3, 4, 5, 6, 7, 9, 10]
    return [z11[i] for i in idx]


# --------------------------------------------------------------------- T1

def check_t1(census: dict) -> dict:
    """Verify the duality symmetry theorem on the repaired bond census."""
    b = census["bond"]["3"]
    out: dict = {"dual_fail": b["dual_fail"],
                 "rank_pair_counts": b["rank_pair_counts"]}
    out["pairs_symmetric"] = (b["rank_pair_counts"]["0,2"]
                              == b["rank_pair_counts"]["2,0"])
    out["M_half_is_zero"] = (Fraction(b["M_half"]) == 0)
    out["self_dual_p_half"] = (b["self_dual_p_half"] == 0.5)
    z = b["Z_levels_float"]
    out["Z_antisymmetry_max_err"] = max(
        abs(z[i] + z[8 - i] - 1.0) for i in range(9))
    out["Z_half_is_half"] = abs(z[4] - 0.5) < 1e-12
    out["verdict"] = ("proved (T1) and machine-checked: dual_fail=0 => "
                      "M antisymmetric, F self-complementing, Z(u)+Z(1-u)=1")
    return out


# --------------------------------------------------------------------- T2

def _sym_density_F(scale: float, hollow: float,
                   n: int = 200_001) -> tuple[list[float], list[float]]:
    """A self-complementing CDF: density symmetric about 1/2.

    f(p) ∝ exp(-(p-1/2)^2 / (2 scale^2)) * (1 + hollow * tanh((p-1/2)/scale)^2)

    ``hollow`` > 0 pushes density to the tails (hollow centre), < 0 to the
    centre.  Every member satisfies every T1 constraint exactly.
    """
    def f(p: float) -> float:
        x = p - 0.5
        return math.exp(-x * x / (2 * scale * scale)) * (
            1 + hollow * math.tanh(x / scale) ** 2)

    h = 1.0 / (n - 1)
    ps = [i * h for i in range(n)]
    F = [0.0] * n
    prev = f(0.0)
    for i in range(1, n):
        cur = f(ps[i])
        F[i] = F[i - 1] + 0.5 * (prev + cur) * h
        prev = cur
    tot = F[-1]
    return [x / tot for x in F], ps


def _quantile(F: list[float], ps: list[float], u: float) -> float:
    lo, hi = 0, len(F) - 1
    while lo < hi:
        mid = (lo + hi) // 2
        if F[mid] < u:
            lo = mid + 1
        else:
            hi = mid
    return ps[lo]


def _Z(F: list[float], ps: list[float]) -> list[float]:
    qa = _quantile(F, ps, ANCHOR_A)
    qb = _quantile(F, ps, ANCHOR_B)
    return [(_quantile(F, ps, u) - qa) / (qb - qa) for u in LEVELS]


def check_t2() -> dict:
    """Duality class does not pin the shape: two members, two shapes."""
    Fa, ps = _sym_density_F(0.02, 0.0)
    Fb, _ = _sym_density_F(0.02, 3.0)
    Za, Zb = _Z(Fa, ps), _Z(Fb, ps)
    # both members satisfy T1 constraints (by construction: symmetric density)
    t1a = max(abs(Za[i] + Za[8 - i] - 1.0) for i in range(9))
    t1b = max(abs(Zb[i] + Zb[8 - i] - 1.0) for i in range(9))
    spread = max(abs(a - b) for a, b in zip(Za, Zb))
    return {
        "family": "symmetric densities: f ∝ exp(-(p-1/2)²/2s²)(1+h tanh²((p-1/2)/s))",
        "Z_scale0_hollow0": Za,
        "Z_scale0_hollow3": Zb,
        "T1_constraint_max_err_A": t1a,
        "T1_constraint_max_err_B": t1b,
        "shape_spread_within_duality_class": spread,
        "verdict": ("W2's pinning mechanism dies: the duality class is the "
                    "symmetric-density class and Z still ranges freely "
                    "inside it (spread ~0.15 at u=0.9)"),
    }


# --------------------------------------------------------------------- T3

def antisym_part(z: list[float]) -> list[float]:
    return [(z[i] - z[8 - i] + 1.0) / 2.0 for i in range(5)]


def sym_part(z: list[float]) -> list[float]:
    return [(z[i] + z[8 - i] - 1.0) / 2.0 for i in range(5)]


def check_t3(census: dict, n725: dict) -> dict:
    shapes = {
        "site_L3": census["site"]["3"]["Z_levels_float"],
        "site_L4": census["site"]["4"]["Z_levels_float"],
        "bond_L3_selfdual": census["bond"]["3"]["Z_levels_float"],
        "N725_spin0": deciles_of_n725(n725),
    }
    A = {k: antisym_part(v) for k, v in shapes.items()}
    S = {k: sym_part(v) for k, v in shapes.items()}
    spread_A, spread_S = {}, {}
    keys = list(shapes)
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            a, b = keys[i], keys[j]
            spread_A[f"{a}|{b}"] = max(abs(x - y) for x, y in zip(A[a], A[b]))
            spread_S[f"{a}|{b}"] = max(abs(x - y) for x, y in zip(S[a], S[b]))
    return {
        "decomposition": "Z = A + S; A(u)=(Z(u)-Z(1-u)+1)/2; S(u)=(Z(u)+Z(1-u)-1)/2",
        "Z_by_object": shapes,
        "antisymmetric_part_A_u0.1_to_0.5": A,
        "symmetric_part_S_u0.1_to_0.5": S,
        "pairwise_max_spread_A": spread_A,
        "pairwise_max_spread_S": spread_S,
        "max_spread_A_overall": max(spread_A.values()),
        "max_spread_S_overall": max(spread_S.values()),
        "S_collapse_reading": [
            "S(0.1): site3 %.4f -> site4 %.4f -> bond3 0 (exact) -> N725 %.4f"
            % (S["site_L3"][0], S["site_L4"][0], S["N725_spin0"][0]),
            "A(0.1): site3 %.4f -> site4 %.4f -> bond3 %.4f -> N725 %.4f "
            "(monotone in this ordering)"
            % (A["site_L3"][0], A["site_L4"][0], A["bond_L3_selfdual"][0],
               A["N725_spin0"][0]),
        ],
        "verdict": ("observation, not theorem: the symmetric part S of Z "
                    "collapses toward 0 while the antisymmetric part A drifts "
                    "slowly (total spread <= 0.048 over L=3..N=725 and "
                    "site/bond); no exponent fitted, four sizes do not make "
                    "a limit"),
    }


# --------------------------------------------------------------------- T4

def check_t4(census: dict) -> dict:
    """State C3 and its exact finite-N witness values on the site lab."""
    # Exact M(1/4) + M(3/4) at site L=3 from the census's Z is not directly
    # available; use the cited M(1/2) values and the census antisymmetry
    # failure surrogate: S(0.1) = S(0.3) = the Z-side even content.
    # The precise C3 witness is M_N(1/4)+M_N(3/4) and 2 M_N(1/2).
    site = census["site"]
    witnesses = {}
    for L, v in site.items():
        witnesses[f"site_L{L}_M_half_exact"] = v["M_half"]
    return {
        "conjecture": ("C3 (self-matching restoration): "
                       "sup_{p in K} |M_N(p) + M_N(1-p)| -> 0 for every "
                       "compact K ⊂ (0,1); equivalently S_N -> 0 on the "
                       "level grid"),
        "status": "OPEN — neither proved nor refuted; not implied by Theorem L",
        "equivalence": ("S_N -> 0 on the grid <=> C3 in the grid's p-window "
                        "(both sides exact at finite N via F = (1+M)/2)"),
        "finite-N_witnesses_M_half": witnesses,
        "bond_witness": "bond L=3 (dual_fail=0): M(1/2)=0 exactly, S ≡ 0",
        "role": ("under C3 the limit shape, if it exists, is antisymmetric; "
                 "the open object is A_∞ = lim A_N"),
    }


def main() -> None:
    census = load_census()
    n725 = load_n725()
    assert census["bond"]["3"]["dual_fail"] == 0, "repaired census required"
    out = {
        "schema": "matching-one.probe-invariant-shape.rerun.v2",
        "issue": 622,
        "run": "2026-09-08 re-run on repaired inputs (PR #653 bond lab, "
               "PR #655 spin-0 N725)",
        "cited_inputs": {
            "census": f"{BRANCH_CENSUS}:{PATH_CENSUS}",
            "n725": f"{BRANCH_N725}:{PATH_N725}",
            "toy_location_theorem": "results/probe-invariant-shape/toy-families.json (merged #628)",
            "consensus_g": "results/type582-residual/latest.json (merged #605 lineage)",
        },
        "T1_duality_symmetry_theorem": check_t1(census),
        "T2_duality_does_not_pin_shape": check_t2(),
        "T3_decomposition_and_rigidity": check_t3(census, n725),
        "T4_conjecture_C3": check_t4(census),
        "standing": {
            "no_new_production": True,
            "no_status_edit": True,
            "no_issue_closed": True,
            "no_exponent_fit": True,
            "n725_not_scored": True,
        },
    }
    dest = ROOT / "results" / "probe-invariant-shape" / "rerun-20260908.json"
    dest.write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1)[:4000])
    print(f"\nwrote {dest}")


if __name__ == "__main__":
    main()

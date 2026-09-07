# C12: mission update — what twelve turns of evidence say about the north-star question

**Date:** 2026-09-06.  This closes the planned C1–C12 cycle (C6 deferred on
data access).  The question was: *is "state dimension" really a family of
quantities indexed by a declared experiment language and its symmetry, rather
than one intrinsic integer of a stochastic process?*

## 1. Verdict on the north-star proposition

**Supported, with sharpened boundaries.**  On every repository-native object
examined, the apparent "state complexity" is a functional of a declared
object — and the probe produced exact statements about *which* functional:

* P398: 14..1430 (microscopic) > 10..750 (C2 orbit / exact positive
  quotient, Gate-1) > up to 209 (static D2 readout classes) > 3–4
  (finite-horizon balanced order).  Each number is a different well-defined
  functional of a *declared* dictionary/task/inner product (C4, C8).
* #549: k+1 predictive classes coexist with response rank ≤ 2 for the
  single-root language; rank grows only with language depth (d roots ⇔ degree
  d ⇔ rank ≤ d+1).  Class counting and linear independence are provably
  different coordinates (C2).
* The quantities that *do* survive observer changes are exactly the structural
  ones: selection zeros (representation-theoretic), the quotient-factorisation
  identity, per-sector inertness, and the language protocol algebra (C8).

**Corrections the evidence forced on the original framing.**

1. The naive list "many numbers, pick your favourite" is not the point:
   the numbers are *ordered by a factorisation* — symmetry quotient first,
   then task compression (Kalman for C2; isotypic decomposition in general).
2. "Even = orbit-constant" is a C2 accident; the general first factor is the
   isotypic decomposition, and exposing a nontrivial fibre requires a
   character-marked channel (round 2, C5).
3. The probe's own Phase-C result was (correctly) absorbed as a special case
   of Kalman/Wonham, and its orbit-count arithmetic as Callan–Smiley
   combinatorics — a useful discipline check (C1).
4. The strongest *negative* results are the no-go certificates: an unbranched
   finite-horizon order cannot constrain a fork-different mechanism (#435/
   #549), and an even task cannot constrain any hidden-sector property,
   including the location of its ergodicity loss (C7).

## 2. What the probe leaves that is durable

| deliverable | file |
|---|---|
| Complexity ladder v1 (every arrow P/F/C/U + counterexamples X1–X9) | `probe-complexity-ladder-v1-20260906.md` |
| Experiment-language formalisation (kappa_L, r_L, minimal distinguishing depth) | `probe-experiment-languages-formal-20260906.md` |
| Markov-generator selection rule (formal, fills #601 Q3 gap) | `probe-markov-selection-rule-and-l2pi-fibre-20260906.md` |
| No-go certificates (unbranched-vs-fork; sector surgery) | `probe-threshold-no-go-theorem-20260906.md` |
| Invariant catalogue | `probe-invariant-catalogue-20260906.md` |
| Ranked bridge conjectures with cheap tests executed | `probe-bridge-conjectures-ranked-20260906.md` |
| Reproducible scripts (deterministic, no sampling) | `scripts/probe/*.py` |
| One-branch submission pack | `matching-one-probe-full-20260906.bundle` + `submit-pack-v2/` |

## 3. Remaining open problems, ranked (short version of C10)

1. exact `r_lin(w8)` — blocked on #593's definition (a naming/definition
   issue, not a compute issue);
2. `r_mem(w9)` saturation — compute is cheap (C3 cost probe); gate is porting
   the frozen convention;
3. two-root fork degree on the real N16 network (needs N16 site data, no
   sampling);
4. whether any current declared task sees the sector that `E_p[X] = 0`
   targets (positive threshold bridge, burden set by C7);
5. compositional depth ⇔ algebraic degree for general tree languages.

## 4. Recommendation for the repository

* **Science to keep visible:** the language-first complexity hierarchy (C4),
  the Markov selection rule (round 2/C5, which #601 Q3 declares a literature
  gap), and the two no-go certificates (C7).  These are the only items the
  sibling agents (#600/#601/#603) did not already cover.
* **Submission form:** one open PR with the bundle contents (no auto-merge)
  plus comments on #598/#599/#593/#594 where the probe directly resolves or
  sharpens an open question.  A live GitHub token is required; the pack is
  ready in `submit-pack-v2/` and as `matching-one-probe-full-20260906.bundle`.
* **Process notes:** sibling overlap on Phase C and orbit counts was handled
  by re-labelling (cross-confirmation / known combinatorics) rather than by
  competing claims — keep that practice.

## 5. Boundaries retained

All statements are exact finite statements on P398, the published #549/#435
protocol algebra, the declared C3 synthetic object, or the two-copy no-go
family.  Nothing here is a site-percolation threshold theorem; C2 parity is
not continuum spin; the balanced/memory numbers remain task- and
convention-relative (C5/C8).  The probe's value is in drawing the boundaries,
not in erasing them.

## 6. Continuation

If the probe continues, the highest-information next cycle is: (i) obtain the
#593 `r_linear` definition and close w8; (ii) port the frozen memory
construction and run w9 (C3 decision table); (iii) two-root N16 enumeration;
(iv) — only with archived data access — the #582 discriminator run (C6).
Items (i)–(iii) need no sampling and no new compute budget beyond a normal
session.

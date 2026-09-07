# Probe MID — After `Q_N(u) → p_c`: which rates are *theorems from already-named inputs*?

**Assign to:** one Agent with solid mathematical writing. Ordinary CPU. **No production, no N=725 rescoring, no optimizer, no arXiv search.**

**Literature is out of scope.** A fourth document owns retrieval. You work **only** from sources already named in `notes/p613-quantile-convergence-20260907.md`, PR #606’s ledger, PR #615, and the in-repo Russo formula. If a step needs a paper that is not already named, write `BLOCKED_ON_LITERATURE_PROBE: <what you need>` and continue the next program. Do not go looking.

---

## Standing

Long-horizon *proof-skeleton* probe. You turn the scoped corollary that already exists into a wall of what it can and cannot imply.

- Does **not** enter `docs/STATUS.md`.
- Does **not** close #276, #321, #613, #606, #615, #275.
- Negative answers are the main product.
- Mark T / C / F / X on every displayed statement:
  - **T** = theorem from named inputs already in the #613 note or the #606 ledger;
  - **C** = conditional on a named open hypothesis (F1, F2, RSW-on-square, …);
  - **F** = would be a fit (you will not perform it);
  - **X** = not well-posed.
- Do not retarget charts (HIGH probe). Do not re-enumerate L=3,4 (exact-controls). Do not reopen #275. Do not fit `ν`.

Read:

```text
PR #614   notes/p613-quantile-convergence-20260907.md     (Theorem L, H1–H3, DKS parent, thin-torus failures)
PR #606   notes/homological-balance-root-ledger-20260906.md  (F1, F2, conditional F, Outcome J)
PR #615   notes/bounded-task-rank-threshold-no-go-20260907.md (CSC no-go; escape D = signed observable)
in-repo   M'(p) = pivotal_primal(p) + pivotal_matching(1-p)   (finite Russo; STATUS / #43)
frontier  claude/matching-one-workspace-pwr5pv @ 8b5f9d1a
```

PR against **frontier**, not `main`.

---

## Why this exists

#613 / PR #614 buys **location**:

```text
Q_N(u) → p_c     uniformly on compact subsets of (0,1)
```

for honest periodic square-cell tori with `ℓ_N / log N → ∞`, under H1–H3 as written in that note. Primitive Gaussian squares have `ℓ_N = √N`. Mixed-sign bounded weights inherit the limit.

It buys **no rate**. The note is silent on `N^{-ω}`, `L^{-4}`, and operators. Your job is to make that silence *provable* from the argument as written, and to write the conditional statements that *would* follow if F1 and/or F2 were granted — without granting them, and without searching for extra papers.

Working suspicion to prove, falsify, or replace:

```text
No polynomial rate for Q_N(u) − p_c follows from H1–H3 plus
subcritical exponential decay. The empirical L^{-4} root shift
is a product of two independent open inputs (F1, F2), not a
third mechanism. #582's 9-vector has no right to a single
exponent from location alone. An Aff(1)-invariant shape Z_N
is invisible to Theorem L.
```

---

## Named inputs you may use (closed list)

If it is not on this list, it is blocked.

From the #613 note (primary quotes already there — *re-use them*, do not re-summarise):

- H1: `r_G + r_Ĝ = 2` (repo digital Alexander).
- H2: DCT, CMP 343 (2016), arXiv:1502.03050, Thm 1.1 item 3 + site remark §1.2; Aizenman–Barsky; Menshikov.
- H3: Grimmett–Li RSA 65 (2024) prove `p_u(G)+p_c(G*)=1`; (H3) on amenable Z² uses `p_c=p_u` (Burton–Keane). van den Berg JMP 22 (1981) counterexamples outside a restricted class.
- Union bound with `R_N = ⌊ℓ_N/4⌋`, the NN and NN+NNN distance conversions written in the note.
- Strict monotonicity via Margulis–Russo plus one rank-increasing step.
- DKS arXiv:2011.11903 / AIHP 2025 as *parent*, square-site = same argument + H3, **no novelty claim**.
- Failures already listed: thin tori `ℓ_N = O(log N)`, `ℓ_N ≤ √2` not honest, unbounded extrapolation weights, `u_N → 0` or `1`.

From the #606 ledger:

- F1: `M_L'(p_c) ≍ L^{3/4}` ⇔ `α_4 ≍ L^{-5/4}`, `ν=4/3`; triangular proved (Smirnov/LSW); square open.
- F2: `M_L(p_c) ≍ L^{-13/4}` (Q4, `x=21/4`) — conjecture.
- Conditional F: F1+F2 ⇒ `p_L^H − p_c ~ L^{-4}` under the ledger’s amplitude and error hypotheses.
- Unconditional: self-matching lattices have `p_L^H = p_c = 1/2`.
- Outcome J: F1 and F2 not currently in the same rigorous model.
- Census numbers: cite exact-controls JSON; do not rerun.

From #615:

- Bounded task rank has no threshold content without CSC.
- Escape A = CSC (necessary). Escape D = signed observable with a sign-change theorem (independent of rank). Route D *is* this homological-balance route.

From the repo:

- `F_N(p) = [1+M_N(p)]/2`.
- `M'(p) = pivotal_primal(p) + pivotal_matching(1-p)` exactly.

---

## Programs (all required)

### R0 — Quote sheet, no paraphrasing

One page. For H1, H2 (Thm 1.1 item 3 + site sentence), H3 (what Grimmett–Li actually prove vs what the note uses), DKS (what they prove in d=2, i=1), van den Berg’s counterexample sentence: **copy the quotes already in the #613 note**, with the note’s page/theorem pointers. If the note is missing a quote it *claims* to have, that is a defect in #613 and you record it. You do not go to the PDF unless the literature probe has already deposited it; if you cannot verify a quote from the note itself, mark `UNVERIFIED_IN_NOTE`.

This program exists so later programs do not drift.

### R1 — No polynomial rate from the union bound (main theorem-shaped negative)

Write a complete proof, in your note, of:

```text
Theorem R1 (tag T, from H1–H3 + DCT as used in #613).
The argument of Theorem L yields, for each fixed p < p_c,
    P_p(r_G > 0) ≤ N · A(p) e^{c(p)} exp(−c(p) ℓ_N / 4) → 0
whenever ℓ_N / log N → ∞. The constant c(p) > 0 is allowed to
tend to 0 as p ↑ p_c. Therefore this bound does not imply
    |Q_N(u) − p_c| ≤ C_u N^{−ω}
for any ω>0, any u ∈ (0,1), or any C_u independent of N.
```

The proof must exhibit *why* `c(p)→0` blocks a power: the N at which the bound becomes small may recede faster than any polynomial as `p_N = p_c − N^{-ω}` approaches `p_c`. Write that receding-window calculation with `c(p)` kept as an unspecified positive function with `c(p)→0`. Do not insert `ν` or RSW to rescue it — that would be tag C, program R4.

If you believe a polynomial *does* follow, you must derive it using only the closed list. Unexpected. Write it as a claimed T and expect it to be killed.

### R2 — What Theorem L *does* give, stated sharply

A list of T-tagged corollaries, each ≤ 10 lines, no extras:

- interior quantile location, uniform on `[ε,1−ε]`;
- mixed-sign bounded weights inherit location (prove it from “combine quantile functions, not CDFs” — this is in the note; write the three-line proof);
- unbounded weights fail (give the amplification argument already in the note, fully);
- `u_N → 0` or `1` not claimed;
- thin tori fail, with the quasi-1D sentence;
- `ℓ_N ≤ √2` not honest, with the identified-cell sentence;
- **no** statement about `F_N(p_c)`, **no** rate, **no** operator.

If you want a corollary not on this list, it needs a proof from the closed list or it is blocked.

### R3 — `p_L^H` is not a matching-odd vanishing of `Q(1/2) − 1/2`

This is an identity program, numbers cited from exact-controls.

Write, as T from the declared contract `F = (1+M)/2` and strict increase of `M` (census: `M(0)=-1`, `M(1)=+1`, monotone):

```text
If M is continuous and strictly increasing [−1,+1], then
p_L^H := unique root of M  satisfies  F(p_L^H)=1/2,
hence Q_L(1/2) = p_L^H.
```

So the slogan “`p_L^H` is not the median” is **false** under the contract, if “median” means `Q(1/2)`. What *is* true, and what F1/F2 use:

```text
p_L^H ≠ 1/2,  equivalently  M(1/2) ≠ 0.
Q_L(1/2) − 1/2  is not a self-dual vanishing.
Self-duality / self-matching would force M(1/2)=0 and p_L^H=1/2.
M(1/2) is the matching-odd content at the self-dual *parameter*,
a different coordinate from Q_L(u) − p_c.
```

Write the first-order expansion that the #606 ledger uses:

```text
M(p_c) + M'(p_c) (p_L^H − p_c) + remainder = 0
p_L^H − p_c = − M(p_c)/M'(p_c)   +  higher order.
```

Then, *conditional* (tag C) on the ledger’s F1 and F2:

```text
M'(p_c) ≍ L^{3/4},   M(p_c) ≍ L^{-13/4}
⇒  p_L^H − p_c ≍ L^{-4}.
```

Write the remainder hypothesis exactly as the ledger states it. Do not improve it. Tag Outcome J: you may **not** treat this C as a T.

Cite L=3,4 numbers (`M(1/2)`, `p_L^H`) from exact-controls; if that JSON is not yet on frontier, cite PR #606’s JSON and leave a hole.

### R4 — Conditional statements, locked to named opens

A table, one row per extra hypothesis **already named** in #606/#613. You do not add rows by searching.

| extra (already named) | implied rate | applies to | tag | still open because |
|---|---|---|---|---|
| F1 only | `M' ≍ L^{3/4}` | denominator of root shift | C | square-site conformal invariance |
| F2 only | `M(p_c) ≍ L^{-13/4}` | numerator | C | operator identification |
| F1+F2 + ledger remainder | `p_L^H − p_c ~ L^{-4}` | root, not the 9-vector | C | Outcome J |
| RSW / box-crossing as *named* in the 606/613 notes (triangular / Ising, not square) | whatever those notes already say | **not** square-site Q_N | C or X | square site not in the named theorems |
| Russo formula (in-repo, exact) | `M'` is a pivotal mass, no rate | derivative, not Q | T | — |

If you are tempted to write “RSW on square site ⇒ quantile window `N^{-3/4}`”, that extra is **not named as a square-site theorem** in the closed list. Row = `BLOCKED_ON_LITERATURE_PROBE`. Do not invent the implication.

### R5 — One exponent for the 9-vector is not implied by location

Tag T or X, from R1+R2 only:

```text
Theorem L + no polynomial rate  ⇒  the 9-vector
(Q_N(0.1),…,Q_N(0.9)) has no theorem-given common ω.
Uniformity on compact u is a statement about limits, not about
a common power. Even if each u had a rate, the rates need not
agree. #582’s “one transferable direction” is therefore not a
corollary of location.
```

Write this so that the HIGH probe can cite it when refusing a single-power 7-vector fit as a *location* argument. (They may still fit as F — that is their stop rule, not yours.)

### R6 — `Z_N` is invisible to Theorem L

Let `Z_N(u) = (Q_N(u)−Q_N(0.5))/(Q_N(0.8)−Q_N(0.2))` (or any Aff(1)-invariant of Q). Prove:

```text
Theorem L constrains Q_N(u) → p_c, hence both numerator and
denominator of Z_N tend to 0. The limit of Z_N is an 0/0 form
and is not determined by Theorem L.
```

This is the interface lemma with the HIGH probe. Tag T (elementary analysis + Theorem L). Do not compute Z_N on production data.

### R7 — Scaling function as a *definition*, if F1 were true

Grant F1 hypothetically (tag C). Define

```text
τ = L^{1/ν} (p − p_c) = L^{3/4} (p − p_c)
Φ_L(u, τ) := F_L( p_c + L^{-3/4} τ )     (or the inverse for Q)
```

and write, as a definition not a fit, what “a single scaling function `Φ(u,τ)` exists” would mean for:

- adjacent-size first differences of Q at fixed u;
- three-size second differences;
- why a *chart* (HIGH probe) can still fake a curvature even if `Φ` exists;
- why F2 (matching-odd at `p_c`) is a *different* coordinate (`τ=0` value of M, not a derivative).

Do not estimate `ν`. Do not use N=725. Do not use #582’s amplitudes except as “this is the object that would have to be a vector field on `Φ`”.

### R8 — Russo / pivotal, what it cannot do

From the exact formula `M'(p) = pivotal_primal(p) + pivotal_matching(1-p)`:

- T: `M'` is a sum of two positive masses; sign of `M'` is plus; used in R3’s strict increase.
- T: this does not give `M'(p_c) ≍ L^{3/4}` (that is F1).
- T: integrating `M'` from `1/2` to `p_c` recovers `M(p_c)−M(1/2)`, which is not F2.
- X: using this formula to “identify Q4” is category error.

No Monte Carlo of pivotals.

### R9 — #615 escape D is this route

One page, tag T as a scope statement:

```text
The homological-balance root is #615’s escape D (signed
observable + sign-change), not escape A (CSC) and not a
bounded-rank statement. Theorem L is the location half of D
(quantiles squeezed to p_c). F1/F2 are the rate half of D and
are open. P398 balanced order remains irrelevant to p_c.
```

Do not reopen #615’s no-go.

### R10 — Punch list for the literature probe

Every `BLOCKED_ON_LITERATURE_PROBE` you wrote, collected as a numbered list of *questions*, not papers. This is the only handover. Typical expected holes (delete any you actually closed from the named list):

- does any *already-proved* planar RSW imply a quantile-window power, even off square site?
- is there a theorem that inverse-CDFs inherit sharp-threshold windows?
- what finite-size *correction* (not location) is proved on triangular site / FK-Ising for crossing probabilities?

You do not answer these. You make them precise enough that the literature probe cannot waffle.

---

## Stop rules

- If R1 fails (you think a power follows), that becomes the paper; stop adding conditionals until it is killed or accepted.
- If you need a source not on the closed list, block and move on — do not search.
- No STATUS, no close, no MC, no `ν` fit.

---

## Deliverables

```text
notes/probe-quantile-rate-after-613-YYYYMMDD.md
  R0 quote sheet
  R1 full proof (no-power)
  R2 corollary list
  R3 identity + C expansion
  R4 table
  R5 9-vector negative
  R6 Z_N interface lemma
  R7 Φ definition
  R8 Russo wall
  R9 #615 scope
  R10 literature-probe questions
```

Scripts only if they check an algebraic identity (R3 expansion as a formal power series in sympy is allowed; no data). Ordinary CPU.

## Interface

- HIGH: they may cite R5 and R6. They do not ask you for a rate of `g`.
- Exact-controls: you cite their `M(1/2)`, `p_L^H`, `Q(u)+Q(1-u)`.
- Literature probe: your R10 is their input. You do not duplicate it.

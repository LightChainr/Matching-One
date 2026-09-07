# Probe HIGH-2 — Claim-survival atlas, the Alexander pair, and F_N in function space

**Assign to:** one high-reasoning Agent. Ordinary-to-serious CPU is enough (exact L=3,4 already enumerated; re-reading committed histograms at extra quantile grids). **Not** the #617 affine-gauge scorer: do not redo λ-sweeps, holonomy, or the 4% jackknife. Cite #617’s freeze if that PR exists; otherwise use the default chart written below.

**Not for:** literature retrieval (#620), rate-skeleton proofs from named inputs (#618), or the exact-controls checklist (#619). Import their JSON when present; do not duplicate their programs.

**Literature is out of scope.** No arXiv. No new production. No STATUS. No ticket close.

---

## Standing

Long-horizon probe, same epistemic level as #599 / #610 / #617.

- Does **not** enter `docs/STATUS.md`.
- Does **not** close #582, #584, #609, #610, #612, #275, #276, #321, #581, #608, #599.
- Does **not** reopen #275. You will *write* the identification contract the next candidate would have to meet; you will not propose a candidate.
- Does **not** fit `A1 N^{-ω1}+A2 N^{-ω2}`. **1.55 is retired.**
- Does **not** dump results into PR #616. New PR against frontier, comment on the issue.

Frontier: `claude/matching-one-workspace-pwr5pv` @ `8b5f9d1a`. Packet of briefs: PR #616. Sister HIGH probe: #617. Read PR #614, #611, #605, #606, #615.

Default declared chart, unless #617 has overwritten it in a merged or reviewable protocol card:

```text
declared chart  = span{1, Q_base, g_frozen}
transport       = #612 identity, applied before any residual
weightings      = spin0 AND equal, always together
reconstruction  = production nine-decile inverse-CDF, named
```

If you need a number that only #617 can produce (λ-sweep, holonomy scalar, 4% 7-vector), write `BLOCKED_ON_617` and continue. Do not secretly become #617.

---

## Why this exists

#617 asks: what is the gauge groupoid of the 9-vector pipeline, and what (if anything) is invariant?

That is necessary and not sufficient. The repository currently talks as if the scientific object *were* that 9-vector on a single NN model with one weighting. The last week showed four silent choices that can manufacture “physics”:

```text
(i)   affine-chart attachment (base vs middle)     → fake 1.55
(ii)  orientation weighting (spin0 vs equal)       → 12% / 25σ
(iii) 9-decile readout of F_N                      → unquantified truncation
(iv)  single model G instead of the Alexander pair
      (G, Ĝ) with r_G + r_Ĝ = 2                    → fake matching-even g
```

Astra already killed (iv) as an interpretation of fitted `g`. Nobody has written the *positive* object: even/odd coordinates on `G ⊕ Ĝ`, built from the **same** configurations (black 4-connect / white 8-connect). Exact L=3,4 have every config. Production histograms may only have stored `r_G`. That difference is itself a result.

Separately, `Q_N(u)` at nine deciles is a chart on the space of CDFs, not the CDF. Theorem L (#613) is a statement about the function `F_N` (or all interior quantiles uniformly). #582’s “one transferable direction” is a statement about a 9-vector in a covariance metric. Those are different categories. A high-value negative is: the 9-vector `g` does not lift to a vector field on `F_N` modulo `Aff(1)`.

#275 remains P0 and closed. The reason it is hard is the same reason 1.55 was fake: a low-dimensional direction, a pretty exponent, and a cross-lineage ratio can all move when (i)–(iv) move. The reopening rule should be rewritten in this language, not relaxed.

---

## North-star question

For every live scientific sentence about square-site threshold law in this repository, decide:

```text
survives (i) chart transport?
survives (ii) orientation weighting?
survives (iii) replacing the 9-vector by F_N as a function
              (or by the law of r)?
survives (iv) replacing the single NN model by the Alexander pair?
```

Then **construct**, without fitting and without production, the two objects the pipeline never formed:

```text
Pair object:     even/odd coordinates on G ⊕ Ĝ,
                 from the existing Alexander pair on each config.
Function object: Aff(1)-quotient of F_N in function space,
                 of which the 9-vector Z_N / r_N is a discretisation.
```

Working suspicion, to prove, falsify, or replace:

```text
Most sentences currently spoken about "g", "N^{-0.970}", and
"the 4%" are sentences about a 9-vector in one chart with one
weighting on one model. The process-level sentences that survive
(i)–(iv) are at most:

  • r_G + r_Ĝ = 2                         (exact, every config)
  • Q_N(u) → p_c on compact u             (Theorem L, H1–H3)
  • M(1/2) ≠ 0                            (exact finite, not self-matching)
  • P398-style selection is representation parity, and
    uniform p is even                     (#615, #610)

Everything else is observer geometry until a survival proof
is written. #275 may not reopen until a candidate specifies
(i)–(iv) in advance and produces original-U forward observables
under that freeze.
```

A clean impossibility — “`g` cannot be an irrep of the matching involution on the pair, and also cannot be a function-space geodesic of `F_N` modulo Aff(1)” — is a successful outcome.

---

## What is already settled (cite, do not re-derive)

1. Chart identity, 97% of 1.55, leftover ~4%, N=725 5.56σ-inside-band — PR #611/#614. **1.55 retired.**
2. Gate 3 — PR #605: no pre-existing discrete label indexes `r_N`.
3. Theorem L — `notes/p613-quantile-convergence-20260907.md`. Location, not rate.
4. `M_3(1/2) = -21/64`, `M_4(1/2) = -13757/32768` — PR #606. Self-symmetry of a single NN model is false.
5. Matching involution `M_Ĝ(p) = -M_G(1-p)` is a **two-model** relation (#610 / #614).
6. #608 wrap/`X` identity is algebra + Harris, typed observer-separation control, not a mechanism. Exact-controls (#619) owns the L=3 bond check.
7. #615: bounded task rank has no threshold content without CSC. P398 selection needs an odd perturbation; uniform `p` is even.
8. #599 asked for an experiment object `E = (dynamics, source language, readout language, horizon, symmetry)`. This probe **adds chart, weighting, discretisation, and single-vs-pair to `E`**, and asks which claims survive the enlarged `E`. You do not rebuild #599’s P398/cut-network atlas.

---

## Programs (all required)

A program is done only with a written survival verdict or a constructed object plus a killing test. “Interesting” is not done.

### S0 — Inventory of live sentences

List, as numbered claims C1, C2, …, every sentence the repo currently treats as live about square-site (or about the pair). Minimum set; add if you find more in STATUS / #582 / #610 / #614, but do not hunt old closed tickets:

| id | sentence (paraphrase allowed here; quote in the note) |
|---|---|
| C1 | `F_N = [1+M_N]/2` is the declared contract |
| C2 | interior `Q_N(u) → p_c` (Theorem L, H1–H3) |
| C3 | `M(1/2) ≠ 0` on honest square-site tori |
| C4 | one transferable non-affine direction `g` carries the 9-vector flow |
| C5 | `A(N) ~ N^{-0.970}` is a finite-size forecast in the declared chart |
| C6 | leftover ~4% / 7.3σ curvature after chart transport |
| C7 | spin0 vs equal disagree at ~12% |
| C8 | Gate 3 remainder `r_N` has no pre-existing discrete label |
| C9 | fitted `g` is matching-even (this sentence is **already supposed to be dead** — confirm) |
| C10 | `p_L^H − p_c ~ L^{-4}` (conditional F of #606) |
| C11 | #608 coupling is a mechanism fingerprint (already supposed dead) |
| C12 | P398 balanced order locates `p_c` (already supposed dead, #615) |
| C13 | original Matching-One `U` is identified with `g` or with H4/H8 (#275) |

You may split C’s if a sentence is actually two. You may not add a C14 after seeing S2 numbers.

### S1 — Survival table under (i)–(iv)

For each Ck, four cells: survives chart / weighting / function-space readout / Alexander pair? Allowed entries:

```text
YES     — proof or exact identity
NO      — counterexample you compute or cite
OPEN    — needs #617 / #618 / #619 / #620
DEAD    — already retired; do not revive
N/A     — the sentence is not about that axis
```

C4 and C5 are the ones people will want to be YES. They are not allowed to be YES without a proof that names the freeze of (i)–(iv). Default: C4/C5 are **NO** for (i) and (ii) unless #617 produced an invariant; **OPEN** for (iii) until S2; **NO** for (iv) by Astra’s involution note unless you construct a pair-irrep that matches `g` (S3).

Deliverable: a markdown table and a one-paragraph verdict per Ck. This table **is** the atlas.

### S2 — Nine deciles are a chart on the space of CDFs

Object: `F_N : [0,1] → [0,1]`, monotone, `F_N = (1+M_N)/2`. The pipeline’s `Q_N ∈ R^9` is the inverse on `{0.1,…,0.9}`.

On exact L=3 and L=4 (import #606 enumerator / #619 polynomials; if #619 is absent, invert `M` by bisection yourself — 512 / 65536 is in budget):

- compute `Q(u)` on a 17-point grid `{0.05,0.10,…,0.95}` and on the 9-point production grid;
- define a function-space remainder after projecting the 17-vector onto `span{1, Q_9-interpolated, g_9-lifted}` (lift `g` by piecewise-linear interpolation of the 9-vector direction; **name the lift**, it is a choice);
- report whether the extra eight points are explained by that lift (truncation small) or not (9-decile `g` is grid-dependent).

On **one** committed production histogram (prefer N=290 or N=325, raw hist on PR #614 or Gate 3 assets): re-read the inverse-CDF on the same 17-point grid. Same test. If raw hist is missing, `BLOCKED` that half of S2; still deliver exact L=3,4.

Kill: if the 17-point shape is not in the span of the 9-point lift within jackknife/exact error, then C4 is a sentence about a grid, and the function object is not `g`. Do not then fit a new `g_17`.

Also write, as a lemma (tag T, elementary):

```text
Theorem L is a statement about F_N (all interior quantiles).
#582's χ² is a statement about a 9-vector in a batch-covariance
metric. Neither implies the other. Aff(1)-invariants of F_N
are 0/0 forms under Theorem L (#618 R6); that does not give
them a finite-N law.
```

### S3 — Even/odd on `G ⊕ Ĝ`, from the same configs

On every exact L=3,4 site config the enumerator already has `(r_b, r_w)` with `r_b+r_w=2`. Identify `G` = black 4-connect, `Ĝ` = white 8-connect (state the convention; if #619 C8 coded the involution polynomial, import it).

Define, and then accept or kill:

```text
M_G(p)      = E[r_b] - 1
M_Ĝ(p)      = E[r_w] - 1           (at the same p, not at 1-p)
M_even(p)   = (M_G(p) - M_G(1-p)) / 2
M_odd(p)    = (M_G(p) + M_G(1-p)) / 2
```

Check exact polynomial identities (L=3,4):

```text
M_G(p) + M_Ĝ(1-p)  ?=  0          (involution)
M_odd(1/2)         ?=  M_G(1/2)   (the −21/64)
M_even(1/2)        ?=  0
```

Then form quantile objects of the pair. Candidates, each with a killing test:

| id | object | kill if |
|---|---|---|
| P1 | `Q` of `F = (1+M_G)/2` (current pipeline) | — baseline |
| P2 | `Q` of `F_even = (1+M_even)/2` | undefined because `M_even` may not run `[-1,1]` |
| P3 | law of `r_b − r_w` (pure odd, values in `{−2,0,+2}`) | independent of `p` (then it is not a threshold coordinate) |
| P4 | `P(r_b=2) − P(r_b=0)` vs `M` | they are linearly equivalent (then you have not found a new observable) |

Verdict required: **fitted production `g` is not P1–P4 of the pair**, because productions never formed the pair. Write the one-sentence no-go:

```text
An irrep of the matching involution is a property of a function
on G ⊕ Ĝ. The fitted g is a property of nine quantiles of G.
These live in different categories. No numerical closeness of
g to a discretisation of M_even/odd can identify them
(GOVERNANCE §2E: derived views of one histogram).
```

If production JSON stored only `r_G` histograms, record `PRODUCTION_PAIR_NOT_RECOVERABLE` and specify the extra histogram a future production would have to write (`n_black` is not enough; you need `r_w` or the white wrapping). **Do not request that production.**

### S4 — Function-space Aff(1) quotient, definition only

Define the object #582 *would* have been if it had been about `F_N` not `Q_N ∈ R^9`:

```text
Let P be the space of CDFs on [0,1] (or monotone maps
[0,1]→[0,1] with F(0)=0, F(1)=1, after the usual cutoff).
Aff(1) acts by F ↦ F ∘ ϕ^{-1} with ϕ(p) = α p + β, restricted
to the range where this stays a CDF — **or** state why that
action is the wrong one and the right action is on the
quantile function Q(u) ↦ α Q(u) + β.
Write the action you actually believe, as a formula.
Then the shape is the orbit space P / Aff(1).
A finite-size "direction g" would be a tangent vector
to that orbit space along N.
```

Kill the wrong action: location-scale on `p` (the occupation parameter) is **not** the same as location-scale on `Q(u)` (the inverse-CDF). The pipeline used the latter. Theorem L is about the former’s limit point `p_c`. Write that distinction as a lemma. If you cannot, S4 is not done.

You may compute, on L=3,4 exact `F`, a discretised distance to the Aff(1)-orbit of a two-parameter family (e.g. beta CDFs, or the affine-warped `F_*((p-p_c)/w)`). This is a **toy**, tagged as such. If a two-parameter warp already fits L=3,4 to machine precision, the exact census has no extra shape; say so. If not, report the residual function, unnamed.

Forbidden: fitting that toy to N=65…725.

### S5 — #275 identification contract, rewritten, not reopened

#275 stays closed. Write a one-page contract that a future candidate **must** fill *before* any profile-rank:

```text
1. declared chart (i), weighting (ii), readout grid (iii),
   single vs pair (iv)   — frozen in the candidate note, SHA-dated
2. original-U forward map: candidate → at least two raw
   coordinates already in the #275 covariance
3. allowed amplitude class (signed-real / complex / …)
   as already frozen in #275; no new class
4. nuisance: anything in the Aff(1) × weighting span is
   nuisance, not signal
5. one frozen scoring, then stop
```

If a candidate cannot fill (1)–(4), the verdict is `UNIDENTIFIABLE_WITH_CURRENT_ASSETS` under #275’s own rule. Do not name H4, H8, or `g` as that candidate. Do not score.

### S6 — Antisymmetric perturbation, design only (P398 scope)

Uniform `p` is even (#610 / #615). Write a **design**, not a run, of a site-dependent tilt `p_x = p + ε σ(x)` with `σ(Rx) = −σ(x)` for a chosen spatial reflection R of the torus (or of P398’s interval).

Deliver:

- the exact P398 lemma’s hypotheses matched line by line (`RG_0=G_0 R`, `RHR=−H`, endpoints invariant) vs which of them a percolation tilt can satisfy;
- why this would *not* locate `p_c` (it is a selection rule for a response, not a threshold);
- why it must not be used to rescue C12.

No P398 width-9 diagonalisation unless the matrices are already in-tree and a 50-line script can check the lemma on `w=4`. If not, design-only is enough.

### S7 — Rewrite the scientific map in eight lines

After S0–S6, eight lines, no seventh mechanism:

```text
exact:     r + r̂ = 2
location:  Q_N(u) → p_c                 (H1–H3)
not self-matching: M(1/2) ≠ 0
observer:  chart × weighting × grid × pair   (#617 + this probe)
forecast:  g in a declared chart, few-percent, not a law
remainder: unnamed, only after transport
P398:      task observability, not threshold
#275:      closed until the contract in S5 is filled
```

If your atlas disagrees with these eight lines, change the atlas or change the lines, and say which. Do not add a ninth line about a second scale.

---

## Compute you may do

- L=3,4 exact (512, 65536) — yes.
- Re-read **one** committed histogram at 17 quantile points — yes.
- Import #606 / #619 polynomials — yes.
- λ-sweep, N=725 equal/spin0 rescoring, 1e4 label permutations — **no**, that is #617.
- New N, extra batches, height-4 — no.
- arXiv — no.

---

## Stop rules

- If S1’s table has no YES outside C1–C3, that is a result: write `ALMOST_NOTHING_SURVIVES_OBSERVER_CHANGE` and still do S3–S5.
- If S2 shows 9-decile `g` lifts to 17 points, do not declare victory for C4; that is only (iii) on two exact sizes.
- If S3’s involution polynomial fails, stop S3 and report; do not change connectivity conventions silently.
- No STATUS, no close, no production, no second exponent, no #275 scoring.

---

## Deliverables

New PR against `claude/matching-one-workspace-pwr5pv`, **not** #616:

```text
notes/probe-claim-survival-pair-function-YYYYMMDD.md    (S0–S7 atlas)
scripts/probe/survival_s2_quantile_grid.py
scripts/probe/survival_s3_pair_even_odd.py
scripts/probe/survival_s4_aff_toy_L3L4.py               (toy, tagged)
results/probe-claim-survival/latest.json
```

S5 and S6 may be note-only. Comment the PR URL on the issue. Leave the issue open.

## Interface

- #617 owns the groupoid, 4% vector, weighting plane. You cite; you do not recompute.
- #618 owns no-power-from-Theorem-L. You cite R6 for S2’s lemma.
- #619 owns `M` polynomials, wrap/`X`. You import.
- #620 owns papers. You do not search.
- #275: contract only.
- #599: you add four axes to `E`; you do not rebuild the complexity ladder.

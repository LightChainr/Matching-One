# Probe HIGH — Does the affine-invariant shape of the threshold law have a limit?

This is a research probe. It is not a reproduction ticket, not a census checklist, and not a request to confirm that 1.55 was chart transport.

A negative theorem, a strict separation, or a proof that two analogies are unrelated is a successful outcome. Confirming a JSON you already have is not.

---

## Standing

- Does **not** enter `docs/STATUS.md`.
- Does **not** close #582, #584, #609, #610, #612, #613, #275, #276, #321, #617, #618, #619, #620, #621.
- Does **not** start a production.
- Does **not** fit `A1 N^{-ω1}+A2 N^{-ω2}`. The ratio 1.55 is retired; do not revive it as a second scale.
- Do not reopen #275. Do not score N=725. Do not become the #617 scorer (no λ-sweep, no 4% jackknife as the job).
- Literature search is #620’s job. You may use any theorem you actually know, with a quote if you lean on it; you may not open a retrieval side-quest.
- Cite vs claim. Mark conjecture. A new object you define is a claim and needs a transformation law.

Frontier: `claude/matching-one-workspace-pwr5pv` @ `8b5f9d1a`.
Briefs packet: PR #616. Sister tickets: #617 (chart groupoid — cite, do not redo), #618 (no polynomial rate from Theorem L’s union bound), #619 (exact polynomials — import if present), #621 (claim-survival atlas — orthogonal).

Submit a **new** PR against frontier. Comment the URL on the issue. Leave the issue open.

---

## The gap

Theorem L (#613 / PR #614) is a location theorem:

```text
Q_N(u) → p_c     uniformly on compact subsets of (0,1),
```

under H1–H3, for honest periodic square-cell tori with `ℓ_N / log N → ∞`.

If that is true, every interior quantile goes to the same point. Any location-and-scale invariant of `Q_N` is then a **0/0 form**:

```text
Z_N(u)  =  ( Q_N(u) − Q_N(a) ) / ( Q_N(b) − Q_N(a) )
```

for fixed anchors `0 < a < b < 1`. Numerator → 0, denominator → 0. Theorem L is silent on `Z_N`. The #618 skeleton was asked to record that silence, not to fill it.

The pipeline’s scientific sentences (`g`, `N^{-0.970}`, leftover 4%) all live *below* this quotient, in a particular affine chart of a particular 9-vector. That is observer geometry. The question this probe owns is one level up:

```text
Does the shape have a limit?
If it does, is the limit topological (Alexander + torus + duality),
or is it extra physics (a scaling function at τ=0, Cardy-like),
or does it fail to exist (no limit, grid-dependent, pair-dependent)?
```

Nobody in this repository has a theorem either way. Exact L=3 and L=4 exist and have never been looked at as shapes. Self-dual models (square bond, triangular site) have `p_c = 1/2` exactly, so Theorem L is trivial there and **shape is the only remaining object**. That is the clean laboratory.

---

## North-star question

Is there an object `S_N`, built from the threshold law, such that

```text
(1)  S_N is invariant under Aff(1) acting on the quantile function
     Q ↦ αQ + β, and you have written the action;
(2)  S_N is defined for the function F_N (or its inverse), not only
     for nine deciles;
(3)  the sequence (S_N), or (S_L) on honest tori, either
       (A) converges to a limit S_∞ you can name,
       (B) is proved not to converge,
       (C) converges only after extra hypotheses you state
           (RSW, F1, self-duality, …), and those hypotheses
           are not silently the thing being claimed;
(4)  you can say whether S_∞, if it exists, is determined by
     H1 (Alexander) plus duality plus torus topology,
     or is independent information.
```

The working suspicions below are to be **killed or promoted**, not recited.

---

## Working suspicions (adversarial)

Treat each as guilty until it has a proof or a counterexample.

**W1.** *The shape has no topological limit.*
H1–H3 plus exponential decay fix where the mass sits (`p_c`) and do not fix how it sits. `Z_L` at L=3 and L=4 already disagree by an amount that will not go to zero along a topological mechanism. If you can turn this into a theorem — even a theorem about a toy monotone family with a sharp threshold but unconstrained shape — do that. A continuum of shapes compatible with Theorem L is a high-value negative for anyone who wanted `g` to be “the” correction.

**W2.** *On self-dual lattices the shape *does* have a limit, and it is the scaling function at τ=0.*
Square-bond / triangular-site: `p_L^H = 1/2`, location is free. The only finite-size object is shape. If W2 is true, square-site `Z_L` is that same limit plus a matching-odd contamination visible as `M(1/2) ≠ 0`. Then the pair `(Z_L, M_L(1/2))` splits “Cardy-like shape” from “failure of self-matching”. This would be a new scientific sentence. It is not a fit of `ν`.

**W3.** *Nine-decile `g` is not a discretisation of any function-space tangent to `S_N`.*
Even if `S_N` exists, the pipeline’s `g` (extracted in `span{1,Q,g}` in a covariance metric on nine points, five transitions) does not lift. Then #582 is about a grid-and-metric artefact, and the leftover 4% is not something to interpret in that basis. Do not “improve `g`” by going to 17 points and fitting again. Either exhibit a lift with a transformation law, or separate the two objects.

**W4.** *The Alexander pair determines an odd shape that square-site NN cannot see in Q_G alone, and that odd shape is the leading finite-L correction to Z on square site.*
`M(1/2) = -21/64` at L=3 is the odd content at the self-dual parameter. The question is whether that number *governs* the shape of `Q_L`, or is a single moment orthogonal to shape. If you can write an exact identity on L=3,4 that expresses `Z_L` in terms of the pair `(r_b, r_w)` and then shows the even part of `Z` is small, that is a discovery. If `Z_L` is not a function of `M(1/2)`, W4 dies.

**W5.** *There is no Aff(1)-invariant S_N worth forming, because the physically rigid object is a near-critical coordinate τ = L^{1/ν}(p−p_c) acting on p, which is a different Aff(1) from the pipeline’s action on Q.*
Location-scale on occupation `p` is not location-scale on the inverse-CDF. Theorem L is about the former’s limit point. The pipeline quotiented the latter. If these two group actions are not compatible, “the invariant shape of the threshold law” is not well-posed until you choose one, and choosing the pipeline’s action is not innocent. Write the two actions as formulas and either exhibit an intertwiner or prove there is none.

Kill W1–W5 with constructions and counterexamples. Do not “evaluate” them rhetorically.

---

## Research directions (open; not a sequence of gates)

There is no mandatory order. There is no P0 reproduction gate. Use exact L=3,4 as a laboratory for *new* quantities. Import `exact_torus_enum.py` (PR #606) as a library; do not re-verify its JSON as a deliverable.

### 1. Define S, then compute it where it is exact

Write `S` as a formula. Candidates you may accept, modify, or throw away:

- a function `Z(u)` with declared anchors (and a proof that a change of anchors is a reparametrisation, not a new observable);
- an `Aff(1)`-quotient of `F` in a function-space metric you name (Wasserstein after an affine warp of the axis, Kolmogorov after warp, Fisher–Rao, …);
- a projective coordinate on the 9-vector (#579 / Fieller) only if you then lift it off the 9-grid;
- something of the pair `(G,Ĝ)` that is not a quantile of `G`.

What you may not do: take `g`’s amplitude, declare it invariant, and stop.

On L=3 and L=4, compute `S` as an exact (or bisection-to-1e-14) object. Two sizes do not make a limit. They are enough to **kill** “S is already constant” and enough to **feed** W4. Square-bond L=3 (`2^{18}` configs) is in budget if you need a self-dual comparison; square-bond L=4 (`2^{32}`) is not.

### 2. A toy monotone family with location but no shape

Construct, mathematically, a sequence of CDFs `F_N` that satisfy Theorem L’s conclusion (`Q_N(u)→p_*` uniformly on compact u) and whose `Z_N` either fails to converge or converges to an arbitrary prescribed shape. If you can, W1 becomes a theorem about location theorems in general, and percolation has to supply extra input to pin `S_∞`. That extra input is then named (RSW, arm events, self-duality, …) without pretending Theorem L contained it.

This direction is analysis. It does not need a percolation engine.

### 3. Self-dual laboratory

On square-bond L=3 exact (and L=4 only if an enumerator already exists in-tree): `p_c=1/2` is free. Compute `S`. Compare to square-site L=3 `S`. The difference is a candidate for “matching-odd contamination of shape”. It may be large, small, or incomparable because the two models are not the same experiment. If they are incomparable, say so — that is W2 dying for a good reason, not a failed plot.

Do not cite Cardy’s formula as a finite-L shape. Cardy is a limit of crossing probabilities at `p_c`, not an inverse-CDF of a homological rank.

### 4. Two Aff(1)s

Write, as formulas:

```text
action on p:     p  ↦  α p + β
action on Q:     Q(u) ↦ α Q(u) + β
action on F:     the pushforward that makes either of the above true
```

Prove they do not coincide (they do not: one warps the occupation axis, one warps the quantile axis). Then answer: which action is Theorem L invariant under? Which action did #582 quotient by? Is there a natural intertwiner along the family `N ↦ F_N`? If not, the phrase “finite-size scaling shape” has been used for two different quotients, and `g` cannot be “the” scaling correction.

A proof that there is no intertwiner is a high-value negative.

### 5. Pair-odd and Q_G

Using `(r_b, r_w)` on every L=3,4 site config, form even/odd combinations that are not `M`. Ask whether `Z_L` of `G` is a function of the odd sector. An exact identity, or an exact counterexample (two configs/measures with the same `M(1/2)` and different `Z`), kills or promotes W4. Do not summarise `M(1/2)=-21/64` as if that were the answer; that number is one moment.

### 6. Lift, or separate, `g`

Take the published 9-vector `g` as a *cited* object (PR #605 / #614). Ask only: is there an `S` from (1) whose finite-N tangent, discretised to nine deciles, is proportional to `g` on L=3,4 or on one committed histogram you re-read? Yes with a transformation law, or no with a separation. “Approximately, on five transitions, in one chart” is not a lift.

If no: #582 and the shape-limit question are different projects. Write that sentence. It is allowed to be the main result.

---

## What this is not

- Not “reproduce the #612 identity”. That identity is a cited fact.
- Not “confirm M(1/2)=-21/64”. That is a cited fact.
- Not “rescore the 4%”. If your `S` happens to make the 4% look like truncation of a function, that is a remark in direction 6, not a new exponent.
- Not a literature review. If you need a named theorem you do not have, write `BLOCKED_ON_#620: <the question>` and keep going.
- Not #275. A candidate that cannot state which Aff(1) it quotients by is not a candidate; you may record that as a contract sentence, then stop.

---

## High-value outcomes (any one is enough)

- A theorem: location theorems do not pin shape (direction 2), with the extra input percolation would need named.
- A constructed `S` with a transformation law, computed on exact L=3,4 and square-bond L=3, and a verdict on W1–W4 that is not rhetoric.
- A proof that the two Aff(1) actions do not intertwine, so “the” invariant shape was not well-posed.
- A strict separation: pipeline `g` is not a discretisation of any `S` you consider admissible.
- An exact identity or counterexample on the pair `(r_b, r_w)` versus `Z_L`.

Low-value outcomes: a survival table of old claims; a second exponent; a protocol card that restates #617; a recap of Theorem L.

---

## Claim discipline

- New objects need a definition, an action of Aff(1), and one killing test.
- Two exact sizes are not a continuum limit. They kill constancy; they do not prove `N^{-ω}`.
- Self-dual and non-self-dual models are different experiments until you have a coupling.
- GOVERNANCE §2E: two readings of the same configs are one experiment.
- Do not treat triangular-site / Cardy as square-site finite-L shape.

---

## Deliverables

New PR against `claude/matching-one-workspace-pwr5pv`:

```text
notes/probe-invariant-shape-limit-YYYYMMDD.md
    north-star verdict (A/B/C, or “not well-posed”)
    W1–W5: killed / still live / promoted, each with the argument
    the formula for S, or a proof there is none worth forming
scripts/  only for new quantities (Z, pair-odd vs Z, bond L=3 S, toy CDF family)
results/probe-invariant-shape/latest.json
```

Empty rhetoric with no formula for `S` and no theorem in direction 2 is a failed probe, not a partial atlas.

## Interface

#617: groupoid of the 9-vector pipeline. You may use its declared chart as a *negative* control (“this is what we are quotienting *out*”). You do not recompute it.
#618: no polynomial rate from the union bound. You may use that as a lemma that Theorem L does not give `|Q−p_c|`.
#619: polynomials for `M_L`. Import.
#620: papers. Do not search.
#621: survival of old sentences. Different job.
#275: closed.

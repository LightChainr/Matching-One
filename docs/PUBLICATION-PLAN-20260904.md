# Publication plan — 2026-09-04

Status: planning / manuscript coordination only. This document changes no scientific claim level, experiment authorization, Issue priority, frozen protocol, or integration status.

## 1. Decision

The repository should no longer treat `P1–P5 / R1–R3` as the only publication ontology. That portfolio remains useful as a historical research map, but several paper-sized objects now cut across those track boundaries.

The next publication program should be organized around **paper objects with bounded claims and explicit entry gates**, not around the historical order in which the research branches were created.

The immediate sequence is:

1. **P1 revision audit** — determine whether the submitted manuscript already contains the post-portfolio strengthening in PRs #549 and #550. If not, decide immediately whether those exact theorems belong in the current submission/revision rather than becoming a separate paper.
2. **P2 submission hardening** — PR #558 has produced the complete draft and closed the historical low-complexity radical gap. Do not expand the search class by default. Finish the remaining literature/provenance checks and submit.
3. **P3 manuscript assembly** — build the dependency-aware evidence table first, then perform canonical rescoring, then draft the empirical finite-size paper from already-completed production. No new size/angle/harmonic/free-exponent acquisition by default.

In parallel, run only cheap theorem/novelty gates for the next wave:

- homological-balance convergence (#276);
- R1 finite-algebra novelty check (#552);
- R2 effective-model scope freeze (#554 / PR #533);
- R3 proof-carrying elimination gate audit (#555);
- #537 three-proof-obligation asymptotic gate.

Do **not** draft P4 or the flagship P5 yet unless their stated theorem/identifiability gates close.

## 2. Why the publication map needs revision

The 2026-09-01 publication portfolio was a useful snapshot, but the repository acquired stronger theorem assets immediately afterward.

### 2.1 P1 is stronger than the portfolio snapshot

PR #549 upgrades the branching story from a finite counterexample to a growing-family theorem: one complete-unbranched-survival equivalence class contains at least `k+1` distinct branching-predictive classes on `8k` future vertices.

PR #550 then proves that even a bounded summary containing the complete safe-subset polynomial, `H2/b2`, and the radius-1 terminal neighbourhood can fail for a frozen depth-2 compositional language. It also records the cut-law equivalence boundary and a general realization lemma for the finite plane two-terminal category.

Therefore the current P1 theorem chain is:

```text
rank-one cut-network representation
  -> finite branching no-go
  -> growing predictive-class lower bound
  -> bounded-summary non-compression
  -> realization / cut-gauge clarification
```

This is a revision-level change in paper strength, not another exploratory side result.

### 2.2 Exact topology creates a separate probability-theory paper object

The digital-Alexander rank bridge gives, on honest periodic square-cell tori,

```text
r_black + r_white = 2,
2 q = r_black - r_white.
```

Consequently the finite matching root is the exact **homological balance point**

```text
M_N(p_N)=0
<=> E_{p_N}[r]=1
<=> P_2(p_N)=P_0(p_N).
```

Issue #276 already isolates the clean theorem target `p_N -> p_c` using monotonicity, sharpness, and torus homology rather than finite-size-scaling assumptions. If that convergence theorem closes, it should become its own probability/topology paper rather than being absorbed into P3 or P5.

### 2.3 Original-U computational statistics is not merely an R3 appendix

PRs #530/#531 contain a coherent methods object:

- the full influence function for the pooled-moving-root normalized observable;
- exact nuisance/common-root handling;
- optimal rank-sector importance allocation and a robust dyadic bound;
- a concrete lower-bound diagnosis showing that ordinary proposal reweighting does not remove the signed-cancellation barrier;
- thermal-gauge-invariant rank-sector quotient coordinates;
- analytic-subtraction targets that change the integrand rather than only the proposal.

This work should be kept distinct from the proof-carrying elimination framework. The two methods programs may share reproducibility standards, but they solve different mathematical problems.

### 2.4 #537 is an asymptotic theorem program, not merely a P5 subgate

The finite contact/thermal program has already produced exact finite obstructions, including rejection of a pure-thermal rank-one cancellation route. The live issue now allows only three proof obligations:

1. contractible-collar quotient identity;
2. bounded normalized pivotal domination;
3. near-critical uniform transport from exact `p_c` to the pooled root.

That is a self-contained asymptotic theorem program. If the three obligations close, it deserves an independent structural/asymptotic manuscript. If any fails, the resulting named obstruction can itself sharply delimit the flagship P5 mechanism program.

## 3. Paper-object map

### A. Immediate manuscript/submission objects

#### A1. P1 — continuation representations and predictive non-compression

**Current status:** submitted manuscript exists; post-submission theorem strengthening is present in open Draft PRs #549/#550.

**Central claim:** unbranched survival information does not close delayed branching; the failure amplifies to an unbounded number of exact branching-predictive classes, and natural bounded summaries still fail inside the finite two-terminal continuation category.

**Next action:** perform a submission-version audit. Compare the submitted manuscript against #549/#550 and classify each new result as:

- already included;
- revision-worthy central theorem;
- supplement/appendix only;
- separate future result.

**Stop rule:** no larger parallel products, no new descriptor ladder, no broader series-parallel census merely to increase witness counts.

#### A2. P2 — certified bounded algebraic exclusion

**Current status:** complete draft merged by PR #558.

**Central claim:** complete, certified exclusion within scientifically motivated finite algebraic-complexity classes; not non-algebraicity or transcendence.

**Next action:** submission hardening only:

1. verify the remaining primary references or weaken the historical-catalog wording;
2. if inexpensive, add an independent exact-filter implementation as referee insurance;
3. regenerate manuscript evidence tables from committed artifacts;
4. submit to the selected venue.

**Stop rule:** do not raise degree/height by default. A larger search requires a new declared scientific hypothesis and an approach-resolution argument showing the search remains informative.

#### A3. P3 — matching-odd finite-size response beyond scalar closure

**Current status:** evidence nucleus complete; manuscript ticket #553 open.

**Central claim:** the matching-odd finite-size signal is nonzero and orientation-sensitive; frozen H4-like transfer outperforms declared H8/H12 aliases in specific prospective designs; stronger scalar closure models fail; the remaining correction state is multicomponent/observer-dependent rather than one fitted scalar amplitude.

**Mandatory order:**

1. dependency-aware evidence table;
2. canonical nullspace-safe rescoring;
3. figure set;
4. results prose;
5. literature bridge and discussion.

**Stop rule:** no new size, angle, harmonic, or free exponent unless a specific reviewer-grade alternative cannot be tested with existing assets.

## 4. Next-wave theorem-gated paper objects

### B1. Homological balance points for site percolation

**Source:** PR #271 + issue #276.

**Paper gate:** prove, for the declared growing torus sequence,

```text
p_N -> p_c
```

from monotonicity/sharpness/homology, with a clean scope for nondegenerate square-cell tori.

**Minimum publishable theorem stack:**

1. exact rank/matching identity;
2. strict monotonicity and uniqueness of `p_N`;
3. subcritical `P_0 -> 1` and supercritical `P_2 -> 1` localization;
4. convergence of the balance point;
5. Russo/rank-birth influence interpretation as a corollary.

**Do not require for first paper:** a finite-size shift exponent, CFT operator identification, or optimal-estimator claims.

### B2. #537 thermal/contact asymptotics

**Paper gate:** close all three named proof obligations, or obtain a single structural counterexample that permanently kills the proposed asymptotic closure.

**Positive paper object if gate closes:** a theorem connecting a finite exact contact/rank decomposition to a controlled near-critical pooled-root asymptotic observable.

**Negative paper object if a gate fails:** a precise obstruction theorem identifying why the natural local/thermal closure cannot propagate to original `U`.

**Stop rule:** no third-size fitting, no new descriptor grid, no top-up of N145.

### B3. R1 — 15-state ordered serial algebra

**Paper gate:** external novelty memo against finite semigroup / diagram / partition / planar algebra literature.

If novel enough, draft the pure finite-algebra paper immediately and keep the probability bridge as optional future strengthening. If already known up to isomorphism, stop the pure-algebra manuscript path and return to the probability-comparison/obstruction route.

### B4. R2 — closed-source directed capillary model

**Paper gate:** define one explicit effective model under which the two-cloud, interface, Bessel/Karlin-McGregor, and Catalan-Toeplitz statements are unconditional.

The title and theorem statements must name the effective model, not the unrestricted square-lattice original-U problem.

**Stop rule:** do not hold this paper for unresolved physical beta-cloud, endpoint-transfer, rank-sector-normalization, second-thermal/root, or full original-U bridges.

### B5. R3a — proof-carrying elimination

**Paper gate:** at least three scientifically distinct instances implementing the same reusable fail-closed contract, plus a convincing fourth-case/API thought experiment or demonstration.

Keep this manuscript about machine-verifiable bounded hypothesis elimination. Do not force the moving-root influence-function program into the same paper.

### B6. R3b candidate — inference for rare signed moving-root observables

**Sources:** PRs #530/#531 and related exact statistical infrastructure.

**Current status:** strong methods nucleus, not yet a paper by default.

**Paper gate:** at least one of:

- a successful conditional-integration / Rao-Blackwell implementation showing end-to-end variance reduction on the same target with root/nuisance uncertainty retained;
- a second scientifically different moving-root signed observable demonstrating that the influence-function and variance-barrier framework is not Matching-One-specific;
- a theorem characterizing when proposal-only importance sampling cannot overcome the signed-integrand barrier, together with a constructive integrand-changing remedy.

Until then, keep the material as a methods dossier rather than expanding ordinary MC.

## 5. Flagship / incubation tracks

### C1. P4 — operator-structure non-identifiability

Do not draft yet. First extract a bounded realization/non-identifiability proposition with typed definitions of microscopic exact state dimension, observer-visible predictive/Hankel dimension, and continuum-operator representation dimension. Existing width examples are evidence for theorem design, not a manuscript substitute.

### C2. P5 — identifiability of the original homology-marked observable

Keep as flagship and do not force an early progress-report paper.

Draft only after at least one of these occurs:

1. two actual candidate forward columns are written in the same original-U normalized data space;
2. a rank calculation proves the candidates observationally equivalent under current assets;
3. the rank calculation identifies one unique missing coordinate and a frozen information/acquisition plan;
4. #537 closes a theorem or obstruction that materially changes the candidate space.

Negative identifiability is an acceptable final scientific outcome. The gate is coherence, not positivity.

## 6. Concrete execution schedule

This is an ordering rule, not a wall-clock promise.

### Pass 1 — current manuscript closure

**Primary:** P1 revision audit.

Deliverable: one comparison memo between submitted P1 and PRs #549/#550, with an explicit `INCLUDE_IN_REVISION / SUPPLEMENT / SEPARATE` decision per theorem.

**Then:** P2 submission hardening.

Deliverable: literature/provenance checklist closed; final evidence regeneration; submission-ready source.

**Then:** P3 dependency table and canonical rescore.

Deliverable: one table that makes dependency/chronology impossible to misread; one frozen score table from the current scorer; only then manuscript prose.

### Pass 2 — cheap parallel gates

These may proceed without new stochastic production:

- #276 proof memo: isolate exactly which sharpness/wrapping lemma is missing for `p_N -> p_c`;
- #552 novelty memo;
- #554 effective-model definition;
- #555 three-case audit;
- #537 proof-obligation status table.

Each gate must return either `READY_FOR_MANUSCRIPT`, `ONE_NAMED_THEOREM_MISSING`, or `STOP/RESCOPE`.

### Pass 3 — choose exactly one next theorem paper

After Pass 2, choose the paper with the largest manuscript delta per missing theorem, not the paper with the most accumulated artifacts.

Default tie-break:

1. homological balance if the convergence proof is one standard-input lemma away;
2. R1 if novelty clears;
3. R2 if the effective model is already closed;
4. #537 if two of three obligations are closed and the third has a precise target;
5. R3a if the reusable three-case gate is genuinely satisfied.

## 7. Manuscript hygiene rules

Every manuscript dossier should contain a one-page claim ledger with:

- central theorem/falsification claim;
- exact hypothesis/observer/model class;
- evidence dependency groups;
- prospective / held-out / post-reveal status where statistical;
- strongest explicit nonclaim;
- one stop rule;
- one sentence describing what new result would materially change the paper.

For empirical papers, derived projections of one batch stream are not independent evidence. For exact papers, search-class completeness and implementation verification must be separated. For effective-model papers, every theorem statement carries its model qualifier. For identifiability papers, model compatibility is not model identity.

## 8. Repository actions proposed by this plan

This planning change itself should remain documentation-only.

After review, follow-up repository actions should be separate and small:

1. add the P1 submitted-version audit as a manuscript note or Issue comment;
2. update #553 only after the dependency table exists;
3. open a dedicated homological-balance manuscript ticket only when #276's convergence proof gate is classified as near-closed;
4. keep R3a and the signed-observable statistics dossier separate in future portfolio refreshes;
5. do not change scientific priority labels merely because a track is publication-ready.

The publication program succeeds when each new unit of work changes a paper claim, closes a manuscript gate, or removes a real referee vulnerability. Additional census, sizes, descriptors, or parameter scans that do none of those things are not publication progress.

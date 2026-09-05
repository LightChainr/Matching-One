# Matching One publication portfolio

**Snapshot basis:** 2026-09-01.  
**Purpose:** convert the repository from an ever-growing research DAG into a small set of paper claims with explicit theorem, falsification, validation, and stop gates.

This document is **not** the live execution queue and does not override Issue labels, `docs/STATUS.md`, frozen protocols, or production authorization. Publication readiness, scientific importance, integration state, and execution priority are separate coordinates.

The unit of progress here is not a PR, Issue, simulation, or additional fit. It is a change in a paper claim:

```text
conjecture -> theorem / counterexample
compatible -> rejected / identified / unidentifiable
analogy -> explicit observable map
branch result -> independently reproduced / manuscript lemma
many results -> one bounded central claim
```

A negative or impossibility result is publication progress when it removes a natural mechanism class or establishes a genuine information boundary.

---

## Portfolio at a glance

| ID | Candidate paper | Readiness | Potential impact | Current publication gate |
|---|---|---:|---:|---|
| P1 | **Continuation representations and nonclosure for topological birth processes** | high | high | add a genuine minimality / no-compression result, or sharply delimit why one is unavailable |
| P2 | **Certified exclusion of bounded algebraic relations for a critical threshold** | high | medium-high | canonical provenance + search-class motivation + manuscript-level completeness statement |
| P3 | **Matching-odd finite-size response beyond scalar closure** | medium-high | high | dependency-aware evidence synthesis + canonical rescoring / protocol audit + literature positioning |
| P4 | **Finite positive dynamics and non-identifiability of operator structure** | medium-high | high | promote finite examples into a general realization / de-identification proposition |
| P5 | **Identifiability of a homology-marked critical observable** | medium | highest | actual candidate forward columns in the same original-`U` data space; then a rank decision |
| R1 | **A finite terminal algebra for connectivity reductions** | high assets, low paper closure | medium | connect the algebra to a probability comparison theorem or a precise impossibility theorem |
| R2 | **Closed-source two-cloud and capillary effective theory** | mixed; overall physical closure low | potentially high | freeze a reduced effective model and publish only statements proved in that model |
| R3 | **Proof-carrying computational elimination / topological-observable statistics** | medium | method/companion | generalize beyond project-specific utilities and benchmark on multiple scientifically distinct cases |

Two different orderings should be kept explicit.

**By probability of producing a defensible manuscript soon:**

```text
P1 cut-network/nonclosure
> P2 certified algebraic exclusion
> P3 matching-odd/scalar-closure paper
> P4 positive-dynamics de-identification
> P5 original-U flagship
> R1 terminal algebra
> R2 capillary full-model claim
```

**By ultimate scientific upside:**

```text
P5 original-U flagship
> P1 cut-network/nonclosure ~ R2 capillary if the physical bridge closes
> P3 matching-odd/scalar-closure
> P4 positive-dynamics de-identification
> P2 certified exclusion
> R1 terminal algebra without a probability theorem
```

The project should therefore work on more than one publication horizon: finish one or two mature papers while preserving the highest-impact flagship as a theorem/identifiability program.

---

# P1 — Continuation representations and nonclosure for topological birth processes

## Proposed central claim

A rank-one embedded torus continuation problem admits an exact cut-network representation, while complete unbranched survival information does not in general close branching continuation. The natural finite continuation state is therefore richer than a scalar survival law or a short list of hazard summaries.

A manuscript title could be:

> **Continuation representations and nonclosure for topological birth processes on embedded graphs**

## Existing theorem / evidence nucleus

The line already contains unusually clean components:

- PR #491: in the stated embedded rank-one scope, cut along a simple occupied essential cycle and contract old occupied components; future ambient rank reaches two exactly when the two cut boundaries connect in the resulting planar **vertex** network;
- pair triggers are bipartite in that representation and arise from direct L/R vacancy edges or neutral-component bicliques;
- genuine minimal triples survive as longer switchable terminal structures;
- keeping the cut fixed gives an update-closed network-valued continuation representation;
- PR #435: two finite states can have exactly the same complete unbranched survival law yet different delayed-fork continuation, with an exact branching gap; the recomputed survival-signature process also fails a strong Markov closure test in the declared finite setting;
- PR #492: complementary dual-cycle blocker certificates on the real N425 witnesses;
- #334 contains a prospective population intervention that rejected two frozen low-dimensional residual projections; that experiment is useful motivation, but it is not part of the exact theorem proof.

## Highest-value tasks

### P1.1 Freeze the paper theorem scope

Write one theorem section independent of the Matching-One repository history:

```text
embedded graph assumptions
rank-one checkpoint assumptions
cut construction
network state
future-rank equivalence
pair-trigger corollary
branch/update rule
```

Every assumption should be necessary, justified, or accompanied by a counterexample when removed.

### P1.2 Seek a real no-compression theorem

The best possible strengthening is not another descriptor. Try to prove one of the following, in descending value:

1. for a declared family of embedded graphs, no continuation state built from any fixed number `k` of specified scalar summaries can be update-closed;
2. a lower bound on the number of distinguishable continuation states / predictive equivalence classes that grows with graph size;
3. a narrower theorem showing that complete unbranched survival laws, all bounded-horizon survival moments, or a named finite scalar hierarchy still fail branching closure.

Do **not** claim a continuum-dimension lower bound from a finite combinatorial state count.

If no clean growing lower bound can be proved, record the obstruction and stop rather than manufacturing additional scalar descriptors.

### P1.3 Make cut dependence mathematically explicit

Determine whether different admissible occupied cuts give:

- canonically isomorphic network states;
- only equivalent future event laws;
- or genuinely different state presentations requiring a quotient / covariant rule.

This is important if any later paper uses the network as a physical mark rather than only as a proof device.

### P1.4 Literature positioning

Connect the result precisely to:

- finite-state lumpability and predictive-state representations;
- network reliability / two-terminal vertex connectivity;
- topological percolation / homology birth processes;
- graph cut and transfer-state methods.

The literature section must state what theorem is imported and what is genuinely new. Similar vocabulary is not enough.

### P1.5 Manuscript engineering

Build a paper-sized artifact set:

- one minimal branching counterexample figure;
- one cut-network construction figure;
- one real N425 mechanism figure/table;
- a theorem dependency map;
- a machine-verifiable exact supplement.

## Stop rule

Do not launch new population Monte Carlo for this paper unless a theorem derived from the cut state produces a **predeclared population prediction not reconstructible from existing archives**. The paper can succeed as an exact representation + nonclosure paper.

---

# P2 — Certified exclusion of bounded algebraic threshold relations

## Proposed central claim

For explicitly frozen polynomial/constant families and method-specific threshold intervals, complete exact searches exclude all low-complexity relations in the declared class. This is a bounded certified exclusion result, not evidence of transcendence.

Possible title:

> **Certified exclusion of low-degree algebraic relations for a critical threshold**

## Existing nucleus

`main` already contains a rare level of computational completeness:

- sign-normalized primitive degree-1 through degree-4 search classes at coefficient height 100;
- complete cubic and quartic interval treatment across multiple frozen method intervals;
- quartic space size `157,309,446,881`;
- exact near-candidate certification and Sturm/root filtering;
- standard-constant and lattice-native candidate families;
- an exact Kagome positive control;
- look-elsewhere counts and precision/endpoint stability audits.

This is substantially stronger than an informal PSLQ search.

## Highest-value tasks

### P2.1 Reconcile provenance before writing the claim

Issue #1's historical body still encodes a provenance dependency that does not visually match the large completed census now present on `main`. Build one canonical table:

```text
method / source
point estimate
quoted or constructed interval
boundary conditions / observable
provenance location
why it is kept separate from other methods
```

If source intervals disagree at their stated precision, preserve the disagreement. Do not manufacture a single preferred ultra-narrow interval.

### P2.2 State the completeness theorem for the search class

The paper needs a theorem/algorithm statement of the form:

> every primitive sign-normalized polynomial in class C(degree, height) is either certified away by the screening bound or enters an exact root test; therefore the census is exhaustive for C.

Separate:

- counting completeness;
- numerical/interval screening completeness;
- exact root certification;
- implementation verification.

### P2.3 Motivate the search bounds scientifically

The degree and coefficient-height bounds cannot look arbitrary. Survey exact-threshold traditions and previously proposed/simple algebraic forms, then explain what natural conjecture class height 100 is meant to cover.

The manuscript should say explicitly what scientifically motivated family is being ruled out.

### P2.4 Strengthen calibration only where it changes credibility

The existing false-positive and positive-control machinery is useful. Additional calibration is justified only if it addresses a likely reviewer objection, for example:

- interval perturbation sensitivity for higher-degree near hits;
- an independent implementation of the final exact filter;
- deterministic synthetic constants constructed to have known relations at comparable height.

Do not expand the constant library merely to create more negative counts.

### P2.5 Write the negative result as a bounded theorem

The conclusion must remain:

```text
excluded within declared finite complexity class
!= transcendental
!= non-algebraic
!= no exact representation of another type
```

## Stop rule

Once provenance, completeness, search-class motivation, and independent verification are manuscript-ready, **stop increasing degree/height by default**. A larger census without a new scientific hypothesis does not materially strengthen the paper.

---

# P3 — Matching-odd finite-size response beyond scalar closure

## Why this should be a separate paper

The repository already has enough empirical structure to support a paper whose conclusion is deliberately weaker than operator identification:

> a robust orientation-sensitive matching-odd finite-size response is observed across independent designs, while several natural one-scalar finite-size closures fail prospectively or exactly.

This paper should not wait for the original-`U` continuum identity problem.

Possible title:

> **Matching-odd finite-size response beyond scalar correction closure in square-lattice site percolation**

## Existing nucleus

Potential core components include:

- independent P43/P57 global-zero versus fixed-H4 evidence;
- prospective N185/N265 new-geometry block;
- norm-5 H4/H12/H8 discriminator;
- N145->N290 held-out failure of a one-multiplier curve law;
- derivative-channel survival/failure pattern, including the pure `P4[S']` failure;
- norm-5 thermal-jet rejection of scalar width/rank-gap closures;
- N100/N400/N900 shape-redistribution and width results as later multicomponent context;
- observable-channel errata and covariance-nullspace QA as explicit evidence-discipline cases.

## Highest-value tasks

### P3.1 Build one dependency-aware evidence table

For every displayed score, record:

```text
raw dependency block
prospective / held-out / post-reveal
primary or derived
observable semantics
geometry
hypothesis class actually tested
current interpretation after errata
```

The paper must not turn multiple derived views of one batch stream into independent replication.

### P3.2 Canonically rescore only the claim-bearing statistics

Run the current nullspace-safe / typed scorer on the minimum set of historical sufficient statistics needed for the manuscript. The goal is a stable canonical table, not another broad reanalysis campaign.

Any result whose interpretation is cutoff-sensitive should be described that way.

### P3.3 Define the central falsification claim

Do not write "H4 is the field." A stronger defensible structure is:

1. global zero is strongly disfavored by independent primary evidence;
2. a frozen H4-like orientation transfer succeeds against declared aliases in specific designs;
3. several stronger scalar closures fail;
4. therefore the unresolved object is multicomponent / observer-dependent rather than a single fitted correction amplitude.

### P3.4 Literature bridge

Position the result against finite-size scaling, irrelevant operators, torus/wrapping observables, anisotropic corrections, and percolation topology literature. Explicitly separate:

- scaling exponent evidence;
- modular/angular fingerprint;
- lattice-observable selection;
- continuum operator identity.

### P3.5 Prebuild manuscript figures from existing data

Prefer figures that communicate decisions rather than model proliferation:

- global-zero vs H4-like primary evidence;
- prospective harmonic alias discriminator;
- one-multiplier closure failure;
- thermal-jet residual mode;
- chronology/erratum panel showing why typed semantics matter.

## Stop rule

Do not add another size, angle, harmonic, free exponent, or modulus unless it is required by a **specific pre-paper reviewer-level alternative** that existing data cannot address.

---

# P4 — Finite positive dynamics and non-identifiability of operator structure

## Proposed central claim

Finite positive/reversible dynamics can reproduce signatures that are often interpreted as evidence for non-semisimple/Jordan structure; observer-visible realization dimension and microscopic state dimension are also generator-dependent. Therefore several finite correlation fingerprints are non-identifying without an explicit microscopic-to-observer map.

Possible title:

> **Non-identifiability of operator structure from finite positive correlation dynamics**

## Existing nucleus

P398 supplies the strongest core:

- exact finite positive transfer/dynamics constructions;
- current deletion: removing stationary probability current does not remove the key inversion/crossing phenomenon;
- reversible positive mixtures can therefore produce part of the phenomenology;
- retaining instantaneous current direction still fails to recover propagation because hidden reversible-force geometry remains;
- controlled rate interventions separate common time reparametrization from genuine cross-ray response.

P250 provides a useful complementary warning:

- a finite-window low-rank description can coexist with an exact >=100-mode lower bound for the complete raw spatial endpoint series;
- compressed commutators can be produced by projection leakage even when underlying operators commute.

The two lines should be combined only at the level of a general **identifiability / realization** principle, not as independent votes about the same physical mechanism.

## Highest-value tasks

### P4.1 Prove a general realization proposition

The most important missing step is to move beyond width-4/5/8 examples. Seek a theorem such as:

- a nontrivial class of crossing / fast-start-slow-tail correlation fingerprints is realizable by ordinary reversible positive mixtures;
- or a class of finite Hankel/correlation signatures cannot distinguish semisimple from near-Jordan realizations under finite tolerance;
- or a quantitative approximation theorem near a Jordan collision with positive semisimple models.

Even a sharply bounded finite-dimensional proposition would significantly strengthen the paper.

### P4.2 Formalize the three dimensions

Use precise definitions for:

```text
microscopic exact state dimension
observer-visible predictive / Hankel dimension
continuum/operator representation dimension
```

Then prove which implications fail. Avoid using "state dimension" without a typed observer/generator/context.

### P4.3 Literature positioning

Connect to:

- positive realization theory;
- hidden Markov / lumpability / system identification;
- Hankel rank and minimal realization;
- nonnormal dynamics and pseudospectral ambiguity;
- LCFT/Jordan finite-size diagnostics only where the observable map is explicit.

### P4.4 Use existing interventions as examples, not more votes

Do not extend width or scan rates by default. The current finite systems should become worked counterexamples / benchmark constructions supporting the general proposition.

## Stop rule

No new width-10, mark, source, or rate grid unless a formal proposition has a concrete missing case that one bounded computation can decide.

---

# P5 — Identifiability of a homology-marked critical observable

## Flagship question

This is the highest-impact line but should not be forced into a premature "operator identification" paper.

The correct paper question is:

> Which continuum/mechanism hypotheses make distinct predictions for the **same normalized homology-marked observable**, and what can the existing data identify once normalization, moving root, nuisance amplitudes, and covariance support are treated correctly?

Possible title:

> **Identifiability of a homology-marked critical observable**

## Current state

The repository has already learned an important negative fact: observer labels and spin labels do not automatically provide a map into original square-site `U`. The current #275 work has reached a current-assets identifiability boundary for named candidates because candidate-specific same-source restricted-trace thermal jets are not yet supplied in the required original-`U` data space.

The flagship therefore needs **theory columns before new samples**.

## Highest-value tasks

### P5.1 Construct actual candidate forward columns

For each named candidate family to be compared, provide in the same convention:

```text
named microscopic / continuum source
restricted rank-sector traces or q/E coordinates
thermal 0/1-jets
physical partition normalizer
rank-1 denominator
pooled moving root and root counterterm
cross-geometry amplitude/phase nuisance class
map to original U
```

A common spin value or modular label is not a forward prediction.

### P5.2 Produce an identifiability certificate before acquisition

Using existing covariance, compute the nuisance-profiled column spaces.

Required outputs:

- numerical and symbolic rank where possible;
- principal angles / intersection dimension;
- treatment of deterministic covariance-null directions;
- explicit missing coordinate if spaces coincide on current assets.

Decision:

```text
same profiled column space
    -> UNIDENTIFIABLE_WITH_CURRENT_ASSETS
    -> no new broad acquisition

distinct full-rank directions
    -> one frozen score on existing covariance/data

unique missing coordinate proven
    -> design exactly one acquisition that measures it
```

### P5.3 Keep #537 as the theorem side, not another fit side

The full original-`U` thermal/contact route is now proof-facing. Progress means proving or refuting the named obligations, especially the exact/local quotient, normalized pivotal domination, and near-critical pooled-root transport required for the desired asymptotic conclusion.

A counterexample is valuable: it changes the flagship from a closure theorem into a no-go/extra-state theorem.

### P5.4 Integrate computational reachability only after the target is fixed

PRs #530/#531/#532 show that ordinary sampling can be a poor match to the signed original-`U` estimating problem and that analytic subtraction / quotient coordinates may matter greatly. Use these results to support an estimator only after the physical observable and candidate columns are fixed.

Do not change the observable to obtain lower variance.

### P5.5 Publication outcomes are allowed to be negative

Any of the following can anchor a strong paper:

- two major candidate mechanisms become genuinely distinguishable and one fails;
- the candidates are provably observationally equivalent under the current observable;
- a natural thermal/pivotal closure is disproved, requiring an additional state variable;
- a full forward map is established but the data are underpowered in a quantifiable direction.

## Stop rule

No third candidate, new angle, new source, descriptor, modulus, or free amplitude class may be added post-reveal to rescue a failed comparison. New acquisition is allowed only after the rank calculation identifies one unique missing coordinate and a frozen power/information calculation shows that measuring it is useful.

---

# R1 — Terminal algebra: reserve until it touches probability

The 15-state ordered serial algebra and its subsemigroups, ideals, Green relations, centralizers, actions, congruences, and operator semigroups are mature exact assets. That maturity alone does not guarantee a compelling paper.

## Only high-value next tasks

1. identify a structural theorem linking algebraic classes to concrete connectivity semantics;
2. connect #13/#14/W5 to a law-preserving local replacement, stochastic domination, Strassen-type coupling, reliability inequality, or rigorous threshold comparison;
3. if the finite algebra cannot support such a comparison, prove a precise obstruction theorem;
4. compare the resulting semigroup with known diagram/partition/planar algebras and finite semigroup classes.

## Explicit stop

Do not continue adjacent census tasks merely because another finite classification is cheap. No new catalogue is publication progress unless it resolves one of the structural questions above.

---

# R2 — Closed-source two-cloud / capillary effective theory

PR #533 contains a large amount of elegant asymptotic and combinatorial structure, but its own maintenance boundary still distinguishes a small set of C1 statements from unresolved physical interfaces. The correct publication strategy is therefore **scope reduction**, not promotion of the full model.

## Viable task direction

Define a closed-source / directed / axis effective model whose assumptions are explicit and then isolate the statements actually proved there, for example:

- two-cloud asymptotics and root chart;
- interface generating functions;
- Bessel / Karlin-McGregor determinant structure;
- Catalan-Toeplitz resummation in the precise relaxed class;
- the strongest rigorously justified sign regime;
- the directed positive-subfamily threshold.

## Required strengthening

- a dependency table separating unconditional exact identities, conditional lemmas, and unresolved physical endpoint transfer;
- independent verification of the key determinant / resummation identities;
- a literature bridge to directed interfaces, nonintersecting paths, Toeplitz/Bessel asymptotics, and cluster/contour expansions;
- a precise theorem statement that never silently upgrades the effective model to the full square-lattice original-`U` problem.

## Stop rule

Do not market the full Matching-One capillary theorem until the physical beta-cloud, endpoint transfer, rank-sector normalization, thermal/root and original-`U` bridges required by that claim are actually closed.

---

# R3 — Proof-carrying elimination and topological-observable statistics

This is a **companion/method** track, not yet a default standalone paper.

A methods paper becomes plausible only if the project-specific utilities are generalized into a reusable framework such as:

```text
scientific hypothesis
-> frozen bounded parameter class
-> exact / interval / ideal / SOS / covariance-support certificate
-> fail-closed machine-verifiable verdict
```

or, on the statistics side:

```text
moving-root estimating equation
-> influence function with shared nuisance/root
-> structured rare-event / signed-integrand variance analysis
-> analytic subtraction / conditional integration
-> auditable finite-sample or asymptotic guarantee
```

## Publication gate

Demonstrate the framework on at least three scientifically different cases, for example:

- threshold relation exclusion;
- finite thermal/rank elimination;
- terminal/polynomial identity or covariance-nullspace decision.

The method must provide a reusable theorem/API/verification contract, not only scripts that happened to work in this repository.

---

# Publication DAG

The portfolio should be treated as a dependency graph rather than one mega-paper:

```text
P1 cut-network/nonclosure --------------------> standalone manuscript

P2 certified algebraic exclusion ------------> standalone manuscript

P3 global matching-odd + scalar failures ----> standalone empirical/statistical-physics manuscript
                                      \
                                       \ contextual motivation
                                        v
P5 original-U identifiability <---- #275 forward columns
            ^                      + #537 theorem/counterexample
            |                      + P530-532 estimator reachability
            |
P4 realization/de-identification -- conceptual guardrail, not evidence vote

R1 terminal algebra --(probability theorem required)--> possible standalone manuscript
R2 capillary --(effective-model scope closure)-------> possible standalone manuscript
R3 proof/statistics --(generalization required)------> methods companion
```

P1, P2, P3, and P4 should not be held hostage by the final outcome of P5.

---

# How an Agent should choose work from this portfolio

A recurring Agent should not simply pick the highest-impact paper. For each candidate task, estimate:

```text
paper_delta      = how much the task changes manuscript readiness
information_gain = how many plausible claims/mechanisms it distinguishes
reusability      = whether the result supports more than one section/paper
cost             = compute + implementation + review complexity
leakage_risk     = post-selection / semantic / dependency risk
```

Prefer high `paper_delta * information_gain / cost` with low leakage risk.

Default ordering of work types:

```text
exact theorem / existing-data discrimination
> deterministic counterexample
> identifiability / rank calculation
> bounded reproduction / provenance repair
> literature-to-object theorem bridge
> small frozen pilot
> large new acquisition
```

Every run should select **one PRIMARY paper task** and at most two reserve tasks.

---

# Manuscript-readiness checklist

A track should transition from research expansion to manuscript writing when it has:

- one sentence beginning `We show that ...` that is stable under current evidence;
- a finite list of theorem/claim dependencies;
- the strongest natural alternatives either rejected or explicitly outside scope;
- clear independent/dependent evidence accounting;
- a reproducible artifact set;
- a literature position based on actual imported theorems, not analogy;
- explicit nonclaims;
- no obvious need for another unconstrained size/angle/descriptor scan.

At that point the highest-value task is usually writing, figure construction, independent reproduction, or external review — not adding another result.

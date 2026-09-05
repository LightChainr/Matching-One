# Matching One

> **A computational research program on square-lattice site percolation, exact matching topology, finite-size response, and the problem of identifying the state behind the signal.**

Matching One started from the square-site matching identity and a numerical threshold problem. It has since become a broader program: exact topology and finite algebra are used to define observables; frozen finite-size experiments test predictive laws; failures are used to enlarge or reject state descriptions; only after those steps do we attempt continuum/operator identification.

The exact structural anchor is

\[
p_c^{\mathrm{site}}(\mathbb Z^2)
+
p_c^{\mathrm{site}}(\mathrm{NN+NNN})
=1.
\]

The main empirical fact is also simple to state:

> **There is strong evidence for a nonzero, orientation-sensitive matching-odd finite-size structure. There is not yet a unique microscopic or continuum identification of that structure.**

That distinction is the organizing principle of the repository.

---

## The observable is topological

For a site configuration \(\omega\) on a torus, let

\[
r(\omega)=\operatorname{rank}\operatorname{im}
\left[H_1(K(\omega))\to H_1(T^2)\right]\in\{0,1,2\}.
\]

Define the centered topological variable

\[
X=r-1\in\{-1,0,+1\}.
\]

Then the two canonical unmarked source coordinates are

```text
A_top = <X>   = P2 - P0      # matching-odd coordinate
E_top = <X^2> = P2 + P0      # topology-even partner
```

and configurationwise

```text
X^3 = X.
```

So the unmarked ambient-rank source has an exact two-dimensional nonconstant algebra. Any additional state needed by the full thermal/finite-size response must enter through something that this coarse rank variable does not retain: thermal/history dependence, projective line or subgroup data, local/defect structure, geometry, or another bulk response direction.

This exact source algebra is **not** itself a two-state CFT/Jordan module. The RG/continuum action on these observables is a separate problem.

---

## Scientific stack

```text
exact topology / observable semantics
              ↓
typed sufficient statistics + covariance
              ↓
prospective finite-size prediction
              ↓
mechanism elimination / predictive-state tests
              ↓
observable identifiability
              ↓
continuum / operator naming
```

A failure high in this stack does not erase a lower-level result. Conversely, a successful finite fit does not automatically identify a continuum field.

---

## What is established, and what is not

| Question | Current answer | Important boundary |
|---|---|---|
| Is there a global matching-odd orientation signal? | **Yes, strongly supported.** Independent primary blocks disfavor global zero; new geometries reproduce the signal. | This establishes a finite-size structure, not a unique field. |
| Does the frozen H4-like transfer beat the tested H12/H8 aliases? | **Yes.** The norm-5 prospective discriminator strongly favors H4 over those frozen aliases. | The N325/N425 child block alone does not reject zero. |
| Does one scalar finite-size correction close the full curve? | **No.** Pure `S'`, one-multiplier curve transfer, scalar rank-gap/width and related simple closures fail. | Some center/slope/root relations survive after the scalar curve law fails. |
| Is q=2 the surviving norm-4 mechanism? | **Specific scalar/common-generator q=2 versions have been rejected.** | Rejection of those models is not a theorem excluding every semisimple realization. |
| Is Jordan identified? | **No.** Jordan/log families repeatedly survive some scores, but finite-noise ordinary models can approach a Jordan collision. | Compatibility is not identity. An orthogonal structural/physical fingerprint is required. |
| Is the relevant state low-dimensional in every sense? | **No.** State dimension is observer-, generator- and context-dependent. | Endpoint/spatial rank, dilation transfer rank and branching state are different realization questions. |
| Is a simple closed form for square-site `p_c` known? | **No.** Large bounded low-complexity relation families have instead been exactly excluded. | Bounded exclusions do not imply transcendence. |

---

## Headline predictive evidence

### 1. New-geometry matching-odd test

The prospective N185/N265 block gives

```text
x = 21/4 H4-like:  chi2 =  3.04598 / 2
zero:               chi2 = 29.40938 / 2
x = 17/4 adversary: chi2 = 30.24613 / 2
```

This is one of the cleanest pieces of evidence that the global matching-odd signal is not a small-geometry artifact.

The same block also produced an important protocol correction in the matching-even sector. The registered source amplitude was `either/even`, while the target scorer measured `cross/even`. The exact finite-torus map is

```text
DeltaS_cross = -DeltaS_either.
```

The literal registered score fails badly; the no-refit post-reveal exact-map repair gives

```text
corrected cross/even: chi2 = 0.57003 / 2
```

The lesson is methodological as well as physical: **observable typing is part of the experiment.**

### 2. Norm-5 harmonic discriminator

The N325/N425 prospective child block gives

```text
H4:   chi2 =  0.4163 / 2
H12:  chi2 = 35.1931 / 2
H8:   chi2 = 16.0120 / 2
zero: chi2 =  1.7764 / 2
```

This resolves the frozen H4/H12/H8 transfer alias in favor of H4. It is a harmonic/transfer discrimination result, not an independent nonzero-effect detection.

### 3. The full curve does not collapse to one multiplier

The held-out N145→N290 production block rejects the three-level one-multiplier thermal-even curve law:

```text
thermal-even DeltaM transfer: chi2 = 9.3520 / 2, p = 0.009316
```

The rejection is driven by a resolved shape mode. At the same time:

```text
bare slope 2^(3/8):          z = -22.690
frozen scalar+H4 correction: z =  -0.666
root-ratio tests:            compatible
P4[D']:                       z ≈ -0.009
P4[S']:                       z ≈  2.695
```

So the useful conclusion is not “the H4 picture failed.” The stronger conclusion is that **the global signal survives while a one-scalar finite-size state does not.**

---

## What changed the project most

### Scalar corrections repeatedly fail

The repository has accumulated several independent ways of saying the same thing:

- `P4[S'] ~ N^-5/4` fails prospectively on N185/N265;
- the N145→N290 curve is not one scalar multiplier;
- a scalar rank-gap width does not close the higher Krawtchouk/Hermite thermal jet;
- a constant rank-gap correction fails strongly;
- N100 three-modulus `A_top/E_top` shape cannot be represented by one common affine scalar shape;
- even two translated copies of a common symmetric positive kernel are obstructed by the measured higher moments.

The live finite-size object is therefore better thought of as a **multicomponent response/state problem**, not another free correction exponent.

### State dimension depends on the question being asked

A major conceptual result of the later work is that there is no useful context-free number called “the state dimension.”

Under different observer/generator languages, the same microscopic system can exhibit very different realization complexity:

```text
spatial translation / endpoint spectrum
Gaussian or cover composition
occupation-growth continuation
ordered intervention / branching
thermal/source response
```

P250 supplies an especially sharp warning: the complete raw spatial endpoint series has an exact spatial-rank lower bound of at least 100 in the studied setting. A compatible low finite-window Hankel model therefore cannot be read as a continuum field count.

### Full future traces need not close branching

The continuation program around P334/#429 produced a second sharp state-space result: two states can have the same complete unbranched survival law and still have different delayed-fork/branching continuation.

The finite object that naturally retains the missing information is not another scalar hazard. After cutting a rank-one torus along an occupied essential cycle, continuation becomes an exact two-terminal **vertex-connectivity network**. The pair-trigger layer is provably bipartite in the supported embedded-NN scope, while genuine three-site cooperation survives. Scalar counts such as `H2`, pair count, `W2`, `c3`, ... are projections of this network-valued state.

### Jordan-like behavior does not require Jordan

Finite positive-process controls around P398 show that ray inversion / slow-tail effects can occur in ordinary reversible positive mixtures. Removing stationary current changes long-time contrast but does not remove the inversion. Retaining the full instantaneous current direction still does not recover propagation; hidden reversible-force geometry remains.

This is why the project now separates

```text
Jordan algebra exists
Jordan/log fits a finite observable
Matching One physically realizes that Jordan module
```

as three different statements.

### Statistical correctness can change scientific interpretation

Two examples are now part of the scientific history, not merely software history:

1. the N185/N265 matching-even source/target channel mismatch described above;
2. generalized chi-square nullspace QA: a residual can lie in a discarded zero-variance direction and be incorrectly assigned `chi2=0` by a naive pseudoinverse/cutoff recipe.

Historical rescoring largely preserved the displayed numbers, but at least one result became explicitly cutoff-sensitive in interpretation. The repository therefore treats covariance support, chronology and observable semantics as part of the evidence contract.

---

## Current frontier

**Live priority source of truth is the Issue label/state, not this README.** The snapshot below is current as of 2026-09-01.

### P0 — original observable identifiability

[#275](https://github.com/LightChainr/Matching-One/issues/275) is the current P0. The task is to close one chain without changing observable semantics midstream:

```text
raw q/E or trace coordinates
        ↓
physical normalizer
        ↓
pooled moving-root U
        ↓
two candidate forward-prediction vectors
        ↓
covariance-weighted profile rank
```

The key point is that a common spin label is not a map between observables.

A primitive-C3 finite subgate has already produced a deliberately narrow negative result under the signed-real contract:

```text
pure H4:  chi2 = 73.641 / 1
H8 alias: chi2 =  1.112 / 1
```

This rejects pure H4 **for that finite observer/contract** and does not overturn the global-channel H4 evidence. The P0 question is precisely how these observer-specific statements map, or fail to map, to the original normalized `U`.

### P1 — proof-level thermal/contact asymptotics

[#537](https://github.com/LightChainr/Matching-One/issues/537) is now proof-driven rather than “fit another size.” Finite N25/N65/N145 gates are complete; N145 is formally unresolved and is not being topped up.

The remaining route requires proof or counterexample for:

```text
contractible-collar quotient identity
bounded normalized pivotal domination
near-critical uniform transport: exact p_c -> pooled root
```

A third-size fit is not a substitute for those statements.

Other P1 programs include the norm-4/source question (#154), projective-birth/global-transmission question (#334), and intrinsic homology-source program (#337). Their completed finite results remain useful, but they are no longer automatic production queues.

---

## Research landscape

### Global matching / finite-size response

The best-established empirical lane: matching-odd orientation response, Gaussian lineages, norm-5 harmonic discrimination, full-curve derivatives, root/slope structure and shape transport.

**Key boundary:** strong finite-size evidence does not uniquely name the continuum operator.

### Original topological source and `U`

The exact ambient-rank source gives canonical `A_top/E_top` coordinates. Later work studies normalizers, pooled roots, source/thermal derivatives, influence functions, sector quotients and the actual statistical reachability of the original `U`.

**Key boundary:** improving an estimator is different from changing the physical observable.

### Predictive state / continuation geometry

P334/#401/#403/#429 develop exact birth clocks, trigger incidence, branching, cut networks, site-collision observables and network-valued continuation states.

**Key boundary:** finite network state is not automatically a bounded-dimensional continuum memory field.

### Positive finite transfer / hidden geometry

P398 studies positive finite generators, rooted/charged retained modules, current deletion, reversible controls and hidden propagation geometry.

**Key boundary:** nonnormal/Jordan-like phenomenology is not sufficient to identify a Jordan block.

### Primitive square-bond / C3 sector

A distinct topology/character program with exact phase arithmetic, reflection nulls, norm-2 production and later multi-character behavior.

**Key boundary:** do not transport its finite observer label directly to square-site `U` without an explicit observable map.

### Rigorous threshold / exact finite algebra

Issue #1 performs bounded exact relation exclusion. Issues #13/#14 build finite terminal-partition/gadget algebra and explicit periodic primal/dual objects.

**Key boundary:** mature finite algebra is not yet a stochastic comparison theorem or a new rigorous square-site threshold bound.

### Computational statistics and experiment design

Covariance-aware scoring, threshold-rank sufficient statistics, synthetic red-team experiments, influence-function analysis, importance-sampling bounds, analytic subtraction and sequential/prequential controls are first-class research components.

**Key boundary:** smaller Monte Carlo variance is useful only if the estimator still targets the same typed physical quantity.

---

<details>
<summary><strong>Bounded exact search for simple threshold relations</strong></summary>

The repository has run a large, explicitly bounded negative search rather than promoting decimal coincidences.

At coefficient height 100, the frozen primitive polynomial counts are

```text
degree 1:          12,175
degree 2:       3,355,121
degree 3:     749,507,743
degree 4: 157,309,446,881
```

The degree-3 family is excluded on all four frozen method intervals.

For degree 4, the complete family has been screened exactly:

```text
Jacobsen interval:      0 surviving root-containing quartics
Mertens p-med interval: 1 surviving quartic
Mertens p-cell interval:15 surviving quartics
Yang-Zhou interval:     0 surviving root-containing quartics
```

The wider-interval survivors are **not** promoted as formulas: the narrower intervals exclude them. Six frozen standard-constant relation families and several lattice-native candidates have also been checked with exact interval/Sturm controls, along with look-elsewhere counts, positive controls and precision-stability audits.

These are bounded exclusions. They do not establish transcendence or exclude more complicated exact representations.

</details>

<details>
<summary><strong>Why the finite terminal algebra matters — and why it is not the answer yet</strong></summary>

The #13/#14 program has built much more than a gadget sketch:

- canonical terminal partitions and exact connectivity corpora;
- port-aware gluing and a 15-state ordered serial monoid;
- exact proof that the seven D4 orbit labels are not a deterministic serial quotient;
- submonoids, subsemigroups, ideals, congruences, Green relations, automorphisms, centralizers and two-sided operator actions;
- an explicit W5 relative-dual/periodic primal-dual object.

This is a mature exact finite-algebra asset. The missing theorem-facing step is a probability comparison, local transformation, stochastic domination/Strassen-type relation, or a precise obstruction showing that the current finite class cannot deliver such a comparison.

No threshold formula follows merely from the existence of the algebra.

</details>

---

## How to read this repository

**Do not infer the scientific frontier from `main` alone.** Integration state and scientific maturity are separate coordinates. Some important results are intentionally open-PR, branch-only or closed-unmerged assets; some merged PRs are exact controls rather than physical evidence.

Recommended entry points:

1. [`docs/RESEARCH-ATLAS.md`](docs/RESEARCH-ATLAS.md) — the broad visibility map, including underexposed negative results, branch-only science, state-space work, estimator/reachability analysis and mature side programs.
2. [`docs/STATUS.md`](docs/STATUS.md) — claim ledger for its stated snapshot.
3. [`docs/RESEARCH-MAP.md`](docs/RESEARCH-MAP.md) — compact scientific track map.
4. [`docs/ROADMAP.md`](docs/ROADMAP.md) — information-gain priorities for its stated snapshot.
5. [`analysis/research_ledger.yaml`](analysis/research_ledger.yaml) — machine-readable work/evidence state.
6. [`analysis/artifact_registry.yaml`](analysis/artifact_registry.yaml) — artifact/navigation registry.
7. [`results/evidence-ledger/latest.md`](results/evidence-ledger/latest.md) — primary-only predictive evidence view.

When reconstructing history, read the relevant Issue/PR comments and result artifacts. An open Issue can contain completed research; a closed PR can have been superseded rather than scientifically abandoned; a derived score on the same random block is not a new independent experiment.

---

## Evidence discipline

The project uses a few rules aggressively because they have already mattered in practice:

1. **Do not rewrite frozen predictions or committed result history.** Errata are append-only.
2. **Do not silently mix observable semantics.** `either`, `cross`, rank, homology line, local contact and source coordinates are not interchangeable labels.
3. **Do not count correlated views of one raw block as independent primary evidence.**
4. **Keep chronology explicit.** Prospective, held-out, post-reveal and exploratory results have different evidential roles.
5. **Check covariance support, not only a pseudoinverse chi-square.** Deterministic/null directions must be respected.
6. **Separate compatibility from identification.** A model can survive because the experiment is underpowered or the nuisance spaces are indistinguishable.
7. **Type every rank/state claim by observer, generator and context.**

Failures, corrections and null results are retained because they define the current model space.

---

## Reproducibility

Production archives preserve sufficient statistics, metadata, batch structure and covariance whenever practical, rather than only final decimals. The same threshold-rank data can support roots, slopes, derivative channels, Krawtchouk/Hermite coordinates, rank-gap observables and selected source/continuation analyses without pretending that those derived views are independent samples.

Run what you changed:

```bash
python3 -m unittest tests.test_<the_thing_you_changed>
```

The whole suite exists and passes, but running it is not a step in the workflow. Run
it when you have reason to think you broke something far away — not before every push.

This project verifies late and deliberately. `GOVERNANCE.md` §0 and §2 give the whole
rule set for exploratory work; the full apparatus — digests, provenance chains,
preregistration, independent implementations — lives in
[`docs/PUBLICATION-CHECKLIST.md`](docs/PUBLICATION-CHECKLIST.md) and applies when
there is a paper. Time spent on assurance is time not spent exploring, and that trade
has gone the wrong way here before.

---

## Nonclaims

Matching One currently does **not** claim:

- an exact or closed-form value of square-site `p_c`;
- transcendence of `p_c`;
- a unique H4/Q4/Jordan continuum identification;
- that every H4-labelled observer is the same physical state;
- that finite Hankel rank is a continuum field count;
- that the P334 continuation network has a bounded-dimensional continuum limit;
- that branch-only/open-PR results are already integrated into `main`;
- that another larger simulation can substitute for an unresolved identifiability or proof problem.

The strongest current statement is narrower and more useful:

> **A robust global matching-odd finite-size signal exists. Simple scalar explanations have repeatedly failed. The frontier is to identify, with typed observables and covariance-aware forward predictions, what state actually carries that signal.**

---

## License

MIT. See [`LICENSE`](LICENSE).

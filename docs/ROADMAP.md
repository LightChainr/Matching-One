# Roadmap

**Updated:** 2026-08-29

This roadmap ranks **research attention by expected information gain**. It is not a permission system: priorities do not lock tasks, prevent parallel analysis, require branch consolidation, or imply that a lower-ranked hypothesis is false.

The repository keeps only its existing evidence-integrity rules: preserve frozen result history, compare compatible observables, and do not count correlated views of one raw block as independent primary evidence.

For the shortest decision view, see [`docs/NEXT-TARGETS.md`](NEXT-TARGETS.md).

## Priority A — highest expected model-space reduction

### A1. Resolve the lower spin-4 competitor and the global selection rule — #257/#262/#250/#264

The exact continuum preflight now contains two serious spin-4 candidates:

```text
V_(2,+/-2):  x=17/4, |spin|=4
Q4 epsilon:   x=21/4, |spin|=4
```

The lower field is more relevant, so the highest-information question is not another exponent fit but whether the physical `Q=1` global matching observable can couple to it.

Useful next outputs include:

- Potts representation / categorical-projector status at `Q=1`;
- exact or controlled zero/nonzero coupling to the global singlet matching observable;
- charged/four-leg positive-control insertion;
- `Q`-velocity or structure-constant derivative only after the generic-`Q` field definition is typed.

A result in either direction is valuable: an exact exclusion strengthens Q4 identification; an allowed overlap forces a more honest mixture/crossover model.

### A2. Infer the smallest predictive state before naming fields — #249/#180/#253/#255

Several current frameworks describe the same residual structure from different coordinates. Use the already archived norm-2, norm-5, N290, score-mode and compatible local blocks to ask one basis-independent question first:

```text
What is the smallest state that predicts held-out cover/radial contexts?
```

Prefer covariance-aware rank, composition and minimal-polynomial tests over adding another scalar correction model.

Treat diagonalizable rank-2, Jordan rank-2, cover-enriched rank-3 and no-low-rank closure as equally live outcomes until the data separate them.

### A3. Calibrate the identification machinery on independent controls — #246/#234/#42/#44/#155/#106

The target model should not be the only place where the method is judged.

High-value controls include:

- the direct triangular-site energy/log-pair program from #246/#234;
- exact square-bond `p_c=1/2` duality control #42;
- C4 self-matching site family and explicit microscopic tangent #44/#155;
- exactly-critical tunable anisotropy / improved-action family #106.

A deliberately tunable `A4(lambda)` or an externally known logarithmic pair is more identifying than another target-only power law.

### A4. Buy an orthogonal prospective discriminator — #205/#154/#159

Current high-specificity production choices probe different axes:

- **#205 same-N norm-5 coalescence:** H4 interpolation, conjugation and Smith/quotient sensitivity without a radial exponent;
- **#154 norm-4 dyadic closure:** exact scale composition and noncyclic deck structure;
- **#159 modulus/shape fingerprint:** continuum-shape information orthogonal to pure radial scaling once the H4-isolating observable is operational.

Run whichever has the best current discrimination-per-cost and a frozen typed target. There is no need to serialize them if compute and analysis capacity allow parallel work.

## Priority B — build the observable/operator bridge

### B1. Matching action on RG / continuum observables — #61/#114/#233

The project still lacks a first-principles statement of what the lattice matching observable becomes in continuum language.

Useful partial results are enough to advance the program:

- an explicit RG tangent map or intertwiner;
- an FK/Potts `Q -> 1` topological-sector formula for the matching observable;
- a bounded-locality defect/interface construction;
- a precise obstruction showing why the available local representation fails.

This line explains selection rules; it should coexist with numerical identification rather than block it.

### B2. Transfer matrix as operator spectroscopy — #120

Use small-width transfer operators to resolve symmetry/topological/Potts sectors and matrix elements. A clean selection-rule zero can be more informative than pushing width solely for another threshold decimal.

The first valuable result need not be a frontier-size computation; a validated small-width sector decomposition is already useful.

### B3. Cross-microscopic dimensionless ratios — #118

Derive metric-free amplitude ratios before measuring them across square-site, exact self-matching, square-bond and tunable exactly-critical controls.

One cross-model invariant is more identifying than several same-model exponent refinements.

## Priority C — threshold-value origin

Operator spectroscopy and the numerical value of `p_c` are distinct research questions. Keep a visible line aimed directly at the microscopic threshold.

### C1. Correlated-hyperedge/self-dual embedding — #123

Search for a local enlarged state space in which the Bernoulli site rule is a slice of a self-dual or integrable critical manifold. A symbolic closure and a rigorous obstruction are both useful outcomes.

### C2. Exact finite-polynomial / critical-manifold structure — #3/#13/#17/#29

Use exact controls to determine whether finite matching polynomials, decorated cells or topological polynomial structure expose a persistent mechanism rather than a decimal coincidence.

### C3. Post-leading annihilator structure — #47

The historical accelerated-root correction remains useful when interpreted as a mechanism discriminator, not merely an improved estimate of `p_c`.

The question is which correction sector survives after the leading H4 contribution is cancelled.

## Priority D — useful parallel analysis and exact side programs

These may continue whenever they are cheap or generate a new information axis:

- #119 multi-`u` / intrinsic-coordinate analysis;
- #227/#245 Boolean noise-semigroup work;
- primitive square-bond KdV/shape work after the already-resolved sign-transfer question;
- exact finite polynomial/Galois, Betti/Euler and reliability programs;
- Q-score / tangent-crossing calculations that share a common typed generic-`Q` observable definition;
- higher-point/rank-3 work when it is tied to a specific observable rather than used as a free rescue of one-insertion data.

Priority D does not mean “stop”; it means the default next unit of attention should first go to work that can eliminate more live mechanisms.

## Completed information that should guide, not constrain, new work

- **N145->290 full curve:** corrected center/slope/root structure survives, but a one-scalar multiplier does not describe the full curve.
- **Norm-5 N325/425:** H4 transfer strongly outperforms H12/H8 aliases; the child block alone does not independently reject zero.
- **Matching-odd primary synthesis:** global zero is strongly disfavored while fixed H4 remains compatible.
- **Current N130/N170 local tangent rows:** a second direction is not resolved by buying more identical samples.
- **Two primitive norm-2 generations:** the sign/phase question is already well tested; further work should target shape/correction identity instead.
- **Deck-character selection rule:** unmarked deck-invariant global observables have exact zero linear response to nontrivial deck characters, so charged information requires covariant/twisted/marked or nonlinear channels.
- **Exact Q4/Jordan module construction:** establishes a real candidate representation but not lattice overlap.
- **Exact `x=17/4` spin-4 competitor and Q-velocity preflight:** establishes a serious alternative that must be handled by selection/overlap, not ignored.

## Low default priority because the next output is mostly repetitive

Still allowed when cheap or strategically useful:

- adding precision to already completed N290 scalar scores;
- more samples of the same N130/N170 tangent rows;
- a third primitive norm-2 generation whose only target is another sign flip;
- another unconstrained scalar width/free-exponent fit on the same P57 block;
- large production that keeps only final scalar summaries rather than reusable sufficient statistics.

A task should move back up immediately if it acquires a new independent observable, a new exact prediction, a much cheaper implementation, or a stronger adversarial comparison.

## Default choice rule

Choose the next work by the ambiguity it can most cheaply change:

```text
x=17/4 vs Q4 / global overlap        -> representation/projector/charged control
Jordan vs ordinary low-rank mixing    -> minimal realization + composition
continuum field vs quotient memory    -> coalescence/deck/modulus controls
method validity                        -> known/tunable exact-critical controls
matching selection mechanism           -> RG/FK/defect/transfer bridge
why p_c has its microscopic value       -> hyperedge/self-dual/exact-threshold track
```

The objective is not to minimize the number of active ideas. It is to maximize how quickly active ideas make different predictions.
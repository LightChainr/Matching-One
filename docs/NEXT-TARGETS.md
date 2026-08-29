# Next Targets

**Updated:** 2026-08-29

This page is the fast decision board for the project. It ranks attention; it does **not** lock tasks, forbid parallel work, require branch consolidation, or turn a lower priority into a rejection.

The project should preserve competing hypotheses until a discriminating observable, exact selection rule, or independent control separates them.

## Three questions that organize the repository

### Q1. Which continuum sector produces the global square-site matching H4 signal?

Current live possibilities include:

- thermal `Q4` descendant / inherited rank-2 Jordan sector at `x=21/4`, spin 4;
- lower four-leg `V_(2,+/-2)` sector at `x=17/4`, spin 4, if the physical `Q=1` observable has nonzero overlap;
- ordinary low-rank multi-field mixing that imitates part of the Jordan phenomenology;
- a cover/defect/topological or combinatorial-map sector not reducible to one local bulk field;
- higher-rank structure if a genuinely basis-invariant rank-2 description fails.

### Q2. Why does the lattice matching observable select that sector?

The missing bridge is between the exact finite matching/complement identity and a continuum/RG/FK/defect observable. `matching-even/odd` is a useful empirical pair-exchange label; a continuum operator grading still needs a derived map, intertwiner, projector, or selection rule.

### Q3. Can any of this explain or structurally constrain the value of square-site `p_c` itself?

Operator identification and threshold-value explanation are related but distinct. The finite-size CFT program may succeed without deriving `p_c`; the exact/self-dual/hyperedge/critical-manifold program should remain visible as a separate line.

## Default attention order

Priority here means **expected reduction of the live model space per unit effort**, not permission.

| Rank | Target | Main issues | Why it is high information | Useful next output |
|---:|---|---|---|---|
| 1 | Resolve the lower spin-4 competitor and its selection rule | #257, #262, #250, #264 | `x=17/4` is lower than `x=21/4`; if it is allowed in the global observable, the Q4 story must account for it | Potts/projector/multiplicity result, exact zero/nonzero coupling, or controlled charged insertion |
| 2 | Infer the smallest predictive state before naming it | #249, #180, #253, #255 | Several frameworks currently describe the same residual structure; existing data can test rank/composition without new production | covariance-aware held-out rank, composition and minimal-polynomial residuals |
| 3 | Obtain an independent positive control | #246, #234, #42, #44, #155, #106 | A known or tunable system can validate the identification machinery without relying on the square-site target | reproduce a known log pair, tune an H4 amplitude, or recover a symmetry-forced sector split |
| 4 | Buy one orthogonal prospective discriminator | #205, #154, #159 | Same-N coalescence, norm-4 composition and modulus response test different axes from another large-`N` radial fit | one preregistered residual vector that sharply separates surviving models |
| 5 | Derive the matching observable in continuum language | #61, #114, #233, #120 | This explains *why* a field is visible or absent, rather than only fitting what is visible | explicit RG/FK/defect/transfer-matrix map or a precise obstruction |
| 6 | Keep the threshold-origin line alive | #123, #3, #13, #17, #29, #47 | This is the line that can address why `p_c` is that microscopic number rather than only its critical corrections | symbolic closure/no-go, self-dual embedding, or bounded exact threshold mechanism |

These ranks can change whenever a cheap exact result, a newly revealed target, or a new control materially changes model separation.

## Immediate next moves that can run in parallel

### Existing-data / exact lane

1. Use #257/#262/#264 to determine the physical `Q=1` status of `V_(2,+/-2)` and whether the global singlet matching observable can couple to it.
2. Run the #249 minimal-realization pass on the already archived norm-2, norm-5, N290 and compatible local/score-mode blocks. Treat #180/#253/#255 as candidate structures inside that common rank question rather than independent evidence streams.
3. Continue #246 as the known-logarithmic-pair calibration. It is especially valuable because it supplies an external answer against which Jordan diagnostics can be tested.
4. Advance #61/#114/#233 in small exact steps: a partial observable formula, a transfer/intertwiner construction, or a clean obstruction is already useful.
5. Keep #123 as a separate high-risk symbolic line for the threshold value itself; it does not need to wait for operator identification.

### Prospective-compute lane

Compute may proceed whenever a frozen target and useful information-per-cost estimate exist. The current high-specificity choices are:

- #205 same-N norm-5 coalescence: angular law + Smith/quotient sensitivity without a radial exponent;
- #154 norm-4 dyadic closure: exact scale composition with noncyclic deck structure;
- #159 modulus/shape fingerprint: orthogonal to pure radial scaling once the observable projector is operational;
- independent exact-critical controls from #42/#44/#106 when they provide a tunable or symmetry-forced comparison.

There is no requirement to wait for one lane to finish before another runs.

## Candidate-discrimination matrix

The project should prefer measurements that make different mechanisms predict qualitatively different outcomes.

| Test axis | Q4/Jordan | `x=17/4` four-leg | ordinary low-rank mixing | cover/topological memory |
|---|---|---|---|---|
| Potts representation / projector | singlet thermal-family expectation | charged/four-leg representation is central | depends on chosen fields | may require twisted/defect labels |
| Q-velocity / confluent tangent | thermal-family velocity; collision/Jordan structure possible | distinct fixed velocity | mixture of fixed velocities | no simple bulk-field velocity required |
| Modulus / torus shape | specific Q4/Jordan shape fingerprint | different primary/leg-field shape | linear combination of field shapes | may retain quotient/map dependence at fixed modulus |
| Marked/charged local insertion | overlap depends on lattice dictionary | natural positive-control channel | can mix both | character/seam observables may be required |
| Norm-4 / cover composition | common nilpotent generator should compose | ordinary scaling character if present | diagonalizable matrix composition | rank may increase only after Smith/deck-sensitive contexts |
| Known exact-critical control | should reproduce a calibrated logarithmic/descendant mechanism | should reproduce a deliberate four-leg insertion | should remain ordinary mixing | arithmetic memory should separate from continuum controls |

A future run is especially valuable when at least two columns predict a sign, zero, rank, representation, or shape that cannot be removed by an amplitude refit.

## How to choose the next expensive target

Use this default scoring logic:

```text
priority ~
    distinct model predictions
  x independence from already-used raw blocks
  x robustness to nuisance amplitudes
  x reusable sufficient statistics
  / compute and engineering cost.
```

Prefer, in order:

1. exact selection rules or zero/nonzero matrix elements;
2. parameter-free ratios, phases, ranks or shape fingerprints;
3. new independent controls;
4. new prospective raw blocks;
5. additional precision on an already-resolved scalar effect.

This is a ranking heuristic, not a prohibition. A low-cost exploratory calculation is still useful when it can reveal a new information axis.

## Repository task organization

Use task state to describe **what kind of information it can add**, not whether it is allowed to run.

Recommended navigation classes:

- **frontier-identification** — directly separates live mechanisms for Q1;
- **observable-bridge** — addresses Q2;
- **positive-control** — validates the method on known/tunable systems;
- **prospective-discriminator** — buys a new independent target;
- **existing-data-synthesis** — reorganizes one or more already-revealed blocks without creating new evidence;
- **threshold-origin** — addresses Q3;
- **exact-side-program** — cheap algebraic/finite-volume work that may delete a mechanism or create a new tool.

No class implies locking, closure, or branch consolidation.

## Priority update rule

Raise a task when it gains a new ability to distinguish live mechanisms, becomes much cheaper, or acquires a strong independent control. Lower its default priority when its next output mainly repeats an already resolved sign/exponent/mean or is a reparameterization of the same raw block.

Do not interpret a priority decrease as scientific rejection. The purpose is to keep attention pointed at the next piece of information that can change the research picture.

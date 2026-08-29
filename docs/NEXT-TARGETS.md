# Next Targets

**Updated:** 2026-08-29

This is the project’s fast decision board. It answers a narrow question:

> **What new piece of information would change the research picture most, given what is already implemented and already known?**

It is not a permission system. Nothing here locks a task, forbids parallel work, requires branch consolidation, or turns lower priority into scientific rejection.

## Read this by lifecycle, not by issue number

One scientific line may simultaneously have an open issue, an old task entry, a mature note, a frozen prediction, tests, a result capsule, an open PR and a later main-branch correction. Therefore:

- an open issue is not evidence that the work is unstarted;
- an open PR is not automatically the research frontier;
- a theory note is not automatically an implemented analysis;
- a frozen prediction plus tests can be much closer to execution than a newer P0 issue;
- a post-reveal reanalysis is useful for mechanism design but is not a new evidence block.

The useful navigation fields are **actual maturity**, **last durable result**, **next missing input**, and **cheapest informative move**.

## Three questions that organize the project

### Q1. Which continuum sector produces the global square-site matching H4 signal?

Live explanations include the thermal `Q4` descendant / inherited rank-2 Jordan sector at `x=21/4`, the lower four-leg `V_(2,+/-2)` sector at `x=17/4` if the physical `Q=1` observable overlaps it, ordinary low-rank multi-field mixing, and cover/defect/topological or combinatorial-map structure.

### Q2. Why does the lattice matching observable select that sector?

The exact finite matching/complement identity is established. What remains incomplete is the bridge from that finite identity to an RG/FK/defect/transfer-matrix object whose selection rules can explain why a continuum field is visible or absent.

### Q3. Can the structure constrain or explain the microscopic value of square-site `p_c` itself?

Operator identification and threshold-value explanation are related but distinct. The finite-size CFT program can succeed without deriving `p_c`; exact/self-dual/hyperedge/critical-manifold work should therefore remain a visible parallel line.

## Current work state: what is actually missing next

| Research line | Actual maturity now | Last durable result | Next missing input | Cheapest informative move |
|---|---|---|---|---|
| **Lower spin-4 representation / selection** — #257/#262/#250/#264 | Exact continuum preflight complete; lattice overlap unresolved | `V_(2,±2)` has `x=17/4`, four legs, `|s|=4`; its `Q` velocity is distinct from thermal Q4 | Physical `Q=1` Potts projector/multiplicity and zero/nonzero global matching matrix element | Do exact representation/projector work or a controlled charged four-leg insertion before another radial fit |
| **Minimal predictive state** — #249/#180/#253/#255 | Matrix/semigroup assets exist, but the context-Hankel minimal-realization program itself is not implemented on `main` | One scalar width and one scalar full-curve multiplier are already insufficient | First covariance-aware partial realization and held-out rank/composition test | Implement the no-new-production pass on archived norm-2/norm-5/N290 blocks |
| **Triangular known-log pair control** — #246/#234 | Exact protocol, tiny oracle and production code ready; scaling control not yet completed | A direct externally defined energy/log-pair observable is available instead of inferring Jordan from square-site data | Small scaling production that demonstrates the diagnostic on the known pair | Run the smallest block that can distinguish rank-2/Jordan behavior from an ordinary control |
| **Same-N norm-5 coalescence** — #205 | Frozen target and regression tests ready; target C nodes unrevealed | Exact amplitude-free H4/H8/H12 interpolation also changes Smith class | Fresh common-field C-node measurement | Variance-only pilot, then the smallest decisive A/B/C block |
| **Norm-4 dyadic closure** — #154 | Variance pilot complete and scorer ready; comparatively expensive | Current budget gives about three-sigma expected q2/Jordan separation; source covariance dominates | Decide whether this is still the best information/CPU purchase under the enlarged model set | Rerun acquisition ranking before committing the large production budget |
| **Modular-scalar shape spectroscopy** — #103/#159/#220 | Exact modular/Q4 assets are strong; new modular-scalar channel classification/typing makes a clean route possible | `cross/either` are modular scalars; old primitive C3/Pell simple-H4 bridge failed and is a separate sector | Typed square-site modular-scalar observable with a parameter-free shape score | Use one cross/either scalar on one informative new modulus; do not recycle the primitive C3 result as Q4 evidence |
| **Matching observable → continuum bridge** — #61/#114/#233/#120 | Substantial theory plus an exact microscopic self-matching tangent; continuum identification still missing | The UV checkerboard family gives an exact tangent involution, but finite matching exchange is not yet a local CFT/OPE involution | RG/FK/defect/TM intertwiner or a precise obstruction that predicts selection | Small exact or transfer-matrix positive control before large square-site calculations |
| **Metric-free cross-lattice ratios** — #118 | Derivation substantially complete | Declared ratios can cancel the thermal metric factor; raw H4 amplitude ratios are not automatically universal | Same typed observable on a second microscopic realization at matched modulus | Reuse exact-critical control data when available |
| **Threshold-origin line** — #123/#3/#13/#17/#29/#47 | Several bounded exact no-go / finite-polynomial results already complete | Ordinary finite independent-bond gadgets cannot exactly reproduce a Bernoulli four-terminal site in the declared class; simple persistent finite factors are disfavored | Correlated-hyperedge/self-dual closure or obstruction; independent post-annihilator correction structure | Continue symbolic/exact work in parallel; no need to wait for Q1/Q2 |
| **Acquisition optimization** — #102/#126 | The first information/CPU optimizer and sequential-stopping calibration are already implemented | The old optimizer correctly selected norm-5 for the then-live H4/H12 question; sequential stopping showed that question could often finish far below fixed budget | Utility model updated to today’s representation/modulus/rank/control questions | Re-score candidate experiments by expected model-space reduction and stopping cost |
| **Lifecycle synchronization** | Navigation debt, not a new scientific program | Canonical views disagree on some lifecycle states; e.g. N145→290 is completed in STATUS while the primary evidence ledger still says `PENDING_REVEAL` | Reconciled generated status without rewriting historical preregistration | Regenerate status from committed scored artifacts and keep frozen history immutable |

## What I would look at first today

This is a default attention order, not a queue:

1. **Cheap exact selection work around `x=17/4`.** It can remove an entire competitor without Monte Carlo, or force the Q4 story to explain a lower field.
2. **A known positive control (#246).** Before arguing about Jordan geometry in an unknown square-site observable, make the same diagnostics work on a logarithmic pair whose lattice definition is externally known.
3. **#205 coalescence if its variance pilot is favorable.** It is unusually clean because it removes radial exponent, parent amplitude and `p_c` while simultaneously probing H4 and quotient/Smith sensitivity.
4. **Implement the #249 minimal-state pass on existing data.** This should decide whether several current descriptions are one rank question before buying another state variable.
5. **Use the new modular-scalar channel route for one genuinely typed shape test.** Keep it separate from the primitive C3/Pell line that already falsified the simple H4 bridge.
6. **Run norm-4 only after an updated acquisition calculation says the expected information justifies its high source/target cost.** It is ready, not mandatory.

The threshold-origin line can run in parallel with all six.

## Candidate-discrimination view

| New information axis | Q4 / Jordan | `x=17/4` four-leg | ordinary low-rank mixing | cover / topological memory |
|---|---|---|---|---|
| Potts representation / projector | thermal/singlet-family expectation | charged/four-leg representation is central | depends on chosen fields | may require twisted/defect labels |
| `Q` velocity after a valid field continuation | thermal-family velocity; collision geometry can be non-semisimple | distinct exact velocity | mixture of fixed velocities | no single bulk-field velocity required |
| Modulus / torus shape | specific Q4/Jordan fingerprint | different primary/leg-field shape | linear combination of shapes | may retain quotient/map dependence at fixed modulus |
| Marked/charged insertion | overlap is a lattice-dictionary question | natural positive channel | can excite several fields | character/seam observables may be essential |
| Cover composition / norm 4 | one common nilpotent generator should compose | ordinary scaling character if present | diagonalizable matrix composition | rank may rise only after Smith/deck-sensitive contexts |
| Known exact-critical control | should reproduce the claimed logarithmic/descendant diagnostic | can deliberately insert a four-leg sector | should look like ordinary mixing | arithmetic memory should separate from continuum control |

Prefer measurements where two live columns predict a different **zero, sign, representation, rank, phase or shape**, not merely slightly different fitted exponents.

## Reuse the optimizer, do not rebuild it

The Gaussian information-per-CPU machinery already exists and successfully selected the norm-5 experiment for the earlier model set. The next version should change its utility function, not restart the methodology.

A useful current objective is roughly

```text
expected reduction of live mechanism classes
× independence from already-used raw blocks
× robustness to nuisance amplitudes
× reuse value of saved sufficient statistics
÷ expected compute under a valid stopping rule.
```

The design space should now include representation/projector tests, positive controls, modulus tests and rank-revealing measurements—not only Gaussian radial/harmonic alternatives.

## Synchronization rule

Do not hand-edit a historical frozen artifact to make the dashboards agree. When a scored block is complete but a generated ledger still says pending, fix or rerun the lifecycle-generation path and preserve the old artifact in history.

The aim is a trustworthy answer to:

> **What is already known, what is merely implemented, what is genuinely unrevealed, and what single datum would move each line next?**

That answer should be easier to obtain than opening dozens of issues or reconstructing chronology from commit history.

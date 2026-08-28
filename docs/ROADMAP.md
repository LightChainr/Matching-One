# Roadmap

This roadmap optimizes for information gain, not issue-count completion. The integrated research view is `notes/SYNTHESIS-20260828.md`.

## Tier A — active discriminators

### A1. Norm-5 H4 versus H12 — #57

This is now the highest-information dedicated experiment.

The current odd central law survives the prospective N=185/265 test, but `1+i` cannot distinguish H4 from H12/H20. Norm 5 does: the frozen N=325/425 children give opposite raw signs for H4 and H12.

Execution:

1. run a small threshold-rank pilot (start at 1M per child, >=100 batches);
2. estimate actual SE and throughput;
3. choose final production count from measured H4-vs-H12/zero power;
4. reuse the same full curves for derivative and root/correction scores.

Do not assume billions of replicas before the pilot.

### A2. Third full-curve lineage — #50

Score `145 -> 290` with the already-frozen correction structure.

Report in order:

```text
raw DeltaM semigroup law
bare 2^(3/8) slope baseline
frozen finite-size slope correction
raw/corrected root-ratio targets
derivative channels
```

This tests whether the small slope correction resolved in the first two clean norm-2 lineages predicts a third lineage.

### A3. Even/derivative correction analysis — #48

The prospective N=185/265 result falsified the old positive `P4[S] ~ N^-1` sign assignment and the pure `P4[S'] ~ N^-5/4` law.

Before new compute, use existing N=65..265 curves to:

- reconstruct signed normalized S/D derivative sequences under the corrected P4 convention;
- verify the S sign reversal directly from raw primal/matching curves;
- test simple finite-size crossing/correction forms for S;
- retain q=2 and log/Jordan as the first S-prime correction models;
- examine cross-channel covariance.

Reuse A1/A2 as fresh data before commissioning a dedicated new run.

## Completed major prospective gate

### N=185/265 — #43

Completed and closed.

```text
DeltaM x=21/4 H4: 3.046 / 2
DeltaM zero:       29.409 / 2
DeltaM x=17/4:     30.246 / 2

DeltaS positive N^-1: 240.247 / 2
```

Interpretation: the odd law survives and strengthens; the simple even companion law fails in sign. Do not rerun these sizes merely to rescue the old conjunction.

## Tier B — useful parallel work

### B1. Axis annihilator / q=3 scalar test — #47 / PR #98

Use the axis-specific threshold-rank engine and exact adjacent-size common-random-number coupling to test fixed post-leading corrections `q=2,3,4,6,...` without consuming orientation-program compute.

Pilot first; production only if the coupling actually improves information per wall time.

### B2. Exact complex-zero frontier — #124 / #84

Keep this bounded and exact. The tiny zero map made prospective L=5 predictions; the L=5 result already falsifies the simple imaginary-RMS extrapolation. Treat this as an exact side result, not a new CFT narrative unless a stronger invariant emerges.

### B3. Paired same-N motif controls — #40

Keep only if they materially improve the actual orientation-difference statistic.

### B4. Literature completion — #4

Finish missing primary-source transcription when convenient. It does not block current numerical work.

## Tier C — theory after the next discriminators

### C1. LCFT/operator identification — #37/#61

The `x=21/4` spin-4 candidate remains the leading interpretation of the robust odd central sector. Unique H4 content and matching/OPE parity must be established before upgrading the operator claim.

### C2. Post-leading correction spectrum — #47/#58

Use actual residual/annihilator evidence to distinguish q=2,3,4,6. Do not let a named operator determine the exponent before the data.

### C3. N=1105 four-angle decomposition

The exact H0/H4/H8/H12 projector and minimality result are already on `main`. Run the expensive four-angle campaign only after norm 5 and the third full-curve lineage have clarified the cheaper sectors.

### C4. kappa3 / continuum bridge — #25/#54

Independent universality track; useful but not on the critical path of the orientation mechanism.

## Engineering policy

### Merge useful research quickly

C0-C2 notes/tools/results may live on `main` when they are understandable and tests run. Evidence labels, not branch exile, control claim strength.

### Test proportionally

- important CLI: smoke/import coverage;
- central transforms: regression tests;
- topology/RNG/rank engines: exact or independent oracle when practical.

Do not pursue production-service coverage targets.

### Spend compute on leverage

The synthetic red-team shows that new multipliers/sizes are more useful for mechanism selection than simply shrinking old five-size SEs. CPU first; use GPU only for a measured end-to-end information gain.

## Deprioritized

- more N=185/265 replicas for the failed conjunction;
- another free-exponent fit on N=65..170;
- N=1105 before norm 5;
- broad PSLQ;
- large Pell scans;
- generalized infrastructure work without an active experiment;
- theory notes that do not produce a distinct exact test or frozen prediction.

## Milestone after A1/A2

Rewrite the synthesis around whichever branch survives:

1. norm-5 H4 passes and third-lineage correction predicts -> move toward paper-level odd-sector mechanism analysis;
2. H12 sign wins -> rebuild the angular interpretation while retaining the empirical odd semigroup structure;
3. odd law survives but slope/root correction fails -> expand the finite-size thermal-metric model;
4. even-sector sign/correction becomes coherent -> rebuild a corrected multi-sector theory;
5. no compact even-sector model survives -> present the robust odd sector separately rather than forcing a unified operator story.

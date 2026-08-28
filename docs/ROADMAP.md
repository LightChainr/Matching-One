# Roadmap

This roadmap optimizes for information gain, not issue-count completion. The integrated research view is `notes/SYNTHESIS-20260828.md`.

## Tier A — active discriminators

### A1. Norm-5 H4 versus H12 — #57

This is the highest-information dedicated experiment.

The central odd law survives the prospective N=185/265 test, but `1+i` cannot distinguish H4 from H12/H20. Norm 5 does: the frozen N=325/425 children give opposite raw signs for H4 and H12.

Execution:

1. run a small threshold-rank pilot, starting at 1M per child with >=100 batches;
2. estimate actual SE and throughput;
3. choose final production count from measured H4-vs-H12/zero power;
4. reuse the same full curves for derivative, root, and S-prime correction scores.

Do not assume billions of replicas before the pilot.

### A2. Third full-curve lineage — #50

Score `145 -> 290` with the already-frozen correction structure.

Report in order:

```text
raw DeltaM semigroup law
bare 2^(3/8) slope baseline
frozen finite-size slope correction
raw/corrected root-ratio targets
P48 derivative channels
```

This tests whether the small slope correction resolved in the first two clean norm-2 lineages predicts a third lineage.

### A3. S-prime correction and coordinate/channel analysis — #48

The prospective N=185/265 intrinsic-center P48 score is now:

```text
P4[S]   ~ N^-1:     1.139 / 2
P4[D]   ~ N^-13/8:  0.281 / 2
P4[D']  ~ N^-5/8:   0.088 / 2
P4[S']  ~ N^-5/4:  52.716 / 2
```

So `S`, `D`, and `D'` pure laws survive. The unresolved pure-law failure is `S'`.

Before new compute:

- retain q=2 and log/Jordan as the first S-prime correction models;
- assemble the signed S-prime sequence and covariance across available sizes;
- compare fixed `p_ref`, intrinsic-center, and derivative transport;
- preserve the Issue #43 cross/either erratum as a separate wrapping-channel contract problem;
- reuse A1/A2 for new multiplier/lineage leverage.

Do not schedule a dedicated #48 production run yet.

### A4. Use existing full curves harder

Parallel zero/low-compute analyses with direct model-reduction value:

- prequential evidence ledger #95;
- pivotal/Russo bridge #100;
- intrinsic quantile-center spectroscopy #101;
- multi-u thermal response #119;
- joint operator-mixing matrix #125.

These may proceed alongside A1/A2 and should inform any later expensive geometry.

### A5. Information-optimal next geometry — #102

After A1/A2 and the zero-compute analyses, choose any later Gaussian target by expected discrimination gain rather than by size alone.

## Completed major prospective gate

### N=185/265 — #43

Completed and closed.

Matching-odd result:

```text
DeltaM x=21/4 H4-like: 3.046 / 2
DeltaM zero:            29.409 / 2
DeltaM x=17/4:          30.246 / 2
```

Registered matching-even artifact and erratum:

```text
frozen source channel: either/even
target channel:        cross/even
original registered score: 240.247 / 2
exact map: DeltaS_cross = -DeltaS_either
channel-corrected score: 0.570 / 2
```

The original failed score remains historical provenance. The corrected score uses no target refit and restores compatibility after the exact source-to-target channel map.

The same new geometries independently give the prospective intrinsic-center P48 result above: `S`, `D`, and `D'` survive, while pure `S'` fails.

Do not rerun N=185/265 for current questions.

## Tier B — useful parallel work

### B1. Axis annihilator / q=3 scalar test — #47

The axis-specific threshold-rank engine and exact adjacent-size CRN coupling are canonical. Pilot first; run production only if the coupling improves information per wall time.

### B2. Exact complex-zero frontier — #124 / #84

Keep this bounded and exact. The tiny zero map made prospective L=5 predictions; the L=5 result falsifies the simple imaginary-RMS extrapolation. Treat this as an exact side result, not a new CFT narrative unless a stronger invariant emerges.

### B3. Paired same-N motif controls — #40

Keep only if they materially improve the actual orientation-difference statistic.

### B4. Literature completion — #4

Finish missing primary-source transcription when convenient. It does not block current numerical work.

## Tier C — theory and orthogonal controls

### C1. LCFT/operator identification — #37/#61

The `x=21/4` spin-4 candidate remains the leading interpretation of the robust odd central sector. Unique H4 content and matching/OPE parity must be established before upgrading the operator claim.

The prospective survival of intrinsic-center `S`, `D`, and `D'` strengthens the empirical parity decomposition but still does not prove a local OPE/interchiral automorphism.

### C2. Post-leading correction spectrum — #47/#58

Use actual residual/annihilator evidence to distinguish q=2,3,4,6. Do not let a named operator determine the exponent before the data.

### C3. Orthogonal mechanism probes

Useful routes already tracked in the frontier program include:

- torus-modulus spectroscopy #103;
- isoradial/star-triangle anisotropy controls #106;
- Euler/Betti configuration-level matching identities #111;
- FK/Potts torus-sector derivation #114;
- exact N=10 Beta-family falsification at N=26 #115;
- universal amplitude ratios #118;
- four-arm/pivotal anisotropy #121.

### C4. N=1105 four-angle decomposition

The exact H0/H4/H8/H12 projector and minimality result are on `main`. Run the expensive four-angle campaign only after norm 5 and the third full-curve lineage have clarified the cheaper sectors.

### C5. kappa3 / continuum bridge — #25/#54

Independent universality track; useful but not on the critical path of the orientation mechanism.

## Engineering policy

### Merge useful research quickly

C0-C2 notes/tools/results may live on `main` when they are understandable and tests run. Evidence labels, not branch exile, control claim strength.

### Test proportionally

- important CLI: smoke/import coverage;
- central transforms: regression tests;
- topology/RNG/rank engines: exact or independent oracle when practical;
- frozen source/target contracts: explicit observable/channel checks so a scorer cannot silently change the statistic.

Do not pursue production-service coverage targets.

### Spend compute on leverage

The synthetic red-team shows that new multipliers/sizes are more useful for mechanism selection than simply shrinking old five-size SEs. CPU first; use GPU only for measured end-to-end information gain.

## Deprioritized

- more N=185/265 replicas;
- another free-exponent fit on N=65..170;
- treating the Issue #43 `240/2` channel-contract failure as a physical even-sector result;
- treating the S-prime correction as failure of the entire P48 parity spectrum;
- N=1105 before norm 5;
- broad PSLQ;
- large Pell scans;
- generalized infrastructure work without an active experiment;
- theory notes that do not produce a distinct exact test or frozen prediction.

## Milestone after A1/A2

Rewrite the synthesis around whichever branch survives:

1. norm-5 H4 passes and third-lineage correction predicts -> move toward paper-level mechanism analysis;
2. H12 sign wins -> rebuild the angular interpretation while retaining empirical Gaussian-semigroup structure;
3. odd law survives but slope/root correction fails -> expand the finite-size thermal-metric model;
4. norm-5/third-lineage data distinguish q=2 versus log in S-prime -> tighten the derivative/operator interpretation;
5. correction models remain degenerate -> add one geometry only after #102 quantifies which design gives the largest discrimination gain.

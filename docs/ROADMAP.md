# Roadmap

This roadmap optimizes for information gain, not issue-count completion. The integrated research view is `notes/SYNTHESIS-20260828.md`; the relationship among long-lived tracks is in `docs/RESEARCH-MAP.md`.

## Tier A — active discriminators

### A1. Norm-5 H4 versus H12 — #57

This remains the highest-information dedicated new-compute experiment.

The current odd central law survives the prospective N=185/265 test, but `1+i` cannot distinguish H4 from H12/H20. Norm 5 does: the frozen N=325/425 children give opposite raw signs for H4 and H12.

The score path is now protocol-ready on `main`:

- frozen 2026-08-28 fixed-p H4/H12/H8/zero scorer kernel is preserved;
- canonical typed wrapper validates the exact `either/odd -> cross/odd` identity `D_either=D_cross` before scoring;
- the frozen intrinsic full-curve q=2/Jordan cocycle scorer is also preserved;
- its typed wrapper requires the same cross-channel primitive descriptors and size-local P4 normalization at every N.

The numerical freezes were not changed by these semantic gates.

Execution:

1. use a disjoint target production block with the already-frozen geometries and sample plan;
2. score the raw H4 target first, then H12, H8, zero in the registered order;
3. reuse the same full curves for intrinsic-center derivative/correction scores;
4. preserve all rank histograms and covariance groups regardless of which model wins.

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

### A3. S-prime correction and existing full-curve structure — #48/#101/#119/#125

The prospective N=185/265 block now has a canonical no-refit four-channel scorer:

```text
P4[S]   ~ N^-1:     1.13878 / 2
P4[D]   ~ N^-13/8:  0.28085 / 2
P4[D']  ~ N^-5/8:   0.08761 / 2
P4[S']  ~ N^-5/4:  52.71634 / 2
```

So the active problem is not a broad parity failure. S, D, and D-prime survive their frozen pure laws; S-prime is nonzero but requires a subleading correction.

Before commissioning another dedicated geometry:

- use intrinsic quantile-center spectroscopy to isolate nonlinear thermal-coordinate effects;
- score the frozen multi-u functional response as one correlated observable;
- compare q=2 and Jordan through the norm-5 functional cocycle;
- model the four primary channels jointly before adding more basis fields.

## Completed high-value gates and exact controls

### N=185/265 prospective gate — #43

Completed and closed.

```text
DeltaM x=21/4 H4:             3.046 / 2
DeltaM zero:                  29.409 / 2
DeltaM x=17/4:                30.246 / 2

DeltaS unconverted mismatch: 240.247 / 2  (historical invalid channel comparison)
DeltaS corrected cross score: 0.5700315436 / 2
```

`DeltaS_cross = -DeltaS_either` is now an executable typed channel map rather than an informal correction. Do not rerun N=185/265 merely to repeat this gate.

### Finite Russo/pivotal oracle — #100 first stage

The exact finite-volume identity

```text
M'(p) = pivotal_mass_primal(p) + pivotal_mass_matching(1-p)
```

is canonical and independently regressed against exact polynomial derivatives on tiny axis/diamond/Gaussian tori; cross also agrees with threshold ranks. The remaining #100 question is the orientation-resolved/continuum pivotal or four-arm interpretation, not the scalar derivative identity.

### C4 self-matching microscopic tangent — #44/#155 foundation

The two-sublattice family has exact complement action `(t,lambda)->(-t,-lambda)` and N=10 response

```text
[[0,    0],
 [15/8, 5/4]]
```

for rows `(Rplus,Rminus)` and columns `(t,lambda)`. Independent likelihood-score derivatives reproduce it exactly. The open task is the large-N RG decomposition of this microscopic odd tangent, not construction of the finite involution.

### N=26 self-matching Beta test — #115

Completed exact negative result. Exhaustive `2^26` enumeration falsifies both frozen laws without a rescue fit:

```text
Beta(5,5): first k=5 coefficient difference = -96
Beta(7,7): first k=5 coefficient difference = +156
```

The associated reliability-signature/minimal-degree theorem explains why the N=10 Beta(3,3) control is algebraically special. Do not fit a generalized Beta family to N=26 as a replacement hypothesis; use the committed activation signature/profile if more finite information is needed.

## Tier B — useful parallel work

### B1. Prequential evidence ledger — #95

Canonical on `main` through #184. Derived views from the same histograms are one correlated evidence block, not independent votes. Pairwise model comparisons should use only common predicted endpoints and include source+target covariance and predictive log score. Remaining work is to add later prospective blocks to the same ledger, not to rebuild it.

### B2. Information-optimal Gaussian design — #102

Use existing variance/throughput and exact harmonic arithmetic to choose any post-#57 expensive geometry by expected model-space reduction per CPU-second. New multiplier/modulus leverage is more valuable than automatically shrinking old five-size error bars.

### B3. Axis annihilator / post-leading correction spectrum — #47/#58

Use the axis-specific threshold-rank engine and exact adjacent-size coupling only after a pilot demonstrates useful information per wall time. Keep q=2,3,4,6 as mechanism-labeled fixed alternatives before a free exponent.

### B4. Exact finite algebra / complex zeros — #84/#104/#113/#144

Keep these bounded and exact. The global complex-zero imaginary-RMS forecast already failed at exact L=5. Continue only with local-root structure, Galois/reliability certificates, or a surface-topological polynomial representation that produces a new exact discriminator.

### B5. Literature completion — #4/#6

Finish the missing primary-source transcription before revealing the frozen widths 22–24 challenge. This is valuable threshold methodology but does not block the orientation/operator program.

## Tier C — continuum/operator bridges

### C1. LCFT/operator identification — #37/#61

The `x=21/4`, spin-4 thermal level-4 quasiprimary remains the leading local candidate. The ordinary Virasoro state and its `g2/E4` torus one-point fingerprint are exact within that module, but lattice overlap and matching/RG parity are not proved.

The exact self-matching tangent supplies a better finite starting point for #61: derive or measure how a microscopic pair-exchange-odd tangent decomposes into thermal and irrelevant continuum directions rather than assuming a scalar parity sign for every field.

### C2. Pivotal/four-arm and FK/Potts bridges — #100/#114/#121

The scalar Russo identity is now exact. The next useful step is a genuinely orientation-resolved pivotal/four-arm H4 observable or an FK/Potts torus-sector derivation that predicts a continuum channel independently of exponent matching.

### C3. Torus modulus / CM-isogeny spectroscopy — #103/#138/#145/#159

Use multiplier composition and modulus dependence as representation-level fingerprints, but keep Pell production gated until the lattice observable is an explicit typed/projected object. Do not infer a lattice Q4 coupling solely from the exact ordinary-CFT E4 zero.

### C4. kappa3 / universal threshold profile — #16/#25/#54/#122

The new reliability-signature interpretation gives exact finite baselines and shows why low-order central derivatives can miss a categorical finite-profile failure. Treat the full standardized threshold profile/signature as the higher-level object; do not promote a rational-looking `kappa3` without cross-realization and continuum support.

## Engineering policy

### Preserve history; reduce active entry points

- frozen predictions are never rewritten after reveal;
- immutable result archives remain in place even after an interpretation changes;
- superseded PRs/queues are closed or indexed, not force-rewritten;
- canonical current state lives in `STATUS`, `RESEARCH-MAP`, `SYNTHESIS`, `ROADMAP`, and the two analysis ledgers.

### Channel semantics are executable

No scorer should compare channel-bearing quantities unless descriptors are identical or a registered exact map exists. The #43 even correction and both #57 target scorer families now use this contract. Broader direction/both audit remains under #146.

### Spend compute on leverage

CPU first. GPU or very large campaigns require a frozen information target, exact CPU agreement, and a measured variance/throughput reason.

## Deprioritized

- more N=185/265 replicas merely to repeat the corrected no-fit gate;
- another free-exponent fit on N=65..170;
- generalized Beta fitting after the exact N=26 failure;
- N=1105 before norm 5 and the third full-curve lineage;
- broad PSLQ;
- large Pell scans without an operational H4 projector;
- transfer-matrix width solely for more threshold digits;
- theory notes that do not create an independent derivation, control, or frozen prediction.

## Milestone after A1/A2

Rewrite the main synthesis around whichever branch survives:

1. norm-5 H4 passes and the third-lineage correction predicts -> move toward paper-level odd-sector mechanism analysis;
2. H12 sign wins -> retain the empirical odd semigroup structure but rebuild the angular/operator interpretation;
3. H4 survives but the radial/slope correction fails -> expand only the finite-size thermal-metric model;
4. q=2/Jordan becomes distinguishable through the linked full-curve score -> update the derivative mechanism without changing the central evidence;
5. no compact derivative-correction model survives -> present the robust central sectors separately rather than forcing a unified operator story.

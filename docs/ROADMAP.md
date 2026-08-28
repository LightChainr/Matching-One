# Roadmap

This roadmap is intentionally execution-oriented. The project should optimize for information gain, not for completing every open issue in dependency order.

The current integrated synthesis is `notes/SYNTHESIS-20260828.md`.

## Tier A — run next

### A1. Full-curve Gaussian semigroup triptych — #49/#50

Use threshold-rank data on

```text
65 -> 130
85 -> 170
145 -> 290
```

to test, in exact lineage order,

```text
DeltaM_(2N) / DeltaM_N       = -2^(-13/8)
mean(M')_(2N) / mean(M')_N   =  2^(3/8)
DeltaRoot_(2N) / DeltaRoot_N = -1/4
```

This is the highest-value extension of the three successful fixed-`p` Gaussian lineages because it closes the mechanism at the root level.

Use the available covariance diagnostics; do not delay the experiment to perfect a general covariance library.

### A2. Norm-5 H4 versus H12 — #57

Run the frozen moderate-size Gaussian multiplier design that gives sharply different H4/H12 predictions.

The primary question is simple: does the child signal have the H4-predicted sign/magnitude or the H12 alternative? This is more informative than another free-exponent fit.

### A3. Exact C4 self-matching parity control — #44

Implement/run the self-matching C4 site triangulation control. Test whether generic square anisotropy survives while the matching-odd central residual vanishes at the exact self-matching center.

If this works, continue with square-bond #42 and derivative-parity #48.

## Tier B — useful parallel work

### B1. Prospective unused sizes — #43

Score frozen `N=185,265` predictions. Useful confirmation, but less discriminating than A1/A2 because it mostly adds more points to the existing model.

### B2. Paired same-N motif controls — #40

Try exact-zero-mean paired controls on the actual orientation-difference statistic. Keep them only if they materially improve information per wall time.

### B3. Literature completion — #4

Finish primary-source transcription of the 2024 Yang–Zhou tables and Jacobsen Reply when the full text is conveniently available. The current literature layer is already usable; do not block numerical work on the missing table transcription.

### B4. Derivative parity — #48

Once clean threshold-rank curves are already being produced, reconstruct the linked `S,D,S',D'` channels. This can become a strong operator discriminator at low marginal simulation cost.

## Tier C — theory after the next discriminator

### C1. LCFT operator identification — #37

Keep the `x=21/4, spin-4` quasiprimary as the leading theory candidate, but do not spend major effort on module taxonomy until A2/A3 tell us whether H4 and matching parity survive direct tests.

### C2. Post-leading correction spectrum — #47/#58

Use annihilator/semigroup residuals to distinguish candidate relative corrections `q=2,3,4,6`. Worth doing once the leading sector is secure enough that a subleading fit is meaningful.

### C3. kappa3 / continuum bridge — #25/#54

Continue as an independent universality track, but not at the expense of resolving the orientation mechanism.

## Engineering policy

### Merge fast

Exploratory scripts, notes, and result archives should enter `main` once they run and are understandable. C0–C2 material is allowed on `main`; the claim level is the warning label.

### Test proportionally

- every important CLI: smoke test or at least compile/import coverage;
- central numerical transforms: regression tests;
- topology/RNG/rank machinery: exact or independent oracle where practical.

Do not pursue product-style exhaustive coverage.

### Reproduce when it matters

#39 clean replay and #46 covariance hardening are precision-quality tasks. Attach them to the quantitative result that needs them instead of globally blocking research.

### Hardware

CPU first. Use GPU/high-memory machines only when they unlock an experiment that is otherwise impractical or provide a measured end-to-end advantage for the actual statistic.

## Deprioritized for now

- broad PSLQ searches;
- N=1105 four-angle runs before the cheaper norm-5 discriminator;
- large Pell scans;
- rewriting already-usable analysis infrastructure for generalized edge cases;
- new theory notes that do not produce a distinct numerical or exact test.

## Milestone after Tier A

After A1/A2/A3, rewrite the synthesis around whichever of these outcomes occurs:

1. H4 + root semigroup + parity all pass — move aggressively toward a paper-level operator/mechanism analysis;
2. odd harmonic survives but H4 fails — rebuild the angular theory around H12/higher harmonics;
3. full-curve root law fails while fixed-`p` law survives — investigate thermal-coordinate/slope contamination;
4. parity control fails — revise the matching-odd operator interpretation while retaining the empirical orientation sector.

Any of those outcomes would substantially clarify the mathematics, which is the point of the next phase.

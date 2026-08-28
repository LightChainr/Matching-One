# Roadmap

This roadmap orders work by information value, reproducibility risk, and dependency. It is deliberately narrower than the full issue list.

## Phase 0 — Repository control and integration

### G0. Governance baseline

- Merge README, governance, contribution, reproducibility, security, templates, and CI.
- Require pull requests and passing CI for `main` through a repository ruleset.
- Keep `docs/STATUS.md` as the canonical claim ledger.

### G1. Unroll the foundational stack

Integrate in ancestry-preserving order:

1. PR #15 into `main`;
2. retarget PR #18 to `main`, verify the reduced diff, then merge;
3. retarget the active server/research layer only after its base is canonical.

Do not squash stacked scientific commits when their identities appear in predictions or result provenance.

### G2. Separate executable logic from bulk evidence

For the PR #21 layer, prefer reviewable boundaries:

- source, tests, protocols, and small exact vectors;
- immutable raw result archive and manifests;
- narrative conclusions and claim-ledger update.

If one historical PR must retain all files, add a compact manifest and generated/raw/source classification before merging.

## Phase 1 — Reproducibility closure

### R0. Cross-size covariance audit

Complete PR #46 and commit the generated audit outputs. Report both full-covariance and diagonal scores when they differ. Do not silently replace the historical analysis.

### R1. Freeze the RNG policy

Resolve #39 before new threshold-rank production:

- choose explicit `N` domain separation, or
- retain deliberate coupling and propagate full cross-size covariance.

The choice must be frozen before model scoring.

### R2. Clean-checkout replay

Re-run the fixed-`p` and threshold-rank confirmation matrices from clean committed source. Record full commit, source hash, executable hash, compiler flags, dirty-tree state, RNG domains, and checksums. Preserve the original runs as historical evidence.

### R3. Canonical literature/data manifest

Complete #4 before further PSLQ or last-digit claims. Represent each published sequence and threshold estimate separately with exact provenance and transcription checks.

## Phase 2 — Prospective empirical tests

### E0. Angular-normalized root amplitude

Run #45 on the primary `N=65,85` sizes using fresh threshold-rank statistics. Score the frozen amplitude before fitting correction models.

### E1. Unused held-out sizes

Execute #43 at `N=185,265` with pilot-powered sample counts and fresh domains. Score frozen `Delta M` and matching-even predictions before using secondary sizes or flexible alternatives.

### E2. Exact parity controls

Prioritize low-cost, discriminating controls:

- square-bond self-dual control (#42);
- C4 self-matching site triangulation (#44).

These test the proposed parity mechanism more directly than simply adding target-lattice sizes.

### E3. Paired motif controls

Test #40 on the actual same-`N` orientation difference. Promote controls only if fresh-sample variance per wall time improves at multiple declared sizes.

## Phase 3 — Theory and universality

### T0. Operator-sector discrimination

Only after Phase 1 and the prospective Phase 2 tests:

- project matching parity, thermal parity, and spin-4 components from full curves;
- compare pure power, logarithmic companion, free exponent, and competing sectors;
- require an independent amplitude/sign/geometry prediction.

This is the gate for upgrading #37 from a candidate to a C4 interpretation.

### T1. Scaling-function controls

Complete same-modulus exact-threshold controls (#25) before interpreting `kappa_3`, higher derivatives, or an analytic profile. Preserve the full standardized profile, not only one ratio.

### T2. Algebraic and rigorous routes

Keep bounded exact-polynomial/GCD work and rigorous gadget bounds as independent tracks. They must state finite search bounds and may not imply transcendence from finite data.

## Resource policy

CPU is the default for current confirmation and analysis work. GPU or high-memory rental requires:

- a validated CPU oracle;
- measured bottlenecks;
- a power/information model;
- a frozen output contract;
- an end-to-end gain, not only kernel throughput.

## Work explicitly deferred

Until the above gates pass, do not prioritize:

- broad PSLQ searches against an unsettled decimal interval;
- a large four-angle campaign solely to rescue the spin-4 model;
- large-Pell fixed-parameter scans;
- GPU production for work already cheap on CPU;
- transfer-matrix frontier extension without published state/memory profiling;
- new conjecture branches that duplicate an existing issue without a sharper discriminator.

## Completion criterion

The project is ready for a paper-oriented release when:

1. canonical source and result archives are integrated to `main`;
2. CI and repository protection are active;
3. literature data have a provenance-complete manifest;
4. clean replay and covariance policy are complete;
5. at least one prospective unused-size or exact-control test has been scored;
6. the claim ledger clearly separates finite-size evidence from asymptotic interpretation.

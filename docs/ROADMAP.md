# Roadmap

This roadmap orders work by information value, reproducibility risk, and dependency. It is deliberately narrower than the full issue list. `docs/STATUS.md` is the claim ledger; this file is the execution order.

## Phase 0 — Repository control and canonicalization

### G0. Governance baseline — integrated

Completed on `main`:

- governance, contribution, reproducibility, security, citation, templates, and CI;
- foundational research program from PR #15;
- exact matching/orientation reference layer from PR #18.

Remaining hosting task: enable the `main` ruleset and require the four CI jobs described in #52. Until then the repository has policy-level protection but not GitHub-enforced branch protection.

### G1. Curate the server archive — #59

PR #21 is a provenance archive, not a branch to merge wholesale into `main`.

Canonicalize it in reviewable layers:

1. production source and tests;
2. compact machine-readable archive manifest;
3. bounded immutable result families with checksums/reproduction contracts;
4. small claim/decision PRs only after their evidence is canonical.

Preserve historical commit/source hashes, failed/null results, and old working-tree runs. Do not rewrite them to match current conclusions.

### G2. Keep stacked research gates explicit

- PR #46 is a draft until covariance contracts are hardened.
- PR #56 is a draft frozen decision layer blocked by #46/#39 and later canonical archive import.
- Obsolete stacked PRs should be closed as superseded only after their retained evidence is reachable from canonical manifests.

## Phase 1 — Reproducibility closure

### R0. Harden cross-size covariance — PR #46

Before merge:

- enforce equal batch/sample weights in the core audit or implement a mathematically correct weighted estimator;
- validate covariance SPD/rank/conditioning in the scorer itself;
- define finite-batch score semantics/calibration rather than treating an estimated-covariance quadratic form as automatically asymptotic chi-square;
- add exact diagonal/correlated/near-singular/unequal-sample synthetic regressions;
- replay the archive from the current server head and regenerate checksums.

Historical full-covariance and diagonal diagnostics remain preserved regardless of the outcome.

### R1. Freeze the RNG and provenance policy — #39

For new threshold-rank production:

- domain-separate distinct `N` by default, or deliberately couple them under a preregistered covariance plan;
- record full git commit, source SHA-256, executable SHA-256, compiler flags, dirty-tree state, seed/counter domains, and batch layout;
- preserve deterministic test vectors and thread-count invariance.

Do not choose the RNG/coupling policy after inspecting the model score it favors.

### R2. Clean-checkout replay — #39

Replay the required fixed-`p` and threshold-rank confirmation matrices from clean committed source. Preserve historical runs as separate artifacts. Compatibility analysis must distinguish scientific changes from provenance-only corrections.

### R3. Canonical literature/data manifest — #4

Complete #4 before further last-digit claims or integer-relation searches. Represent each published sequence and estimate separately with exact provenance and transcription checks. Issue #1 is P2 and blocked by this work.

## Phase 2 — Prospective finite-size tests

### E0. Angular-normalized root amplitude — completed branch evidence

Issue #45 is completed on the server archive: the frozen `A_p=-N^2 DeltaRoot/DeltaCos4` target was tested at `N=65,85`. Treat it as provisional branch evidence until imported under #59; do not rerun merely to improve significance.

### E1. Full-curve Gaussian semigroup/root tests — #49/#50

Use clean threshold-rank sufficient statistics to test, in exact lineage order:

```text
DeltaM_(2N) / DeltaM_N       = -2^(-13/8)
mean(M')_(2N) / mean(M')_N   =  2^(3/8)
DeltaRoot_(2N) / DeltaRoot_N = -1/4
```

Include the third `145 -> 290` lineage. Preserve full covariance and the local closure diagnostic. Pure-power residuals are primary; logarithmic/Jordan alternatives are secondary and preregistered.

### E2. Unused held-out sizes — #43

Score the frozen `N=185,265` predictions before using secondary sizes or fitting extra radial terms. Sample count may be pilot-powered, but target signs/models and stopping rules must remain frozen.

### E3. Paired same-N motif controls — #40

Evaluate exact zero-mean paired controls on the actual orientation-difference target. Promote them only if fresh-sample variance per wall time improves at multiple declared sizes. Single-geometry gains do not satisfy this gate.

### E4. H4 versus higher odd harmonics — #55/#57

Before expensive four-angle `N=1105` work, use the cheaper prospective Gaussian designs that distinguish H4 from H12 by sign/leverage. Score pure H4 first, then declared H12/H8/mixed alternatives.

### E5. Exact parity controls — #42/#44/#48

Use square-bond self-duality, the C4 self-matching site triangulation, and the derivative-parity spectrum to test matching/duality parity rather than merely accumulating target-lattice size points.

These are prerequisites for any C4 operator-level promotion.

## Phase 3 — Theory and correction-spectrum discrimination

### T0. LCFT/operator identification — #37

Only after Phase 1 and the prospective Phase 2 discriminators:

- project matching parity, thermal parity, and lattice-spin sectors from full curves;
- compare pure power, logarithmic companion, free exponent, and competing operators;
- require at least one prediction not used to identify the candidate.

An algebraically allowed `x=21/4, s=4` quasiprimary is not itself evidence that the lattice observable couples uniquely to it.

### T1. Post-leading correction spectrum — #47/#58

Interpret held-out correction exponents structurally rather than by numerical proximity:

- `q=2` / accelerated `w=6`: conditional even-scalar mixing;
- `q=4` / `w=8`: nonlinear `T4*I4^2`, with an H12 sideband expectation;
- `q=6` / `w=10`: next ordinary thermal spin-4 quasiprimary;
- stable `q=3` / `w=7`: evidence for additional/logarithmic/nonminimal structure, not an automatic “next descendant.”

### T2. Scaling-function controls — #25/#54

Complete same-modulus exact-threshold controls before interpreting `kappa_3`, higher derivatives, or an analytic profile. Preserve the full standardized profile and covariance, not one near-rational number.

### T3. Algebraic and rigorous tracks

Bounded exact-polynomial/GCD work and rigorous gadget bounds remain independent tracks. They must state finite search/certification bounds and may not imply transcendence from finite data.

## Resource policy

CPU is the default for current confirmation and analysis work. GPU or high-memory rental requires:

- a validated CPU/exact oracle;
- measured bottlenecks;
- a frozen sufficient-statistic output contract;
- a power/information model for the actual scientific statistic;
- CPU/GPU equality on deterministic vectors where applicable;
- a meaningful end-to-end gain, not only kernel throughput.

## Work explicitly deferred

Until upstream gates pass, do not prioritize:

- broad PSLQ searches against an unsettled decimal interval;
- N=1105 multi-angle production before cheaper H4/H12 discriminators;
- large-Pell fixed-parameter scans;
- GPU production for work already cheap on CPU;
- transfer-matrix frontier extension without measured state/memory profiling;
- new conjecture branches that duplicate an existing issue without a sharper discriminator.

## Paper-oriented release criterion

The project is ready for a paper-oriented research release when:

1. canonical production source and the required result families are integrated to `main` under #59;
2. CI and repository protection are active;
3. literature data have a provenance-complete manifest;
4. clean replay and covariance/RNG policy are complete;
5. at least one prospective unused-size/harmonic/root test and one exact parity control have been scored;
6. the claim ledger clearly separates finite-size C3 evidence from asymptotic/operator C4 interpretation;
7. release manifests make all reported figures/tables reconstructible from canonical artifacts.
# Project Status and Claim Ledger

**Status date:** 2026-08-28

This file is the canonical repository-level summary of what Matching One currently treats as established, observed, provisional, or open. Detailed evidence remains in linked issues, pull requests, reports, and raw result directories. Claim levels are defined in `GOVERNANCE.md`.

## Background and exact structure

| Statement | Level | Status |
|---|---:|---|
| Square-site percolation on `Z^2` has no known closed form for `p_c` | background | Current literature status; the project makes no closed-form claim |
| The square-site and NN+NNN matching-site thresholds satisfy `p_c + p_c_hat = 1` | C5 | Exact structural constraint used throughout the project |
| Square-bond and triangular-site thresholds are `1/2` | C5 | Exact controls |
| A rounded numerical estimate is not a definition of `p_c` | governance | Required project policy |

Published high-precision sequences and estimates differ in their last digits and methods. Issue #4 is the required provenance reconciliation. Until that work is canonical, `constants/pc.yaml` is a historical/reference value file rather than an adjudication of the literature. Issue #1 is explicitly blocked by #4 and may not use a single rounded estimate as the definition of the threshold.

## Canonical finite-size methods on `main`

The following foundations are now integrated and reviewed on `main`:

- leakage-aware finite-size audit and preregistration tooling from PR #15;
- exact finite-torus matching references, tiny exact polynomials, Gaussian-orientation research tooling, correction filters, and associated tests from PR #18;
- project governance, reproducibility policy, contribution templates, security/citation files, and multi-version CI.

These integrations establish methods and exact/reference contracts. They do **not** automatically make every result that exists on a server/research branch canonical.

## Empirical results currently retained as branch/archive evidence

PR #21 and its `server/huawei-analysis-20260828` branch contain the provenance archive for the current numerical campaigns. The evidence is retained and auditable, but the bulk archive has not yet been curated onto `main`; issue #59 governs that import. Therefore the following claims remain **provisional at repository level** even when the underlying finite-size experiment is strong.

| Claim | Level | Evidence | Qualification |
|---|---:|---|---|
| Primitive same-`N` Gaussian tori exhibit a nonzero orientation difference in the matching observable | C3, provisional | PR #21; #22 | Independent-seed finite-size evidence exists; clean-source/canonical replay requirements remain under #39/#59 |
| The sign agrees with `Delta cos(4 theta)` for the tested pairs | C3, provisional | PR #21; #22 | Finite tested geometry set; not a unique harmonic identification |
| A fixed `Delta cos(4 theta) N^(-13/8)` model predicts the declared held-out sizes better than zero effect | C3, provisional | PR #21; #36 | Finite-size held-out result; not proof of the unique asymptotic exponent |
| Tested `cos(8 theta)`, logarithmic, simple power-correction, and free-exponent alternatives do not outperform the fixed H4 model on the current held-out set | C3, provisional | PR #21; #36 | H12 and other odd harmonics require additional prospective discrimination |
| Fresh Gaussian `1+i` doubling data are compatible with the no-fit ratio `DeltaM(2N)/DeltaM(N)=-2^(-13/8)` on two frozen lineages | C3, provisional | PR #21; P37; regression PR #53 | Prospective finite-coordinate test; does not uniquely identify H4 versus higher odd harmonics |
| Finite reconstructed curves satisfy local linear closure between `Delta M`, mean `M'`, and the root gap | C3, provisional | PR #21; #35 | Establishes local conversion for tested finite curves, not the asymptotic root exponent |
| Angular-normalized root amplitude is compatible with the frozen finite-size target at `N=65,85` in the P45 replay | C3, provisional | PR #21; #45 | Needs canonical archive import and broader clean threshold-rank production before asymptotic promotion |
| Wrapping-only GLS provides useful variance reduction | C1 negative result | PR #21 | False for the tested matching-difference channels: configuration-level identities make this route structurally ineffective |
| Single-geometry Euler/local-motif controls pass the multiple-size `>=2x` production gate | C2 negative/mixed result | PR #21; #34 | They do not; paired same-`N` exact controls are the next route (#40) |

The word “provisional” here is a repository-integration/provenance qualification. It is not a statement that all observed finite-size signals are statistically weak.

## Interpretations that remain open

| Interpretation | Current level | Required discriminator |
|---|---:|---|
| The unique asymptotic radial exponent is `13/8` in `N` (`13/4` in length) | C2/C0 | Clean replay, additional frozen held-out sizes, full covariance, competing correction models |
| The root gap follows a stable angular-normalized `N^-2` law | C2 | High-stat full-curve threshold-rank doubling/root-ratio tests (#49/#50) |
| The leading angular signal is uniquely H4 rather than H12 or another odd harmonic | C2/C0 | Prospective norm-5 / orthogonal harmonic tests (#55/#57) |
| The signal is generated by an `x=21/4`, spin-4 thermal-family LCFT operator | C0/C2 | Thermal/matching parity projection, competing-operator exclusion, log/no-log test, independent predictions (#37/#48) |
| A two-spin-4 parity decomposition explains matching-even and matching-odd sectors | C0/C2 | Prospective held-out sizes and exact self-dual/self-matching controls (#42–#44) |
| The post-leading correction is associated with `q=2,3,4`, or `6` | C0 | Held-out annihilator/semigroup residual spectroscopy (#47/#58) |
| `kappa_3 = -5/3` is universal or exact | C0 | Same-modulus exact-threshold controls and held-out finite-size analysis (#25/#54) |
| A simple algebraic/PSLQ expression gives `p_c` | C0 | Provenance-complete interval and bounded exclusion/search protocol (#1/#4) |

None of these interpretations should be described as established in the README, releases, abstracts, or issue titles without a dedicated claim-ledger upgrade.

## Engineering and provenance status

### Canonical branch

`main` is now the reviewed integration line and includes the governance baseline plus PR #15 and PR #18. The previous “minimal main / alternate research history” state has therefore been partially resolved.

The latest integrated `main` CI is configured for:

- Python 3.9;
- Python 3.11;
- Python 3.13;
- C++17 compilation and self-tests.

The current checks pass on the latest integrated main commit. However, GitHub branch protection/ruleset enforcement is still not active; issue #52 records the required hosting-side configuration.

### Server archive and stacked research work

- PR #21 is the large historical/provenance server archive. It must **not** be bulk-merged into `main`; #59 defines a curated code/tests → manifest → bounded raw-result import strategy.
- PR #46 is a draft statistical-infrastructure PR. Its empirical covariance archive is retained, but merge is blocked on equal-batch-weight enforcement in the core audit, covariance numerical-stability/rank contracts, finite-batch score semantics/calibration, stronger synthetic regressions, and replay onto the current server head.
- PR #56 is a draft frozen research-decision layer. It remains blocked by #46 and #39; its prediction artifacts must not be retuned using later outcomes.
- PR #53 has been merged into the server archive branch and regression-tests the fresh P37 Gaussian-doubling score. This does not by itself make the server archive canonical on `main`.

### Known P0 defects or limitations

- #39: historical production runs include source/provenance and RNG-domain issues that require clean-checkout replay and an explicit cross-size RNG policy.
- #46: cross-size covariance infrastructure must fail closed on violated batch-weight/conditioning contracts before it is trusted as a generic scorer.
- #52: `main` has policy-level protection but no active hosting ruleset; direct pushes/force-push/deletion and required CI checks are not yet enforced by GitHub settings.
- #59: the large Huawei result archive is not yet represented by canonical, reviewable manifests/import PRs on `main`.
- #4: published threshold data are not yet represented by a provenance-complete canonical manifest.

## Current decision gates

Work should proceed in this order unless a governance PR explicitly changes it:

1. enable the `main` ruleset and required CI checks described in #52;
2. harden PR #46 and freeze/complete the #39 clean-source and RNG-domain policy;
3. curate production source/tests from PR #21 onto `main` under #59, without bulk-merging the archive branch;
4. import bounded raw result families with manifests, checksums, and reproduction contracts under #59;
5. rebase/retarget the small frozen decision layer (#56) only after its dependencies are canonical;
6. run/score full-curve Gaussian doubling and root-ratio tests (#49/#50);
7. score prospective unused sizes (#43) before adding more flexible radial models;
8. run paired exact motif controls (#40) and the cheaper H4-versus-H12 discriminators (#55/#57) before expensive N=1105 work;
9. run derivative-parity and self-dual/self-matching controls (#48/#42/#44) before promoting the LCFT operator interpretation;
10. pursue annihilator/LCFT/frontier/GPU programs only after their explicit upstream gates pass.

## Explicit non-claims

The project currently does not claim:

- an exact decimal value beyond cited numerical methods;
- a closed form for square-site `p_c`;
- proof of asymptotic exponent `13/8` or root exponent `2` in `N`;
- unique identification of H4 over all higher odd harmonics;
- proof of an `x=21/4` LCFT operator identification;
- universality of a rational derivative invariant;
- a successful GPU production gate;
- a new rigorous percolation bound.

## Updating this ledger

A claim upgrade, downgrade, or retraction requires a pull request that links the exact canonical evidence, identifies the old and new level, and states which alternative explanations were tested. Raw evidence remains preserved even when the ledger changes.
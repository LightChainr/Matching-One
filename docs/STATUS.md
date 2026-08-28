# Project Status and Claim Ledger

**Status date:** 2026-08-28

`main` now contains the governance/reference layer **and** the Huawei numerical research archive. Claim strength is determined by evidence, not by whether a result lives on a side branch.

Claim levels are defined in `GOVERNANCE.md`. The execution-facing summary is `notes/SYNTHESIS-20260828.md`.

## Exact/background facts

| Statement | Level | Status |
|---|---:|---|
| Square-site `p_c` has no known closed form | background | Current project/literature position |
| Square-site and NN+NNN matching-site thresholds satisfy `p_c + p_c_hat = 1` | C5 | Exact structural constraint |
| Square-bond and triangular-site thresholds are `1/2` | C5 | Exact controls |
| A rounded decimal is not a definition of `p_c` | governance | Enforced in the literature provenance layer |

Published threshold estimates are stored as method-specific values in `data/literature_threshold_sources.json`. The 2015, Mertens 2022, and Yang–Zhou 2024 values are intentionally not collapsed into one synthetic confidence interval.

## Strong current finite-size evidence

| Claim | Level | Evidence | Current interpretation |
|---|---:|---|---|
| Primitive same-`N` Gaussian tori have a nonzero orientation-dependent matching signal | C3 | P31 | Reproduced with independent seed at five frozen sizes |
| The tested signs agree with `Delta cos(4 theta)` | C3 | P31 | Strong evidence for an odd square-harmonic orientation sector |
| The current five-size data are compatible with `DeltaM ~ DeltaCos4 N^-13/8` | C3 | P31/P32 | Finite-range support; not unique asymptotic proof |
| Gaussian `1+i` doubling gives the frozen sign/radial transformation on two lineages | C3 | P37 | Parameter-free prospective test passes |
| A third prospective `145 -> 290` Gaussian-doubling lineage passes | C3 | P50-A | Strengthens the semigroup/orientation mechanism |
| Local residual-to-root conversion satisfies `-DeltaRoot*mean(M')/DeltaM ~= 1` | C3 | P35 | Root movement is explained by the measured residual/slope locally |
| Angular-normalized root amplitude passes the frozen primary N=65/85 test | C3 | P45 | Independent clean-source support for the root-moving sector |
| Wrapping-only GLS helps variance | C1 negative | server archive | No: tested matching-difference channels are configuration-identical |
| Single-geometry motif controls clear a robust multi-size production gate | C2 mixed/negative | P34 | No; paired same-N controls remain more relevant (#40) |

### Key numbers

P31 independent confirmation:

```text
N=65   z=16.03
N=85   z=11.23
N=130  z=5.22
N=145  z=5.27
N=170  z=2.58
A4 pooled = 0.7885 +/- 0.0352
```

P37 fresh Gaussian doubling:

```text
frozen target = -0.3242098887...
65 -> 130     = -0.31382 +/- 0.0908
85 -> 170     = -0.34095 +/- 0.1118
```

P50 third prospective lineage:

```text
DeltaM_290 observed = -0.000160648 +/- 0.000040542
frozen target       = -0.0001376564 +/- 0.000024997
residual z           = -0.483
```

P45 root amplitude:

```text
A_p(65) = 0.42034 +/- 0.02157
A_p(85) = 0.39495 +/- 0.03078
frozen  = 0.45101 +/- 0.02013
```

## Open interpretations

| Interpretation | Level | What would move it forward |
|---|---:|---|
| Leading radial law is asymptotically `N^-13/8` | C2/C3 | full-curve semigroup/root tests and additional prospective sizes |
| Leading angular harmonic is uniquely H4 | C2 | norm-5 H4-vs-H12 test (#57) and/or orthogonal angular design (#55) |
| Root gap follows a clean `N^-2` semigroup law | C2 | full-curve `DeltaRoot(2N)/DeltaRoot(N)=-1/4` on multiple lineages (#49/#50) |
| `x=21/4`, spin-4 thermal-family LCFT field is the mechanism | C0/C2 | parity controls, derivative spectrum, log/no-log discrimination (#37/#44/#48) |
| Two spin-4 parity sectors explain even/odd matching structure | C0/C2 | exact self-matching/self-dual controls (#42/#44) |
| `kappa3=-5/3` is universal/exact | C0 | same-modulus controls and continuum bridge (#25/#54) |
| A simple algebraic/PSLQ formula gives `p_c` | C0 | bounded post-provenance search (#1); low priority |

## Current execution priorities

If compute/attention is scarce, do these first:

1. **#49/#50 full-curve Gaussian triptych** — residual, slope, root ratio on three lineages.
2. **#57 norm-5 H4 versus H12** — highest-value harmonic discriminator.
3. **#44 exact C4 self-matching parity control** — cheapest direct test of the matching-parity mechanism.

Run #43 prospective `N=185,265`, #40 paired controls, and #48 derivative spectrum in parallel when convenient.

## Reproducibility and covariance scope

#39 is now a P1 precision/reproducibility task, not a global blocker. Historical provenance/RNG coupling matters most for exact pooled amplitudes and covariance-sensitive scores; it does not erase the independently reproduced orientation signs or the clean later P37/P45/P50 runs.

#46 is treated as a research covariance diagnostic. Full numerical-library hardening can follow if a paper-facing conclusion becomes sensitive to it.

## Engineering status

- Huawei result archive and production tools are integrated into `main` via PR #21.
- Literature provenance layer is on `main` via PR #62.
- CI runs Python 3.9/3.11/3.13 plus C++17 build/self-tests.
- `main` is still not hosting-side protected (#52), but this is repository hygiene rather than a scientific blocker.
- Governance is intentionally optimized for a solo exploratory project: ordinary exploratory merges do not require an external reviewer.

## Explicit non-claims

The project does **not** currently claim:

- a closed form or exact new value for square-site `p_c`;
- proof that `13/8` is the unique asymptotic exponent;
- proof that H4 is the unique angular harmonic;
- proof of the `x=21/4` LCFT operator identification;
- universality of `kappa3=-5/3`;
- a rigorous new percolation bound.

These are targets for discriminating experiments, not reasons to keep useful exploratory results off `main`.

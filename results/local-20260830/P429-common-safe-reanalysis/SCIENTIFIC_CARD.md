# Scientific card: P429 common-safe residual allocation

## Status

`POST-REVEAL / ZERO-NEW-SAMPLE / DESCRIPTIVE`

This analysis reuses the locked P429 production block.  It is a model- and
order-dependent predictive allocation of an already observed secondary remainder, not an independent
confirmation.  The initial contract was committed as `fa69c92`; the
outcome-blind support/collinearity amendment and runnable analyzer were
committed as `e696fc0` before the mechanism score was run.

## Question and unit

Question: how much of the clone1/clone2 survival dependence remaining after an
orientation-specific safe common update is predictable from checkpoint H2,
local multicontact proxies, low-redundancy geometry, age, and physical line?

One analysis row is one pre-common checkpoint whose common update was safe,
carrying the two subsequent clone outcomes.  Clone outcomes are not expanded
into separate rows.  Five folds hold out complete `batch % 5` groups in both
orientations and both sizes at once.  Direction covariance is retained inside
each size; no all-size meta-estimate or p-value is formed.

## Fixed nested result

| state | N325 residual (conditional dispersion) | N325 absorbed point estimate | N425 residual (conditional dispersion) | N425 absorbed point estimate | N325 / N425 log-loss gain |
|---|---:|---:|---:|---:|---:|
| intercept only | 0.0016112 (0.0002257) | 0.0% | 0.0009500 (0.0002111) | 0.0% | 0 / 0 |
| H2 | 0.0004764 (0.0002214) | 70.4% | 0.0000202 (0.0001997) | 97.9% | 0.0102253 / 0.0095809 |
| H2 + local multicontact | 0.0004657 (0.0002216) | 71.1% | 0.0000160 (0.0001999) | 98.3% | 0.0102795 / 0.0096198 |
| H2 + local multicontact + geometry | 0.0004311 (0.0002215) | 73.2% | -0.0000102 (0.0001993) | 101.1% | 0.0104707 / 0.0099448 |
| all above + age + physical line | 0.0004315 (0.0002215) | 73.2% | -0.0000073 (0.0001992) | 100.8% | 0.0104594 / 0.0099321 |

The main predictive allocation is therefore H2: its absorbed-fraction point
estimates are about 70% at N325 and 98% at N425.  These ratios are not precise
estimates: they inherit numerator, denominator, and nuisance-fit uncertainty
that the conditional score dispersion does not fully capture.  The local multicontact
block adds only `0.0000107` (conditional dispersion `0.0000043`) at N325 and `0.0000042`
(conditional dispersion `0.0000031`) at N425.  The full geometry block adds another `0.0000345`
(conditional dispersion `0.0000071`) and `0.0000262` (conditional dispersion `0.0000056`).  Age and physical-line
blocks do not improve held-out residual or log loss in the fixed sequence.

The final N325 point estimate is `0.0004315`; this does not justify either
“fully absorbed” or “persistent residual” language without refit-aware
uncertainty.  The final N425 equal-direction point estimate is near zero, but
it averages a positive first-direction residual (`0.0002342`) and a negative
second-direction residual (`-0.0002488`).  Absorption slightly above 100% at
N425 is a noise-scale/model overcorrection, not evidence for negative clone
dependence or exact full absorption.

H2 is the exact pre-common count of one-step rank-two hazard sites, so its
dominant absorption is mainly risk-set accounting rather than a newly
identified memory mechanism.  The shared common update can change the
successor H2, which was not saved.

## Predictive reading

- The cross-fitted `h2_rate` coefficient is negative and stable across all
  folds: more one-step rank-two hazards at the pre-common checkpoint predict
  lower suffix survival.
- Local multicontact/contact-pair proxies add little beyond H2.  They are not
  the exact cooperative two-safe-site count.
- The low-redundancy geometry state adds a small but repeatable held-out
  improvement beyond H2 and local contacts.
- Age and physical `chi4(P*ell)` add no transportable held-out improvement in
  this fixed hierarchy.

These are predictive associations with the pre-common checkpoint.  The
production files do not contain successor geometry, successor H2, or a full
microscopic state identifier.

## Hard boundaries

- `branch_common_safe=1` is applied within each orientation.  Requiring both
  orientations of the same replica to be safe defines a different post-hoc
  estimand and is not substituted here.
- `H2_theta=H2`; `H2_figure8=H2_separate=0`; and
  `checkpoint_b1_safe_count=(N-k0)-H2`.  These are not separate evidence.
- Cross-fitting limits same-row fit bias; it does not create fresh data,
  causal mediation, full-state sufficiency, Markov closure/nonclosure at
  scale, a continuum memory field, or a scale exponent.
- Brier gain equals absorbed residual algebraically under the symmetry-locked
  clone prediction and is only an audit.  Held-out log loss is the separate
  predictive-loss diagnostic.
- Stored batch quantities are conditional-on-fitted-fold-model score
  dispersions using 100 size-batches and the observed within-batch direction
  covariance.  Shared slopes couple the nuisance fits across sizes; zero
  cross-size blocks are conditional bookkeeping, not an estimate of full
  unconditional covariance.  These are not sampling SEs, confidence
  intervals, or multiplicity-adjusted inference.

## Provenance

- Baseline commit: `751f8b384883b3ce92e5efa77c35f45a86afa84d`
- N325 CSV SHA-256:
  `04d06fc6d5eafaaa9ce0c9fcdebb970cf7cd5eb5422855c4e468b426a6b527e9`
- N425 CSV SHA-256:
  `f11af42bd8c61f9b34170645bbe453365d926a4d6d4db255c1c56b17b7a70848`
- Safe rows: N325 `85,430`; N425 `84,730`
- Full machine-readable result: `score.json`
- Human-readable complete report: `REPORT.md`

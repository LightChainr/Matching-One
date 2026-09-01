# P429 common-safe dependence remainder: blocked cross-fit reanalysis

## Outcome first

This is a zero-new-sample, post-reveal predictive allocation on the locked P429
production files.  Five-fold cross-fitting holds out whole batches simultaneously
across both orientations and both sizes.  Each row remains one paired-clone unit.

The tables below report the symmetry-locked held-out residual product, the amount
predictively absorbed relative to the environment-intercept model, and held-out loss.
Direction means use equal weights with the measured within-batch direction covariance;
the two sizes are deliberately not pooled into a single evidence number.

## Analysis rows

| size | at-risk rows | common-safe rows | first / second safe | jointly safe replicas |
|---|---:|---:|---:|---:|
| N325 | 91,182 | 85,430 | 42,623 / 42,807 | 18,270 |
| N425 | 89,374 | 84,730 | 42,427 / 42,303 | 18,005 |

## Parent-secondary reproduction

| size | parent conventional gap (GLS SE) | cross-fit baseline (conditional score dispersion) |
|---|---:|---:|
| N325 | 0.0016098 (0.0002243) | 0.0016112 (0.0002257) |
| N425 | 0.0009283 (0.0002106) | 0.0009500 (0.0002111) |

The environment-level conventional covariances reproduce the locked parent score exactly.
The small size-summary differences above are expected: the parent uses full-sample,
clone-specific means and GLS direction weights, whereas this reanalysis uses fold-held,
symmetry-locked means and equal direction weights.

## N325: fixed nested candidate states

| candidate state | features | residual dependence (conditional dispersion) | predictively absorbed | absorbed fraction point | log-loss gain | Brier gain audit |
|---|---:|---:|---:|---:|---:|---:|
| `intercept_only` | 0 | 0.0016112 (0.0002257) | 0.0000000 | 0.0% | 0.0000000 | 0.0000000 |
| `H2` | 2 | 0.0004764 (0.0002214) | 0.0011348 | 70.4% | 0.0102253 | 0.0011348 |
| `H2_cooperative` | 4 | 0.0004657 (0.0002216) | 0.0011456 | 71.1% | 0.0102795 | 0.0011456 |
| `H2_cooperative_geometry` | 13 | 0.0004311 (0.0002215) | 0.0011801 | 73.2% | 0.0104707 | 0.0011801 |
| `H2_cooperative_geometry_age` | 14 | 0.0004313 (0.0002216) | 0.0011799 | 73.2% | 0.0104689 | 0.0011799 |
| `H2_cooperative_geometry_age_line` | 16 | 0.0004315 (0.0002215) | 0.0011797 | 73.2% | 0.0104594 | 0.0011797 |

### N325: direction audit (point estimates)

| candidate state | first residual | second residual |
|---|---:|---:|
| `intercept_only` | 0.0016229 | 0.0015996 |
| `H2` | 0.0004897 | 0.0004631 |
| `H2_cooperative` | 0.0004826 | 0.0004487 |
| `H2_cooperative_geometry` | 0.0004548 | 0.0004075 |
| `H2_cooperative_geometry_age` | 0.0004543 | 0.0004083 |
| `H2_cooperative_geometry_age_line` | 0.0004539 | 0.0004092 |

### N325: incremental contribution of each added block

| added block endpoint | predictive residual reduction beyond previous state | log-loss gain beyond previous state |
|---|---:|---:|
| `H2` | 0.0011348 (conditional dispersion 0.0000544) | 0.0102253 |
| `H2_cooperative` | 0.0000107 (conditional dispersion 0.0000043) | 0.0000543 |
| `H2_cooperative_geometry` | 0.0000345 (conditional dispersion 0.0000071) | 0.0001912 |
| `H2_cooperative_geometry_age` | -0.0000002 (conditional dispersion 0.0000007) | -0.0000018 |
| `H2_cooperative_geometry_age_line` | -0.0000002 (conditional dispersion 0.0000009) | -0.0000096 |

## N425: fixed nested candidate states

| candidate state | features | residual dependence (conditional dispersion) | predictively absorbed | absorbed fraction point | log-loss gain | Brier gain audit |
|---|---:|---:|---:|---:|---:|---:|
| `intercept_only` | 0 | 0.0009500 (0.0002111) | 0.0000000 | 0.0% | 0.0000000 | 0.0000000 |
| `H2` | 2 | 0.0000202 (0.0001997) | 0.0009299 | 97.9% | 0.0095809 | 0.0009299 |
| `H2_cooperative` | 4 | 0.0000160 (0.0001999) | 0.0009341 | 98.3% | 0.0096198 | 0.0009341 |
| `H2_cooperative_geometry` | 13 | -0.0000102 (0.0001993) | 0.0009603 | 101.1% | 0.0099448 | 0.0009603 |
| `H2_cooperative_geometry_age` | 14 | -0.0000090 (0.0001992) | 0.0009590 | 100.9% | 0.0099356 | 0.0009590 |
| `H2_cooperative_geometry_age_line` | 16 | -0.0000073 (0.0001992) | 0.0009573 | 100.8% | 0.0099321 | 0.0009573 |

### N425: direction audit (point estimates)

| candidate state | first residual | second residual |
|---|---:|---:|
| `intercept_only` | 0.0012556 | 0.0006444 |
| `H2` | 0.0002498 | -0.0002095 |
| `H2_cooperative` | 0.0002452 | -0.0002133 |
| `H2_cooperative_geometry` | 0.0002315 | -0.0002520 |
| `H2_cooperative_geometry_age` | 0.0002333 | -0.0002513 |
| `H2_cooperative_geometry_age_line` | 0.0002342 | -0.0002488 |

### N425: incremental contribution of each added block

| added block endpoint | predictive residual reduction beyond previous state | log-loss gain beyond previous state |
|---|---:|---:|
| `H2` | 0.0009299 (conditional dispersion 0.0000519) | 0.0095809 |
| `H2_cooperative` | 0.0000042 (conditional dispersion 0.0000031) | 0.0000389 |
| `H2_cooperative_geometry` | 0.0000262 (conditional dispersion 0.0000056) | 0.0003250 |
| `H2_cooperative_geometry_age` | -0.0000012 (conditional dispersion 0.0000005) | -0.0000092 |
| `H2_cooperative_geometry_age_line` | -0.0000017 (conditional dispersion 0.0000006) | -0.0000034 |

## Interpretation and boundaries

- The baseline and all candidate scores use one shared prediction for the two suffix
  streams.  The exact audit in `score.json` reports how little this symmetry lock differs
  from the conventional product-of-separate-means covariance.
- `H2_theta` is exactly `H2` in these rows; `H2_figure8` and `H2_separate` are zero.
  `checkpoint_b1_safe_count` is algebraically redundant with `H2`.  None is presented as
  extra evidence.
- `branch_common_safe=1` is filtered within each orientation.  The jointly safe replica
  counts above are overlap diagnostics only; a both-orientations-safe subset is a different
  post-hoc estimand and is not substituted here.
- The cooperative block contains boundary multicontact/contact-pair proxies.  It is not
  the exact microscopic two-step cooperative-pair count.
- H2 is the exact pre-common count of one-step rank-two hazard sites.  Its dominant
  absorption is therefore primarily risk-set accounting, not a newly identified memory
  mechanism; the common update can still change the unobserved successor H2.
- Every covariate is measured before the shared common update.  The analysis therefore
  tests how much pre-update state predicts successor heterogeneity; it does not observe
  the complete successor state.
- Cross-fitting addresses in-sample prediction bias only.  These already revealed rows
  cannot provide independent confirmation, causal mediation, full-state sufficiency,
  Markov closure/nonclosure at scale, a continuum memory field, or a scale exponent.
- Batch quantities are conditional score dispersions, not nuisance-refit sampling SEs or
  confidence intervals.  Shared cross-size slopes induce unconditional coupling that the
  stored zero cross-size blocks do not estimate.  No significance claim is made.
- Absorbed fractions are point estimates without ratio intervals.  Values above 100% and
  negative residuals are model/noise-scale overcorrection, not exact full absorption or
  evidence of negative dependence.  At N425, direction averaging also hides opposite-sign
  first/second residual point estimates, which are shown explicitly above.
- Brier gain is algebraically identical to absorbed dependence under the shared clone
  prediction, so it is only an internal audit.  Log-loss gain is the separate held-out
  predictive-loss diagnostic.

Full feature transforms, fold fits, four-environment covariance matrices, source hashes,
and support audits are in `score.json`.

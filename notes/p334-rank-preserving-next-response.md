# P334: most of the positive next-label Gamma separates safe and birth-changing labels

The positive common next-label response is **not primarily a covariance among
rank-preserving next sites**. A three-way split of the already stored forks
places about 78% of canonical Gamma in quartets with one safe label and one
birth-changing label. Safe-safe pairs have a smaller positive canonical
residual, rather than carrying the main effect.

## Three masks, one original population denominator

A label is called birth-changing if it changes either orientation's checkpoint
rank. Define neither/safe-safe when U,V are both safe, mixed when exactly one
changes a rank, and both when each changes at least one rank. The masks partition
each original quartet and its signed Bhat; no new sites, suffixes or finer
birth-type search is added.

| Full-population canonical Gamma allocation | N325 | N425 |
|---|---:|---:|
| neither / safe-safe | 6.82% +/- 1.81pp | 7.82% +/- 2.47pp |
| mixed | 78.51% +/- 1.67pp | 78.45% +/- 2.15pp |
| both | 14.67% +/- 0.98pp | 13.73% +/- 0.70pp |

These are **signed Gamma shares**, not mean fractions or complete variance
shares. The safe-safe quartet mass is 80.32%/83.15%, so its small Gamma share
is not explained by a small number of safe quartets. Its absolute canonical
Gamma is `9.41865e-5 +/- 2.68289e-5` at N325 and
`7.79850e-5 +/- 2.69450e-5` at N425. The residual is positive at approximately
3.5/2.9 batch SE, appreciably smaller than the full response.

For integrated Gamma the safe-safe contributions are
`2.13931e-6 +/- 1.29037e-6` and `1.57323e-6 +/- 6.59656e-7`, corresponding
to 9.45% +/- 5.20pp and 13.12% +/- 4.73pp. These finer residuals are less
precisely determined. Within the named 01+10 sector, safe-safe supplies about
8% of canonical Gamma, with 3.9/4.5 percentage-point uncertainty; mixed still
supplies approximately 69%/71%. The score preserves all absolute values and
their common errors, not only percentages.

For each original prefix Z, let pi_safe(Z) be the probability that a uniform
next label preserves **both** ranks, and m its conditional response. The
safe-safe target is exactly

```
E_Z[pi_safe(Z)^2 Cov_U(m(U) | Z, safe)].
```

It is not `E[pi_safe Var(m|safe)]`, an unweighted conditional covariance, or
the square of a pooled safety rate. No division by an estimated pi is made.
Writing ps,pb for the two prefix-specific probabilities and Bs,Bb for the
within-class covariances, the mixed target is
`ps*pb*(Bs+Bb+(mu_s-mu_b)(mu_s-mu_b)^T)`. Thus its dominance does **not**
identify a pure between-class term. With exact prefix probabilities the true
within contribution is `B_neither/ps+B_both/pb`, computed before prefix
averaging; only total minus this contribution is the between-class term.
Neither the current mixed share nor a pooled safety-rate normalization is a
measurement of within-safe heterogeneity or a Markov-memory claim.

## Within-orientation response and common-label cross terms

The scorer extends the original paired eleven-vector with the **unscaled**
single-orientation `[F1f,F2f,F1s,F2s]` at p_ref and after integration. It saves
the complete 19x19 B for all prefixes, all nine original cells, every mask,
and the named 01+10 aggregate. In each endpoint's raw 4x4 block,

```
Gamma_pair = (B12+B34-B14-B23)/delta_cos4².
```

The first two terms are within-orientation responses; the latter two are the
common-label cross-orientation (CRN-cross) correction. Full canonical values
in the original H4 normalization are:

| Contribution | N325 mean +/- SE | N425 mean +/- SE |
|---|---:|---:|
| within sum | 0.00144939 +/- 0.00002913 | 0.00101274 +/- 0.00002612 |
| cross correction | -0.00006806 +/- 0.00002902 | -0.00001537 +/- 0.00003207 |
| paired Gamma | 0.00138133 +/- 0.00004462 | 0.00099737 +/- 0.00004001 |

The overall positive direction is consequently not predominantly supplied by
the cross correction. A small total correction does **not** make every
masked cross term small: at N425 the canonical mixed correction is
`-1.06701e-4 +/- 2.10429e-5`, while both contributes
`+9.67337e-5 +/- 4.08821e-6`. Their opposite signs largely cancel in the
unmasked result. The complete covariance retains this dependence.

All raw within-orientation terms above inherit the **paired** mask. They must
not be silently equated to a separate experiment selecting only one
orientation's R0-preserving labels. Furthermore these are next-label
conditional-response covariances, not after-next suffix noise covariances.

## Source, artifact and covariance handoff

Source `e32a85939279b8574278024d647b56d2d1485247` contains the original 1.28M
completed forks. The previous common factor is `e0494fdf`; this scorer appends
the new matrices and Gamma allocations on the same twenty original paired
batches, with the same 1,000-prefix denominator. No new random sample, DP,
path replay, test suite, covariance inversion or PSD clipping is performed.

The machine output is `results/p334-rank-preserving-next-response/score.json`.
Its two `N*.complete_common_factor.json.gz` companions retain the raw per-batch
matrix coefficients and full common covariance as `factor.T @ factor`, rank
at most 19. Thus the forthcoming contact readout can be appended without
reprocessing these forks or assuming independence.

Reproduce the named readout with
`/Users/lc/python-envs/research-py311/bin/python scripts/p334_rank_preserving_next_response.py`.

Scientific card: this result shifts the main Gamma interpretation toward a
safe-versus-birth label contrast, leaves a smaller rank-preserving response
open, and separates within-orientation from common-label cross terms. It does
not establish an asymptotic field, a new state variable or a causal role for
any contact motif; those require their own declared observables.

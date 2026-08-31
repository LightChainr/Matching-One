# Quarter-power width is not a full-profile collapse

**Allowing a free center, width and amplitude does not make the N100 and
N400 D_A rank profiles coincide.** The four fixed standardized rank-step
moments give `chi2=73.7393/4`, nominal `p=3.68e-15`. The profile becomes more
symmetric, and the later peak moves inward and gains relative height.
This is additional finite shape evolution, not simply a wider observation
kernel or a differently chosen thermal exponent.

The orders and peak readouts were fixed at **0e805abb before this comparison
was calculated**. N100/N400 width and unstandardized peaks were already
known; no independent new-data status is claimed for this auxiliary
mechanism analysis. It does not change the rank-variance primary predictions
for the next N900 stream.

## Exact location-scale invariants

For a signed profile S with nonzero area M, define its signed-area mean
mu, positive centered second moment sigma^2, and

\[
\Gamma_r=\frac{\int(p-\mu)^rS(p)\,dp}{M\sigma^r},
\quad r=3,4,5,6.
\]

If two whole profiles differ only by center, positive width and a nonzero
overall amplitude, the four Gamma_r must be identical. No value of a
thermal scaling exponent is needed: each profile supplies its own center
and width. Equality of these four moments would be necessary, not sufficient,
for a common standardized profile. A limiting universal shape with finite-N
corrections is not ruled out by a finite two-area difference.

The rank-step profile is the exact archived `S(p)=f_j` on
`[j/N,(j+1)/N)`. For normalized signed weights `w_j=f_j/sum(f)`,

\[
\mu_S=(E_wJ+1/2)/N,\qquad
\sigma_S^2=(\operatorname{Var}_wJ+1/12)/N^2,
\]

and every centered moment is integrated directly over those bins:

\[
\mu_r=\sum_jw_j\frac{N}{r+1}
\left[\left(\frac{j+1}{N}-\mu_S\right)^{r+1}
-\left(\frac jN-\mu_S\right)^{r+1}\right].
\]

This does not smooth a selected window. Canonical control moments use the
same weights and exact Bernstein/Beta integration; their raw conditional
moments are `(j+1)_r/(N+2)_r`. Signed moments are not assumed to be the
cumulants of a probability distribution. If a signed variance were
nonpositive the declared normalization would be undefined, not repaired
by taking an absolute value.

## The finite rank shape actually changes

| Standardized rank-step moment | N100 | N400 | difference | shared-source SE of difference |
|---|---:|---:|---:|---:|
| Gamma3 | 0.4077231 | 0.2054798 | -0.2022433 | 0.068732 |
| Gamma4 | 1.6271517 | 1.5751699 | -0.0519817 | 0.061799 |
| Gamma5 | 1.1289886 | 0.4849629 | -0.6440257 | 0.21635 |
| Gamma6 | 3.2865657 | 3.2370862 | -0.0494795 | 0.25902 |

The four-dimensional covariance gives **chi2=73.7393/4**. The individual
odd moments visibly decrease, with Gamma3 approximately halved. Marginal
changes in Gamma4 and Gamma6 look small, but they are strongly correlated:
their difference errors have correlation about **0.9955**. The even-moment
relationship changes more precisely than either marginal error suggests.
Similarly the Gamma3/Gamma5 error correlation is about 0.9930.

This is why simply treating the four marginal SEs as independent would
miss the coherent shape change. The primary covariance condition number
is approximately 8702, comfortably computable in double precision; the
full matrix and all LOO vectors are retained. The two odd/even pairs are
interpretations of the already-fixed four-coordinate result, not extra
independent tests or a source-selected new field direction.

Canonical moments give the same qualitative answer, with the corresponding
fixed four-coordinate score `chi2=111.3785/4`. That control shares all data
with the rank-step result and is not added as another independent vote.
The pre-canonical primary establishes that the shape difference is not
solely binomial smoothing.

## Where the shape moves: the later peak rebalances

Peaks of the raw step function are not smooth, stable landmarks. The
frozen auxiliary therefore uses the smooth canonical D_A curve, preserving
its full ordinal sequence of three critical points. For each source it
reports

\[
z_i=(p_i-\mu_C)/\sigma_C,\qquad
h_i=\sigma_C C(p_i)/\int C.
\]

These are also invariant under an exact center/width/amplitude change of
the *canonical* profile. They are not claimed to be literal extrema of the
microcanonical step function.

| Canonical standardized landmark | N100 | N400 | difference +/- SE |
|---|---:|---:|---:|
| first peak position | -0.6488184 | -0.6933478 | -0.0445294 +/- 0.040230 |
| valley position | 0.4022829 | 0.3050979 | -0.0971849 +/- 0.069876 |
| **second peak position** | **1.2602201** | **1.0644501** | **-0.1957700 +/- 0.047476** |
| first peak unit-area height | 0.5965623 | 0.5749856 | -0.0215767 +/- 0.033692 |
| valley unit-area height | 0.0305229 | 0.0144800 | -0.0160429 +/- 0.040506 |
| second peak unit-area height | 0.4323116 | 0.4958209 | +0.0635094 +/- 0.033888 |

The relative second/first peak height rises from `0.7246712` to `0.8623188`,
with difference `0.1376476 +/- 0.0598845`. Thus the clearest geometric
description is **the later peak moves toward the center and becomes more
competitive with the first peak**, consistent with reduced signed skewness.
The valley-height change alone is not resolved. We do not infer a new
number of peaks or claim monotone tail-weight growth from marginal kurtosis.

The six nonredundant standardized peak coordinates give
`chi2=53.0340/6`, nominal `p=1.15e-9`. Height ratios are also saved for
interpretability but excluded from that omnibus statistic, since they
are functions of the three heights.

## Covariance, scope and the next scale

Each of the 200/400 common-batch LOO calculations jointly recomputes signed
normalization, center, width, all four moments and all three extrema.
Integer-Bernstein certificates find exactly three critical points in both
empirical mean canonical profiles, and every LOO retains their curvature
types. Exact empirical root counts are not population root-count proofs.

Sources are independent archives at `7b30648` and `3e01b49`; all operations
within an archive share one dependency block. The output keeps the entire
16-coordinate covariance, not just separate error bars. No window, moment
order or competing peak was selected after looking for significance.

The quarter-power width from `fb1a944e` remains a useful finite-regime
fingerprint, but **width compatibility does not mean profile collapse**.
At N900, the already-named rank variance remains the primary question.
The unchanged Gamma3..Gamma6 and ordinal-peak coordinates can additionally
ask whether the skewness/late-peak redistribution continues, saturates, or
reverses. These are auxiliary shape readouts, not a replacement primary.

## Artifacts / scientific card

Run `python3 scripts/p267_standardized_rank_shape.py`; outputs are
`results/p267-standardized-rank-shape/{score.json,REPORT.md}`. The protocol
is `experiments/p267_standardized_rank_shape_20260831.json`. Focused tests
cover a uniform rank bin, its symmetric Beta kernel, and exact rank-step
invariance under center/width/amplitude changes.

- **Changed mechanism space:** a common standardized finite rank profile
  is too small even after free location, scale and amplitude. Reduced
  asymmetry and later-peak rebalancing are explicit new shape coordinates.
- **Not established:** no new critical exponent, full asymptotic profile,
  field identity, or state count; no independent replication from this reuse.
- **Observer / geometry:** D_A ordinary P4 rank-step response; canonical
  landmarks as explicitly smoothed auxiliaries; homothetic N100/N400 pairs.
- **Dependency:** the same two existing blocks; no new prefix, MC stream,
  server run, or modification of the N900 main target.
- **Next discriminant:** the same auxiliary moment and ordinal shape vector
  on the independently revealed next scale, alongside—not instead of—the
  frozen rank-width predictions.

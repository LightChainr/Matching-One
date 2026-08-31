# A maximal Gaussian core with three positive moment centers

This is a post-reveal exploratory explanation of the existing `N100`/`N400`
`D_A` rank-step profiles. It does not change the independent `N900` width target.
The fixed construction and unused readouts are in
`experiments/p267_max_gaussian_three_center_20260831.json`.

## Definition before numerical evaluation

Normalize the signed rank-step profile by its area and write its centered,
unit-variance coordinate as `u`. Let `m_r` denote its moments, so that
`m_0=1,m_1=0,m_2=1`. This normalization does not assume that every empirical
rank bin is positive. A *positive moment realization* is an additional claim
to be checked, not an assumption about the source.

For a candidate common Gaussian variance `t`, define

\[
q_r(t)=\sum_{j=0}^{\lfloor r/2\rfloor}
\frac{(-t)^j r!}{2^j j!(r-2j)!}m_{r-2j}.
\]

These are the formal moments after inverse Gaussian convolution. Form
`H_3(t)=[q_{i+j}(t)]_{i,j=0}^3`. Seek the first boundary `t_*` of its PSD
region starting at zero. Feasibility is downward closed: if the functional
at `t` is positive on squares of cubics, then the functional at `s<t` is its
Gaussian convolution with variance `t-s`; averaging translated polynomial
squares preserves positivity. Thus this is one interval endpoint, not a scan
over disconnected mixtures. Its variance condition gives `t<=1`.

If `H_2(t_*)` is positive definite and `H_3(t_*)` is flat of rank three, the
truncated moments have a unique positive three-atom realization. Its monic
orthogonal cubic is

\[
P_3(x)=x^3+c_2x^2+c_1x+c_0,\qquad
H_2(c_0,c_1,c_2)^T=-(q_3,q_4,q_5)^T.
\]

The ordered real roots are the centers, and positive weights reproduce
`q_0,q_1,q_2`. Flatness ensures reconstruction through `q_6`. Convolving these
centers with `Normal(0,t_*)` gives the declared maximal-common-variance
three-center moment model through `m_6`. A failure of positivity, flatness,
or moment reproduction stops this construction; it does not trigger another
component count or family.

The moments `m_7,m_8` are then new algebraic predictions of this construction,
not inputs. They are read directly from the same original histograms and
compared using the full common-batch leave-one-out covariance, with the
entire construction rerun inside each replicate. Their status is **unused
moments in a post-reveal analysis**, not independent held-out data.

Three atoms are not three fields. Even exact agreement of finitely many moments
would not identify physical Gaussian components, and a finite-support
rank-step function cannot literally equal an unbounded Gaussian mixture with
positive variance. The purpose is a minimal positive explanation of the
observed changing shape, with a definite next failure direction.

## Numerical result

Both empirical profiles admit the declared realization. The boundary Hankel
matrices have numerical rank three; their lower `H_2` minimum eigenvalues are
`0.08423` and `0.07250`. All 200 `N100` and 400 `N400` leave-one-common-batch-out
replicates retain three positive centers and reproduce moments zero through
six to absolute error below `7.2e-15`. This is a numerical statement about the
empirical moments, not an exact population positivity proof. The normalized
negative-bin masses of the original signed profiles are `0.0002965` and
`0.0053123`, respectively: no bin clipping was performed.

Here `t_*` is in unit-variance `u` coordinates. Uncertainties are common-batch
jackknife SEs, including normalization, boundary finding, roots and weights.

| Descriptor | N100 | N400 |
|---|---:|---:|
| `t_*`, common Gaussian fraction of total variance | 0.0691064 ± 0.0030121 | 0.1023710 ± 0.0105077 |
| first center `u_0` | −1.102085 ± 0.011406 | −1.374032 ± 0.064735 |
| middle center `u_1` | −0.547823 ± 0.015898 | −0.714255 ± 0.045902 |
| last center `u_2` | 1.266493 ± 0.020072 | 1.106942 ± 0.045597 |
| first weight | 0.180589 ± 0.010761 | 0.065362 ± 0.020212 |
| middle weight | 0.462297 ± 0.013205 | 0.518769 ± 0.022566 |
| last weight | 0.357113 ± 0.007235 | 0.415869 ± 0.019993 |

In the original `p` coordinate the centers are
`(0.288475, 0.406446, 0.792611)` at `N100` and
`(0.353953, 0.454922, 0.733631)` at `N400`. The common Gaussian variances are
`0.00313068 ± 0.00011298` and `0.00239753 ± 0.00017723`. Thus its absolute
width decreases, but more slowly than the total profile width, so its share
increases. The first-center weight changes by `−0.115227 ± 0.022898`, while
the middle/last weights change by `+0.056472 ± 0.026146` and
`+0.058755 ± 0.021262`. This gives an interpretable low-moment picture of
**early-shoulder depletion and rebalancing toward the central/late portion**,
rather than just a fitted width exponent.

There is a useful decomposition in the previously discussed quarter-power
coordinate `x=N^(1/4)(p-p_ref)`:

| `x` variance | N100 | N400 |
|---|---:|---:|
| total profile | 0.453023 ± 0.005199 | 0.468401 ± 0.016714 |
| common Gaussian part | 0.031307 ± 0.001130 | 0.047951 ± 0.003545 |
| between-center part | 0.421716 ± 0.005909 | 0.420450 ± 0.019424 |

Within this representation, the between-center variance is almost unchanged
at the point estimates (`−0.001266 ± 0.020303`), while the common Gaussian
part grows by `0.016644 ± 0.003720`. This is an algebraic decomposition with
full covariance, not an independent test of collapse. It explains how an
approximately stable total quarter-coordinate width can coexist with changing
shape and substantial weight redistribution. The positions, weights and the
unused higher moments still change; nothing here asserts full-profile scaling.

## The unused seventh and eighth moments

The three-center predictions do **not** get to refit these moments.

| Scale/order | Observed standardized moment | Prediction from `m_0..m_6` | Observed minus predicted ± SE |
|---|---:|---:|---:|
| N100, 7 | 2.90434349 | 2.93700512 | −0.03266163 ± 0.01016892 |
| N100, 8 | 7.58863277 | 7.62445423 | −0.03582146 ± 0.01281973 |
| N400, 7 | 1.04890902 | 1.00979429 | +0.03911474 ± 0.04229860 |
| N400, 8 | 8.09319784 | 8.05425636 | +0.03894148 ± 0.04281641 |

The joint, covariance-aware reference scores are `chi2=10.36365/2`, nominal
`p=0.005618` at `N100`, and `chi2=3.25593/2`, nominal `p=0.1963` at `N400`.
The much larger uncertainties of the separate observed/predicted moments
must not be treated as independent: their shared-batch covariance leaves
the smaller residual SEs above. The complete covariance and every LOO vector
are retained in `results/p267-max-gaussian-three-center/score.json`.

So this is a **positive explanatory structure through degree six with a small,
resolved higher-degree remainder at N100**, not an exact Gaussian-mixture law.
The N100 model slightly overpredicts both unused moments. N400 does not resolve
that remainder at comparable precision; its nominal nonrejection is not proof
that the Gaussian representation becomes exact. Existence at the truncated
Hankel boundary is itself a generic moment-geometry construction, not evidence
for three microscopic components. “Maximal” refers only to the retained
degree-six constraints, not to Gaussian deconvolvability of the full profile.

## Scientific card

- **Changed mechanism space:** the two rank-step profiles have a unique
  maximal-common-variance, positive three-center degree-six realization.
  Its scale evolution resolves shoulder depletion and separates a growing
  relative Gaussian core from the between-center part of the width.
- **Not proved:** three fields, microscopic Gaussian components, complete
  profile collapse, or a new critical exponent. The N100 unused moments already
  show that the finite-moment approximation is not exact.
- **Observer / sector / geometry:** the same ordinary odd rank-step
  `D_A=P4[A_top](4i)-P4[A_top](2i)` at N100/N400; no new observer selected.
- **Source / dependency:** N100 `7b30648be558df0652a7ff22143cc87ed399d042`,
  N400 `3e01b495b5b637b0070705e37b4137a9a0ef0d8b`; common shape batches within
  each source and independent source scales. Same data as the width/shape notes.
- **Next discriminating quantity:** the unused-moment remainder and whether
  the early-center weight/core-variance evolution persists at another scale.
  These are auxiliary post-reveal descriptors, not changes to the N900 primary.

Construction freeze: `191c20e2`. Reproduce with
`python3 scripts/p267_max_gaussian_three_center.py`.
The focused two-case synthetic check covers an exactly known common-Gaussian
mixture and exact uniform-bin integration through degree eight; it does not
rerun any production simulation. No new Monte Carlo was generated and N900
was not read.

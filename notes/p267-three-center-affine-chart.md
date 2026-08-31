# A source-defined affine shape chart, not three physical components

This continues the fixed degree-six realization in `2a824e96`. It uses only
its already reported N100/N400 values and paper algebra; no old significance
score is recalculated and no partial N900 artifact is read.

## Exact covariance and a nonredundant chart

Let a normalized source profile have mean `mu`, variance `s^2`, maximal
degree-six Gaussian variance `T`, ordered physical centers `c0<c1<c2` and
positive weights `w0,w1,w2`. Under an orientation-preserving affine change
`y=a*p+b`, `a>0`,

\[
T_y=a^2T_p,\quad c_{j,y}=a c_{j,p}+b,\quad w_{j,y}=w_{j,p},
\qquad s_y^2=a^2s_p^2.
\]

This follows at the level of the whole feasible interval, not just the fitted
point: polynomial translation/scaling gives a congruence of the Hankel forms,
and Gaussian deconvolution variance changes by `a^2`. PSD boundaries and the
unique flat three-atom measure therefore transform together. A nonzero
overall profile amplitude cancels when normalizing its signed area.

Four nonredundant, dimensionless coordinates are

\[
\alpha=T/s^2,\qquad w_0,\qquad w_2,\qquad
\rho=\frac{c_1-c_0}{c_2-c_0}\in(0,1).
\]

They reconstruct the entire standardized realization. Put
`w1=1-w0-w2`, `a=(0,rho,1)`, `abar=sum(w*a)`, and
`v_a=sum(w*(a-abar)^2)`. Then

\[
u_j=\sqrt{\frac{1-\alpha}{v_a}}(a_j-\bar a).
\]

Thus these are coordinates for the same four-dimensional standardized
`m3..m6` information on the positive flat-rank3 branch, not extra observables
or extra fields. The normalized between-center variance is exactly
`1-alpha`: it is **not** a fifth independent shape coordinate. Reflection
`a<0` reverses center order, swaps `w0,w2`, and sends `rho` to `1-rho`; “early
weight” is anchored to the declared increasing physical `p` orientation.

The existing point estimates are:

| Coordinate | N100 | N400 |
|---|---:|---:|
| alpha | 0.0691064 | 0.1023710 |
| early weight w0 | 0.1805891 | 0.0653620 |
| late weight w2 | 0.3571135 | 0.4158688 |
| relative middle gap rho | 0.2340064 | 0.2659347 |

These are an exact re-expression of the previously reported realization, not
a new fit or a new significance test.

## Two physically distinguishable changes

Pure shift/dilation of the source leaves all four coordinates unchanged.
Adding a common independent Gaussian blur after an affine map is a different
mechanism. If its variance in the new coordinate is `eta>=0`, then

\[
T_{\rm new}=a^2T+\eta,\quad s_{\rm new}^2=a^2s^2+\eta,
\quad
\alpha_{\rm new}=\frac{a^2s^2\alpha+\eta}{a^2s^2+\eta},
\]

but **`w0,w2,rho` remain unchanged**. All standardized centers contract by
the same factor `sqrt((1-alpha_new)/(1-alpha))`. The exact covariance of the
maximal degree-six boundary makes these statements true of the descriptor,
not merely a heuristic for fitted Gaussians.

Consequently a change of the ordered weights or relative gap identifies
atom-geometry redistribution that cannot be attributed to a common blur and
an affine coordinate change alone. The already reported early-weight shift
therefore has a distinct interpretation from the increase of `alpha`. This
does not elevate the three atoms to physical components: the N100 unused
`m7,m8` already show a remainder beyond this moment representation.

## Quarter-width decomposition

In `x=N^(1/4)(p-p_ref)`, write `Vx=N^(1/2)s^2`. Then exactly within this model,

\[
V_x=G_x+B_x,\qquad G_x=\alpha V_x,\qquad B_x=(1-\alpha)V_x.
\]

`Vx` is an area-dependent physical width coordinate; `alpha,w0,w2,rho` are
affine shape coordinates. One should not call `B_x` affine invariant or count
`B_x` and `alpha` as independent shape information when `Vx` is already known.
The old values `G_x:0.031307->0.047951` and
`B_x:0.421716->0.420450` give a compact explanation of how near-constant total
quarter-width can coexist with changing internal shape. They do not predict
that either part is exactly constant at the next scale.

## One-shot completed-source scorer

Once the root has supplied the **completed, committed** N900 artifact, run
this command in this branch's worktree, replacing `FINAL_SHA` with that commit:

```sh
python3 scripts/p267_max_gaussian_three_center.py \
  --source-commit FINAL_SHA \
  --source-directory results/etop-n900-rank-width \
  --output results/p267-max-gaussian-three-center-n900
```

This path is supplied by the root; it has not been inspected during this
follow-up. The scorer reads the declared `score.json` contract and the existing
`raw/tau_2i.hist.csv`, `raw/tau_4i.hist.csv` archive through the original loader,
then reconstructs all standardized rank-step moments through order eight.
The source's `rank_profile.batch_raw_moments` containing only orders zero
through two is **not enough** and is not substituted. The first two shape
entries must be the same named `2i,4i` order defining `D_A` in the original
schema. Only that one source is scored; N100/N400 are not recomputed.

The callable moment-level entry is already available for a caller holding
the full moments and their LOO arrays:

```python
vector, certificate = features_from_moments(m0_to_m8, mean_p, variance_p, N)
# Apply the same call to every original common-batch LOO moment vector.
covariance = jackknife_covariance(loo_vectors)
```

Here `m0_to_m8` must be the **standardized integrated rank-step moments**,
starting `[1,0,1]`, not canonical Bernstein moments or the raw moments of `K2`.
The function returns the same labels as `LABELS`; it does not change the model.

Output fields are `sources["900"].estimate`, `se`, `full_covariance`,
`leave_one_common_batch_out_vectors`, `realization`, `LOO_construction`,
`unused_moment_readout`, `affine_shape_chart`, and
`quarter_coordinate_width_decomposition`. Failed PSD/flatness/positive-weight
gates are reported without trying another component count. The seventh/eighth
moments remain unused by the construction. This is an auxiliary completed-source
application; the N900 primary width prediction remains unchanged.

No N900 production was started, stopped, read, or modified here.

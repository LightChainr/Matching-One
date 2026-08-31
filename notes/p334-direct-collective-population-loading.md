# Original single-site gates and collective birth load the population differently

The completed paired R1 archive has a concrete source decomposition. At N325,
the original-checkpoint single-site and collective sources have opposite
orientation contrasts; the collective term reverses the integrated contrast.
At N425, their classified point estimates have the same sign. This is a new
readout of the full fixed 20k-counter populations, not another selected-prefix
clock calculation.

## Exact competing-source identity

At the rank-one checkpoint let d be the remaining labels and D the H2=h
original singleton triggers. Let f_j count j-subsets that remain safe and
`S(j)=f_j/binomial(d,j)`. Under a uniform remaining permutation, T is the
first subsequent rank-two insertion. Monotonicity guarantees that every
site in D stays a trigger. The event `T>j-1` therefore leaves all h labels
in D uninserted. The next label is uniform among the d-j+1 remaining labels:

\[
 P(T=j,V_{\rm final}\in D)
 =S(j-1)\frac{h}{d-j+1},\qquad
 P(T=j)=S(j-1)-S(j).
\]

The collective mass is their difference. It includes final sites outside
the original D, even when earlier safe insertions have made those sites
new singleton triggers. In particular, this is not a comparison of the
numbers of single-site and multi-site minimal triggers. Both source terms
depend on the entire physical survival curve, not only on H2.

The implementation divides nonnegative exact integer counts:

\[
 p_D(j)=\frac{h f_{j-1}}{\binom d{j-1}(d-j+1)},\quad
 p_G(j)=\frac{(d-j+1-h)f_{j-1}-j f_j}
 {\binom d{j-1}(d-j+1)}.
\]

For each channel, apply the existing kernels
`P(Binom(N,p_ref)>=k0+j)` and `(N-k0-j+1)/(N+1)`, with unchanged
`p_ref=0.59274605079`. The second is the exact p-integrated canonical kernel.
This does not delete any site or assign a unique knockout-causal effect.

## Preserve the complete paired population

The input is `0d1e586dafbade5e7d1f9bfc598170d0c881e337`,
`results/p334-paired-clock-loading/batches/`. All 20 original batches of
1000 paired counters are retained at each size. Checkpoints stay fixed at
N325:k0=193 and N425:k0=252; no age, line or solvability stratum is selected.

- `exact_pair`: apply the identity to every rank-one orientation's saved
  coefficients. An orientation outside the specified R1 contribution is zero.
- `whole_pair_fallback`: both original Y orientations go entirely into an
  unclassified source U, even when one exact clock was saved before failure.
- `outside_rank_one`: the defined contribution remains zero, without calling
  the absent state rank zero or rank two.

The counts remain 13907/13803 exact pairs, 47/164 whole-pair fallbacks and
6046/6033 both-outside pairs. Thus `D+G+U=Y` for every counter/orientation/readout
as an exact identity. Numerical evaluation differs by at most `7.22e-16` per
row; the saved hybrid batch means are reproduced exactly. No network or DP
was rerun.

## Source contrasts

The table divides first-minus-second by the original difference in cos(4theta).
These are archive point estimates, not sign-established population effects.

| Size/readout | Original H2 direct D | Collective G | Unclassified U | Total |
|---|---:|---:|---:|---:|
| N325 canonical | +0.0008485573 | -0.0002661861 | +0.0000276068 | +0.0006099781 |
| N325 integrated | +0.0010160879 | -0.0015632820 | -0.0000409825 | -0.0005881767 |
| N425 canonical | +0.0008318802 | +0.0002531131 | +0.0000508864 | +0.0011358797 |
| N425 integrated | +0.0016788882 | +0.0017213305 | +0.0000909628 | +0.0034911815 |

Among the classified channels, `1-|D+G|/(|D|+|G|)` is 47.76% for N325's
canonical contrast and 78.79% for its integrated contrast. The integrated
collective direction outweighs the positive direct direction. At N425, the
classified integrated contributions are almost equal and reinforce each
other. No cross-size sign-switch significance is inferred here.

For the positive orientation loadings before subtraction, original direct
gates account for approximately 82% of the classified canonical mass and
66% of the integrated mass. Consequently, which source dominates a small
orientation difference cannot be read from which source dominates the
positive loading itself. Direct gates may carry most mass while collective
loading determines the sign of the contrast, as in N325's integrated result.

Unclassified U is not silently negligible. For either classified source,
its unnormalized first-minus-second contrast has allocation range
`[D_f-D_s-U_s, D_f-D_s+U_f]`; divide by delta_cos4 and sort. This allows any
nonnegative allocation of U in the fixed hybrid estimator, not population
uncertainty. The two source allocations are coupled.

At N325 the direct-positive/collective-negative pattern survives every such
allocation for both readouts. At N425 the canonical direct range remains
positive, but the collective range includes zero; both integrated source
ranges include zero. A small *net* U contrast does not imply a small possible
source-specific allocation uncertainty. Full numerical ranges are in the
score and report; none are confidence intervals.

## Shared covariance handoff and scientific card

- Changes: finite-time physical birth is decomposed by a named checkpoint
  site class across the whole paired archive. N325's two sources oppose each
  other in both kernels, more strongly after integration; N425's classified
  point means reinforce.
- Observer/sector: rank-one checkpoint contribution to canonical F2 and its
  p-integral, original single-site versus remaining collective final birth.
  The H4 normalization is a geometry contrast, not global A_top completion
  or identification of a field.
- Source/dependency: the same N325/N425 20k paired streams and original
  20 batches as the root hybrid and prevalence/conditional-clock analyses.
  All D/G/U columns are jointly dependent with those quantities.
- Lifecycle: saved complete safe coefficients -> exact source partition ->
  original 20-batch loading vectors -> common covariance coordinator. Zero
  new MC, zero new DP and no source substitution or sampling filter.
- Next discriminant: jointly propagate this source partition with the
  prevalence/clock decomposition, instead of combining separate marginal
  errors or treating source signs as already established.

Output `results/p334-direct-collective-population-loading/score.json` has
`labels` ordered as orientation(first,second) x readout(canonical,integrated)
x source(original_H2_direct,collective,unclassified_original_Y).
Each size's `joint_20_batch_means_orientation_readout_source` is the aligned
20-by-12 matrix. `batch_ids` are 0..19, with denominator 1000 per row and
20000 for the mean. The artifact additionally retains the original four
hybrid batch columns, all source hashes, six projected source batch columns
and allocation envelopes. The common covariance coordinator can therefore
join the same batch IDs with its existing risk/C/L columns; no independent
source error bars are manufactured here.

```sh
/Users/lc/python-envs/research-py311/bin/python scripts/p334_direct_collective_population_loading.py
```

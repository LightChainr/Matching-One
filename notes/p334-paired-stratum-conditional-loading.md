# P334: the complete rank-one stratum, with paired conditional loading

The next object is an unconditional, stratum-weighted contribution over the
two saved20k counter populations, not a selected-prefix survival law. The
N325 checkpoint is fixed at k0=193 and N425 at252. In either orientation,
the geometry archive emits a row exactly when the checkpoint is rank one.
Missing rows are not all rank two: they include rank zero.

## Declared observable and denominator

For a fixed size and a uniform insertion order, put

\[
 X_i(p)=1\{R_i(k_0)=1\}\,P\{\mathrm{Bin}(N,p)\ge K_{2,i}\},
 \quad i\in\{\text{first},\text{second}\}.
\]

The sample denominator is all20000 original counters. This is the rank-one
stratum's contribution to canonical F2, **not** the complete F2 or A_top.
Outside this stratum the defined X is zero without classifying the omitted
state. Its integrated counterpart is

\[
 \int_0^1 X_i(p)\,dp
 =1\{R_i(k_0)=1\}\frac{N-K_{2,i}+1}{N+1}.
\]

Both orientations share the same HNF-label permutation and the same k0 at
a given size. Consequently they share the same occupied label set at the
checkpoint, even though their embedded adjacency differs. Conditional on
this common prefix set, the remaining permutation is uniform. All
occupation/age/line strata are retained; there is no narrow age10 or line
filter and no assumption that the147 earlier prefixes represent this pool.

## Complete conditional suffix replacement

For each rank-one orientation, the generalized existing two-terminal
construction gives exact safe counts f_j over the d=N-k0 remaining labels.
With S(j)=f_j/binomial(d,j),

\[
 P(T=j\mid\text{prefix})=S(j-1)-S(j),\qquad K_2=k_0+T.
\]

Substituting this entire law into the two kernels above gives the exact
conditional mean. The only geometry generalization is to the four declared
HNF matrices `[[N,shear],[0,1]]` and their saved primitive line. In physical
coordinates that line is `(N*ell_u+shear*ell_v,ell_v)`; its primitive
transverse covector replaces the old hard-coded `(19,8)`. The connectivity
DP and exact integer polynomial arithmetic are unchanged.

Some prefixes may be too expensive under the declared one-second whole-pair
and12000-state computation budget. The entire paired vector is replaced
only if **all** its rank-one orientations are solved. Otherwise both
orientations retain their original suffix observations. Solved coefficients
from a partially completed pair are retained in the artifact, but are not
substituted into the paired estimator.

## Why whole-pair replacement matters

Let F be the common prefix information and let R select the pairs whose
required conditional expectations were obtained. R uses prefix geometry
and the declared computation budget, not observed K2, the suffix value, or
an effect-size gate. Permutation decoding and baseline scoring occur
outside that budget. With exogenous execution timing included in F, define

\[
 Z=R\,E[X\mid F]+(1-R)X.
\]

Then E[Z]=E[X], and the residual X-Z is orthogonal to Z. In particular,

\[
 \boxed{\operatorname{Cov}(X)=\operatorname{Cov}(Z)
 +E[(X-Z)(X-Z)^\mathsf T].}
\]

This identity holds for the joint canonical/integrated orientation vector,
so every fixed linear contrast inherits non-increasing population variance.
Replacing just one member of a pair when the other solve fails would not
have this common-conditional-expectation form and could lose that guarantee.
Exogenous timing is not a reason to tune budgets after observing responses.

The empirical residual second-moment matrix estimates the removed suffix
noise. It is not forced to equal a finite-sample covariance subtraction;
both the original and hybrid batch covariances remain visible. It is also
not a universal computational speedup or the full conditional noise of
prefixes left in fallback.

## Geometry contrast and the next decomposition

The reported orientation contrast divides `mean(first)-mean(second)` by
the original representatives' difference in cos(4theta). This merely states
the existing geometric normalization. A nonzero rank-one contribution does
not identify a global H4 field, because the other checkpoint strata and
the complementary P0 observable have not been included here.

Let r_i=P(rank_i(k0)=1) and m_i=E[X_i|rank_i(k0)=1]. The contribution can
subsequently be separated, with shared batch covariance, by the identity

\[
 r_fm_f-r_sm_s
 =\tfrac12(r_f-r_s)(m_f+m_s)
 +\tfrac12(r_f+r_s)(m_f-m_s).
\]

This distinguishes stratum prevalence from within-stratum clock loading;
it is a symmetric descriptive decomposition, not causal attribution.

## Execution and source

- Original source: e81dd59ff6be69056e504e0e81cfeccf73dc5e97, N325/N425,
  all20000 counters each and their original20 batches.
- Original shared permutation decoding and exact network solver are reused.
- Policy/code:657e29e5, with decoding outside the prefix budget at5b81e2f9.
- No new random samples or remote job. Four local workers process independent
  original batches; no new widths, windows, or fitting families are selected.
- Per-batch compressed artifacts preserve source rows, original observations,
  exact coefficients when available, pair status and the substituted vector.
- The declared full pool has completed:40 original batches and40000 paired
  counters. All outputs are under `results/p334-paired-clock-loading/`.

## First population-level result

For the H4-normalized orientation difference, the estimated fraction of
individual-observation variance removed by suffix replacement is49.15% at
N325 and50.03% atN425 for the canonical readout. For the integrated readout
it is only0.816% and0.681%, respectively. The corresponding hybrid means
and original-batch standard errors are:

| Size | Canonical contribution | Integrated contribution |
|---|---:|---:|
| N325 | 0.00060998 +/- 0.00156001 | -0.00058818 +/- 0.00315435 |
| N425 | 0.00113588 +/- 0.00118417 | 0.00349118 +/- 0.00220881 |

The full-population noise allocation therefore depends strongly on the
readout kernel. Averaging the suffix can remove about half the canonical
noise without materially reducing the integrated noise. The earlier
equal-weight147-prefix mixture was a different, selected estimand; its
83.95% integrated conditional-noise fraction cannot be transferred here.
The prevalence/conditional-clock decomposition above is the next direct
way to understand what dominates the stratum-weighted integrated signal.

The pair policy produced13907/13803 exact replacements,6046/6033
both-outside-stratum zero vectors, and47/164 whole-pair fallbacks for
N325/N425. No unresolved pair was partially substituted. These estimates
remain archive-derived R1 contributions and do not complete global A_top.
